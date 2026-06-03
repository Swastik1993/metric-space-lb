// Subsample.h
//
// Node-subsample with connectivity check and retry-on-failure (Reviewer
// R1.W3 / R3.D5 + v16 multi-seed variance estimation).
//
// Purpose:
//   Given a full graph G(V, E) loaded by --sparse-graph mode (typically the
//   n=100k HNSW edge file), induce a random node-subgraph of `target_n`
//   vertices. Used to obtain multiple independent subsamples of the same
//   high-dim dataset for variance estimation. Each subsample seed yields a
//   reproducibly-different subset; the experiment runs identically on each.
//
// Algorithm:
//   1. Shuffle [0..full_n) using mt19937(seed).
//   2. Take the first `target_n` IDs as the survivor set; sort ascending for
//      determinism.
//   3. Build the induced subgraph (keep edges where both endpoints survive),
//      remapping IDs to {0..target_n-1}.
//   4. BFS-count connected components.
//   5. If 1 component: accept. Otherwise increment seed by 1 and retry, up
//      to `max_retries` times.
//   6. On persistent failure: hard error (no silent fallback to largest
//      component, since that would change `n` from what the slurm requested
//      and invalidate variance-across-seeds reporting).
//
// Why retry-on-failure: at the proposed (n=100k base, target=50k, M=24)
// configuration the induced graph has avg degree ~8.3 (=33.3 * 0.25), so a
// random subsample is almost always connected. But for safety / smaller M /
// edge-density combinations the retry path protects against rare bad draws.
//
// Reproducibility note: the seed actually used (post-retries) is logged on
// every run, so an aggregator can join across seeds cleanly even when the
// requested seed differs from the accepted seed.

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <list>
#include <queue>
#include <random>
#include <utility>
#include <vector>

namespace subsample {

using AdjList = std::vector<std::list<std::pair<unsigned int, double>>*>;

struct InduceResult {
    AdjList* new_adj_lst = nullptr;  // owned by caller
    unsigned int n_new = 0;          // post-LCC count if LCC fallback fires
    std::size_t m_new = 0;
    unsigned int seed_used = 0;
    unsigned int components_pre = 0; // components BEFORE LCC pruning (== 1 in strict mode)
    unsigned int components = 0;     // components in the returned graph (always 1 on success)
    unsigned int retries_consumed = 0;
    unsigned int lcc_dropped = 0;    // nodes dropped by LCC fallback (0 in strict mode)
    unsigned int requested_n = 0;    // what the caller asked for
    bool success = false;
};

// Count connected components of an undirected adj-list graph (BFS).
inline unsigned int count_components(const AdjList& adj_lst, unsigned int n) {
    if (n == 0) return 0;
    std::vector<char> visited(n, 0);
    unsigned int comps = 0;
    for (unsigned int start = 0; start < n; ++start) {
        if (visited[start]) continue;
        ++comps;
        std::queue<unsigned int> q;
        q.push(start);
        visited[start] = 1;
        while (!q.empty()) {
            unsigned int u = q.front();
            q.pop();
            for (const auto& nb : *adj_lst[u]) {
                unsigned int v = nb.first;
                if (!visited[v]) {
                    visited[v] = 1;
                    q.push(v);
                }
            }
        }
    }
    return comps;
}

// Sum of all adjacency-list sizes is the directed-edge count; undirected is
// half that since every (u,v) edge appears as both u->v and v->u.
inline std::size_t count_undirected_edges(const AdjList& adj_lst,
                                          unsigned int n) {
    std::size_t directed = 0;
    for (unsigned int i = 0; i < n; ++i) directed += adj_lst[i]->size();
    return directed / 2;
}

// Sample `target_n` distinct node IDs from {0..full_n-1} using mt19937(seed).
// Returns the sorted-ascending survivor set so the new_id ordering is
// deterministic given (full_n, target_n, seed).
inline std::vector<unsigned int>
sample_node_ids(unsigned int full_n, unsigned int target_n, unsigned int seed) {
    std::vector<unsigned int> all(full_n);
    for (unsigned int i = 0; i < full_n; ++i) all[i] = i;
    std::mt19937 rng(seed);
    std::shuffle(all.begin(), all.end(), rng);
    all.resize(target_n);
    std::sort(all.begin(), all.end());
    return all;
}

// Build induced subgraph from `full` keeping only edges with both endpoints
// in `survivors`. Remaps IDs to {0..|survivors|-1} in survivor-sorted order.
// Returns a heap-allocated AdjList* compatible with clean_up_adj_list().
inline AdjList* build_induced_subgraph(
    const AdjList* full,
    unsigned int full_n,
    const std::vector<unsigned int>& survivors)
{
    const unsigned int target_n = (unsigned int)survivors.size();

    // old_id -> new_id, or full_n (sentinel: dropped)
    std::vector<unsigned int> remap(full_n, full_n);
    for (unsigned int new_id = 0; new_id < target_n; ++new_id) {
        remap[survivors[new_id]] = new_id;
    }

    auto* out = new AdjList(target_n);
    for (unsigned int new_id = 0; new_id < target_n; ++new_id) {
        (*out)[new_id] = new std::list<std::pair<unsigned int, double>>();
    }

    for (unsigned int new_u = 0; new_u < target_n; ++new_u) {
        unsigned int old_u = survivors[new_u];
        for (const auto& nb : *full->at(old_u)) {
            unsigned int old_v = nb.first;
            double w = nb.second;
            unsigned int new_v = remap[old_v];
            if (new_v == full_n) continue;  // neighbor dropped
            (*out)[new_u]->emplace_back(new_v, w);
        }
    }
    return out;
}

// Free an AdjList* produced by build_induced_subgraph(). Identical to the
// project's clean_up_adj_list(); duplicated here to avoid the include cycle.
inline void free_adj_list(AdjList* adj_lst) {
    if (!adj_lst) return;
    for (auto* l : *adj_lst) delete l;
    delete adj_lst;
}

// ---------------------------------------------------------------------------
// Largest-connected-component (LCC) extraction
// ---------------------------------------------------------------------------
//
// HNSW layer-0 graphs subsampled at 50% retention consistently produce 10-25
// connected components: one giant component covering ~99.97% of nodes plus a
// handful of fringe singletons/pairs/triples from heavily-pruned high-degree
// outliers. Retrying with new seeds doesn't help (the fringe structure is a
// property of the graph, not the sample). The LCC fallback drops those
// fringe nodes and accepts the giant component as the working graph for
// downstream metrics. n drifts by ~20 per seed on the production
// configuration; the aggregator picks this up via the post-LCC `n` reported
// in the ACCEPTED line.

// BFS-label every node's component, return component id per vertex and the
// total component count.
inline void label_components(const AdjList& adj_lst, unsigned int n,
                              std::vector<int>& comp_id,
                              unsigned int& num_components_out) {
    comp_id.assign(n, -1);
    unsigned int comps = 0;
    for (unsigned int start = 0; start < n; ++start) {
        if (comp_id[start] >= 0) continue;
        std::queue<unsigned int> q;
        q.push(start);
        comp_id[start] = (int)comps;
        while (!q.empty()) {
            unsigned int u = q.front();
            q.pop();
            for (const auto& nb : *adj_lst[u]) {
                unsigned int v = nb.first;
                if (comp_id[v] < 0) {
                    comp_id[v] = (int)comps;
                    q.push(v);
                }
            }
        }
        ++comps;
    }
    num_components_out = comps;
}

// Build a new AdjList containing only nodes in the largest connected
// component, remapped to {0..lcc_n-1}. Returns the new graph; sets
// `lcc_n_out` and `num_components_out`.
inline AdjList* extract_lcc(const AdjList* induced, unsigned int n,
                             unsigned int& lcc_n_out,
                             unsigned int& num_components_out) {
    std::vector<int> comp_id;
    label_components(*induced, n, comp_id, num_components_out);

    // Find largest component by counting members.
    std::vector<unsigned int> comp_size(num_components_out, 0);
    for (unsigned int u = 0; u < n; ++u) {
        if (comp_id[u] >= 0) comp_size[(unsigned)comp_id[u]]++;
    }
    unsigned int max_id = 0;
    for (unsigned int c = 1; c < num_components_out; ++c) {
        if (comp_size[c] > comp_size[max_id]) max_id = c;
    }

    // Collect LCC members in original-id order (stable remap).
    std::vector<unsigned int> members;
    members.reserve(comp_size[max_id]);
    for (unsigned int u = 0; u < n; ++u) {
        if (comp_id[u] == (int)max_id) members.push_back(u);
    }
    lcc_n_out = (unsigned int)members.size();

    // old_id -> new_id; sentinel n for "not in LCC".
    std::vector<unsigned int> remap(n, n);
    for (unsigned int new_id = 0; new_id < lcc_n_out; ++new_id) {
        remap[members[new_id]] = new_id;
    }

    auto* out = new AdjList(lcc_n_out);
    for (unsigned int i = 0; i < lcc_n_out; ++i) {
        (*out)[i] = new std::list<std::pair<unsigned int, double>>();
    }
    for (unsigned int old_u = 0; old_u < n; ++old_u) {
        unsigned int new_u = remap[old_u];
        if (new_u == n) continue;
        for (const auto& nb : *induced->at(old_u)) {
            unsigned int new_v = remap[nb.first];
            if (new_v == n) continue;  // safety: should not occur for LCC members
            (*out)[new_u]->emplace_back(new_v, nb.second);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------------
//
// Two modes:
//
//   strict_connectivity = false  (DEFAULT, v16.1)
//     Build one induced subgraph; if it has >1 component, prune to LCC.
//     n_new may be slightly less than target_n (typically ~20 less on the
//     production SIFT config). No retries needed.
//
//   strict_connectivity = true   (legacy v16 behavior, opt-in via
//                                  --subsample-strict-connectivity)
//     Retry with seed+1, seed+2, ... until a single-component sample is
//     drawn; hard-fail after max_retries. Use this if the experimental
//     setup requires the exact requested n.

inline InduceResult induce_node_subgraph(
    const AdjList* full,
    unsigned int full_n,
    unsigned int target_n,
    unsigned int requested_seed,
    unsigned int max_retries,
    bool strict_connectivity = false)
{
    InduceResult res;
    res.seed_used = requested_seed;
    res.requested_n = target_n;

    if (target_n == 0 || target_n > full_n) {
        std::cerr << "[Subsample] ERROR: target_n=" << target_n
                  << " must be in (0, " << full_n << "]" << std::endl;
        return res;
    }
    if (target_n == full_n) {
        std::cerr << "[Subsample] WARNING: target_n == full_n; "
                  << "no subsampling needed." << std::endl;
        return res;
    }

    std::cout << "[Subsample] Requested target_n=" << target_n
              << " from full_n=" << full_n
              << "  seed=" << requested_seed
              << "  strict_connectivity="
              << (strict_connectivity ? "true" : "false");
    if (strict_connectivity) {
        std::cout << "  max_retries=" << max_retries;
    }
    std::cout << std::endl;

    if (strict_connectivity) {
        // ---- Strict mode: retry until single-component or hard fail ----
        for (unsigned int attempt = 0; attempt <= max_retries; ++attempt) {
            unsigned int this_seed = requested_seed + attempt;
            auto survivors = sample_node_ids(full_n, target_n, this_seed);
            AdjList* induced = build_induced_subgraph(full, full_n, survivors);
            std::size_t m = count_undirected_edges(*induced, target_n);
            unsigned int comps = count_components(*induced, target_n);

            if (comps == 1) {
                std::cout << "[Subsample] ACCEPTED  seed_used=" << this_seed
                          << " (attempt " << (attempt + 1) << "/"
                          << (max_retries + 1) << ")"
                          << "  n=" << target_n
                          << "  m=" << m
                          << "  components=" << comps
                          << "  avg_deg="
                          << (target_n ? (2.0 * m / target_n) : 0.0)
                          << std::endl;
                res.new_adj_lst = induced;
                res.n_new = target_n;
                res.m_new = m;
                res.seed_used = this_seed;
                res.components = comps;
                res.components_pre = comps;
                res.retries_consumed = attempt;
                res.success = true;
                return res;
            }

            std::cout << "[Subsample] REJECTED  seed_used=" << this_seed
                      << " (attempt " << (attempt + 1) << "/"
                      << (max_retries + 1) << ")"
                      << "  n=" << target_n
                      << "  m=" << m
                      << "  components=" << comps
                      << "  (strict mode -- retrying with seed+1)"
                      << std::endl;
            free_adj_list(induced);
        }

        std::cerr << "[Subsample] HARD ERROR (strict mode): failed to find "
                  << "a connected " << target_n << "-node subsample within "
                  << (max_retries + 1) << " attempts starting from seed "
                  << requested_seed << ". Drop --subsample-strict-connectivity "
                  << "to use LCC-fallback mode instead." << std::endl;
        return res;
    }

    // ---- LCC fallback mode (default) -------------------------------------
    unsigned int this_seed = requested_seed;
    auto survivors = sample_node_ids(full_n, target_n, this_seed);
    AdjList* induced = build_induced_subgraph(full, full_n, survivors);

    unsigned int comps_pre = 0;
    unsigned int lcc_n = 0;
    AdjList* lcc = extract_lcc(induced, target_n, lcc_n, comps_pre);
    free_adj_list(induced);

    unsigned int dropped = target_n - lcc_n;
    std::size_t m_lcc = count_undirected_edges(*lcc, lcc_n);
    unsigned int comps_post = count_components(*lcc, lcc_n);

    std::cout << "[Subsample] LCC-fallback  seed_used=" << this_seed
              << "  requested_n=" << target_n
              << "  components_pre=" << comps_pre
              << "  lcc_n=" << lcc_n
              << "  dropped=" << dropped
              << "  (giant_frac="
              << (target_n ? (double)lcc_n / (double)target_n : 0.0)
              << ")" << std::endl;

    std::cout << "[Subsample] ACCEPTED  seed_used=" << this_seed
              << "  n=" << lcc_n
              << "  m=" << m_lcc
              << "  components=" << comps_post
              << "  avg_deg=" << (lcc_n ? (2.0 * m_lcc / lcc_n) : 0.0)
              << std::endl;

    res.new_adj_lst = lcc;
    res.n_new = lcc_n;
    res.m_new = m_lcc;
    res.seed_used = this_seed;
    res.components_pre = comps_pre;
    res.components = comps_post;
    res.retries_consumed = 0;
    res.lcc_dropped = dropped;
    res.success = true;
    return res;
}

}  // namespace subsample
