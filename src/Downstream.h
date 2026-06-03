// Downstream.h
//
// Downstream-task evaluation (Reviewer R3.D2): construct a Minimum Spanning
// Tree via Prim's algorithm and a PAM k-medoid clustering using each method's
// lower-bound distances, then compare to the true-distance reference.
//
// Inputs:
//   - The full true-distance matrix `distance` (n*n) from the input dataset.
//     This is what `get_adj_matrix_file` loads in legacy mode. In
//     sparse-graph mode `distance` is nullptr and the entire downstream
//     pass is skipped (which is what we want for SIFT n=50000 -- the
//     downstream cost there is infeasible anyway and the user explicitly
//     scoped MST/PAM to non-sparse runs).
//
//   - One `lookup(u,v)` function per method (TriSearch + each CurveModel).
//     Each method's full n*n LB matrix is materialized once via the lookup
//     fn, then reused for that method's MST + every PAM cluster_k.
//
// Outputs (per method):
//   MST line:
//     overlap  = | E_true intersect E_method | / (n - 1)
//     jaccard  = | intersect | / | union |
//     plus the MST weight under TRUE distances for both edge sets, so a
//     reader can see the geometric quality of the approximate MST as a
//     real tree (not just edge-set similarity).
//
//   PAM line per cluster_k:
//     cost_method   = sum of nearest-medoid distances under the method's
//                     LB matrix (what PAM was minimizing).
//     cost_true     = sum of nearest-medoid distances under the TRUE matrix
//                     using the SAME final medoid set (per user's Q2:
//                     fair-comparison evaluator).
//     cost_true_ratio = cost_true / pam_on_truth.cost_true. Equals 1.0 when
//                     the method recovers the same clustering as PAM-on-truth.
//
// PAM initialization: identical random medoid set across every method, seeded
// from `seed` (the same RNG seed that drives the rest of the binary). Per
// user's Q3: "ensuring initialization is consistent for fair comparison."
//
// Algorithmic notes:
//   - Prim's: O(n^2) dense version. n <= 4900 in scope, so 24M ops per MST.
//   - PAM: classical greedy global-best swap. Swap cost evaluation uses
//     cached nearest/second-nearest distances per point, so each swap eval
//     is O(n) and a full swap iteration is O(k * (n - k) * n). After an
//     accepted swap we recompute assignments from scratch (O(n * k)) which
//     is cheap relative to the search; this avoids the second-nearest
//     bookkeeping bugs that bite incremental updates.

#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "TriSearch.h"

namespace downstream {

using LookupFn = std::function<double(unsigned, unsigned)>;

struct MstEdge {
    unsigned u;
    unsigned v;
    double w;
};

static inline bool ds_isfinite(double x) { return std::isfinite(x); }

// ---------------------------------------------------------------------------
// Matrix materialization
// ---------------------------------------------------------------------------

// Materialize an n*n symmetric LB matrix by calling `lookup` for every i<j.
// For TriSearch this is O(n^2). For ELM/SPT methods each lookup is O(k) so
// this is O(n^2 * k). At n=4900, k=256: ~6e9 ops per method. Acceptable
// once; we reuse it for MST + every cluster_k in PAM.
inline std::vector<std::vector<double>>
materialize_lookup_matrix(const LookupFn& lookup, unsigned n) {
    std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = i + 1; j < n; ++j) {
            double w = lookup(i, j);
            if (!ds_isfinite(w) || w < 0) w = 0.0;
            M[i][j] = w;
            M[j][i] = w;
        }
    }
    return M;
}

// Copy the full true-distance matrix into a local n*n vector. Cheap; called
// once per downstream pass.
inline std::vector<std::vector<double>>
materialize_true_matrix(std::vector<std::vector<double>*>* distance, unsigned n) {
    std::vector<std::vector<double>> T(n, std::vector<double>(n, 0.0));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            double w = distance->at(i)->at(j);
            if (!ds_isfinite(w) || w < 0) w = 0.0;
            T[i][j] = w;
        }
    }
    return T;
}

// ---------------------------------------------------------------------------
// Prim's MST (dense O(n^2))
// ---------------------------------------------------------------------------

inline std::vector<MstEdge>
prims_mst(const std::vector<std::vector<double>>& M, unsigned n) {
    if (n < 2) return {};
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> dist(n, INF);
    std::vector<int> parent(n, -1);
    std::vector<char> in_tree(n, 0);
    dist[0] = 0.0;
    std::vector<MstEdge> mst;
    mst.reserve(n - 1);

    for (unsigned iter = 0; iter < n; ++iter) {
        int u = -1;
        double best = INF;
        for (unsigned v = 0; v < n; ++v) {
            if (!in_tree[v] && dist[v] < best) {
                best = dist[v];
                u = (int)v;
            }
        }
        if (u < 0) break;
        in_tree[u] = 1;
        if (parent[u] >= 0) {
            unsigned a = std::min((unsigned)u, (unsigned)parent[u]);
            unsigned b = std::max((unsigned)u, (unsigned)parent[u]);
            mst.push_back({a, b, M[(unsigned)u][(unsigned)parent[u]]});
        }
        for (unsigned v = 0; v < n; ++v) {
            if (!in_tree[v]) {
                double w = M[(unsigned)u][v];
                if (w < dist[v]) {
                    dist[v] = w;
                    parent[v] = u;
                }
            }
        }
    }
    return mst;
}

// ---------------------------------------------------------------------------
// MST comparison
// ---------------------------------------------------------------------------

struct MstComparison {
    unsigned true_edges;         // n - 1
    unsigned method_edges;       // n - 1 (also)
    unsigned overlap;            // |intersect|
    unsigned union_size;         // |union|
    double jaccard;              // overlap / union
    double overlap_frac;         // overlap / true_edges
    double true_weight;          // sum of true MST edge weights (TRUE matrix)
    double method_weight_lb;     // sum of method MST edge weights (LB matrix)
    double method_weight_true;   // sum of method MST edges UNDER true matrix
};

inline MstComparison compare_msts(const std::vector<MstEdge>& true_mst,
                                  const std::vector<MstEdge>& method_mst,
                                  const std::vector<std::vector<double>>& T) {
    std::set<std::pair<unsigned, unsigned>> true_set;
    std::set<std::pair<unsigned, unsigned>> method_set;
    for (const auto& e : true_mst) true_set.insert({e.u, e.v});
    for (const auto& e : method_mst) method_set.insert({e.u, e.v});

    unsigned overlap = 0;
    for (const auto& e : method_set) {
        if (true_set.count(e)) ++overlap;
    }
    unsigned union_size = (unsigned)(true_set.size() + method_set.size() - overlap);

    double jaccard = union_size > 0
        ? (double)overlap / (double)union_size
        : 0.0;
    double overlap_frac = !true_set.empty()
        ? (double)overlap / (double)true_set.size()
        : 0.0;

    double true_w = 0.0;
    for (const auto& e : true_mst) true_w += e.w;

    double method_w_lb = 0.0;
    for (const auto& e : method_mst) method_w_lb += e.w;

    double method_w_true = 0.0;
    for (const auto& e : method_mst) method_w_true += T[e.u][e.v];

    return {(unsigned)true_set.size(), (unsigned)method_set.size(),
            overlap, union_size, jaccard, overlap_frac,
            true_w, method_w_lb, method_w_true};
}

// ---------------------------------------------------------------------------
// PAM K-Medoids
// ---------------------------------------------------------------------------

// Internal state used to amortize the swap-cost evaluation.
struct PamState {
    std::vector<unsigned> medoids;     // size = cluster_k
    std::vector<unsigned> nearest_mi;  // for each point i: which medoid index is nearest
    std::vector<double>   nearest_d;
    std::vector<double>   second_d;
    double total_cost = 0.0;
};

// Recompute nearest/second-nearest medoid distances and total cost from
// scratch. O(n * cluster_k).
inline void pam_recompute(const std::vector<std::vector<double>>& D,
                          PamState& s, unsigned n) {
    const double INF = std::numeric_limits<double>::infinity();
    const unsigned cluster_k = (unsigned)s.medoids.size();
    s.nearest_mi.assign(n, 0);
    s.nearest_d.assign(n, INF);
    s.second_d.assign(n, INF);
    double total = 0.0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned mi = 0; mi < cluster_k; ++mi) {
            double d = D[i][s.medoids[mi]];
            if (!ds_isfinite(d) || d < 0) d = 0.0;
            if (d < s.nearest_d[i]) {
                s.second_d[i] = s.nearest_d[i];
                s.nearest_d[i] = d;
                s.nearest_mi[i] = mi;
            } else if (d < s.second_d[i]) {
                s.second_d[i] = d;
            }
        }
        if (!ds_isfinite(s.nearest_d[i])) s.nearest_d[i] = 0.0;
        total += s.nearest_d[i];
    }
    s.total_cost = total;
}

// Cost delta for swapping medoid at index mi_swap with non-medoid h.
// Uses cached nearest_d/second_d -> O(n) per evaluation.
inline double pam_swap_delta(const std::vector<std::vector<double>>& D,
                             const PamState& s,
                             unsigned mi_swap,
                             unsigned h,
                             unsigned n) {
    double delta = 0.0;
    for (unsigned i = 0; i < n; ++i) {
        double d_h = D[i][h];
        if (!ds_isfinite(d_h) || d_h < 0) d_h = 0.0;
        if (s.nearest_mi[i] == mi_swap) {
            // Removing i's current nearest. New nearest is min(h, second).
            double new_d = std::min(d_h, s.second_d[i]);
            if (!ds_isfinite(new_d)) new_d = 0.0;
            delta += (new_d - s.nearest_d[i]);
        } else {
            // i's current nearest is some other medoid (unchanged).
            // The new medoid h is a candidate; if closer, i migrates to h.
            if (d_h < s.nearest_d[i]) {
                delta += (d_h - s.nearest_d[i]);
            }
        }
    }
    return delta;
}

struct PamResult {
    std::vector<unsigned> medoids;
    double cost_method = 0.0;   // total cost under the matrix PAM optimized on
    double cost_true   = 0.0;   // total cost under the TRUE matrix (eval-only)
    unsigned iterations = 0;
    bool converged = false;
};

// Total nearest-medoid cost under matrix D given a medoid set. O(n * k).
inline double pam_eval_cost(const std::vector<std::vector<double>>& D,
                            const std::vector<unsigned>& medoids,
                            unsigned n) {
    const double INF = std::numeric_limits<double>::infinity();
    double total = 0.0;
    for (unsigned i = 0; i < n; ++i) {
        double best = INF;
        for (unsigned m : medoids) {
            double d = D[i][m];
            if (!ds_isfinite(d) || d < 0) d = 0.0;
            if (d < best) best = d;
        }
        if (!ds_isfinite(best)) best = 0.0;
        total += best;
    }
    return total;
}

// Greedy global-best swap. PAM original variant.
// M = distance matrix used for swap decisions (LB matrix for method runs,
//     TRUE matrix for the reference "PAM on truth" run).
// T = TRUE distance matrix used to score the final medoid set
//     (cost_true). For the reference run, pass T == M so cost_method == cost_true.
inline PamResult run_pam(const std::vector<std::vector<double>>& M,
                         const std::vector<std::vector<double>>& T,
                         unsigned n,
                         unsigned cluster_k,
                         const std::vector<unsigned>& initial_medoids,
                         unsigned max_iters = 100) {
    PamResult res;
    res.medoids = initial_medoids;
    if (cluster_k == 0 || n == 0 || cluster_k > n) {
        res.cost_method = 0.0;
        res.cost_true = 0.0;
        res.converged = true;
        return res;
    }

    PamState state;
    state.medoids = initial_medoids;
    pam_recompute(M, state, n);

    unsigned iter = 0;
    bool changed = true;
    while (changed && iter < max_iters) {
        changed = false;
        std::set<unsigned> medoid_set(state.medoids.begin(), state.medoids.end());

        int best_mi = -1;
        int best_h = -1;
        double best_delta = 0.0;  // only accept strict improvement
        for (unsigned mi = 0; mi < cluster_k; ++mi) {
            for (unsigned h = 0; h < n; ++h) {
                if (medoid_set.count(h)) continue;
                double delta = pam_swap_delta(M, state, mi, h, n);
                if (delta < best_delta - 1e-12) {
                    best_delta = delta;
                    best_mi = (int)mi;
                    best_h  = (int)h;
                }
            }
        }
        if (best_mi >= 0) {
            state.medoids[best_mi] = (unsigned)best_h;
            pam_recompute(M, state, n);
            changed = true;
        }
        ++iter;
    }

    res.medoids = state.medoids;
    res.cost_method = state.total_cost;
    res.cost_true   = pam_eval_cost(T, state.medoids, n);
    res.iterations = iter;
    res.converged = !changed;
    return res;
}

// Same initial medoid set for every method, seeded reproducibly.
inline std::vector<unsigned>
sample_initial_medoids(unsigned n, unsigned cluster_k, unsigned seed) {
    if (cluster_k > n) cluster_k = n;
    std::vector<unsigned> all(n);
    std::iota(all.begin(), all.end(), 0u);
    std::mt19937 rng(seed);
    std::shuffle(all.begin(), all.end(), rng);
    all.resize(cluster_k);
    std::sort(all.begin(), all.end());
    return all;
}

// ---------------------------------------------------------------------------
// Augustine et al. SIGMOD 2021 [4] protocol: count oracle (= true-distance)
// calls when the method's LB is used to skip comparisons that would otherwise
// require a distance lookup. The OUTPUT of the algorithm is identical to
// running on true distances; the METRIC is the count of oracle calls saved.
//
// This is reviewer R3.D2's literal request. Distinct from our existing
// quality metrics (MST overlap, PAM cost_true_ratio), which characterize how
// the output degrades when LBs *replace* true distances. Both blocks are
// emitted side by side in the downstream output.
//
// Implementation notes
// ====================
// Augustine's protocol wraps every distance comparison `dist(u,v) < theta` as:
//   - if LB(u,v) >= theta : skip without oracle call         (savings)
//   - if UB(u,v) <  theta : the inequality is true; we still
//                           need the exact value to update,
//                           so an oracle call is forced.
//   - otherwise           : inconclusive; oracle call.
// UB therefore yields no Prim's/PAM savings on its own because the algorithms
// must propagate the resolved value forward. The savings come entirely from
// the LB side: rejecting non-improving edges and far-away medoid candidates
// without invoking the oracle. We use only LB here. The TRUE matrix T is
// passed in so we can call it the "oracle"; a `resolved` bitmap memoizes
// edges so we count each unique oracle call once. Known edges (adj_mat >= 0)
// are pre-resolved and never charged.

struct OracleStats {
    unsigned long long oracle_calls = 0;   // unique true-distance lookups
    unsigned long long comparisons = 0;    // total bound-vs-threshold checks
    unsigned long long saved_by_bounds = 0;// comparisons resolved without oracle
};

// Build an n*n known-mask: 1 if adj_mat[u][v] >= 0 (i.e. (u,v) was a known
// input edge), else 0. adj_mat uses -1 as the unknown sentinel; the diagonal
// has 0 and never contributes oracle calls.
inline std::vector<std::vector<char>>
build_known_mask(std::vector<std::vector<double>*>* adj_mat, unsigned n) {
    std::vector<std::vector<char>> mask(n, std::vector<char>(n, 0));
    if (adj_mat == nullptr) return mask;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (i == j) continue;
            double d = adj_mat->at(i)->at(j);
            if (d >= -0.1 && ds_isfinite(d)) mask[i][j] = 1;
        }
    }
    return mask;
}

// Resolve true_dist(u,v): charges an oracle call if not pre-known and not
// yet memoized. Diagonal returns 0 with no charge.
inline double oracle_resolve(unsigned u, unsigned v,
                             const std::vector<std::vector<double>>& T,
                             const std::vector<std::vector<char>>& known,
                             std::vector<std::vector<char>>& resolved,
                             OracleStats& stats) {
    if (u == v) return 0.0;
    if (!known[u][v] && !resolved[u][v]) {
        resolved[u][v] = resolved[v][u] = 1;
        stats.oracle_calls++;
    }
    return T[u][v];
}

// ---------------------------------------------------------------------------
// Prim's MST under Augustine's oracle-saving protocol.
// Returns the same MST as prims_mst(T, n); reports oracle_calls saved.
// ---------------------------------------------------------------------------

inline std::pair<std::vector<MstEdge>, OracleStats>
prims_mst_augustine(const std::vector<std::vector<double>>& LB,
                    const std::vector<std::vector<double>>& T,
                    const std::vector<std::vector<char>>& known,
                    unsigned n) {
    OracleStats stats;
    std::vector<MstEdge> mst;
    if (n < 2) return {mst, stats};
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> dist(n, INF);
    std::vector<int> parent(n, -1);
    std::vector<char> in_tree(n, 0);
    std::vector<std::vector<char>> resolved(n, std::vector<char>(n, 0));
    dist[0] = 0.0;
    mst.reserve(n - 1);

    for (unsigned iter = 0; iter < n; ++iter) {
        int u = -1;
        double best = INF;
        for (unsigned v = 0; v < n; ++v) {
            if (!in_tree[v] && dist[v] < best) {
                best = dist[v];
                u = (int)v;
            }
        }
        if (u < 0) break;
        in_tree[u] = 1;
        if (parent[u] >= 0) {
            unsigned a = std::min((unsigned)u, (unsigned)parent[u]);
            unsigned b = std::max((unsigned)u, (unsigned)parent[u]);
            // Edge (parent[u], u) is the chosen MST edge for this iteration;
            // we already paid the oracle to learn its weight when it became
            // dist[u]. Record true weight.
            mst.push_back({a, b, T[(unsigned)u][(unsigned)parent[u]]});
        }
        // Relaxation. For each v not in tree, decide whether edge (u,v)
        // improves dist[v]. This is the wrapped comparison.
        for (unsigned v = 0; v < n; ++v) {
            if (in_tree[v] || (unsigned)u == v) continue;
            stats.comparisons++;
            const double lb_uv = LB[(unsigned)u][v];
            // If LB(u,v) >= dist[v], the edge cannot improve dist[v] and we
            // skip without an oracle call. This is the only savings branch
            // available without UB; see header note above.
            if (lb_uv >= dist[v]) {
                stats.saved_by_bounds++;
                continue;
            }
            // Otherwise we must resolve to compare and possibly update.
            double w = oracle_resolve((unsigned)u, v, T, known, resolved, stats);
            if (w < dist[v]) {
                dist[v] = w;
                parent[v] = u;
            }
        }
    }
    return {mst, stats};
}

// ---------------------------------------------------------------------------
// PAM under Augustine's oracle-saving protocol.
// Output medoids identical to run_pam(T, T, ...); reports oracle_calls saved.
//
// The Augustine wrap is applied to two operations:
//   (a) pam_recompute: for each point i, find nearest medoid. Comparison
//       D[i][m_j] < nearest_d[i]: skip if LB[i][m_j] >= nearest_d[i],
//       else oracle. Also need second-nearest for fast swap deltas, so
//       a second pass tracks the second smallest using the same protocol.
//   (b) pam_swap_delta: for each i, examine the new candidate h. Same
//       LB-vs-threshold check, oracle only when inconclusive.
// ---------------------------------------------------------------------------

inline void pam_recompute_augustine(
    const std::vector<std::vector<double>>& LB,
    const std::vector<std::vector<double>>& T,
    const std::vector<std::vector<char>>& known,
    std::vector<std::vector<char>>& resolved,
    OracleStats& stats,
    PamState& s, unsigned n)
{
    const double INF = std::numeric_limits<double>::infinity();
    const unsigned cluster_k = (unsigned)s.medoids.size();
    s.nearest_mi.assign(n, 0);
    s.nearest_d.assign(n, INF);
    s.second_d.assign(n, INF);
    double total = 0.0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned mi = 0; mi < cluster_k; ++mi) {
            unsigned m = s.medoids[mi];
            stats.comparisons++;
            // Two-tier check: skip if LB cannot beat current second_d
            // (no chance of becoming top-2). This is the conservative skip
            // since we maintain BOTH nearest and second-nearest.
            const double lb_im = LB[i][m];
            if (lb_im >= s.second_d[i] && lb_im >= s.nearest_d[i]) {
                stats.saved_by_bounds++;
                continue;
            }
            double d = oracle_resolve(i, m, T, known, resolved, stats);
            if (!ds_isfinite(d) || d < 0) d = 0.0;
            if (d < s.nearest_d[i]) {
                s.second_d[i] = s.nearest_d[i];
                s.nearest_d[i] = d;
                s.nearest_mi[i] = mi;
            } else if (d < s.second_d[i]) {
                s.second_d[i] = d;
            }
        }
        if (!ds_isfinite(s.nearest_d[i])) s.nearest_d[i] = 0.0;
        total += s.nearest_d[i];
    }
    s.total_cost = total;
}

inline double pam_swap_delta_augustine(
    const std::vector<std::vector<double>>& LB,
    const std::vector<std::vector<double>>& T,
    const std::vector<std::vector<char>>& known,
    std::vector<std::vector<char>>& resolved,
    OracleStats& stats,
    const PamState& s, unsigned mi_swap, unsigned h, unsigned n)
{
    double delta = 0.0;
    for (unsigned i = 0; i < n; ++i) {
        stats.comparisons++;
        // Decision: does i benefit from h being the new medoid?
        // i's current contribution is either nearest_d (if its nearest is not
        // mi_swap) or second_d (if mi_swap is being removed). Define the
        // threshold once and apply the LB check.
        double threshold;
        if (s.nearest_mi[i] == mi_swap) {
            threshold = s.second_d[i];
        } else {
            threshold = s.nearest_d[i];
        }
        const double lb_ih = LB[i][h];
        if (lb_ih >= threshold) {
            // h cannot improve i's contribution. Compute the no-h delta:
            //   if i's nearest is mi_swap: i migrates to second_d
            //   else: no change
            stats.saved_by_bounds++;
            if (s.nearest_mi[i] == mi_swap) {
                double new_d = s.second_d[i];
                if (!ds_isfinite(new_d)) new_d = 0.0;
                delta += (new_d - s.nearest_d[i]);
            }
            continue;
        }
        // Otherwise resolve d(i,h) via oracle.
        double d_h = oracle_resolve(i, h, T, known, resolved, stats);
        if (!ds_isfinite(d_h) || d_h < 0) d_h = 0.0;
        if (s.nearest_mi[i] == mi_swap) {
            double new_d = std::min(d_h, s.second_d[i]);
            if (!ds_isfinite(new_d)) new_d = 0.0;
            delta += (new_d - s.nearest_d[i]);
        } else {
            if (d_h < s.nearest_d[i]) {
                delta += (d_h - s.nearest_d[i]);
            }
        }
    }
    return delta;
}

inline std::pair<PamResult, OracleStats>
run_pam_augustine(const std::vector<std::vector<double>>& LB,
                  const std::vector<std::vector<double>>& T,
                  const std::vector<std::vector<char>>& known,
                  unsigned n, unsigned cluster_k,
                  const std::vector<unsigned>& initial_medoids,
                  unsigned max_iters = 100)
{
    OracleStats stats;
    PamResult res;
    res.medoids = initial_medoids;
    if (cluster_k == 0 || n == 0 || cluster_k > n) {
        res.converged = true;
        return {res, stats};
    }

    std::vector<std::vector<char>> resolved(n, std::vector<char>(n, 0));
    PamState state;
    state.medoids = initial_medoids;
    pam_recompute_augustine(LB, T, known, resolved, stats, state, n);

    unsigned iter = 0;
    bool changed = true;
    while (changed && iter < max_iters) {
        changed = false;
        std::set<unsigned> medoid_set(state.medoids.begin(),
                                      state.medoids.end());
        int best_mi = -1;
        int best_h = -1;
        double best_delta = 0.0;
        for (unsigned mi = 0; mi < cluster_k; ++mi) {
            for (unsigned h = 0; h < n; ++h) {
                if (medoid_set.count(h)) continue;
                double delta = pam_swap_delta_augustine(
                    LB, T, known, resolved, stats, state, mi, h, n);
                if (delta < best_delta - 1e-12) {
                    best_delta = delta;
                    best_mi = (int)mi;
                    best_h  = (int)h;
                }
            }
        }
        if (best_mi >= 0) {
            state.medoids[best_mi] = (unsigned)best_h;
            pam_recompute_augustine(LB, T, known, resolved, stats, state, n);
            changed = true;
        }
        ++iter;
    }

    res.medoids = state.medoids;
    res.iterations = iter;
    res.converged = !changed;
    // Final cost evaluations under TRUE matrix (these read T but are post-hoc
    // measurements not part of the algorithm; they do not count as oracle
    // calls in Augustine's framing).
    res.cost_method = pam_eval_cost(T, res.medoids, n);
    res.cost_true   = res.cost_method;
    return {res, stats};
}

// Naive baseline (Augustine's "Without Plug" style): count distance-oracle
// calls with no LB skip gate, while preserving the same graph-update / unique
// oracle-call convention used by oracle_resolve(). These helpers are
// reporting-only and do not change the MST/PAM execution.

// Count input-known and input-unknown unordered pairs. These are reporting
// baselines only; they do not affect the Augustine MST/PAM execution.
inline unsigned long long count_known_unordered_edges(
    const std::vector<std::vector<char>>& known, unsigned n)
{
    unsigned long long calls = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = i + 1; j < n; ++j) {
            if (known[i][j]) calls++;
        }
    }
    return calls;
}

inline unsigned long long count_unknown_unordered_edges(
    const std::vector<std::vector<char>>& known, unsigned n)
{
    unsigned long long calls = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = i + 1; j < n; ++j) {
            if (!known[i][j]) calls++;
        }
    }
    return calls;
}

// Naive baseline (Augustine's "Without Plug" style) for Prim's under the
// same graph-update / memoization convention used by oracle_resolve(): if no
// lower-bound gate is available, dense Prim eventually considers every
// unordered pair once. Pre-known input edges are free, so the charged oracle
// calls are exactly the initially unknown unordered pairs.
inline unsigned long long count_naive_prims_calls(
    const std::vector<std::vector<char>>& known, unsigned n)
{
    return count_unknown_unordered_edges(known, n);
}

// Naive baseline for PAM under the same graph-update / memoization convention:
// a full PAM swap scan touches every pair at least once (initial medoid
// assignment covers medoid-medoid and point-medoid pairs; swap scans cover
// non-medoid candidates). Therefore the charged unique oracle calls are again
// the initially unknown unordered pairs. This is the denominator used for the
// paper-style exact-output oracle-saving percentage.
inline unsigned long long count_naive_pam_unique_calls(
    const std::vector<std::vector<char>>& known, unsigned n, unsigned cluster_k)
{
    if (n < 2 || cluster_k == 0 || cluster_k > n) return 0;
    return count_unknown_unordered_edges(known, n);
}

// Optional diagnostic only: an uncached repeated-lookup count for PAM's inner
// loops. This can be useful for debugging CPU/oracle-pressure, but it is not
// used for the paper-style saved_vs_naive_pct because oracle_resolve() follows
// the graph-update convention and charges each unique unknown pair once.
inline unsigned long long estimate_pam_uncached_repeated_calls(
    const std::vector<std::vector<char>>& known,
    unsigned n,
    unsigned cluster_k,
    unsigned iterations,
    const std::vector<unsigned>& initial_medoids)
{
    if (n == 0 || cluster_k == 0 || cluster_k > n) return 0;

    unsigned long long calls = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned m_v : initial_medoids) {
            if (i != m_v && !known[i][m_v]) calls++;
        }
    }

    unsigned long long unknown_ordered_pairs = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (i != j && !known[i][j]) unknown_ordered_pairs++;
        }
    }
    double pair_count = (double)n * (double)(n - 1);
    double unk_frac = pair_count > 0
        ? (double)unknown_ordered_pairs / pair_count
        : 0.0;
    calls += (unsigned long long)(
        (double)iterations *
        (double)cluster_k * (double)(n - cluster_k) * (double)n * unk_frac);
    return calls;
}

// ---------------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------------

inline void run_downstream(
    std::vector<std::vector<double>*>* distance,
    std::vector<std::vector<double>*>* adj_mat,  // -1 sentinel for unknown
    TriSearch& tri,
    const std::vector<LookupFn>& model_lookups,
    const std::vector<std::string>& model_names,
    unsigned n,
    const std::vector<unsigned>& cluster_ks,
    unsigned seed,
    bool lb_as_distance = false)  // OFF by default: only paper-faithful metrics
{
    std::cout << "\n==================== DOWNSTREAM TASKS ====================\n";
    if (distance == nullptr) {
        std::cout << "[Downstream] SKIPPED: distance matrix is null "
                  << "(sparse-graph mode; not in scope per R3.D2)\n";
        std::cout << "==========================================================\n";
        return;
    }
    if (n < 2) {
        std::cout << "[Downstream] SKIPPED: n=" << n << " < 2\n";
        std::cout << "==========================================================\n";
        return;
    }

    std::cout << "[Downstream] n=" << n << "  cluster_ks=[";
    for (size_t i = 0; i < cluster_ks.size(); ++i) {
        std::cout << cluster_ks[i] << (i + 1 < cluster_ks.size() ? "," : "");
    }
    std::cout << "]  seed=" << seed << "\n";

    auto t_start = std::chrono::steady_clock::now();

    // 1. Materialize true matrix and true MST.
    auto T = materialize_true_matrix(distance, n);

    // Augustine protocol: which (u,v) pairs were known at input time.
    auto known = build_known_mask(adj_mat, n);
    unsigned long long known_edge_count = count_known_unordered_edges(known, n);
    unsigned long long unknown_edge_count = count_unknown_unordered_edges(known, n);
    unsigned long long naive_mst_calls = count_naive_prims_calls(known, n);

    auto true_mst = prims_mst(T, n);
    double true_mst_w = 0.0;
    for (auto& e : true_mst) true_mst_w += e.w;
    std::cout << "[MST true] edges=" << true_mst.size()
              << "  total_weight=" << true_mst_w
              << "  known_edges=" << known_edge_count
              << "  unknown_edges=" << unknown_edge_count
              << "  naive_prims_oracle_calls=" << naive_mst_calls << "\n";

    // 2. Sample one shared initial medoid set per cluster_k.
    std::map<unsigned, std::vector<unsigned>> initial_medoids;
    for (unsigned ck : cluster_ks) {
        initial_medoids[ck] = sample_initial_medoids(n, ck, seed);
        std::cout << "[PAM init cluster_k=" << ck << "] medoids={";
        for (size_t i = 0; i < initial_medoids[ck].size(); ++i) {
            std::cout << initial_medoids[ck][i]
                      << (i + 1 < initial_medoids[ck].size() ? "," : "");
        }
        std::cout << "}\n";
    }

    // 3. Reference PAM on TRUE distances. cost_method == cost_true here.
    std::map<unsigned, PamResult> true_pam;
    for (unsigned ck : cluster_ks) {
        auto pr = run_pam(T, T, n, ck, initial_medoids[ck]);
        true_pam[ck] = pr;
        std::cout << "[PAM true cluster_k=" << ck << "]"
                  << "  cost_true=" << pr.cost_true
                  << "  iters=" << pr.iterations
                  << "  converged=" << (pr.converged ? "yes" : "no")
                  << "  medoids={";
        for (size_t i = 0; i < pr.medoids.size(); ++i) {
            std::cout << pr.medoids[i]
                      << (i + 1 < pr.medoids.size() ? "," : "");
        }
        std::cout << "}\n";
    }

    // 4. Per-method MST + PAM. One materialization per method.
    // Collect each method's MST_augustine result so we can print a
    // cross-method comparison (vs SPLUB) at the end: same MST weight for all,
    // differing only in oracle calls.
    struct MstSummary { std::string name; unsigned long long oracle_calls; double weight; bool exact; };
    std::vector<MstSummary> mst_summary;
    auto run_for_method = [&](const std::string& name, const LookupFn& lookup) {
        auto t0 = std::chrono::steady_clock::now();
        std::cout << "\n[Downstream method=" << name << "]\n";

        auto M = materialize_lookup_matrix(lookup, n);
        auto t_mat = std::chrono::steady_clock::now();
        double mat_s = std::chrono::duration<double>(t_mat - t0).count();
        unsigned long long materialize_pairs = ((unsigned long long)n * ((unsigned long long)n - 1) / 2);
        double materialize_avg_us = materialize_pairs > 0
            ? (mat_s * 1000000.0) / (double)materialize_pairs
            : 0.0;
        std::cout << "  materialize_s=" << mat_s
                  << "  materialize_pairs=" << materialize_pairs
                  << "  materialize_avg_us=" << materialize_avg_us
                  << "\n";

        // MST -- quality metric (LB-as-distance). NON-PAPER framing: this
        // builds the MST on lower-bound values AS IF they were the true
        // distance, which produces a different (wrong) tree per method. The
        // paper never does this. OFF unless --downstream-lb-as-distance.
        if (lb_as_distance) {
            auto method_mst = prims_mst(M, n);
            auto cmp = compare_msts(true_mst, method_mst, T);
            std::cout << "  [non-paper] MST_lb_as_distance: overlap=" << cmp.overlap
                      << "/" << cmp.true_edges
                      << " (" << (100.0 * cmp.overlap_frac) << "%)"
                      << "  jaccard=" << cmp.jaccard
                      << "  true_w=" << cmp.true_weight
                      << "  method_w_lb=" << cmp.method_weight_lb
                      << "  method_w_true=" << cmp.method_weight_true << "\n";
        }

        // MST -- Augustine protocol (LB as oracle-skip gate; EXACT output).
        // When LB cannot resolve an edge comparison, the true-distance oracle
        // is called, so the resulting MST -- and its total weight -- is exact
        // and IDENTICAL for every method. Only oracle_calls differ.
        auto aug_mst = prims_mst_augustine(M, T, known, n);
        const auto& aug_mst_edges = aug_mst.first;
        const auto& aug_mst_stats = aug_mst.second;
        // Total weight of the produced MST (uses true edge weights).
        double aug_mst_w = 0.0;
        for (const auto& e : aug_mst_edges) aug_mst_w += e.w;
        // Assertion: Augustine output must equal true MST (edge sets identical
        // up to representation). Sort and compare canonical (min,max) pairs.
        bool mst_identical = (aug_mst_edges.size() == true_mst.size());
        if (mst_identical) {
            std::set<std::pair<unsigned,unsigned>> true_set, aug_set;
            for (const auto& e : true_mst) true_set.insert({e.u, e.v});
            for (const auto& e : aug_mst_edges) aug_set.insert({e.u, e.v});
            mst_identical = (true_set == aug_set);
        }
        // Weight matches true MST weight to numerical tolerance.
        bool weight_matches = (std::fabs(aug_mst_w - true_mst_w) <= 1e-6 * (1.0 + std::fabs(true_mst_w)));
        double mst_save_pct = (naive_mst_calls > 0)
            ? (100.0 * ((double)naive_mst_calls -
                        (double)aug_mst_stats.oracle_calls)
                       / (double)naive_mst_calls)
            : 0.0;
        std::cout << "  MST_augustine: oracle_calls=" << aug_mst_stats.oracle_calls
                  << "  naive_calls=" << naive_mst_calls
                  << "  saved_vs_naive_pct=" << mst_save_pct
                  << "  comparisons=" << aug_mst_stats.comparisons
                  << "  saved_by_lb=" << aug_mst_stats.saved_by_bounds
                  << "  mst_weight=" << aug_mst_w
                  << "  weight_matches_true=" << (weight_matches ? "yes" : "no")
                  << "  output_matches_true=" << (mst_identical ? "yes" : "no")
                  << "\n";
        // Record for the cross-method / vs-SPLUB comparison summary.
        mst_summary.push_back({name, aug_mst_stats.oracle_calls, aug_mst_w,
                               mst_identical && weight_matches});

        // PAM at each cluster_k.
        for (unsigned ck : cluster_ks) {
            // PAM quality metric (LB-as-distance, NON-PAPER). OFF by default.
            if (lb_as_distance) {
                auto pr = run_pam(M, T, n, ck, initial_medoids[ck]);
                double ref_true = true_pam[ck].cost_true;
                double ratio = (ref_true > 0) ? pr.cost_true / ref_true : 1.0;
                std::cout << "  [non-paper] PAM_lb_as_distance cluster_k=" << ck
                          << "  cost_method=" << pr.cost_method
                          << "  cost_true=" << pr.cost_true
                          << "  cost_true_ratio=" << ratio
                          << "  medoids={";
                for (size_t i = 0; i < pr.medoids.size(); ++i) {
                    std::cout << pr.medoids[i]
                              << (i + 1 < pr.medoids.size() ? "," : "");
                }
                std::cout << "}\n";
            }

            // PAM under Augustine's protocol (exact output, oracle savings).
            auto aug_pam_t0 = std::chrono::steady_clock::now();
            auto aug_pr = run_pam_augustine(
                M, T, known, n, ck, initial_medoids[ck]);
            auto aug_pam_t1 = std::chrono::steady_clock::now();
            double aug_pam_s = std::chrono::duration<double>(
                aug_pam_t1 - aug_pam_t0).count();
            const auto& aug_pam = aug_pr.first;
            const auto& aug_pam_stats = aug_pr.second;
            // Reporting-only baselines. The paper-style denominator follows
            // the same graph-update / unique-oracle-call convention as
            // oracle_resolve(). The uncached number is printed separately as
            // a diagnostic and is not used for saved_vs_naive_pct.
            unsigned long long naive_pam_calls =
                count_naive_pam_unique_calls(known, n, ck);
            unsigned long long naive_pam_uncached_calls =
                estimate_pam_uncached_repeated_calls(
                    known, n, ck, aug_pam.iterations, initial_medoids[ck]);
            bool pam_medoids_match =
                (std::set<unsigned>(aug_pam.medoids.begin(),
                                    aug_pam.medoids.end())
                 == std::set<unsigned>(true_pam[ck].medoids.begin(),
                                       true_pam[ck].medoids.end()));
            double pam_save_pct = (naive_pam_calls > 0)
                ? (100.0 * ((double)naive_pam_calls -
                            (double)aug_pam_stats.oracle_calls)
                          / (double)naive_pam_calls)
                : 0.0;
            std::cout << "  PAM_augustine cluster_k=" << ck
                      << "  oracle_calls=" << aug_pam_stats.oracle_calls
                      << "  naive_calls=" << naive_pam_calls
                      << "  saved_vs_naive_pct=" << pam_save_pct
                      << "  naive_uncached_calls~=" << naive_pam_uncached_calls
                      << "  comparisons=" << aug_pam_stats.comparisons
                      << "  saved_by_lb=" << aug_pam_stats.saved_by_bounds
                      << "  iters=" << aug_pam.iterations
                      << "  medoids_match_true=" << (pam_medoids_match ? "yes" : "no")
                      << "  pam_s=" << aug_pam_s
                      << "\n";
        }
    };

    // 4a. TriSearch.
    run_for_method("TriSearch",
                   [&](unsigned u, unsigned v) { return tri.lookup(u, v); });

    // 4b. Each CurveModel.
    for (size_t mi = 0; mi < model_lookups.size(); ++mi) {
        run_for_method(model_names[mi], model_lookups[mi]);
    }

    auto t_end = std::chrono::steady_clock::now();
    double total_s = std::chrono::duration<double>(t_end - t_start).count();

    // -----------------------------------------------------------------------
    // MST cross-method comparison (vs SPLUB). Under the Augustine protocol the
    // produced MST is exact, so the total weight is IDENTICAL across methods;
    // methods differ only in the number of true-distance oracle calls. SPLUB
    // (exact tightest lower bound) is the reference: it skips the most edges,
    // so it is the lower envelope on oracle calls among the bound methods.
    // -----------------------------------------------------------------------
    std::cout << "\n[MST comparison vs SPLUB]  (true MST weight=" << true_mst_w << ")\n";
    unsigned long long splub_calls = 0;
    bool have_splub = false;
    for (const auto& s : mst_summary) {
        if (s.name == "SPLUB") { splub_calls = s.oracle_calls; have_splub = true; break; }
    }
    bool all_same_weight = true;
    for (const auto& s : mst_summary)
        if (std::fabs(s.weight - true_mst_w) > 1e-6 * (1.0 + std::fabs(true_mst_w))) all_same_weight = false;
    std::cout << "  all_methods_same_mst_weight=" << (all_same_weight ? "yes" : "no") << "\n";
    std::cout << "  " << std::left << std::setw(22) << "method"
              << std::right << std::setw(14) << "mst_weight"
              << std::setw(14) << "oracle_calls"
              << std::setw(16) << "vs_SPLUB_calls"
              << std::setw(10) << "exact" << "\n";
    for (const auto& s : mst_summary) {
        std::cout << "  " << std::left << std::setw(22) << s.name
                  << std::right << std::setw(14) << s.weight
                  << std::setw(14) << s.oracle_calls;
        if (have_splub) {
            long long delta = (long long)s.oracle_calls - (long long)splub_calls;
            std::cout << std::setw(15) << (delta >= 0 ? "+" : "")
                      << delta;
        } else {
            std::cout << std::setw(16) << "n/a";
        }
        std::cout << std::setw(10) << (s.exact ? "yes" : "no") << "\n";
    }
    if (have_splub)
        std::cout << "  (SPLUB oracle_calls=" << splub_calls
                  << " is the exact-bound reference; positive vs_SPLUB means a"
                     " looser bound forced more oracle calls.)\n";

    std::cout << "\n[Downstream] total_seconds=" << total_s << "\n";
    std::cout << "==========================================================\n";
}

}  // namespace downstream
