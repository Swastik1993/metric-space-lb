#ifndef ELM_SPTREE_H
#define ELM_SPTREE_H

#include "Dijkstra.h"
#include "GraphDefinitions.h"
#include "GraphUtils.h"
#include "ShortestPathTree.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <list>
#include <map>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <ctime>

using namespace std;

// Per-component microsecond clock for the ablation timing report. Uses
// CLOCK_PROCESS_CPUTIME_ID to match cpu_time_us() in fixed_main.cpp. Used to
// attribute SPT (Dijkstra) construction and RMQ (LCA-table) preprocessing cost
// to separate accumulators inside ELM_SPTree.
static inline uint64_t elm_cpu_us() {
  struct timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  return uint64_t(ts.tv_sec) * 1000000ULL + ts.tv_nsec / 1000ULL;
}

class ELM_SPTree {
private:
  vector<list<pair<unsigned int, double>> *> *adj_list;
  unsigned int no_nodes;
  unsigned int no_landmarks;
  unsigned int no_samples;

  map<pair<unsigned int, unsigned int>, double> *edges_map;
  map<pair<unsigned int, unsigned int>, double> *landmarks;
  vector<pair<pair<unsigned int, unsigned int>, double>> landmark_order;

  // root -> SPT(root) nodes (indexed by vertex id)
  // direct-indexed (size=no_nodes). nullptr means "not built yet".
  vector<vector<shortest_path_tree*>*> sp_tree_vec;

  // Only used for last round of HCA lookups
  vector<unsigned int> temp_HCA;

  // Temporary storage for HCA
  vector<unsigned int> _hca_information;

  // Generation counters for O(1) HCA reset. The original _reset_HCA() refilled
  // both _hca_information and temp_HCA with sentinel values in O(no_nodes)
  // time, which on a 50000-node graph with thousands of (candidate, query)
  // pairs costs trillions of fill operations and dominates the wall clock.
  // The two _gen vectors record the generation at which each cell was last
  // written; a cell is considered "fresh" iff its gen matches _hca_gen.
  // Bumping _hca_gen is O(1) and effectively invalidates every cell.
  std::vector<uint64_t> _hca_information_gen;
  std::vector<uint64_t> _temp_HCA_gen;
  uint64_t _hca_gen = 0;

  // Optional cross-query LB cache. When enabled (via enable_lb_cache()),
  // lookup() consults the cache before walking landmarks; on miss it
  // computes, stores, and returns. The cache is cleared whenever the
  // landmark set changes (e.g. after get_*_landmarks). Off by default to
  // preserve the original lookup() semantics and memory footprint.
  bool lb_cache_enabled_ = false;
  std::unordered_map<uint64_t, double> lb_cache_;

  // Ablation flags (Reviewer R3.D7). Each toggle isolates one Section-4
  // optimization. All off by default => v12 / paper behavior.
  // --ablation-no-hca   : get_hca_landmarks_parallel uses _lookup
  //                       (edge-folding + landmark-tree §4.1 bound) instead
  //                       of _lookup_HCA (which adds the query-endpoint
  //                       §4.2.1 HCA term). Isolates the §4.2.1 contribution.
  // --ablation-no-rmq   : when set, find_LCA / _lookup_vertex skip the
  //                       O(1) LCA/RMQ preprocessing and walk the parent
  //                       chain naively (O(depth) per query). Isolates the
  //                       §4.1.3 amortization benefit.
  // --ablation-no-spt   : disables the SPT/tree-path tightening terms in
  //                       _lookup and _lookup_HCA, leaving only edge-landmark
  //                       folding through the selected landmark edge. This is
  //                       the meaningful quality ablation for the SPT bound.
  bool ablation_no_hca_ = false;
  bool ablation_no_rmq_ = false;
  bool ablation_no_spt_ = false;

  // Component build-phase timers for the ablation report (R3.D7 timing).
  // spt_build_us_ : cumulative time spent constructing shortest-path trees
  //                 (the Dijkstra() calls inside ensure_spt). This is the SPT
  //                 build cost.
  // rmq_build_us_ : cumulative time spent in preprocess_spt_queries() building
  //                 the Euler-tour / sparse-table / Cartesian RMQ structures
  //                 used by O(1) find_LCA. This is the RMQ build cost.
  // Both accumulate only along the sequential build path (ensure_spt and
  // preprocess_landmark_roots), so they are race-free. reset_build_timers()
  // lets a caller measure a single preprocess_landmark_roots() pass in
  // isolation (the ablation report builds a fresh twin and times exactly that).
  uint64_t spt_build_us_ = 0;
  uint64_t rmq_build_us_ = 0;

  inline pair<unsigned int, unsigned int> norm_edge(unsigned int a, unsigned int b) const {
    return make_pair(min(a, b), max(a, b));
  }

  inline uint64_t pair_code(unsigned int a, unsigned int b) const {
    auto e = norm_edge(a, b);
    return (uint64_t(e.first) << 32) | uint64_t(e.second);
  }

  void store_map();

  inline void ensure_spt(unsigned int root) {
    if (root >= no_nodes) return;
    if (sp_tree_vec[root] == nullptr) {
      uint64_t _t0 = elm_cpu_us();
      sp_tree_vec[root] = Dijkstra(adj_list, root, true);
      spt_build_us_ += (elm_cpu_us() - _t0);
    }
  }

  double _lookup_vertex(unsigned int root_node, unsigned int u, unsigned int v);
  double _lookup(unsigned int u, unsigned int v,
                 pair<unsigned int, unsigned int> landmark_edge,
                 double edge_dist);

  // Explicit-ablation variants used by the single-run ablation report. They
  // take per-call no_spt / no_rmq flags instead of reading the instance
  // members, so the SAME selected model can be evaluated with a component ON
  // and OFF over an identical query set without mutating instance state (and
  // without disturbing the selected landmark set). The member-flag versions
  // above simply delegate to these with (ablation_no_spt_, ablation_no_rmq_),
  // so existing behavior is byte-identical.
  double _lookup_vertex_ex(unsigned int root_node, unsigned int u,
                           unsigned int v, bool no_rmq, bool skip_lca = false);
  double _lookup_ex(unsigned int u, unsigned int v,
                    pair<unsigned int, unsigned int> landmark_edge,
                    double edge_dist, bool no_spt, bool no_rmq,
                    bool skip_lca = false);

  double _lookup_vertex_HCA(unsigned int root_node, unsigned int u, unsigned int v);
  double _lookup_HCA(unsigned int u, unsigned int v,
                     pair<unsigned int, unsigned int> landmark_edge,
                     double edge_dist);

  vector<pair<unsigned int, unsigned int>> sample_missing_pairs(unsigned int count,
                                                                std::mt19937 &rng) const;
  vector<pair<pair<unsigned int, unsigned int>, double>> sample_candidate_edges(
      unsigned int count,
      std::mt19937 &rng) const;
  vector<pair<unsigned int, unsigned int>> all_missing_pairs() const;

  // Optimization (#19): hashed mirror of edges_map for O(1) membership tests.
  // Built lazily on first query that asks for it. Keeps the original
  // std::map (ordered iteration is used elsewhere in the class).
  std::unordered_map<uint64_t, double> edges_hash_;
  bool edges_hash_built_ = false;
  inline void ensure_edges_hash() {
    if (edges_hash_built_) return;
    edges_hash_.reserve(edges_map->size() * 2 + 8);
    for (auto &kv : *edges_map) {
      edges_hash_.emplace(pair_code(kv.first.first, kv.first.second), kv.second);
    }
    edges_hash_built_ = true;
  }
  inline bool has_edge_fast(unsigned int u, unsigned int v) {
    ensure_edges_hash();
    return edges_hash_.count(pair_code(min(u, v), max(u, v))) > 0;
  }

public:
  ELM_SPTree(vector<list<pair<unsigned int, double>> *> *adj_list,
             unsigned int nodes,
             unsigned int k,
             unsigned int no_samples);

  ~ELM_SPTree();

  void preprocess_landmark_roots();

  // Persist / restore the SELECTED landmark set (landmark_order) to/from disk.
  // The landmark set is the only expensive output of get_*_landmarks(); after
  // load_landmarks(...) the caller must call preprocess_landmark_roots() to
  // rebuild the SP-trees + RMQ, which fully reconstructs a queryable model
  // (lookup() depends only on those, not on selection-time state).
  // save_landmarks writes a small text file: a metadata header line, then one
  // line "a b weight" per landmark. load_landmarks returns false if the file
  // is missing/unreadable or the metadata header does not match `meta`.
  bool save_landmarks(const std::string& path, const std::string& meta) const;
  bool load_landmarks(const std::string& path, const std::string& meta);

  void get_sampling_landmarks();
  void get_exact_landmarks();
  // Add `cand_cap`: upper bound on the per-iteration candidate sample size
  // produced by the Mirzasoleiman formula (n/k * log(1/eps)). Without a cap,
  // small k on a large graph asks for hundreds of thousands of candidate
  // evaluations per iteration, each one triggering two full SPT+RMQ builds.
  // Pass UINT_MAX (the default) to preserve the original unbounded behavior.
  void get_stochastic_landmarks(double eps,
                                unsigned int query_samples,
                                unsigned int seed = 73,
                                unsigned int cand_cap = std::numeric_limits<unsigned int>::max());

  // ---- New methods (additive; existing API unchanged) ----

  // Lazy/accelerated greedy GEL (Minoux 1978). Same approximation guarantee
  // as get_exact_landmarks() but typically 10-100x fewer marginal-gain
  // computations on real datasets. Maintains a max-heap of stale upper-bound
  // gains and refreshes lazily; submodularity keeps the bound valid.
  void get_exact_landmarks_lazy();

  // Optional OpenMP-parallel inner loop over candidate edges in
  // get_exact_landmarks(). Disabled when not compiled with -fopenmp; in that
  // case, falls back to a single-threaded run that is byte-identical to
  // get_exact_landmarks(). When enabled, each thread accumulates a private
  // best-gain candidate; the main thread reduces them.
  void get_exact_landmarks_parallel();

  // Exact greedy GEL using the Zoom/Suraj ARLB-aligned objective.
  // Same candidate and query sets as get_exact_landmarks_parallel(), but
  // each marginal gain is normalized by 1/OPT(q), where OPT(q) is the
  // exact/SPLUB lower bound for the missing query edge q. Queries with
  // OPT(q)=0 are excluded. This is the exact-greedy counterpart of
  // ELM Sampling RelObj.
  void get_exact_landmarks_relative_parallel();

  // Paper-faithful §4.2 SGEL: samples query pairs only, evaluates EVERY known
  // edge as a candidate per query using the HCA-based bound (one upward walk
  // in each of the two query-side SPTs, memoized). This matches the paper's
  // claimed O(sk(m + n log n)) preprocessing and avoids the per-candidate
  // Dijkstra cost of get_stochastic_landmarks().
  //
  // Difference from get_stochastic_landmarks():
  //   - That one ALSO subsamples candidates (Mirzasoleiman 2015 style) and
  //     pays full Dijkstra+RMQ for each sampled candidate.
  //   - This one keeps the candidate set complete and uses the cheaper HCA
  //     scoring. Better matches the paper's complexity story; recommended
  //     entry point if you want the algorithm exactly as described.
  void get_hca_landmarks(unsigned int query_samples,
                         unsigned int seed = 73,
                         bool relative_objective = false);

  // OpenMP-parallel version of get_hca_landmarks. Same algorithm, same
  // approximation guarantees; differs from the serial version only in
  // tie-breaking order among equal-gain candidates.
  //
  // If relative_objective=true, the greedy marginal gain is weighted by
  // 1/OPT(q), where OPT(q) is the exact/SPLUB lower bound of the sampled
  // training query q. Queries with OPT(q)=0 are ignored. This is the
  // Zoom-discussed ARLB-aligned objective: maximize average relative lower
  // bound instead of raw L1 lower-bound increase.
  //
  // Uses per-thread HCA buffers with generation-counter resets, so the
  // O(n) reset cost per (candidate, query) pair is reduced to O(1). On a
  // 16-core run for n = 50000, k = 256, query_samples = 500 this typically
  // achieves an 8-12x speedup over the serial version (not 16x due to
  // memory bandwidth saturation in the SPT pointer chases).
  //
  // When compiled without -fopenmp, falls back to the serial implementation.
  void get_hca_landmarks_parallel(unsigned int query_samples,
                                  unsigned int seed = 73,
                                  bool relative_objective = false);

  // Batched per-source LB fill for kNN evaluation.
  // Equivalent to looping `LB[v] = lookup(u, v)` over all v but with one walk
  // over landmarks and tighter inner loops. Substantially faster on warm
  // caches.
  void lookup_from_source(unsigned int u, vector<double> &LB);

  // Cross-query LB cache control (item #11 of the audit).
  // When enabled, lookup() returns cached results for repeated (u, v) pairs
  // and stores fresh ones. Useful when downstream code hits the same pair
  // multiple times (e.g., kNN pair-resolution sampling). Default: off.
  void enable_lb_cache(bool enable = true) { lb_cache_enabled_ = enable; if (!enable) lb_cache_.clear(); }
  void clear_lb_cache() { lb_cache_.clear(); }
  size_t lb_cache_size() const { return lb_cache_.size(); }

  // Ablation toggles (R3.D7).
  void set_ablation_no_hca(bool v) { ablation_no_hca_ = v; }
  void set_ablation_no_rmq(bool v) { ablation_no_rmq_ = v; }
  void set_ablation_no_spt(bool v) { ablation_no_spt_ = v; }
  bool ablation_no_hca() const { return ablation_no_hca_; }
  bool ablation_no_rmq() const { return ablation_no_rmq_; }
  bool ablation_no_spt() const { return ablation_no_spt_; }

  // ---- Single-run ablation report support (R3.D7 timing + dual quality) ----
  // Build-phase component timers, populated along the sequential build path.
  uint64_t spt_build_us() const { return spt_build_us_; }
  uint64_t rmq_build_us() const { return rmq_build_us_; }
  void reset_build_timers() { spt_build_us_ = 0; rmq_build_us_ = 0; }

  // Evaluate the lower bound for (u, v) with per-call ablation flags...
  // (see header notes). skip_lca additionally bypasses find_LCA in the SPT
  // tightening terms (using the tree root as a stand-in LCA) so the report can
  // separate find_LCA time from the surrounding tree-path arithmetic by
  // differencing a normal vs skip_lca pass. The returned value is meaningful
  // only when skip_lca == false; with skip_lca == true it is a timing probe.
  double lookup_ablate(unsigned int u, unsigned int v, bool no_spt, bool no_rmq,
                       bool skip_lca = false);

  // Loop the selected landmarks doing negligible work -- measures the
  // per-query landmark-iteration / dispatch overhead (the "other" query cost).
  double lookup_skeleton(unsigned int u, unsigned int v);

  double lookup(unsigned int u, unsigned int v);
  vector<double> *lookup_multiple(unsigned int u, unsigned int v);
  size_t _sizeof() const;

  void clean_unwanted_shortest_paths();
  void _reset_HCA();
};

#endif
