#ifndef EDGELANDMARK_H
#define EDGELANDMARK_H

#include "Dijkstra.h"
#include "GraphDefinitions.h"
#include <algorithm>
#include <boost/heap/pairing_heap.hpp>
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <utility>
#include <vector>


using namespace std;
using namespace boost::heap;

class EdgeLandMark {
private:
  vector<list<pair<unsigned int, double>> *> *adj_list;
  unsigned int no_nodes;
  unsigned int no_landmarks;
  unsigned int no_samples;
  vector<vector<double> *> *sp_vector;
  vector<vector<double> *> *le_vector;
  map<pair<unsigned int, unsigned int>, double> *landmarks;
  map<pair<unsigned int, unsigned int>, double> *edges_map;

  // Reproducible RNG (seedable per instance). Replaces ad-hoc rand() usage in
  // get_landmarks; rand()-based code paths remain available for backward
  // compatibility (set seed via set_seed; default mirrors prior behavior since
  // seeding still happens but the stream is now mt19937).
  std::mt19937 rng_;

  // Lazy SPT cache keyed by vertex id. Used by lookup_lca / lookup_multiple_lca
  // so that the proper LCA-based bound (paper Eq. 7) can be computed without
  // requiring this class to build SPTs eagerly the way ELM_SPTree does. Each
  // entry is a fully preprocessed SPT (RMQ-LCA + DSU Cartesian) for a single
  // landmark endpoint; built on first lookup that needs it.
  std::map<unsigned int, std::vector<shortest_path_tree*>*> lca_spt_cache_;

  // Helper for the new dense-graph-tolerant sampler (paper §4.2's PDF-based
  // alternative). Returns a uniformly random missing edge; falls back to
  // rejection sampling when missing edges are abundant, switches to a
  // PDF-based two-step sampler when the graph approaches dense.
  pair<unsigned int, unsigned int> sample_missing_edge_robust();

public:
  // Constructor.
  // `precompute_all_pairs` (default true for backward compat): runs Dijkstra
  // from every vertex up-front so that sp_vector / le_vector are populated for
  // every (root, target) pair. This is what the original code did. It costs
  // O(n) Dijkstras and 2 * O(n^2) doubles of RAM (about 40 GB for n = 50,000).
  //
  // For workloads where the only required SP rows are those of the landmark
  // endpoints (kLE, kFE, and any post-selection lookup path), pass `false`
  // and call `ensure_sp_for_landmarks()` after landmark selection. That cuts
  // build cost to O(k) Dijkstras and memory to O(k * n) doubles -- e.g.
  // ~200 MB at k = 256 instead of 40 GB.
  EdgeLandMark(vector<list<pair<unsigned int, double>> *> *adj_list,
               unsigned int n_nodes, unsigned int k,
               unsigned int sampling_size,
               bool precompute_all_pairs = true);
  ~EdgeLandMark();
  long _sizeof();
  void store_map();
  void large_edge_heuristic();
  void far_away_heuristic();

  // Populate sp_vector / le_vector rows only for the endpoints currently in
  // `landmarks`. Idempotent: rows already present are not recomputed. Safe to
  // call at any time after a get_*_landmarks() / *_heuristic() call. Required
  // before lookup() / lookup_multiple() / lookup_ub() if the instance was
  // constructed with precompute_all_pairs = false.
  void ensure_sp(unsigned int u);
  void ensure_sp_for_landmarks();
  void get_exact_landmarks();
  void clean_unwanted_shortest_paths();
  void get_landmarks();
  double lookup(unsigned int u, unsigned int v);
  double lookup_ub(unsigned int u, unsigned int v);
  vector<double> *lookup_multiple(unsigned int u, unsigned int v);

  // ---- New methods (additive; no existing call site is affected) ----

  // Seed the internal RNG for reproducible sampling. Existing methods that
  // used rand() pick this up via the new sampler when called via the *_seeded
  // variants below.
  void set_seed(unsigned int seed) { rng_.seed(seed); }

  // Sampling SGEL with seeded mt19937 + dense-graph-safe sampling.
  // Equivalent to get_landmarks() but with reproducible randomness and no
  // pathological behavior near m -> C(n,2). Delegates to the same lookup
  // formula as get_landmarks() so quality is identical when the graph is sparse.
  void get_landmarks_seeded(unsigned int seed);

  // Lazy / accelerated greedy (Minoux 1978) for GEL.
  // Maintains a max-heap of upper-bound marginal gains and recomputes only
  // when the current top is popped. Submodularity guarantees the bound stays
  // valid, so the (1 - 1/e) approximation is preserved while typically
  // doing 10-100x fewer marginal-gain evaluations than the naive scan in
  // get_exact_landmarks(). Output set is not byte-identical to the naive
  // version under ties, but has the same worst-case approximation guarantee.
  void get_exact_landmarks_lazy();

  // Batched: fill LB[v] for all v in [0, n) for a fixed source u.
  // Cheaper than calling lookup(u, v) n-1 times because we walk landmarks
  // once and update all destinations in a single pass per landmark, with
  // better cache locality on the SP / LE arrays.
  void lookup_from_source(unsigned int u, vector<double> &LB);

  // ---- Tighter LCA-based lookup (paper §4.1, equation 7 exactly) ----
  //
  // The class's original lookup() uses
  //   2 * max(LE_x[u], LE_x[v]) - SP_x[u] - SP_x[v]
  // which is the longest edge over the *union* of two root-to-vertex paths
  // in tau_x. This is a valid lower bound but strictly looser than the true
  // vertex-vertex max-edge through the LCA in tau_x.
  //
  // lookup_lca() / lookup_multiple_lca() use the proper LCA-based bound by
  // building (and caching) actual SPTs at landmark endpoints on first call.
  // Same answer as ELM_SPTree::_lookup, so kLE/kFE heuristics now have a
  // fair-quality lookup option without changing the rest of the class.
  double lookup_lca(unsigned int u, unsigned int v);
  vector<double> *lookup_multiple_lca(unsigned int u, unsigned int v);
  // Free cached SPTs (call when changing landmark set or to release memory).
  void clear_lca_cache();
};

#endif // EDGELANDMARK_H
