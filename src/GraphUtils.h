#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include "GraphDefinitions.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/random/linear_congruential.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


using namespace std;

int check_connected(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id);
unsigned long count_better(vector<double> *lower_bounds,
                           vector<double> *upper_bounds);
vector<list<pair<unsigned int, double>> *> *get_adjacency_list(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id,
    vector<vector<double> *> *dist);
vector<vector<double> *> *get_adjacency_matrix(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id,
    vector<vector<double> *> *dist, double default_missing = -1.);
vector<vector<double> *> *distance_matrix(unsigned int nodes, unsigned int dims,
                                          unsigned p = 2);
void clean_up_adj_matrix(vector<vector<double> *> *adj_matrix);
void clean_up_adj_list(vector<list<pair<unsigned int, double>> *> *adj_lst);
map<pair<unsigned int, unsigned int>, double> *convert_adjList_to_knownEdges(
    vector<list<pair<unsigned int, double>> *> *adj_lst);
vector<list<pair<unsigned int, double>> *> *get_adj_list_file(char *filename);
vector<vector<double> *> *get_adj_matrix_file(char *filename,
                                              double default_val = -1.);
Graph waxman_graph_distances(unsigned int n, vector<vector<double> *> *distance,
                             double beta = 0.4, double alpha = 0.1);
pair<Graph, vector<vector<double> *> *>
waxman_graph(unsigned int n, double beta = 0.4, double alpha = 0.1,
             unsigned int p = 1, unsigned int dims = 3);

// ============================================================================
// Optimization (#17): Flat row-major n*n matrix backed by a single allocation.
// ============================================================================
// vector<vector<double>*> requires two pointer dereferences per access and
// scatters its rows across the heap; for n in the thousands this hurts both
// cache locality and allocator pressure. FlatMatrix presents the same logical
// shape (i, j -> double) but stores all n*n cells in one contiguous buffer.
//
// Usage:
//   FlatMatrix M(n, /*default=*/-1.0);
//   M(i, j) = 0.42;
//   double x = M(i, j);
//   double *row_i = M.row(i);    // pointer to the i-th row's contiguous data
//
// This is offered as an opt-in helper; existing vector<vector<double>*> code
// paths are unchanged.
class FlatMatrix {
public:
  FlatMatrix() : n_(0) {}
  FlatMatrix(unsigned int n, double init_val) : n_(n), data_(size_t(n)*n, init_val) {}
  void resize(unsigned int n, double init_val) {
    n_ = n;
    data_.assign(size_t(n)*n, init_val);
  }
  inline double &operator()(unsigned int i, unsigned int j) { return data_[size_t(i)*n_ + j]; }
  inline double  operator()(unsigned int i, unsigned int j) const { return data_[size_t(i)*n_ + j]; }
  inline double *row(unsigned int i) { return data_.data() + size_t(i)*n_; }
  inline const double *row(unsigned int i) const { return data_.data() + size_t(i)*n_; }
  unsigned int size() const { return n_; }
  size_t bytes() const { return data_.size() * sizeof(double); }
private:
  unsigned int n_;
  std::vector<double> data_;
};

// Build a FlatMatrix from a vector<vector<double>*> (one-time copy). For
// callers who want to migrate from the legacy storage to flat storage in a
// hot loop without restructuring upstream code.
FlatMatrix flat_matrix_from_vec_of_vec(vector<vector<double>*> *src,
                                       double default_missing = -1.0);

// ============================================================================
// Cached SPLUB helper (paper §6 baseline: Augustine et al. SPLUB).
// ============================================================================
// SPLUB on a queried missing edge (a, b) requires Dijkstra from BOTH a and b
// followed by a fold-over-all-known-edges scan. Naively, sampling K queries
// costs 2K Dijkstras. With a cache keyed on the source vertex, repeated
// endpoints amortize to one Dijkstra per unique vertex -- substantially
// fewer than 2K when K >> n / E[unique].
//
// Usage:
//   SplubCache cache(adj_lst);
//   double lb = cache.splub(a, b);   // caches DijkstraELM(a) and DijkstraELM(b)
//   cache.clear();                   // free cached SP arrays
//
// The cache holds vectors of doubles keyed by source. Each cached entry is
// O(n) doubles. Memory footprint scales with the number of distinct
// endpoints actually queried, not the total node count.
class SplubCache {
public:
  explicit SplubCache(vector<list<pair<unsigned int, double>> *> *adj_lst);
  ~SplubCache();
  double splub(unsigned int a, unsigned int b);
  // Read cached row for `src` (returns nullptr if not cached). Used by
  // post-build RMSE/ERC code that just wants the exact distance src->v
  // without rerunning the O(m) splub() loop. Public for direct access.
  const vector<double>* ensure(unsigned int src);
  // Check whether `src` is already cached, without triggering a Dijkstra.
  // Used by post-build code that wants to avoid serial Dijkstra cost on
  // freshly-sampled endpoints.
  bool has(unsigned int src) const {
    return sp_cache_.find(src) != sp_cache_.end();
  }
  // Build cache rows for every source in [0, n) in parallel. After this,
  // every subsequent splub()/ensure() call is a guaranteed cache hit and
  // costs O(m) (splub) or O(1) (ensure). Cost: n Dijkstras / OMP_THREADS.
  // Memory: 8 * n * n bytes for the dense cache.
  void prewarm_all_parallel(unsigned int n);
  void clear();
private:
  vector<list<pair<unsigned int, double>> *> *adj_lst_;
  // Per-source cache: dist[s] = vector of n SP distances from s (or nullptr).
  unordered_map<unsigned int, vector<double>*> sp_cache_;
};

#endif // GRAPHUTILS_H
