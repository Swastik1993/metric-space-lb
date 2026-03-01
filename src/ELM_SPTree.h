#ifndef ELM_SPTREE_H
#define ELM_SPTREE_H

#include "Dijkstra.h"
#include "GraphDefinitions.h"
#include "GraphUtils.h"
#include "ShortestPathTree.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <list>
#include <map>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

class ELM_SPTree {
private:
  vector<list<pair<unsigned int, double>> *> *adj_list;
  unsigned int no_nodes;
  unsigned int no_landmarks;
  unsigned int no_samples;

  map<pair<unsigned int, unsigned int>, double> *edges_map;
  map<pair<unsigned int, unsigned int>, double> *landmarks;

  // root -> SPT(root) nodes (indexed by vertex id)
  // OLD: map<unsigned int, vector<shortest_path_tree *> *> sp_tree_map;
  // NEW: direct-indexed (size=no_nodes). nullptr means "not built yet".
  vector<vector<shortest_path_tree*>*> sp_tree_vec;

  // Only used for last round of HCA lookups
  vector<unsigned int> temp_HCA;

  // Temporary storage for HCA
  vector<unsigned int> _hca_information;

  // helper: normalize undirected edge key
  inline pair<unsigned int, unsigned int> norm_edge(unsigned int a, unsigned int b) const {
    return make_pair(min(a, b), max(a, b));
  }

  void store_map();

  // Ensure SPT(root) exists and is preprocessed (RMQ + Cartesian built in Dijkstra(..., lifting=true))
  inline void ensure_spt(unsigned int root) {
    if (root >= no_nodes) return; // defensive
    if (sp_tree_vec[root] == nullptr) {
      sp_tree_vec[root] = Dijkstra(adj_list, root, true);
    }
  }

  // Improved vertex-based bound using SPT structure (calls find_LCA from ShortestPathTree.cpp)
  double _lookup_vertex(unsigned int root_node, unsigned int u, unsigned int v);

  // Improved bound per landmark edge = folding + refined stitching
  double _lookup(unsigned int u, unsigned int v,
                 pair<unsigned int, unsigned int> landmark_edge,
                 double edge_dist);

  double _lookup_vertex_HCA(unsigned int root_node, unsigned int u, unsigned int v);

  double _lookup_HCA(unsigned int u, unsigned int v,
                     pair<unsigned int, unsigned int> landmark_edge,
                     double edge_dist);

public:
  ELM_SPTree(vector<list<pair<unsigned int, double>> *> *adj_list,
             unsigned int nodes,
             unsigned int k,
             unsigned int no_samples);

  ~ELM_SPTree();

  // preprocessing landmarks
  void preprocess_landmark_roots();

  // landmark selection
  void get_sampling_landmarks();
  void get_exact_landmarks();

  // query
  double lookup(unsigned int u, unsigned int v);

  // REQUIRED by main.cpp (per-landmark RMSE curve)
  vector<double> *lookup_multiple(unsigned int u, unsigned int v);

  // Compatibility with old print in main.cpp (proxy, not true heap memory)
  size_t _sizeof() const;

  // optional: cleanup
  void clean_unwanted_shortest_paths();

  void _reset_HCA();
};

#endif



