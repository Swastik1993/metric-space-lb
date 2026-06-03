#ifndef SHORTEST_PATH_TREE_H
#define SHORTEST_PATH_TREE_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stack>
#include <utility>
#include <vector>

using namespace std;

// ----------------------------
// RMQ-based LCA (Sparse Table)
// ----------------------------
struct RMQLCA {
  // Euler tour (indices into some node-array)
  vector<int> euler;
  // Depth for each Euler position
  vector<int> depth;
  // first occurrence per node-index
  vector<int> first;
  // log table and sparse table over Euler positions (stores index in euler[])
  vector<int> lg;
  vector<vector<int>> st;

  RMQLCA() = default;
};

// Forward declaration for Cartesian-tree nodes
struct CTNode;

struct CartesianPreproc {
  // owns all nodes (so we can delete them)
  vector<CTNode *> nodes;
  // leaf index per original vertex id (size = n_vertices)
  vector<int> leaf;
  // RMQ for LCA on the Cartesian tree
  RMQLCA rmq;
  // root node index in nodes[]
  int root = -1;

  ~CartesianPreproc();
};

// Node for Cartesian tree (internal nodes correspond to edges; leaves correspond to vertices)
struct CTNode {
  int idx;              // index in CartesianPreproc::nodes
  double weight;        // edge weight (internal), or 0 for leaves
  CTNode *left = nullptr;
  CTNode *right = nullptr;
};

// ----------------------------
// Shortest path tree (SPT)
// ----------------------------
struct shortest_path_tree {
  unsigned int id, depth;
  double path_length;
  double max_edge; // still used by existing HCA code

  // Stored weight of the edge to this node's parent. For the root this is
  // 0.0. Storing it directly avoids the subtract-and-clamp pattern
  // (path_length - parent->path_length) that is sensitive to floating-point
  // drift on graphs with many short edges. Populated during Dijkstra.
  double parent_edge_w = 0.0;

  shortest_path_tree *parent = nullptr;
  shortest_path_tree *root = nullptr;

  // Existing binary-lifting storage (still used by find_HCA)
  vector<pair<shortest_path_tree *, double>> *jump_pointers = nullptr;
  vector<shortest_path_tree *> *children = nullptr;

  // New preprocessing (stored only at root)
  RMQLCA *rmq_vertex = nullptr;           // LCA over vertices in the SPT
  CartesianPreproc *cartesian = nullptr;  // DSU-built Cartesian tree + RMQ

  // parent edge weight is intentionally NOT stored. It can be recovered as:
  //   w(parent,node) = node.path_length - parent.path_length
  // for shortest-path trees produced by Dijkstra.
  shortest_path_tree(int id_, int depth_, double path_length_, double max_edge_)
      : id(id_), depth(depth_), path_length(path_length_), max_edge(max_edge_) {}

  ~shortest_path_tree();
};

// Existing binary lifting (kept because HCA code depends on jump_pointers)
void binary_lifting(shortest_path_tree *root);

// New: preprocess RMQ-LCA on vertices + DSU Cartesian tree for max-edge queries
void preprocess_spt_queries(vector<shortest_path_tree *> *sp_tree);

// Granular preprocessing: build only the structures the caller will use.
// preprocess_spt_queries above is equivalent to calling both of these.
// Useful when the SGEL/HCA pipeline only needs jump_pointers, or when the
// LCA pipeline only needs the RMQ + Cartesian half.
void preprocess_spt_lca_rmq(vector<shortest_path_tree *> *sp_tree);
void preprocess_spt_cartesian(vector<shortest_path_tree *> *sp_tree);

// ============================================================================
// Optimization (#18): Optional SPT node pool.
// ============================================================================
// Each shortest_path_tree node is an individual `new`'d allocation today.
// For n nodes per tree and k landmarks, that's O(kn) heap ops + fragmentation.
// SPTNodePool owns a single contiguous std::vector<shortest_path_tree> and
// hands out stable pointers into it (vector is reserved up-front to capacity
// = n so no rehashing ever moves the storage). Use is opt-in: pass the pool
// to the new Dijkstra overload below; existing call sites keep using the
// per-node allocator and are unaffected.
struct SPTNodePool {
  std::vector<shortest_path_tree> storage;
  // Reserve capacity exactly once. Critical: never let storage grow past
  // capacity, or the node pointers we hand out get invalidated.
  void reserve(size_t n) { storage.reserve(n); }
  // Allocate a node. Returns a stable pointer (since reserve was called).
  shortest_path_tree *make(unsigned int id, unsigned int depth,
                           double path_length, double max_edge) {
    storage.emplace_back(id, depth, path_length, max_edge);
    return &storage.back();
  }
  size_t size() const { return storage.size(); }
};

// LCA query returning (vertex-LCA node, maximum edge weight along the path between the two vertices)
pair<shortest_path_tree *, double>
find_LCA(vector<shortest_path_tree *> *sp_tree, unsigned int index_i,
         unsigned int index_j);

// Existing HCA query (still uses jump_pointers)
unsigned int find_HCA(vector<shortest_path_tree *> *sp_tree_i,
                      vector<shortest_path_tree *> *sp_tree_j,
                      unsigned int index);

// Does not use jump pointers
unsigned int find_HCA(vector<shortest_path_tree *> *sp_tree_i,
                      vector<shortest_path_tree *> *sp_tree_j,
                      unsigned int index,
                      vector<unsigned int>& memoized);

#endif // SHORTEST_PATH_TREE_H
