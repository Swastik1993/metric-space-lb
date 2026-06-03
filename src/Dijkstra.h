#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "ShortestPathTree.h"
#include <cstdint>
#include <list>
#include <utility>
#include <vector>


using namespace std;

// Flags controlling which auxiliary structures Dijkstra builds on top of the
// raw SPT. Building them is paid up-front per tree; queries pick which ones
// they consume. Defaults preserve original behavior (lifting + LCA + Cartesian
// all built whenever lifting=true).
//
// - need_lifting:    binary-lifting jump pointers, used by find_HCA.
// - need_lca_rmq:    Euler-tour RMQ over the SPT, used by find_LCA's vertex side.
// - need_cartesian:  DSU Cartesian tree + its own RMQ, used by find_LCA's
//                    max-edge side.
//
// If you only call find_LCA, the lifting work is wasted; if you only call
// find_HCA, the LCA RMQ + Cartesian work is wasted. Use the explicit overload
// to skip the unused half.
struct SPTPrepFlags {
  bool need_lifting   = true;
  bool need_lca_rmq   = true;
  bool need_cartesian = true;
};

vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, bool lifting = true);

// Selective overload: builds only the requested structures.
vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, SPTPrepFlags flags);

// Pooled overload (#18). Uses `pool` for all SPT node allocations. The pool
// must outlive the returned vector. The destructor of shortest_path_tree
// (which deletes children/jump_pointers/rmq_vertex/cartesian) is still
// called by ~SPTNodePool when the pool is destroyed; the *node objects*
// themselves live in the pool's contiguous buffer and are not individually
// freed.
//
// Caveat: the returned vector<shortest_path_tree*> aliases the pool's
// internal storage. Do NOT call `delete *t;` on its elements -- that would
// double-free against the pool. Use `delete tree_nodes` (vector) only.
vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, SPTPrepFlags flags, SPTNodePool *pool);

pair<vector<double> *, vector<double> *>
DijkstraELM(vector<list<pair<unsigned int, double>> *> *adj_lst,
            unsigned int source);

// ---- Build-phase timing for the ablation report (R3.D7) --------------------
// Dijkstra(...,lifting=true) builds the SPT tree AND its auxiliary structures
// (binary lifting, Euler/sparse-table LCA-RMQ, Cartesian RMQ) in one call.
// These accumulators attribute cumulative process-CPU microseconds to the
// three phases so the report can separate SPT / LIFT / RMQ construction. Only
// the sequential (ensure_spt) build path is instrumented; safe to read after a
// sequential build. Call dijkstra_timing_reset() to zero them.
void     dijkstra_timing_reset();
uint64_t dijkstra_tree_us();   // core shortest-path-tree construction
uint64_t dijkstra_lift_us();   // binary_lifting (jump pointers, for find_HCA)
uint64_t dijkstra_rmq_us();    // lca_rmq + cartesian (LCA tables, for find_LCA)

#endif // DIJKSTRA_H
