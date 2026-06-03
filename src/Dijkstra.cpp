#include "Dijkstra.h"
#include <boost/heap/pairing_heap.hpp>
#include <ctime>

// Build-phase timing accumulators (see Dijkstra.h). File-scope; updated only on
// the sequential build path used by the ablation report.
namespace {
uint64_t g_dj_tree_us = 0, g_dj_lift_us = 0, g_dj_rmq_us = 0;
inline uint64_t dj_now_us() {
  struct timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  return uint64_t(ts.tv_sec) * 1000000ULL + ts.tv_nsec / 1000ULL;
}
}  // namespace

void     dijkstra_timing_reset() { g_dj_tree_us = g_dj_lift_us = g_dj_rmq_us = 0; }
uint64_t dijkstra_tree_us() { return g_dj_tree_us; }
uint64_t dijkstra_lift_us() { return g_dj_lift_us; }
uint64_t dijkstra_rmq_us()  { return g_dj_rmq_us; }

struct _sps {
  unsigned int node, parent;
  double path_length;
  double max_edge;
  double parent_w;  // weight of the edge from parent to this node

  _sps(unsigned int _node, unsigned int _parent, double _path_len,
       double _max_edge, double _parent_w = 0.0) {
    node = _node;
    parent = _parent;
    path_length = _path_len;
    max_edge = _max_edge;
    parent_w = _parent_w;
  }

  bool operator<(_sps const &sps2) const {
    return path_length > sps2.path_length;
  }
};

// Internal: just builds the raw SPT (no auxiliary structures). Returned
// vector is indexed by vertex id; nullptr entries indicate unreachable nodes.
static vector<shortest_path_tree *> *
DijkstraCore(vector<list<pair<unsigned int, double>> *> *adj_lst,
             unsigned int source) {
  unsigned int nodes = adj_lst->size();
  vector<bool> visited(nodes, false);
  visited.at(source) = true;
  boost::heap::pairing_heap<_sps> H;
  shortest_path_tree *sp_tree_root =
      new shortest_path_tree(source, 0, 0., 0.);
  sp_tree_root->parent = nullptr;
  sp_tree_root->root = sp_tree_root;
  vector<shortest_path_tree *> *tree_nodes =
      new vector<shortest_path_tree *>(nodes, nullptr);

  tree_nodes->at(source) = sp_tree_root;

  vector<boost::heap::pairing_heap<_sps>::handle_type> handles;

  for (unsigned int j = 0; j < nodes; ++j) {
    handles.push_back((boost::heap::pairing_heap<_sps>::handle_type)NULL);
  }

  for (list<pair<unsigned int, double>>::iterator it =
           adj_lst->at(source)->begin();
       it != adj_lst->at(source)->end(); ++it) {
    handles.at(it->first) =
        H.push(_sps(it->first, source, it->second, it->second, it->second));
  }
  while (!H.empty()) {
    unsigned int dest = H.top().node;
    unsigned int parent = H.top().parent;
    double dist = H.top().path_length;
    double max_edge = H.top().max_edge;
    double parent_w = H.top().parent_w;
    H.pop();

    if (visited.at(dest)) {
      continue;
    }

    visited.at(dest) = true;

    shortest_path_tree *sp_tree_temp = new shortest_path_tree(
        dest, 1 + tree_nodes->at(parent)->depth, dist, max_edge);
    sp_tree_temp->parent_edge_w = parent_w;  // store directly; no subtract+clamp
    sp_tree_temp->parent = tree_nodes->at(parent);
    sp_tree_temp->root = sp_tree_root;
    if (tree_nodes->at(parent)->children == nullptr) {
      tree_nodes->at(parent)->children = new vector<shortest_path_tree *>();
    }
    tree_nodes->at(parent)->children->push_back(sp_tree_temp);
    tree_nodes->at(dest) = sp_tree_temp;

    for (list<pair<unsigned int, double>>::iterator it =
             adj_lst->at(dest)->begin();
         it != adj_lst->at(dest)->end(); ++it) {
      unsigned int neighbour = it->first;
      if (visited.at(neighbour)) {
        continue;
      }
      double total = it->second + dist;
      if (handles.at(neighbour) ==
          (boost::heap::pairing_heap<_sps>::handle_type)NULL) {
        handles.at(neighbour) = H.push(_sps(neighbour, dest, total,
                                            max(max_edge, it->second),
                                            it->second));
      } else if ((*handles.at(neighbour)).path_length > total) {
        H.increase(handles.at(neighbour),
                   _sps(neighbour, dest, total,
                        max(max_edge, it->second), it->second));
      }
    }
  }

  return tree_nodes;
}

vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, bool lifting) {
  // Backward-compatible wrapper: lifting=true => build everything (the
  // original behavior); lifting=false => skip all auxiliary structures.
  SPTPrepFlags flags;
  if (!lifting) {
    flags.need_lifting = flags.need_lca_rmq = flags.need_cartesian = false;
  }
  return Dijkstra(adj_lst, source, flags);
}

vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, SPTPrepFlags flags) {
  uint64_t _t0 = dj_now_us();
  vector<shortest_path_tree *> *tree_nodes = DijkstraCore(adj_lst, source);
  // tree_nodes->at(source) is the root by construction
  shortest_path_tree *sp_tree_root = tree_nodes->at(source);
  g_dj_tree_us += (dj_now_us() - _t0);

  if (flags.need_lifting && sp_tree_root) {
    uint64_t _l0 = dj_now_us();
    binary_lifting(sp_tree_root);
    g_dj_lift_us += (dj_now_us() - _l0);
  }
  uint64_t _r0 = dj_now_us();
  if (flags.need_lca_rmq) {
    preprocess_spt_lca_rmq(tree_nodes);
  }
  if (flags.need_cartesian) {
    preprocess_spt_cartesian(tree_nodes);
  }
  g_dj_rmq_us += (dj_now_us() - _r0);
  return tree_nodes;
}

// Pooled SPT builder (#18). Same algorithm as DijkstraCore but every
// shortest_path_tree node is allocated from the supplied pool's contiguous
// storage. The pool must have been reserve()'d to at least `nodes` capacity
// before the call (or have headroom for additional growth).
static vector<shortest_path_tree *> *
DijkstraCorePooled(vector<list<pair<unsigned int, double>> *> *adj_lst,
                   unsigned int source, SPTNodePool *pool) {
  unsigned int nodes = adj_lst->size();
  // Ensure capacity so emplace_back never invalidates pointers.
  if (pool->storage.capacity() < pool->storage.size() + nodes) {
    pool->reserve(pool->storage.size() + nodes);
  }
  vector<bool> visited(nodes, false);
  visited.at(source) = true;
  boost::heap::pairing_heap<_sps> H;
  shortest_path_tree *sp_tree_root = pool->make(source, 0, 0., 0.);
  sp_tree_root->parent = nullptr;
  sp_tree_root->root = sp_tree_root;
  vector<shortest_path_tree *> *tree_nodes =
      new vector<shortest_path_tree *>(nodes, nullptr);
  tree_nodes->at(source) = sp_tree_root;

  vector<boost::heap::pairing_heap<_sps>::handle_type> handles(nodes);
  for (auto it = adj_lst->at(source)->begin(); it != adj_lst->at(source)->end(); ++it) {
    handles.at(it->first) =
        H.push(_sps(it->first, source, it->second, it->second, it->second));
  }
  while (!H.empty()) {
    unsigned int dest = H.top().node;
    unsigned int parent = H.top().parent;
    double dist = H.top().path_length;
    double max_edge = H.top().max_edge;
    double parent_w = H.top().parent_w;
    H.pop();
    if (visited.at(dest)) continue;
    visited.at(dest) = true;
    shortest_path_tree *node = pool->make(
        dest, 1 + tree_nodes->at(parent)->depth, dist, max_edge);
    node->parent_edge_w = parent_w;
    node->parent = tree_nodes->at(parent);
    node->root = sp_tree_root;
    if (tree_nodes->at(parent)->children == nullptr) {
      tree_nodes->at(parent)->children = new vector<shortest_path_tree *>();
    }
    tree_nodes->at(parent)->children->push_back(node);
    tree_nodes->at(dest) = node;
    for (auto it = adj_lst->at(dest)->begin(); it != adj_lst->at(dest)->end(); ++it) {
      unsigned int neighbour = it->first;
      if (visited.at(neighbour)) continue;
      double total = it->second + dist;
      if (handles.at(neighbour) ==
          (boost::heap::pairing_heap<_sps>::handle_type)NULL) {
        handles.at(neighbour) = H.push(_sps(neighbour, dest, total,
                                            max(max_edge, it->second),
                                            it->second));
      } else if ((*handles.at(neighbour)).path_length > total) {
        H.increase(handles.at(neighbour),
                   _sps(neighbour, dest, total,
                        max(max_edge, it->second), it->second));
      }
    }
  }
  return tree_nodes;
}

vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, SPTPrepFlags flags, SPTNodePool *pool) {
  if (!pool) return Dijkstra(adj_lst, source, flags);
  vector<shortest_path_tree *> *tree_nodes = DijkstraCorePooled(adj_lst, source, pool);
  shortest_path_tree *sp_tree_root = tree_nodes->at(source);
  if (flags.need_lifting && sp_tree_root) binary_lifting(sp_tree_root);
  if (flags.need_lca_rmq) preprocess_spt_lca_rmq(tree_nodes);
  if (flags.need_cartesian) preprocess_spt_cartesian(tree_nodes);
  return tree_nodes;
}

pair<vector<double> *, vector<double> *>
DijkstraELM(vector<list<pair<unsigned int, double>> *> *adj_lst,
            unsigned int source) {
  unsigned int nodes = adj_lst->size();
  vector<bool> visited(nodes, false);
  visited.at(source) = true;
  boost::heap::pairing_heap<_sps> H;

  vector<double> *dist_vec = new vector<double>(nodes, 0.0);
  vector<double> *max_edge_vec = new vector<double>(nodes, 0.0);

  vector<boost::heap::pairing_heap<_sps>::handle_type> handles;

  for (unsigned int j = 0; j < nodes; ++j) {
    handles.push_back((boost::heap::pairing_heap<_sps>::handle_type)NULL);
  }

  for (list<pair<unsigned int, double>>::iterator it =
           adj_lst->at(source)->begin();
       it != adj_lst->at(source)->end(); ++it) {
    handles.at(it->first) =
        H.push(_sps(it->first, source, it->second, it->second));
  }
  while (!H.empty()) {
    unsigned int dest = H.top().node;
    double dist = H.top().path_length;
    double max_edge = H.top().max_edge;
    H.pop();

    if (visited.at(dest)) {
      continue;
    }

    visited.at(dest) = true;
    dist_vec->at(dest) = dist;
    max_edge_vec->at(dest) = max_edge;

    for (list<pair<unsigned int, double>>::iterator it =
             adj_lst->at(dest)->begin();
         it != adj_lst->at(dest)->end(); ++it) {
      unsigned int neighbour = it->first;
      if (visited.at(neighbour)) {
        continue;
      }
      double total = it->second + dist;
      if (handles.at(neighbour) ==
          (boost::heap::pairing_heap<_sps>::handle_type)NULL) {
        handles.at(neighbour) =
            H.push(_sps(neighbour, dest, total, max(max_edge, it->second)));
      } else if ((*handles.at(neighbour)).path_length > total) {
        H.increase(handles.at(neighbour),
                   _sps(neighbour, dest, total, max(max_edge, it->second)));
      }
    }
  }

  return make_pair(dist_vec, max_edge_vec);
}
