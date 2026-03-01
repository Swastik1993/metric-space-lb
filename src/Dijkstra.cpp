#include "Dijkstra.h"
#include <boost/heap/pairing_heap.hpp>

struct _sps {
  unsigned int node, parent;
  double path_length;
  double max_edge;

  _sps(unsigned int _node, unsigned int _parent, double _path_len,
       double _max_edge) {
    node = _node;
    parent = _parent;
    path_length = _path_len;
    max_edge = _max_edge;
  }

  bool operator<(_sps const &sps2) const {
    return path_length > sps2.path_length;
  }
};

vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, bool lifting) {
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
        H.push(_sps(it->first, source, it->second, it->second));
  }
  while (!H.empty()) {
    unsigned int dest = H.top().node;
    unsigned int parent = H.top().parent;
    double dist = H.top().path_length;
    double max_edge = H.top().max_edge;
    H.pop();

    if (visited.at(dest)) {
      continue;
    }

    visited.at(dest) = true;

    shortest_path_tree *sp_tree_temp = new shortest_path_tree(
        dest, 1 + tree_nodes->at(parent)->depth, dist, max_edge);
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
                                            max(max_edge, it->second)));
      } else if (handles.at(neighbour).node_->value.path_length > total) {
        H.increase(handles.at(neighbour),
                   _sps(neighbour, dest, total,
                        max(max_edge, it->second)));
      }
    }
  }

  if (lifting) {
    binary_lifting(sp_tree_root);
    // Replace jump-pointer LCA with RMQ-LCA + DSU Cartesian max-edge queries
    // (jump_pointers are still built for HCA queries)
    preprocess_spt_queries(tree_nodes);
  }

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
      } else if (handles.at(neighbour).node_->value.path_length > total) {
        H.increase(handles.at(neighbour),
                   _sps(neighbour, dest, total, max(max_edge, it->second)));
      }
    }
  }

  return make_pair(dist_vec, max_edge_vec);
}
