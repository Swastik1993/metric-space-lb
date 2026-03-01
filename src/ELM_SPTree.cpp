#include "ELM_SPTree.h"

ELM_SPTree::ELM_SPTree(vector<list<pair<unsigned int, double>> *> *adj_list,
                       unsigned int nodes,
                       unsigned int k,
                       unsigned int no_samples)
    : adj_list(adj_list),
      no_nodes(nodes),
      no_landmarks(k),
      no_samples(no_samples),
      edges_map(new map<pair<unsigned int, unsigned int>, double>()),
      landmarks(new map<pair<unsigned int, unsigned int>, double>()),
      sp_tree_vec(nodes, nullptr) {

  // Initially all nodes are invalid
  _hca_information = vector<unsigned int>(no_nodes, no_nodes);

  // Set it to invalid
  temp_HCA = vector<unsigned int>(no_nodes, numeric_limits<unsigned int>::max());

  store_map();
}

ELM_SPTree::~ELM_SPTree() {
  // Free SPT vectors
  for (auto *ptr : sp_tree_vec) {
    delete ptr; // matches existing ownership pattern in your repo
  }
  sp_tree_vec.assign(no_nodes, nullptr);

  delete edges_map;
  delete landmarks;
}

void ELM_SPTree::store_map() {
  edges_map->clear();
  for (unsigned int i = 0; i < adj_list->size(); ++i) {
    for (auto &item : *adj_list->at(i)) {
      unsigned int j = item.first;
      double w = item.second;
      auto e = norm_edge(i, j);
      if (edges_map->find(e) == edges_map->end()) edges_map->insert({e, w});
    }
  }
}

size_t ELM_SPTree::_sizeof() const {
  size_t s = 0;

  s += sizeof(*this);

  s += _hca_information.capacity() * sizeof(unsigned int);
  s += temp_HCA.capacity() * sizeof(unsigned int);

  auto rb_node_bytes = [](size_t value_bytes) -> size_t {
    return value_bytes + 3 * sizeof(void*) + sizeof(char);
  };

  if (edges_map) {
    s += sizeof(*edges_map);
    using KV = std::pair<const std::pair<unsigned int, unsigned int>, double>;
    s += edges_map->size() * rb_node_bytes(sizeof(KV));
  }
  if (landmarks) {
    s += sizeof(*landmarks);
    using KV = std::pair<const std::pair<unsigned int, unsigned int>, double>;
    s += landmarks->size() * rb_node_bytes(sizeof(KV));
  }

  // sp_tree_vec itself
  s += sp_tree_vec.capacity() * sizeof(vector<shortest_path_tree*>*);

  for (auto *vec : sp_tree_vec) {
    if (!vec) continue;

    s += sizeof(*vec);
    s += vec->capacity() * sizeof(shortest_path_tree*);

    for (auto *node : *vec) {
      if (!node) continue;

      s += sizeof(*node);

      if (node->children) {
        s += sizeof(*node->children);
        s += node->children->capacity() * sizeof(shortest_path_tree*);
      }

      if (node->jump_pointers) {
        s += sizeof(*node->jump_pointers);
        s += node->jump_pointers->capacity() *
             sizeof(std::pair<shortest_path_tree*, double>);
      }

      if (node->rmq_vertex) {
        s += sizeof(*node->rmq_vertex);

        s += node->rmq_vertex->euler.capacity() * sizeof(int);
        s += node->rmq_vertex->depth.capacity() * sizeof(int);
        s += node->rmq_vertex->first.capacity() * sizeof(int);
        s += node->rmq_vertex->lg.capacity() * sizeof(int);

        s += node->rmq_vertex->st.capacity() * sizeof(std::vector<int>);
        for (auto const &row : node->rmq_vertex->st) {
          s += row.capacity() * sizeof(int);
        }
      }

      if (node->cartesian) {
        s += sizeof(*node->cartesian);

        s += node->cartesian->leaf.capacity() * sizeof(int);

        s += node->cartesian->nodes.capacity() * sizeof(CTNode*);
        s += node->cartesian->nodes.size() * sizeof(CTNode);

        auto &rmq = node->cartesian->rmq;
        s += rmq.euler.capacity() * sizeof(int);
        s += rmq.depth.capacity() * sizeof(int);
        s += rmq.first.capacity() * sizeof(int);
        s += rmq.lg.capacity() * sizeof(int);
        s += rmq.st.capacity() * sizeof(std::vector<int>);
        for (auto const &row : rmq.st) {
          s += row.capacity() * sizeof(int);
        }
      }
    }
  }

  return s;
}

double ELM_SPTree::_lookup_vertex(unsigned int root_node, unsigned int u, unsigned int v) {
  // ensure_spt(root_node);

  auto *tree = sp_tree_vec[root_node];

  pair<shortest_path_tree *, double> lca_info =
      find_LCA(tree, u, v);

  shortest_path_tree *lca_node = lca_info.first;
  double max_edge_on_path = lca_info.second;

  double dist_root_u   = tree->at(u)->path_length;
  double dist_root_v   = tree->at(v)->path_length;
  double dist_root_lca = lca_node->path_length;

  double tree_dist_uv = dist_root_u + dist_root_v - 2.0 * dist_root_lca;

  return 2.0 * max_edge_on_path - tree_dist_uv;
}

double ELM_SPTree::_lookup(unsigned int u, unsigned int v,
                           pair<unsigned int, unsigned int> landmark_edge,
                           double edge_dist) {
  unsigned int a = landmark_edge.first;
  unsigned int b = landmark_edge.second;

  auto *treeA = sp_tree_vec[a];
  auto *treeB = sp_tree_vec[b];

  double au = treeA->at(u)->path_length;
  double av = treeA->at(v)->path_length;
  double bu = treeB->at(u)->path_length;
  double bv = treeB->at(v)->path_length;

  double best = 0.0;
  best = max(best, edge_dist - au - bv);
  best = max(best, edge_dist - av - bu);

  best = max(best, _lookup_vertex(a, u, v));
  best = max(best, _lookup_vertex(b, u, v));

  return best;
}

double ELM_SPTree::_lookup_vertex_HCA(unsigned int root_node, unsigned int u, unsigned int v) {
  // HCA based computations
  if (_hca_information.at(root_node) == no_nodes) {
    _hca_information.at(root_node) = find_HCA(
        sp_tree_vec[u], sp_tree_vec[v], root_node, temp_HCA
    );
  }

  unsigned int hca_node = _hca_information.at(root_node);

  double lca_longest_edge = max(sp_tree_vec[u]->at(hca_node)->max_edge,
                                sp_tree_vec[v]->at(hca_node)->max_edge);

  double full_path_length = sp_tree_vec[u]->at(hca_node)->path_length
                          + sp_tree_vec[v]->at(hca_node)->path_length;

  return 2 * lca_longest_edge - full_path_length;
}

// Evaluation of possible landmark candidate "landmark"
double ELM_SPTree::_lookup_HCA(unsigned int u, unsigned int v,
                               pair<unsigned int, unsigned int> landmark,
                               double edge_dist) {
  double max_dist = 0.;

  double l1u = sp_tree_vec[u]->at(landmark.first)->path_length;
  double l1v = sp_tree_vec[v]->at(landmark.first)->path_length;
  double l2u = sp_tree_vec[u]->at(landmark.second)->path_length;
  double l2v = sp_tree_vec[v]->at(landmark.second)->path_length;

  max_dist = max(max_dist, max(edge_dist - l1u - l2v, edge_dist - l1v - l2u));

  max_dist = max(max_dist, _lookup_vertex_HCA(landmark.first, u, v));
  max_dist = max(max_dist, _lookup_vertex_HCA(landmark.second, u, v));

  return max_dist;
}

// inside ELM_SPTree
void ELM_SPTree::preprocess_landmark_roots() {
  for (auto &lm : *landmarks) {
    unsigned a = lm.first.first;
    unsigned b = lm.first.second;

    ensure_spt(a);
    ensure_spt(b);

    preprocess_spt_queries(sp_tree_vec[a]);
    preprocess_spt_queries(sp_tree_vec[b]);
  }
}

double ELM_SPTree::lookup(unsigned int u, unsigned int v) {
  auto key = norm_edge(u, v);
  auto it = edges_map->find(key);
  if (it != edges_map->end()) return it->second;

  double best = 0.0;
  for (auto &lm : *landmarks) {
    best = max(best, _lookup(u, v, lm.first, lm.second));
  }
  return best;
}

vector<double> *ELM_SPTree::lookup_multiple(unsigned int u, unsigned int v) {
  vector<double> *vals = new vector<double>();
  vals->reserve(no_landmarks);

  auto key = norm_edge(u, v);
  auto it = edges_map->find(key);
  if (it != edges_map->end()) {
    for (unsigned int i = 0; i < no_landmarks; ++i) vals->push_back(it->second);
    return vals;
  }

  double best = 0.0;
  unsigned int count = 0;
  for (auto &lm : *landmarks) {
    best = max(best, _lookup(u, v, lm.first, lm.second));
    vals->push_back(best);
    ++count;
    if (count == no_landmarks) break;
  }

  while (vals->size() < no_landmarks) vals->push_back(best);

  return vals;
}

void ELM_SPTree::clean_unwanted_shortest_paths() {
  unordered_map<unsigned int, bool> keep;
  keep.reserve(2 * landmarks->size() + 8);

  for (auto &lm : *landmarks) {
    keep[lm.first.first] = true;
    keep[lm.first.second] = true;
  }

  for (unsigned int r = 0; r < no_nodes; ++r) {
    if (sp_tree_vec[r] == nullptr) continue;
    if (keep.find(r) != keep.end()) continue;

    delete sp_tree_vec[r];
    sp_tree_vec[r] = nullptr;
  }
}

void ELM_SPTree::get_exact_landmarks() {
  landmarks->clear();

  // exact version builds SPTs for all nodes
  for (unsigned int i = 0; i < no_nodes; ++i) ensure_spt(i);

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    double best_gain = -1.0;
    pair<unsigned int, unsigned int> best_edge = {0, 0};
    double best_len = 0.0;

    for (auto &cand : *edges_map) {
      if (landmarks->find(cand.first) != landmarks->end()) continue;

      const auto &e = cand.first;
      double e_len = cand.second;

      ensure_spt(e.first);
      ensure_spt(e.second);

      double gain = 0.0;

      for (unsigned int u = 0; u < no_nodes; ++u) {
        for (unsigned int v = u + 1; v < no_nodes; ++v) {
          auto uv = norm_edge(u, v);
          if (edges_map->find(uv) != edges_map->end()) continue;

          double curr = lookup(u, v);
          double cand_lb = _lookup(u, v, e, e_len);
          if (cand_lb > curr) gain += (cand_lb - curr);
        }
      }

      if (gain > best_gain) {
        best_gain = gain;
        best_edge = e;
        best_len = e_len;
      }
    }

    landmarks->insert({best_edge, best_len});
  }

  clean_unwanted_shortest_paths();
}

void ELM_SPTree::get_sampling_landmarks() {
  landmarks->clear();

  std::mt19937 rng(73);
  std::uniform_int_distribution<unsigned int> uni(0, no_nodes - 1);

  cout << "==========>> LANDMARKS:  " << landmarks->size() << endl;
  cout << "==========>> no_samples: " << no_samples << endl;

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    if (landmarks->size() > 0) {
      this->_reset_HCA();
      cout << "Selected " << landmarks->size() << " number of landmarks" << endl;
    }

    vector<pair<unsigned int, unsigned int>> samples;
    samples.reserve(no_samples);

    while (samples.size() < no_samples) {
      unsigned int u = uni(rng);
      unsigned int v = uni(rng);
      if (u == v) continue;
      auto uv = norm_edge(u, v);
      if (edges_map->find(uv) != edges_map->end()) continue;
      samples.push_back(uv);
    }

    double best_gain = -1.0;
    pair<unsigned int, unsigned int> best_edge = {0, 0};
    double best_len = 0.0;

    map<pair<unsigned int, unsigned int>, double> contributions;

    for (auto &uv : samples) {
      this->_reset_HCA();
      unsigned int u = uv.first;
      unsigned int v = uv.second;

      ensure_spt(u);
      ensure_spt(v);

      double curr = lookup(u, v);

      for (auto &cand : *edges_map) {
        if (landmarks->find(cand.first) != landmarks->end()) continue;

        const auto &e = cand.first;
        double e_len = cand.second;

        double cand_lb = _lookup_HCA(u, v, e, e_len);

        if (cand_lb > curr) {
          if (!contributions.contains(uv)) contributions.insert({uv, 0.0});
          contributions[uv] += cand_lb - curr;

          if (best_gain < contributions[uv]) {
            best_gain = contributions[uv];
            best_edge = e;
            best_len = e_len;
          }
        }
      }
    }

    ensure_spt(best_edge.first);
    ensure_spt(best_edge.second);
    landmarks->insert({best_edge, best_len});
  }

  clean_unwanted_shortest_paths();

  vector<unsigned int>().swap(_hca_information);
  vector<unsigned int>().swap(temp_HCA);
}

void ELM_SPTree::_reset_HCA() {
  _hca_information.assign(_hca_information.size(), no_nodes);
  temp_HCA.assign(temp_HCA.size(), numeric_limits<unsigned int>::max());
}

