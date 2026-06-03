#include "ELM_SPTree.h"
#include <queue>
#include <fstream>
#include <sstream>
#include <string>

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
  // Optimization (#20): defer HCA buffer allocation until HCA is actually
  // used. The HCA pipeline currently activates only when get_hca_landmarks()
  // (or _lookup_HCA via that path) is called. Constructing zero-length
  // vectors costs nothing here; the memory bill arrives at first HCA use.
  _hca_information.clear();
  temp_HCA.clear();
  store_map();
}

// Lazy HCA buffer initializer. Idempotent.
static inline void elm_init_hca_buffers(vector<unsigned int> &hca_info,
                                        vector<unsigned int> &temp,
                                        unsigned int no_nodes) {
  if (hca_info.size() == no_nodes) return;
  hca_info.assign(no_nodes, no_nodes);
  temp.assign(no_nodes, std::numeric_limits<unsigned int>::max());
}

// Generation-aware initializer used by the serial HCA fast path. Allocates
// the value buffers and a parallel _gen array sized identically. Once
// allocated, future _reset_HCA() calls only bump the generation counter.
static inline void elm_init_hca_buffers_with_gen(vector<unsigned int> &hca_info,
                                                 vector<unsigned int> &temp,
                                                 std::vector<uint64_t> &hca_info_gen,
                                                 std::vector<uint64_t> &temp_gen,
                                                 unsigned int no_nodes) {
  if (hca_info.size() == no_nodes && hca_info_gen.size() == no_nodes) return;
  hca_info.assign(no_nodes, no_nodes);
  temp.assign(no_nodes, std::numeric_limits<unsigned int>::max());
  hca_info_gen.assign(no_nodes, 0);
  temp_gen.assign(no_nodes, 0);
}

ELM_SPTree::~ELM_SPTree() {
  // Each sp_tree_vec[r] is a vector<shortest_path_tree*> whose elements were
  // individually allocated via `new` in Dijkstra(). The std::vector
  // destructor does NOT call delete on its pointer elements, so we have to
  // do it explicitly here, otherwise every "freed" SPT leaks its ~50000
  // nodes plus their jump_pointers / rmq_vertex / cartesian sub-allocations
  // (~30-100 MB per leaked SPT at n=50000).
  for (auto *ptr : sp_tree_vec) {
    if (!ptr) continue;
    for (auto *node : *ptr) {
      delete node;
    }
    delete ptr;
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
  s += landmark_order.capacity() * sizeof(pair<pair<unsigned int, unsigned int>, double>);

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
        s += node->jump_pointers->capacity() * sizeof(std::pair<shortest_path_tree*, double>);
      }
      if (node->rmq_vertex) {
        s += sizeof(*node->rmq_vertex);
        s += node->rmq_vertex->euler.capacity() * sizeof(int);
        s += node->rmq_vertex->depth.capacity() * sizeof(int);
        s += node->rmq_vertex->first.capacity() * sizeof(int);
        s += node->rmq_vertex->lg.capacity() * sizeof(int);
        s += node->rmq_vertex->st.capacity() * sizeof(std::vector<int>);
        for (auto const &row : node->rmq_vertex->st) s += row.capacity() * sizeof(int);
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
        for (auto const &row : rmq.st) s += row.capacity() * sizeof(int);
      }
    }
  }
  return s;
}

// Naive parent-walk LCA (no RMQ preprocessing). Walks from u up to root,
// marking each visited node id in a temporary set, then walks from v up
// until it hits a marked node. Returns LCA + max edge weight observed on
// either walk before reaching the LCA. Used only by --ablation-no-rmq to
// isolate the §4.1.3 LCA/RMQ amortization benefit. O(depth) per call.
static pair<shortest_path_tree *, double>
find_LCA_naive(vector<shortest_path_tree *> *sp_tree,
               unsigned int index_i, unsigned int index_j) {
  if (!sp_tree) return make_pair((shortest_path_tree *)nullptr, 0.0);
  if (index_i >= sp_tree->size() || index_j >= sp_tree->size()) {
    return make_pair((shortest_path_tree *)nullptr, 0.0);
  }
  shortest_path_tree *node_i = sp_tree->at(index_i);
  shortest_path_tree *node_j = sp_tree->at(index_j);
  if (!node_i || !node_j) return make_pair((shortest_path_tree *)nullptr, 0.0);

  // Mark ancestors of i and track the max edge along the i-walk.
  std::unordered_map<shortest_path_tree *, double> marked;
  shortest_path_tree *cur = node_i;
  double max_i = 0.0;
  while (cur) {
    marked[cur] = max_i;
    if (!cur->parent) break;
    max_i = std::max(max_i, cur->parent_edge_w);
    cur = cur->parent;
  }
  // Walk j up; first node in `marked` is the LCA.
  cur = node_j;
  double max_j = 0.0;
  while (cur) {
    auto it = marked.find(cur);
    if (it != marked.end()) {
      // it->second is max edge from i to this LCA. max_j is max edge from j.
      return make_pair(cur, std::max(it->second, max_j));
    }
    if (!cur->parent) break;
    max_j = std::max(max_j, cur->parent_edge_w);
    cur = cur->parent;
  }
  // Disconnected (shouldn't happen for valid SPTs); fall through to root.
  return make_pair((shortest_path_tree *)nullptr, 0.0);
}

double ELM_SPTree::_lookup_vertex_ex(unsigned int root_node, unsigned int u,
                                     unsigned int v, bool no_rmq, bool skip_lca) {
  auto *tree = sp_tree_vec[root_node];
  shortest_path_tree *lca_node;
  double max_edge_on_path;
  if (skip_lca) {
    // Timing probe: skip the find_LCA call but still execute the surrounding
    // tree-path arithmetic with valid reads (root has path_length 0).
    lca_node = tree->at(root_node);
    max_edge_on_path = 0.0;
  } else {
    pair<shortest_path_tree *, double> lca_info = no_rmq
        ? find_LCA_naive(tree, u, v)
        : find_LCA(tree, u, v);
    lca_node = lca_info.first;
    if (!lca_node) return 0.0;
    max_edge_on_path = lca_info.second;
  }

  double dist_root_u   = tree->at(u)->path_length;
  double dist_root_v   = tree->at(v)->path_length;
  double dist_root_lca = lca_node->path_length;
  double tree_dist_uv = dist_root_u + dist_root_v - 2.0 * dist_root_lca;

  return 2.0 * max_edge_on_path - tree_dist_uv;
}

double ELM_SPTree::_lookup_vertex(unsigned int root_node, unsigned int u, unsigned int v) {
  // Member-flag wrapper; behavior unchanged from the original implementation.
  return _lookup_vertex_ex(root_node, u, v, ablation_no_rmq_);
}

double ELM_SPTree::_lookup_ex(unsigned int u, unsigned int v,
                              pair<unsigned int, unsigned int> landmark_edge,
                              double edge_dist, bool no_spt, bool no_rmq,
                              bool skip_lca) {
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

  // NO-SPT ablation: keep only edge-landmark folding. The SPT/tree-path
  // tightening terms below are the quality contribution being ablated.
  if (!no_spt) {
    best = max(best, _lookup_vertex_ex(a, u, v, no_rmq, skip_lca));
    best = max(best, _lookup_vertex_ex(b, u, v, no_rmq, skip_lca));
  }

  return best;
}

double ELM_SPTree::_lookup(unsigned int u, unsigned int v,
                           pair<unsigned int, unsigned int> landmark_edge,
                           double edge_dist) {
  // Member-flag wrapper; behavior unchanged from the original implementation.
  return _lookup_ex(u, v, landmark_edge, edge_dist,
                    ablation_no_spt_, ablation_no_rmq_);
}

double ELM_SPTree::_lookup_vertex_HCA(unsigned int root_node, unsigned int u, unsigned int v) {
  // Lazy buffer init (#20): only allocate when HCA is first used.
  elm_init_hca_buffers(_hca_information, temp_HCA, no_nodes);
  if (_hca_information.at(root_node) == no_nodes) {
    _hca_information.at(root_node) = find_HCA(sp_tree_vec[u], sp_tree_vec[v], root_node, temp_HCA);
  }

  unsigned int hca_node = _hca_information.at(root_node);
  double lca_longest_edge = max(sp_tree_vec[u]->at(hca_node)->max_edge,
                                sp_tree_vec[v]->at(hca_node)->max_edge);
  double full_path_length = sp_tree_vec[u]->at(hca_node)->path_length
                          + sp_tree_vec[v]->at(hca_node)->path_length;

  return 2 * lca_longest_edge - full_path_length;
}

double ELM_SPTree::_lookup_HCA(unsigned int u, unsigned int v,
                               pair<unsigned int, unsigned int> landmark,
                               double edge_dist) {
  double max_dist = 0.;

  double l1u = sp_tree_vec[u]->at(landmark.first)->path_length;
  double l1v = sp_tree_vec[v]->at(landmark.first)->path_length;
  double l2u = sp_tree_vec[u]->at(landmark.second)->path_length;
  double l2v = sp_tree_vec[v]->at(landmark.second)->path_length;

  max_dist = max(max_dist, max(edge_dist - l1u - l2v, edge_dist - l1v - l2u));

  // True NO-SPT ablation: do not add query-endpoint/landmark-endpoint
  // SPT tree-path tightening terms; keep only edge-landmark folding.
  if (!ablation_no_spt_) {
    max_dist = max(max_dist, _lookup_vertex_HCA(landmark.first, u, v));
    max_dist = max(max_dist, _lookup_vertex_HCA(landmark.second, u, v));
  }

  return max_dist;
}

vector<pair<unsigned int, unsigned int>> ELM_SPTree::all_missing_pairs() const {
  vector<pair<unsigned int, unsigned int>> out;
  if (no_nodes < 2) return out;
  out.reserve((size_t)no_nodes * (size_t)(no_nodes - 1) / 2);
  for (unsigned int u = 0; u < no_nodes; ++u) {
    for (unsigned int v = u + 1; v < no_nodes; ++v) {
      auto uv = norm_edge(u, v);
      if (edges_map->find(uv) == edges_map->end()) out.push_back(uv);
    }
  }
  return out;
}

vector<pair<unsigned int, unsigned int>> ELM_SPTree::sample_missing_pairs(unsigned int count,
                                                                          std::mt19937 &rng) const {
  vector<pair<unsigned int, unsigned int>> out;
  if (no_nodes < 2 || count == 0) return out;

  const uint64_t total_pairs = (uint64_t)no_nodes * (uint64_t)(no_nodes - 1) / 2ULL;
  const uint64_t nonedge_count = total_pairs - (uint64_t)edges_map->size();
  if (nonedge_count == 0) return out;

  if (count >= nonedge_count) return all_missing_pairs();

  out.reserve(count);
  uniform_int_distribution<unsigned int> uni(0, no_nodes - 1);
  unordered_set<uint64_t> used;
  used.reserve(count * 2 + 1);

  while (out.size() < count) {
    unsigned int u = uni(rng), v = uni(rng);
    if (u == v) continue;
    auto uv = norm_edge(u, v);
    if (edges_map->find(uv) != edges_map->end()) continue;
    uint64_t code = pair_code(uv.first, uv.second);
    if (!used.insert(code).second) continue;
    out.push_back(uv);
  }
  return out;
}

vector<pair<pair<unsigned int, unsigned int>, double>> ELM_SPTree::sample_candidate_edges(
    unsigned int count,
    std::mt19937 &rng) const {
  vector<pair<pair<unsigned int, unsigned int>, double>> candidates;
  candidates.reserve(edges_map->size());
  for (auto const &kv : *edges_map) {
    if (landmarks->find(kv.first) == landmarks->end()) candidates.push_back(kv);
  }
  if (count >= candidates.size()) return candidates;
  shuffle(candidates.begin(), candidates.end(), rng);
  candidates.resize(count);
  return candidates;
}

void ELM_SPTree::preprocess_landmark_roots() {
  for (auto &lm : landmark_order) {
    unsigned a = lm.first.first;
    unsigned b = lm.first.second;
    ensure_spt(a);   // self-timed into spt_build_us_
    ensure_spt(b);
    uint64_t _r0 = elm_cpu_us();
    preprocess_spt_queries(sp_tree_vec[a]);
    preprocess_spt_queries(sp_tree_vec[b]);
    rmq_build_us_ += (elm_cpu_us() - _r0);
  }
}

bool ELM_SPTree::save_landmarks(const std::string& path,
                                const std::string& meta) const {
  std::ofstream ofs(path, std::ios::out | std::ios::trunc);
  if (!ofs.is_open()) return false;
  // Metadata header: load_landmarks rejects the file unless this matches the
  // current build config (objective, no-hca, dataset fingerprint, ...), so a
  // cache file selected under a different objective is never silently reused.
  ofs << "#META\t" << meta << "\n";
  ofs << landmark_order.size() << "\n";
  ofs.setf(std::ios::scientific);
  ofs.precision(17);
  for (const auto& lm : landmark_order) {
    ofs << lm.first.first << " " << lm.first.second << " " << lm.second << "\n";
  }
  return ofs.good();
}

bool ELM_SPTree::load_landmarks(const std::string& path,
                                const std::string& meta) {
  std::ifstream ifs(path, std::ios::in);
  if (!ifs.is_open()) return false;
  std::string header;
  if (!std::getline(ifs, header)) return false;
  const std::string expected = std::string("#META\t") + meta;
  if (header != expected) return false;  // config mismatch -> caller rebuilds
  size_t count = 0;
  if (!(ifs >> count)) return false;
  // Reset and repopulate both the map and the ordered list exactly as the
  // get_*_landmarks() routines do.
  landmarks->clear();
  landmark_order.clear();
  landmark_order.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    unsigned int a, b;
    double w;
    if (!(ifs >> a >> b >> w)) return false;
    std::pair<unsigned int, unsigned int> e(a, b);
    landmarks->insert({e, w});
    landmark_order.push_back({e, w});
  }
  return true;
}

double ELM_SPTree::lookup(unsigned int u, unsigned int v) {
  auto key = norm_edge(u, v);
  auto it = edges_map->find(key);
  if (it != edges_map->end()) return it->second;

  // Cross-query cache: a pure read-through cache keyed by the encoded
  // (u, v) pair. Safe to consult anytime the landmark set is unchanged.
  uint64_t code = pair_code(key.first, key.second);
  if (lb_cache_enabled_) {
    auto cit = lb_cache_.find(code);
    if (cit != lb_cache_.end()) return cit->second;
  }

  double best = 0.0;
  for (auto &lm : landmark_order) {
    best = max(best, _lookup(u, v, lm.first, lm.second));
  }
  if (lb_cache_enabled_) lb_cache_.emplace(code, best);
  return best;
}

double ELM_SPTree::lookup_ablate(unsigned int u, unsigned int v,
                                 bool no_spt, bool no_rmq, bool skip_lca) {
  // Same structure as lookup(), but uses the explicit-flag scoring path and
  // intentionally bypasses the LB cache so repeated ON/OFF passes over the
  // same query set do not contaminate each other's timings.
  auto key = norm_edge(u, v);
  auto it = edges_map->find(key);
  if (it != edges_map->end()) return it->second;

  double best = 0.0;
  for (auto &lm : landmark_order) {
    best = max(best, _lookup_ex(u, v, lm.first, lm.second, no_spt, no_rmq, skip_lca));
  }
  return best;
}

double ELM_SPTree::lookup_skeleton(unsigned int u, unsigned int v) {
  // Loop the landmarks doing negligible (non-elidable) work to measure the
  // per-query iteration/dispatch overhead. Touches edge_dist + endpoints so the
  // optimizer cannot remove the loop.
  auto key = norm_edge(u, v);
  auto it = edges_map->find(key);
  if (it != edges_map->end()) return it->second;
  double acc = 0.0;
  for (auto &lm : landmark_order) {
    acc += lm.second + (double)lm.first.first + (double)lm.first.second;
  }
  // Fold in u,v so the signature is used and acc depends on the query.
  return acc * 0.0 + (double)(u ^ v) * 0.0;
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

  // Optimization (#12): the bound from each landmark uses path_length(u) and
  // path_length(v) in the SPT rooted at each endpoint of the landmark. Many
  // landmarks share endpoints across iterations of this loop; even when they
  // don't, the per-tree dereferences incur pointer-chasing. Hoist what we
  // can: cache (treeRoot -> (path_u, path_v)) for the duration of this call.
  std::unordered_map<unsigned int, std::pair<double,double>> spt_uv_cache;
  spt_uv_cache.reserve(no_landmarks * 2 + 4);
  auto get_uv = [&](unsigned int r) -> std::pair<double,double> {
    auto cit = spt_uv_cache.find(r);
    if (cit != spt_uv_cache.end()) return cit->second;
    auto *tree = sp_tree_vec[r];
    std::pair<double,double> pr = {tree->at(u)->path_length, tree->at(v)->path_length};
    spt_uv_cache.emplace(r, pr);
    return pr;
  };

  double best = 0.0;
  unsigned int count = 0;
  for (auto &lm : landmark_order) {
    unsigned int a = lm.first.first, b = lm.first.second;
    double w = lm.second;
    auto au_av = get_uv(a);
    auto bu_bv = get_uv(b);
    double cand = std::max(0.0, std::max(w - au_av.first  - bu_bv.second,
                                          w - au_av.second - bu_bv.first));
    if (!ablation_no_spt_) {
      cand = std::max(cand, _lookup_vertex(a, u, v));
      cand = std::max(cand, _lookup_vertex(b, u, v));
    }
    if (cand > best) best = cand;
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
    // Free every node in the SPT before freeing the holding vector. The
    // vector destructor doesn't recursively delete its element pointers.
    for (auto *node : *sp_tree_vec[r]) {
      delete node;
    }
    delete sp_tree_vec[r];
    sp_tree_vec[r] = nullptr;
  }
}

void ELM_SPTree::get_exact_landmarks() {
  landmarks->clear();
  lb_cache_.clear();
  landmark_order.clear();

  vector<pair<unsigned int, unsigned int>> queries = all_missing_pairs();
  if (queries.empty()) return;

  for (unsigned int i = 0; i < no_nodes; ++i) ensure_spt(i);

  vector<double> current_best(queries.size(), 0.0);

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
      for (size_t qi = 0; qi < queries.size(); ++qi) {
        unsigned int u = queries[qi].first;
        unsigned int v = queries[qi].second;
        double cand_lb = _lookup(u, v, e, e_len);
        if (cand_lb > current_best[qi]) gain += (cand_lb - current_best[qi]);
      }

      if (gain > best_gain) {
        best_gain = gain;
        best_edge = e;
        best_len = e_len;
      }
    }

    if (best_gain < 0.0) break;
    landmarks->insert({best_edge, best_len});
    landmark_order.push_back({best_edge, best_len});

    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first;
      unsigned int v = queries[qi].second;
      double val = _lookup(u, v, best_edge, best_len);
      if (val > current_best[qi]) current_best[qi] = val;
    }
  }

  clean_unwanted_shortest_paths();
}

void ELM_SPTree::get_stochastic_landmarks(double eps,
                                          unsigned int query_samples,
                                          unsigned int seed,
                                          unsigned int cand_cap) {
  // Note (#24): the candidate-side sample size below is the
  // Mirzasoleiman-Badanidiyuru-Vondrak (NeurIPS 2015) "stochastic greedy"
  // formula, giving expected approximation 1 - 1/e - eps for monotone
  // submodular maximization. The paper's Theorem 4.6 derives a different,
  // query-side sample size via Hoeffding's inequality (s >= c^2 ln(2/delta)
  // / (2 eps^2) where c is the graph diameter) for sum-estimation
  // confidence. Both are valid choices; the code uses Mirzasoleiman's for
  // the candidate set (this function) and the user-supplied query_samples
  // for the query set (the paper's Hoeffding bound is the recommended way
  // to pick query_samples, computed externally).
  //
  // `cand_cap` is an upper bound on the per-iteration candidate sample
  // size. Each sampled candidate triggers two SPT+RMQ builds via
  // `ensure_spt`, so on a 50000-node graph the unbounded sample size of
  // n/k * log(1/eps) (~160k at k=8, eps=0.1) leads to runtimes of days.
  // Capping trades a small loss of the (1 - 1/e - eps) guarantee for a
  // bounded, predictable wall clock. Default UINT_MAX preserves the
  // original unbounded behavior; a reasonable cap for large graphs is
  // ~1000-5000.
  landmarks->clear();
  lb_cache_.clear();
  landmark_order.clear();

  if (eps <= 0.0) eps = 1e-3;
  if (eps >= 1.0) eps = 0.999;

  std::mt19937 rng(seed);
  vector<pair<unsigned int, unsigned int>> queries = sample_missing_pairs(query_samples, rng);
  if (queries.empty()) return;

  // IMPORTANT: get_stochastic_landmarks uses _lookup (NOT _lookup_HCA).
  // _lookup reads sp_tree_vec[landmark_endpoint] only -- it never reads
  // sp_tree_vec[u] or sp_tree_vec[v] where (u,v) are the query pair.
  // (Contrast with get_hca_landmarks, where _lookup_HCA reads the (u,v)
  // trees.) So we do NOT need to build query-endpoint trees here. Building
  // them anyway would burn ~95 MB per tree of pure dead memory at n=50000.

  vector<double> current_best(queries.size(), 0.0);

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    const unsigned int remaining_budget = no_landmarks - iter;
    const unsigned int remaining_candidates =
        (unsigned int)edges_map->size() - (unsigned int)landmarks->size();
    if (remaining_candidates == 0) break;

    unsigned int sample_size = (unsigned int)std::ceil(
        (double)remaining_candidates / (double)std::max(1u, remaining_budget) * std::log(1.0 / eps));
    sample_size = std::max(1u, std::min(sample_size, remaining_candidates));
    if (cand_cap > 0 && sample_size > cand_cap) sample_size = cand_cap;

    auto candidates = sample_candidate_edges(sample_size, rng);

    // Pre-build candidate SPTs in parallel. Each candidate has two endpoints;
    // some may already be cached from previous iterations (landmark endpoints
    // are kept; non-landmark candidates from previous iters are not). Build
    // the missing ones with one parallel-for. This is the dominant time cost
    // of the stochastic path; parallelizing it pays back ~10x on 16 cores.
    {
      std::unordered_set<unsigned int> missing;
      missing.reserve(candidates.size() * 2 + 8);
      for (auto &cand : candidates) {
        if (sp_tree_vec[cand.first.first]  == nullptr) missing.insert(cand.first.first);
        if (sp_tree_vec[cand.first.second] == nullptr) missing.insert(cand.first.second);
      }
      std::vector<unsigned int> to_build(missing.begin(), missing.end());
      #pragma omp parallel for schedule(dynamic, 1) if (to_build.size() > 4)
      for (long long i = 0; i < (long long)to_build.size(); ++i) {
        unsigned int r = to_build[i];
        sp_tree_vec[r] = Dijkstra(adj_list, r, true);
      }
    }

    // Marginal-gain evaluation: independent across candidates. Each thread
    // accumulates its own best; one critical section at the end reduces.
    double best_gain = -1.0;
    pair<unsigned int, unsigned int> best_edge = {0, 0};
    double best_len = 0.0;
    {
      #pragma omp parallel
      {
        double local_best_gain = -1.0;
        pair<unsigned int, unsigned int> local_best_edge = {0, 0};
        double local_best_len = 0.0;

        #pragma omp for schedule(dynamic, 16) nowait
        for (long long ci = 0; ci < (long long)candidates.size(); ++ci) {
          const auto &e = candidates[ci].first;
          double e_len = candidates[ci].second;

          double gain = 0.0;
          for (size_t qi = 0; qi < queries.size(); ++qi) {
            unsigned int u = queries[qi].first;
            unsigned int v = queries[qi].second;
            double cand_lb = _lookup(u, v, e, e_len);
            if (cand_lb > current_best[qi]) gain += (cand_lb - current_best[qi]);
          }

          if (gain > local_best_gain) {
            local_best_gain = gain;
            local_best_edge = e;
            local_best_len = e_len;
          }
        }

        #pragma omp critical
        {
          if (local_best_gain > best_gain) {
            best_gain = local_best_gain;
            best_edge = local_best_edge;
            best_len = local_best_len;
          }
        }
      }
    }

    if (best_gain < 0.0) break;

    landmarks->insert({best_edge, best_len});
    landmark_order.push_back({best_edge, best_len});

    ensure_spt(best_edge.first);
    ensure_spt(best_edge.second);
    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first;
      unsigned int v = queries[qi].second;
      double val = _lookup(u, v, best_edge, best_len);
      if (val > current_best[qi]) current_best[qi] = val;
    }

    // INTRA-ITER CLEANUP: free SPTs of non-winning candidates. Keep ONLY
    // selected-landmark endpoint trees (needed for future iters' _lookup
    // updates of current_best, and for post-build query-time lookup()).
    // Query (u,v) endpoint trees are NOT needed by _lookup, so they aren't
    // in the keep set. At n=50000 each SPT is ~95 MB, so this matters.
    //
    // Note the explicit node-deletion loop: each sp_tree_vec[r] holds
    // pointers to individually-`new`-allocated shortest_path_tree nodes.
    // The std::vector destructor does NOT cascade to those pointers; if we
    // only `delete sp_tree_vec[r]`, every "freed" SPT leaks all ~50000 of
    // its nodes (plus their jump_pointers, RMQ tables, Cartesian preproc),
    // and intra-iter cleanup buys us nothing.
    {
      std::unordered_set<unsigned int> keep;
      keep.reserve(landmarks->size() * 2 + 8);
      for (auto &kv : *landmarks) {
        keep.insert(kv.first.first);
        keep.insert(kv.first.second);
      }
      for (unsigned int r = 0; r < no_nodes; ++r) {
        if (sp_tree_vec[r] == nullptr) continue;
        if (keep.find(r) != keep.end()) continue;
        for (auto *node : *sp_tree_vec[r]) {
          delete node;
        }
        delete sp_tree_vec[r];
        sp_tree_vec[r] = nullptr;
      }
    }
  }

  clean_unwanted_shortest_paths();
}

void ELM_SPTree::get_sampling_landmarks() {
  // Backward-compatible entry point: use stochastic greedy with the existing
  // no_samples parameter as the fixed query-set size.
  get_stochastic_landmarks(0.1, no_samples, 73);
}

void ELM_SPTree::_reset_HCA() {
  // Honors the lazy-init contract from #20: if buffers haven't been
  // allocated yet, leave them empty -- callers will init lazily on use.
  if (_hca_information.size() == no_nodes) {
    _hca_information.assign(_hca_information.size(), no_nodes);
  }
  if (temp_HCA.size() == no_nodes) {
    temp_HCA.assign(temp_HCA.size(), numeric_limits<unsigned int>::max());
  }
}

// ============================================================================
// New: paper-faithful §4.2 SGEL using HCA-based candidate scoring.
// ============================================================================
// Per iteration:
//   1. Sample `query_samples` missing edges (training queries).
//   2. For each sampled query (u, v), build SPTs at u and v if not already
//      built. These are the only Dijkstras paid per iteration -- no per-
//      candidate Dijkstra.
//   3. Score every known edge (x, y) against every sampled query using
//      _lookup_HCA, which uses memoized upward walks in tau_u and tau_v.
//   4. Pick the edge with maximum total marginal gain; commit.
//
// This matches the paper's Algorithm 2 + Section 4.2 semantics. The
// current_best vector tracks the running max LB per training query so that
// marginal gain is always computed against the committed landmark set.
void ELM_SPTree::get_hca_landmarks(unsigned int query_samples,
                                   unsigned int seed,
                                   bool relative_objective) {
  landmarks->clear();
  lb_cache_.clear();
  landmark_order.clear();

  // Activate HCA path: allocate buffers now so we don't pay for the lazy
  // check on every iteration.
  elm_init_hca_buffers(_hca_information, temp_HCA, no_nodes);

  std::mt19937 rng(seed);
  vector<pair<unsigned int, unsigned int>> queries =
      sample_missing_pairs(query_samples, rng);
  if (queries.empty()) return;

  // For HCA we need jump_pointers in the query-side SPTs. We force-build
  // SPTs for every query endpoint up-front since they get hit O(k * m) times
  // each round.
  unordered_set<unsigned int> query_endpoints;
  query_endpoints.reserve(queries.size() * 2 + 8);
  for (auto &q : queries) {
    query_endpoints.insert(q.first);
    query_endpoints.insert(q.second);
  }
  for (unsigned int r : query_endpoints) ensure_spt(r);

  // Objective weights for the Zoom-discussed ARLB-aligned selection.
  // Raw/L1 objective: weight = 1 for every sampled query.
  // Relative objective: weight = 1/OPT(q), where OPT(q) is exact/SPLUB LB.
  // Queries with OPT(q)=0 are excluded from the training objective.
  vector<double> query_weights(queries.size(), 1.0);
  unsigned long long zero_opt_training_queries = 0;
  if (relative_objective) {
    const double EPS = 1e-12;
    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first;
      unsigned int v = queries[qi].second;
      auto *tu = sp_tree_vec[u];
      auto *tv = sp_tree_vec[v];
      double opt = 0.0;
      for (auto &kv : *edges_map) {
        unsigned int x = kv.first.first;
        unsigned int y = kv.first.second;
        double w = kv.second;
        if (!tu || !tv || !tu->at(x) || !tu->at(y) || !tv->at(x) || !tv->at(y)) continue;
        double c1 = w - tu->at(x)->path_length - tv->at(y)->path_length;
        double c2 = w - tu->at(y)->path_length - tv->at(x)->path_length;
        opt = std::max(opt, std::max(c1, c2));
      }
      if (opt <= EPS) {
        query_weights[qi] = 0.0;
        ++zero_opt_training_queries;
      } else {
        query_weights[qi] = 1.0 / opt;
      }
    }
  }
  cout << "[SGEL objective] mode=" << (relative_objective ? "relative_ARLB" : "raw_L1")
       << " training_queries=" << queries.size()
       << " zero_opt_training_queries=" << zero_opt_training_queries
       << " valid_training_queries=" << (queries.size() - zero_opt_training_queries)
       << endl;

  vector<double> current_best(queries.size(), 0.0);

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    double best_gain = -1.0;
    pair<unsigned int, unsigned int> best_edge = {0, 0};
    double best_len = 0.0;

    // For each query, the HCA cache is per-(u, v); reset before scoring all
    // candidates against THIS query.
    for (auto &cand : *edges_map) {
      if (landmarks->find(cand.first) != landmarks->end()) continue;
      const auto &e = cand.first;
      double e_len = cand.second;

      double gain = 0.0;
      for (size_t qi = 0; qi < queries.size(); ++qi) {
        unsigned int u = queries[qi].first;
        unsigned int v = queries[qi].second;
        // Reset HCA cache because it depends on the specific (u, v) pair
        // (the cache is keyed by root_node only). This is the correctness
        // fix discussed in the audit.
        _reset_HCA();
        // Ablation R3.D7: if --ablation-no-hca, use the §4.1 edge-folding +
        // landmark-tree bound (_lookup), which omits the §4.2.1 HCA term
        // that involves the (u,v) SPTs. Isolates the contribution of the
        // HCA approach.
        double cand_lb = ablation_no_hca_
                       ? _lookup(u, v, e, e_len)
                       : _lookup_HCA(u, v, e, e_len);
        if (cand_lb > current_best[qi]) gain += query_weights[qi] * (cand_lb - current_best[qi]);
      }

      if (gain > best_gain) {
        best_gain = gain;
        best_edge = e;
        best_len = e_len;
      }
    }

    if (best_gain <= 0.0) break;
    landmarks->insert({best_edge, best_len});
    landmark_order.push_back({best_edge, best_len});

    // For lookup() (LCA-based) we still want full SPT preprocessing at each
    // landmark endpoint. Do that now lazily.
    ensure_spt(best_edge.first);
    ensure_spt(best_edge.second);

    // Update current_best for this committed landmark.
    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first;
      unsigned int v = queries[qi].second;
      _reset_HCA();
      double val = ablation_no_hca_
                 ? _lookup(u, v, best_edge, best_len)
                 : _lookup_HCA(u, v, best_edge, best_len);
      if (val > current_best[qi]) current_best[qi] = val;
    }
  }

  clean_unwanted_shortest_paths();
}

// ============================================================================
// New: lazy/accelerated greedy GEL (Minoux 1978).
// ============================================================================
// Marginal gains are non-increasing under submodular monotone maximization,
// so cached gains from earlier iterations are valid upper bounds. We
// maintain a max-heap of (cached_gain, last_iter, edge); each round we pop
// and refresh until the top is up-to-date, then commit.
namespace {
struct ELMLazyEntry {
  double gain;
  unsigned int last_iter;
  pair<unsigned int, unsigned int> edge;
  double edge_len;
  bool operator<(const ELMLazyEntry &o) const {
    if (gain != o.gain) return gain < o.gain;
    if (edge.first != o.edge.first) return edge.first > o.edge.first;
    return edge.second > o.edge.second;
  }
};
} // namespace

void ELM_SPTree::get_exact_landmarks_lazy() {
  landmarks->clear();
  lb_cache_.clear();
  landmark_order.clear();

  vector<pair<unsigned int, unsigned int>> queries = all_missing_pairs();
  if (queries.empty()) return;

  // Same upfront SPT build as get_exact_landmarks (we'll reuse them).
  for (unsigned int i = 0; i < no_nodes; ++i) ensure_spt(i);

  vector<double> current_best(queries.size(), 0.0);

  auto eval_gain = [&](const pair<unsigned int, unsigned int> &e, double w) -> double {
    double g = 0.0;
    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first, v = queries[qi].second;
      double cand_lb = _lookup(u, v, e, w);
      if (cand_lb > current_best[qi]) g += (cand_lb - current_best[qi]);
    }
    return g;
  };

  std::priority_queue<ELMLazyEntry> heap;
  for (auto &kv : *edges_map) {
    double g = eval_gain(kv.first, kv.second);
    heap.push({g, 0u, kv.first, kv.second});
  }

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    if (heap.empty()) break;
    ELMLazyEntry top = heap.top();
    while (top.last_iter < iter) {
      heap.pop();
      double fresh = eval_gain(top.edge, top.edge_len);
      heap.push({fresh, iter, top.edge, top.edge_len});
      top = heap.top();
    }
    if (top.gain <= 0.0) {
      cout << "[get_exact_landmarks_lazy] No positive gain at iter "
           << iter << "; stopping." << endl;
      break;
    }
    landmarks->insert({top.edge, top.edge_len});
    landmark_order.push_back({top.edge, top.edge_len});
    heap.pop();

    // Refresh current_best for the committed landmark.
    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first, v = queries[qi].second;
      double val = _lookup(u, v, top.edge, top.edge_len);
      if (val > current_best[qi]) current_best[qi] = val;
    }
  }

  clean_unwanted_shortest_paths();
}

// ============================================================================
// New: parallel candidate evaluation for GEL.
// ============================================================================
// Same algorithm as get_exact_landmarks() but the per-iteration scan over
// candidate edges is parallelized via OpenMP if available. Each thread
// keeps a private best-gain triple; main thread reduces.
//
// Compile with -fopenmp to enable; otherwise this is a faithful sequential
// implementation. Output set is byte-identical to get_exact_landmarks() up
// to tie-breaking among equal-gain edges.
void ELM_SPTree::get_exact_landmarks_parallel() {
  landmarks->clear();
  lb_cache_.clear();
  landmark_order.clear();

  vector<pair<unsigned int, unsigned int>> queries = all_missing_pairs();
  if (queries.empty()) return;

  for (unsigned int i = 0; i < no_nodes; ++i) ensure_spt(i);

  vector<double> current_best(queries.size(), 0.0);

  // Snapshot candidates into a vector so OpenMP can index over them.
  vector<pair<pair<unsigned int, unsigned int>, double>> candidates;
  candidates.reserve(edges_map->size());
  for (auto &kv : *edges_map) candidates.push_back(kv);

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    double best_gain = -1.0;
    pair<unsigned int, unsigned int> best_edge = {0, 0};
    double best_len = 0.0;

#ifdef _OPENMP
    #pragma omp parallel
    {
      double local_best_gain = -1.0;
      pair<unsigned int, unsigned int> local_best_edge = {0, 0};
      double local_best_len = 0.0;

      #pragma omp for nowait schedule(dynamic, 32)
      for (long long ci = 0; ci < (long long)candidates.size(); ++ci) {
        const auto &e = candidates[ci].first;
        if (landmarks->find(e) != landmarks->end()) continue;
        double e_len = candidates[ci].second;
        double gain = 0.0;
        for (size_t qi = 0; qi < queries.size(); ++qi) {
          unsigned int u = queries[qi].first, v = queries[qi].second;
          double cand_lb = _lookup(u, v, e, e_len);
          if (cand_lb > current_best[qi]) gain += (cand_lb - current_best[qi]);
        }
        if (gain > local_best_gain) {
          local_best_gain = gain;
          local_best_edge = e;
          local_best_len = e_len;
        }
      }

      #pragma omp critical
      {
        if (local_best_gain > best_gain) {
          best_gain = local_best_gain;
          best_edge = local_best_edge;
          best_len = local_best_len;
        }
      }
    }
#else
    // Sequential fallback (identical to get_exact_landmarks's inner scan).
    for (const auto &cand : candidates) {
      const auto &e = cand.first;
      if (landmarks->find(e) != landmarks->end()) continue;
      double e_len = cand.second;
      double gain = 0.0;
      for (size_t qi = 0; qi < queries.size(); ++qi) {
        unsigned int u = queries[qi].first, v = queries[qi].second;
        double cand_lb = _lookup(u, v, e, e_len);
        if (cand_lb > current_best[qi]) gain += (cand_lb - current_best[qi]);
      }
      if (gain > best_gain) {
        best_gain = gain;
        best_edge = e;
        best_len = e_len;
      }
    }
#endif

    if (best_gain <= 0.0) break;
    landmarks->insert({best_edge, best_len});
    landmark_order.push_back({best_edge, best_len});

    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first, v = queries[qi].second;
      double val = _lookup(u, v, best_edge, best_len);
      if (val > current_best[qi]) current_best[qi] = val;
    }
  }

  clean_unwanted_shortest_paths();
}


// ============================================================================
// New v26: exact greedy GEL with ARLB-aligned relative objective.
// ============================================================================
// This is the exact-greedy counterpart of ELM Sampling RelObj / SGEL RelObj.
// It uses the same full missing-edge query set as get_exact_landmarks_parallel()
// and the same full candidate set of known edges, but weights every marginal
// gain by 1/OPT(q), where OPT(q) is the exact/SPLUB lower bound for that
// missing query q. Queries with OPT(q)=0 are excluded from the objective.
void ELM_SPTree::get_exact_landmarks_relative_parallel() {
  landmarks->clear();
  lb_cache_.clear();
  landmark_order.clear();

  vector<pair<unsigned int, unsigned int>> queries = all_missing_pairs();
  if (queries.empty()) return;

  // Exact greedy path materializes all SPTs once; this is why the method is
  // used only in EXACT/small-graph mode by fixed_main.cpp.
  for (unsigned int i = 0; i < no_nodes; ++i) ensure_spt(i);

  const double EPS = 1e-12;
  vector<double> query_weights(queries.size(), 1.0);
  unsigned long long zero_opt_training_queries = 0;

  for (size_t qi = 0; qi < queries.size(); ++qi) {
    unsigned int u = queries[qi].first;
    unsigned int v = queries[qi].second;
    auto *tu = sp_tree_vec[u];
    auto *tv = sp_tree_vec[v];
    double opt = 0.0;
    for (auto &kv : *edges_map) {
      unsigned int x = kv.first.first;
      unsigned int y = kv.first.second;
      double w = kv.second;
      if (!tu || !tv || !tu->at(x) || !tu->at(y) || !tv->at(x) || !tv->at(y)) continue;
      double c1 = w - tu->at(x)->path_length - tv->at(y)->path_length;
      double c2 = w - tu->at(y)->path_length - tv->at(x)->path_length;
      opt = std::max(opt, std::max(c1, c2));
    }
    if (opt <= EPS) {
      query_weights[qi] = 0.0;
      ++zero_opt_training_queries;
    } else {
      query_weights[qi] = 1.0 / opt;
    }
  }

  cout << "[Exact objective] mode=relative_ARLB"
       << " training_queries=" << queries.size()
       << " zero_opt_training_queries=" << zero_opt_training_queries
       << " valid_training_queries=" << (queries.size() - zero_opt_training_queries)
       << endl;

  vector<double> current_best(queries.size(), 0.0);

  vector<pair<pair<unsigned int, unsigned int>, double>> candidates;
  candidates.reserve(edges_map->size());
  for (auto &kv : *edges_map) candidates.push_back(kv);

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    double best_gain = -1.0;
    pair<unsigned int, unsigned int> best_edge = {0, 0};
    double best_len = 0.0;

#ifdef _OPENMP
    #pragma omp parallel
    {
      double local_best_gain = -1.0;
      pair<unsigned int, unsigned int> local_best_edge = {0, 0};
      double local_best_len = 0.0;

      #pragma omp for nowait schedule(dynamic, 32)
      for (long long ci = 0; ci < (long long)candidates.size(); ++ci) {
        const auto &e = candidates[ci].first;
        if (landmarks->find(e) != landmarks->end()) continue;
        double e_len = candidates[ci].second;
        double gain = 0.0;
        for (size_t qi = 0; qi < queries.size(); ++qi) {
          double wt = query_weights[qi];
          if (wt <= 0.0) continue;
          unsigned int u = queries[qi].first, v = queries[qi].second;
          double cand_lb = _lookup(u, v, e, e_len);
          if (cand_lb > current_best[qi]) gain += wt * (cand_lb - current_best[qi]);
        }
        if (gain > local_best_gain) {
          local_best_gain = gain;
          local_best_edge = e;
          local_best_len = e_len;
        }
      }

      #pragma omp critical
      {
        if (local_best_gain > best_gain) {
          best_gain = local_best_gain;
          best_edge = local_best_edge;
          best_len = local_best_len;
        }
      }
    }
#else
    for (const auto &cand : candidates) {
      const auto &e = cand.first;
      if (landmarks->find(e) != landmarks->end()) continue;
      double e_len = cand.second;
      double gain = 0.0;
      for (size_t qi = 0; qi < queries.size(); ++qi) {
        double wt = query_weights[qi];
        if (wt <= 0.0) continue;
        unsigned int u = queries[qi].first, v = queries[qi].second;
        double cand_lb = _lookup(u, v, e, e_len);
        if (cand_lb > current_best[qi]) gain += wt * (cand_lb - current_best[qi]);
      }
      if (gain > best_gain) {
        best_gain = gain;
        best_edge = e;
        best_len = e_len;
      }
    }
#endif

    if (best_gain <= 0.0) break;
    landmarks->insert({best_edge, best_len});
    landmark_order.push_back({best_edge, best_len});

    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first, v = queries[qi].second;
      double val = _lookup(u, v, best_edge, best_len);
      if (val > current_best[qi]) current_best[qi] = val;
    }
  }

  clean_unwanted_shortest_paths();
}

// ============================================================================
// New: batched per-source LB fill (kNN-friendly).
// ============================================================================
void ELM_SPTree::lookup_from_source(unsigned int u, vector<double> &LB) {
  if (LB.size() != (size_t)no_nodes) LB.assign(no_nodes, 0.0);
  else std::fill(LB.begin(), LB.end(), 0.0);

  // Walk landmarks once; for each, inner-loop over destinations.
  // The lookup formula needs both query endpoints' SPT path lengths through
  // each landmark endpoint's SPT, so this is just a tighter-loop version of
  // calling lookup(u, v) per v -- but with all the per-landmark dispatch
  // hoisted out of the inner loop.
  for (auto &lm : landmark_order) {
    unsigned int a = lm.first.first, b = lm.first.second;
    double w = lm.second;
    auto *treeA = sp_tree_vec[a];
    auto *treeB = sp_tree_vec[b];
    if (!treeA || !treeB) continue;

    double au = treeA->at(u)->path_length;
    double bu = treeB->at(u)->path_length;

    // (u, v) might equal a known edge -- defer that case to lookup() since
    // edges_map fast-path is rare in kNN scans anyway.
    for (unsigned int v = 0; v < no_nodes; ++v) {
      if (v == u) continue;
      double av = treeA->at(v)->path_length;
      double bv = treeB->at(v)->path_length;
      double cand = max(0.0, max(w - au - bv, w - av - bu));
      if (!ablation_no_spt_) {
        cand = max(cand, _lookup_vertex(a, u, v));
        cand = max(cand, _lookup_vertex(b, u, v));
      }
      if (cand > LB[v]) LB[v] = cand;
    }
  }
}


// ============================================================================
// OpenMP-parallel HCA SGEL.
// ============================================================================
// Same algorithm as get_hca_landmarks but with two performance changes:
//   1. The candidate loop is parallelized; each thread maintains a private
//      best-gain entry; the parallel region reduces them.
//   2. Each thread carries its own HCA memoization buffers, gated by a
//      thread-local generation counter. Resetting the HCA cache between
//      queries becomes O(1) (a generation bump) instead of O(no_nodes).
//
// Output set is identical to get_hca_landmarks() up to tie-breaking among
// equal-gain edges (which is also non-deterministic in the serial version
// when ties occur because std::map traversal order is fixed). The (1 - 1/e)
// approximation guarantee is preserved.

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// Inlined find_HCA that consults a generation array. A memoization cell is
// considered valid iff `gen_arr[i] == gen`. Writes update both `vals[i]` and
// `gen_arr[i]`. Identical answer to the global find_HCA; the difference is
// that resets are O(1).
inline unsigned int find_HCA_gen(
    vector<shortest_path_tree*> *sp_tree_u,
    vector<shortest_path_tree*> *sp_tree_v,
    unsigned int index,
    std::vector<unsigned int> &vals,
    std::vector<uint64_t> &gen_arr,
    uint64_t gen)
{
    const size_t N = vals.size();
    if ((size_t)index >= N) return index;
    if (gen_arr[index] == gen) return vals[index];

    unsigned int index_i = index;
    unsigned int index_j = index;

    // Small inline buffer; SPT depth bounded by ~log(n) for sparse graphs.
    unsigned int path_buf[128];
    int path_len = 0;

    while (true) {
        if (index_i != index_j) break;
        if (path_len < 128) path_buf[path_len++] = index_i;

        shortest_path_tree* ni = sp_tree_u->at(index_i);
        shortest_path_tree* nj = sp_tree_v->at(index_j);
        if (!ni || !nj) break;
        if (!ni->jump_pointers || ni->jump_pointers->empty()) break;
        if (!nj->jump_pointers || nj->jump_pointers->empty()) break;

        unsigned int next_i = ni->jump_pointers->at(0).first->id;
        unsigned int next_j = nj->jump_pointers->at(0).first->id;
        if (next_i == index_i || next_j == index_j) break;

        index_i = next_i;
        index_j = next_j;
    }

    unsigned int answer = (path_len == 0) ? index : path_buf[path_len - 1];
    for (int k = 0; k < path_len; ++k) {
        unsigned int v = path_buf[k];
        if ((size_t)v < N) { vals[v] = answer; gen_arr[v] = gen; }
    }
    return answer;
}

// Inlined _lookup_HCA using per-thread buffers + generation counter.
// `gen` is a thread-local generation that must be bumped (++) by the caller
// between distinct (u, v) queries.
inline double lookup_HCA_local(
    vector<vector<shortest_path_tree*>*> &sp_tree_vec,
    unsigned int u, unsigned int v,
    pair<unsigned int, unsigned int> landmark, double edge_dist,
    std::vector<unsigned int> &hca_vals, std::vector<uint64_t> &hca_gen,
    std::vector<unsigned int> &tmp_vals, std::vector<uint64_t> &tmp_gen,
    uint64_t gen, unsigned int no_nodes, bool no_spt = false)
{
    double l1u = sp_tree_vec[u]->at(landmark.first)->path_length;
    double l1v = sp_tree_vec[v]->at(landmark.first)->path_length;
    double l2u = sp_tree_vec[u]->at(landmark.second)->path_length;
    double l2v = sp_tree_vec[v]->at(landmark.second)->path_length;

    double max_dist = max(edge_dist - l1u - l2v, edge_dist - l1v - l2u);
    if (max_dist < 0.0) max_dist = 0.0;

    auto vertex_bound = [&](unsigned int root_node) -> double {
        if ((size_t)root_node >= hca_vals.size()) return 0.0;
        unsigned int hca_node;
        if (hca_gen[root_node] == gen) {
            hca_node = hca_vals[root_node];
        } else {
            hca_node = find_HCA_gen(sp_tree_vec[u], sp_tree_vec[v], root_node,
                                    tmp_vals, tmp_gen, gen);
            hca_vals[root_node] = hca_node;
            hca_gen[root_node] = gen;
        }
        double e_max = max(sp_tree_vec[u]->at(hca_node)->max_edge,
                           sp_tree_vec[v]->at(hca_node)->max_edge);
        double plen = sp_tree_vec[u]->at(hca_node)->path_length
                    + sp_tree_vec[v]->at(hca_node)->path_length;
        return 2.0 * e_max - plen;
    };

    if (!no_spt) {
        max_dist = max(max_dist, vertex_bound(landmark.first));
        max_dist = max(max_dist, vertex_bound(landmark.second));
    }
    return max_dist;
}

} // anonymous namespace

void ELM_SPTree::get_hca_landmarks_parallel(unsigned int query_samples,
                                            unsigned int seed,
                                            bool relative_objective) {
    landmarks->clear();
    lb_cache_.clear();
    landmark_order.clear();

    elm_init_hca_buffers(_hca_information, temp_HCA, no_nodes);

    std::mt19937 rng(seed);
    vector<pair<unsigned int, unsigned int>> queries =
        sample_missing_pairs(query_samples, rng);
    if (queries.empty()) return;

    // Build SPTs for every query endpoint sequentially (these are read by
    // every thread; safer to materialize them before going parallel).
    unordered_set<unsigned int> query_endpoints;
    query_endpoints.reserve(queries.size() * 2 + 8);
    for (auto &q : queries) {
        query_endpoints.insert(q.first);
        query_endpoints.insert(q.second);
    }
    for (unsigned int r : query_endpoints) ensure_spt(r);

    // Objective weights for raw-L1 vs ARLB-aligned relative-gain SGEL.
    vector<double> query_weights(queries.size(), 1.0);
    unsigned long long zero_opt_training_queries = 0;
    if (relative_objective) {
        const double EPS = 1e-12;
        for (size_t qi = 0; qi < queries.size(); ++qi) {
            unsigned int u = queries[qi].first;
            unsigned int v = queries[qi].second;
            auto *tu = sp_tree_vec[u];
            auto *tv = sp_tree_vec[v];
            double opt = 0.0;
            for (auto &kv : *edges_map) {
                unsigned int x = kv.first.first;
                unsigned int y = kv.first.second;
                double w = kv.second;
                if (!tu || !tv || !tu->at(x) || !tu->at(y) || !tv->at(x) || !tv->at(y)) continue;
                double c1 = w - tu->at(x)->path_length - tv->at(y)->path_length;
                double c2 = w - tu->at(y)->path_length - tv->at(x)->path_length;
                opt = std::max(opt, std::max(c1, c2));
            }
            if (opt <= EPS) {
                query_weights[qi] = 0.0;
                ++zero_opt_training_queries;
            } else {
                query_weights[qi] = 1.0 / opt;
            }
        }
    }
    cout << "[SGEL objective] mode=" << (relative_objective ? "relative_ARLB" : "raw_L1")
         << " training_queries=" << queries.size()
         << " zero_opt_training_queries=" << zero_opt_training_queries
         << " valid_training_queries=" << (queries.size() - zero_opt_training_queries)
         << endl;

    // Snapshot candidates for index-based OpenMP iteration.
    vector<pair<pair<unsigned int, unsigned int>, double>> candidates;
    candidates.reserve(edges_map->size());
    for (auto &kv : *edges_map) candidates.push_back(kv);

    // Correctness for --ablation-no-hca: the non-HCA score path calls
    // _lookup(u,v,e,w), which reads SPTs rooted at the candidate edge
    // endpoints e.first/e.second.  The HCA path avoids those candidate-root
    // Dijkstras by scoring from query-endpoint SPTs.  If we are explicitly
    // ablating HCA, materialize the required candidate-root SPTs before the
    // OpenMP scoring loop so the loop only performs read-only lookups.  This is
    // intentionally expensive and is exactly the cost of removing HCA.
    if (ablation_no_hca_) {
      cout << "[ablation-no-hca] materializing candidate endpoint SPTs for non-HCA SGEL scoring" << endl;
      for (auto &cand : candidates) {
        ensure_spt(cand.first.first);
        ensure_spt(cand.first.second);
      }
    }

    vector<double> current_best(queries.size(), 0.0);

    for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
        double best_gain = -1.0;
        pair<unsigned int, unsigned int> best_edge = {0, 0};
        double best_len = 0.0;

#ifdef _OPENMP
        #pragma omp parallel
        {
            // Per-thread HCA buffers (allocated once per parallel region).
            std::vector<unsigned int> tls_hca_vals(no_nodes, no_nodes);
            std::vector<uint64_t>     tls_hca_gen (no_nodes, 0);
            std::vector<unsigned int> tls_tmp_vals(no_nodes,
                std::numeric_limits<unsigned int>::max());
            std::vector<uint64_t>     tls_tmp_gen (no_nodes, 0);
            uint64_t my_gen = 0;

            double local_best_gain = -1.0;
            pair<unsigned int, unsigned int> local_best_edge = {0, 0};
            double local_best_len = 0.0;

            #pragma omp for schedule(dynamic, 64) nowait
            for (long long ci = 0; ci < (long long)candidates.size(); ++ci) {
                const auto &e = candidates[ci].first;
                if (landmarks->find(e) != landmarks->end()) continue;
                double e_len = candidates[ci].second;

                double gain = 0.0;
                for (size_t qi = 0; qi < queries.size(); ++qi) {
                    unsigned int u = queries[qi].first;
                    unsigned int v = queries[qi].second;
                    ++my_gen;  // O(1) reset (logical) between queries
                    double cand_lb = ablation_no_hca_
                        ? _lookup(u, v, e, e_len)
                        : lookup_HCA_local(
                            sp_tree_vec, u, v, e, e_len,
                            tls_hca_vals, tls_hca_gen,
                            tls_tmp_vals, tls_tmp_gen,
                            my_gen, no_nodes, ablation_no_spt_);
                    if (cand_lb > current_best[qi]) gain += query_weights[qi] * (cand_lb - current_best[qi]);
                }

                if (gain > local_best_gain) {
                    local_best_gain = gain;
                    local_best_edge = e;
                    local_best_len = e_len;
                }
            }

            #pragma omp critical
            {
                if (local_best_gain > best_gain) {
                    best_gain = local_best_gain;
                    best_edge = local_best_edge;
                    best_len = local_best_len;
                }
            }
        }
#else
        // Sequential fallback: identical to get_hca_landmarks() body.
        for (auto &cand : *edges_map) {
            if (landmarks->find(cand.first) != landmarks->end()) continue;
            const auto &e = cand.first;
            double e_len = cand.second;

            double gain = 0.0;
            for (size_t qi = 0; qi < queries.size(); ++qi) {
                unsigned int u = queries[qi].first;
                unsigned int v = queries[qi].second;
                _reset_HCA();
                double cand_lb = ablation_no_hca_
                               ? _lookup(u, v, e, e_len)
                               : _lookup_HCA(u, v, e, e_len);
                if (cand_lb > current_best[qi]) gain += query_weights[qi] * (cand_lb - current_best[qi]);
            }

            if (gain > best_gain) {
                best_gain = gain;
                best_edge = e;
                best_len = e_len;
            }
        }
#endif

        if (best_gain <= 0.0) break;
        landmarks->insert({best_edge, best_len});
        landmark_order.push_back({best_edge, best_len});

        ensure_spt(best_edge.first);
        ensure_spt(best_edge.second);

        // Update current_best (sequential is fine; size = query_samples).
        for (size_t qi = 0; qi < queries.size(); ++qi) {
            unsigned int u = queries[qi].first;
            unsigned int v = queries[qi].second;
            _reset_HCA();
            double val = ablation_no_hca_
                       ? _lookup(u, v, best_edge, best_len)
                       : _lookup_HCA(u, v, best_edge, best_len);
            if (val > current_best[qi]) current_best[qi] = val;
        }
    }

    clean_unwanted_shortest_paths();
}