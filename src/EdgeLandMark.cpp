#include "EdgeLandMark.h"
#include <iterator>
#include <queue>

EdgeLandMark::EdgeLandMark(vector<list<pair<unsigned int, double>> *> *adj_list,
                           unsigned int n_nodes, unsigned int k,
                           unsigned int sampling_size,
                           bool precompute_all_pairs)
    : rng_(73 /* default seed; override with set_seed() */) {
  cout << "No nodes : " << n_nodes << endl
       << "No landmarks : " << k << endl
       << "Samples : " << sampling_size << endl;
  no_nodes = n_nodes;
  no_landmarks = k;
  no_samples = sampling_size;
  this->adj_list = adj_list;
  store_map();
  sp_vector = new vector<vector<double> *>(no_nodes, nullptr);
  le_vector = new vector<vector<double> *>(no_nodes, nullptr);
  if (precompute_all_pairs) {
    // Original O(n) Dijkstra eager fill. Required when caller depends on
    // sp_vector[*] being populated regardless of landmarks (e.g. H2 reads
    // last_landmark_*'s rows immediately).
    for (unsigned int i = 0; i < no_nodes; ++i) {
      pair<vector<double> *, vector<double> *> shortest_path_info =
          DijkstraELM(adj_list, i);
      (*sp_vector)[i] = shortest_path_info.first;
      (*le_vector)[i] = shortest_path_info.second;
    }
  }
  // else: callers must call ensure_sp_for_landmarks() before lookup().
  landmarks = new map<pair<unsigned int, unsigned int>, double>();
}

void EdgeLandMark::ensure_sp(unsigned int u) {
  if (u >= no_nodes) return;
  if ((*sp_vector)[u] == nullptr) {
    pair<vector<double> *, vector<double> *> shortest_path_info =
        DijkstraELM(adj_list, u);
    (*sp_vector)[u] = shortest_path_info.first;
    (*le_vector)[u] = shortest_path_info.second;
  }
}

void EdgeLandMark::ensure_sp_for_landmarks() {
  if (!landmarks) return;

  // Collect unique landmark endpoints whose SP/LE rows have not yet been
  // materialized.  The old version ran these DijkstraELM calls serially; for
  // k=128 this can mean up to 256 full SIFT Dijkstras on one core.
  vector<unsigned int> to_build;
  vector<char> seen(no_nodes, 0);
  for (auto &kv : *landmarks) {
    unsigned int endpoints[2] = {kv.first.first, kv.first.second};
    for (unsigned int v : endpoints) {
      if (v >= no_nodes) continue;
      if ((*sp_vector)[v] != nullptr) continue;
      if (seen[v]) continue;
      seen[v] = 1;
      to_build.push_back(v);
    }
  }

  if (to_build.empty()) return;
  cout << "[EdgeLandMark] materializing " << to_build.size()
       << " SP/LE rows for landmark endpoints" << endl;

  vector<pair<vector<double> *, vector<double> *>> built(to_build.size(),
                                                        {nullptr, nullptr});

#pragma omp parallel for schedule(dynamic, 1) if (to_build.size() > 1)
  for (long long i = 0; i < (long long)to_build.size(); ++i) {
    built[(size_t)i] = DijkstraELM(adj_list, to_build[(size_t)i]);
  }

  for (size_t i = 0; i < to_build.size(); ++i) {
    unsigned int v = to_build[i];
    if ((*sp_vector)[v] == nullptr) {
      (*sp_vector)[v] = built[i].first;
      (*le_vector)[v] = built[i].second;
    } else {
      // Defensive cleanup; should not happen in normal single-call usage.
      delete built[i].first;
      delete built[i].second;
    }
  }
}

long EdgeLandMark::_sizeof() {
  // Account for the underlying SP/LE data rows, not just the pointer slots.
  // Original implementation summed only sizeof(pointer) which understated the
  // real footprint by 6+ orders of magnitude on large graphs.
  long size_val = sizeof(no_nodes) + sizeof(no_landmarks) + sizeof(no_samples) +
                  sizeof(*sp_vector) + sizeof(*le_vector) + sizeof(*landmarks) +
                  sizeof(*edges_map);
  for (unsigned int i = 0; i < sp_vector->size(); ++i) {
    size_val += sizeof(void*); // pointer slot
    if ((*sp_vector)[i]) {
      size_val += (long)((*sp_vector)[i]->capacity() * sizeof(double) +
                         sizeof(*(*sp_vector)[i]));
    }
  }
  for (unsigned int i = 0; i < le_vector->size(); ++i) {
    size_val += sizeof(void*);
    if ((*le_vector)[i]) {
      size_val += (long)((*le_vector)[i]->capacity() * sizeof(double) +
                         sizeof(*(*le_vector)[i]));
    }
  }
  return size_val;
}

EdgeLandMark::~EdgeLandMark() {
  if (sp_vector) {
    for (auto it = this->sp_vector->begin(); it != this->sp_vector->end();
         ++it) {
      delete *it;
    }
    delete sp_vector;
  }
  if (le_vector) {
    for (auto it = this->le_vector->begin(); it != this->le_vector->end();
         ++it) {
      delete *it;
    }
    delete le_vector;
  }
  delete landmarks;
  delete edges_map;
  // Free cached SPTs if any were built (lookup_lca path).
  for (auto &kv : lca_spt_cache_) {
    if (!kv.second) continue;
    for (auto *node : *kv.second) delete node;
    delete kv.second;
  }
  lca_spt_cache_.clear();
}

void EdgeLandMark::store_map() {
  edges_map = new map<pair<unsigned int, unsigned int>, double>();
  for (unsigned int i = 0; i < adj_list->size(); ++i) {
    for (auto itr = adj_list->at(i)->begin(); itr != adj_list->at(i)->end();
         ++itr) {
      unsigned int j = itr->first;
      double dist = itr->second;
      unsigned int u = min(i, j), v = max(i, j);
      pair<unsigned int, unsigned int> pr = make_pair(u, v);
      if (edges_map->find(pr) == edges_map->end()) {
        edges_map->insert({pr, dist});
      }
    }
  }
}

void EdgeLandMark::large_edge_heuristic() {
  pairing_heap<min_heap_edge> H;
  for (auto itr = edges_map->begin(); itr != edges_map->end(); ++itr) {
    double dist = itr->second;
    if (H.size() < no_landmarks) {
      H.push(min_heap_edge(itr->first.first, itr->first.second, dist));
    } else if (H.size() >= no_landmarks && H.top().dist < dist) {
      H.pop();
      H.push(min_heap_edge(itr->first.first, itr->first.second, dist));
    }
  }
  while (H.size() > 0) {
    min_heap_edge edge = H.top();
    H.pop();
    landmarks->insert({make_pair(edge.u, edge.v), edge.dist});
  }
}

void EdgeLandMark::far_away_heuristic() {
  cout << "[EdgeLandMark::H2] selecting " << no_landmarks
       << " far-away landmarks from " << edges_map->size() << " edges" << endl;

  if (no_landmarks == 0 || edges_map->empty()) return;

  // Correct farthest-edge heuristic:
  //   1. Start from the first available edge, matching the original seed.
  //   2. Maintain for every candidate edge the accumulated distance to the
  //      already selected landmark edges.
  //   3. At each round, add the distance to the most recently selected
  //      landmark edge and select the candidate with the MAXIMUM accumulated
  //      distance.
  //
  // The previous code had two coupled issues:
  //   - it did not update last_landmark_u / last_landmark_v after inserting a
  //     new H2 landmark, so every later round used only the first edge again;
  //   - it chose the MINIMUM accumulated distance even though this is the
  //     "far-away" heuristic.
  // Both are fixed here. The landmark set will therefore differ from older
  // runs, but it now matches the intended farthest-first behavior.

  auto materialize_rows_for_pair = [&](unsigned int a, unsigned int b) {
    vector<unsigned int> needed;
    if (a < no_nodes && (*sp_vector)[a] == nullptr) needed.push_back(a);
    if (b < no_nodes && b != a && (*sp_vector)[b] == nullptr) needed.push_back(b);
    if (needed.empty()) return;

    vector<pair<vector<double> *, vector<double> *>> built(needed.size(),
                                                          {nullptr, nullptr});
#pragma omp parallel for schedule(dynamic, 1) if (needed.size() > 1)
    for (long long i = 0; i < (long long)needed.size(); ++i) {
      built[(size_t)i] = DijkstraELM(adj_list, needed[(size_t)i]);
    }

    for (size_t i = 0; i < needed.size(); ++i) {
      unsigned int v = needed[i];
      if ((*sp_vector)[v] == nullptr) {
        (*sp_vector)[v] = built[i].first;
        (*le_vector)[v] = built[i].second;
      } else {
        delete built[i].first;
        delete built[i].second;
      }
    }
  };

  const size_t target_landmarks =
      min((size_t)no_landmarks, (size_t)edges_map->size());

  // v18: flatten distance_map into a vector indexed by edge position in the
  // ordered traversal of edges_map. This cuts the per-iter inner loop from
  // O(m log m) (std::map find + sorted iteration) to O(m) plain array access
  // and makes the loop trivially parallel over a contiguous range.
  const size_t M = edges_map->size();
  vector<pair<unsigned int, unsigned int>> edges_vec;
  vector<double> edge_len;
  edges_vec.reserve(M);
  edge_len.reserve(M);
  for (auto &kv : *edges_map) {
    edges_vec.push_back(kv.first);
    edge_len.push_back(kv.second);
  }
  vector<double> accumulated_dist(M, 0.0);
  vector<char> is_landmark(M, 0);

  // Seed with the first edge in std::map order. This matches the seed point
  // chosen by the legacy code, so the first landmark is identical.
  size_t seed_idx = 0;
  is_landmark[seed_idx] = 1;
  landmarks->insert({edges_vec[seed_idx], edge_len[seed_idx]});
  unsigned int last_landmark_u = edges_vec[seed_idx].first;
  unsigned int last_landmark_v = edges_vec[seed_idx].second;
  materialize_rows_for_pair(last_landmark_u, last_landmark_v);
  cout << "[EdgeLandMark::H2] selected 1/" << target_landmarks
       << " landmarks" << endl;

  while (landmarks->size() < target_landmarks) {
    const vector<double> &sp_lu = *sp_vector->at(last_landmark_u);
    const vector<double> &sp_lv = *sp_vector->at(last_landmark_v);

    // Per-thread local reduction; final merge in a critical section.
    // Tie-breaking: smaller index wins. Since edges_vec preserves the
    // std::map traversal order, smaller index == smaller edge key, which
    // matches the legacy tie-breaking and keeps the landmark set
    // bit-identical across OMP thread counts.
    double best_dist = -1.0;
    size_t best_idx = M;

#pragma omp parallel
    {
      double t_best_dist = -1.0;
      size_t t_best_idx = M;

#pragma omp for schedule(static) nowait
      for (long long i = 0; i < (long long)M; ++i) {
        if (is_landmark[(size_t)i]) continue;
        unsigned int x = edges_vec[(size_t)i].first;
        unsigned int y = edges_vec[(size_t)i].second;
        double a = sp_lu[x] + sp_lv[y];
        double b = sp_lv[x] + sp_lu[y];
        accumulated_dist[(size_t)i] += (a < b) ? a : b;
        double d = accumulated_dist[(size_t)i];
        if (d > t_best_dist ||
            (d == t_best_dist && (size_t)i < t_best_idx)) {
          t_best_dist = d;
          t_best_idx = (size_t)i;
        }
      }

#pragma omp critical
      {
        if (t_best_idx < M) {
          if (t_best_dist > best_dist ||
              (t_best_dist == best_dist && t_best_idx < best_idx)) {
            best_dist = t_best_dist;
            best_idx = t_best_idx;
          }
        }
      }
    }

    if (best_idx >= M) break;

    is_landmark[best_idx] = 1;
    landmarks->insert({edges_vec[best_idx], edge_len[best_idx]});

    // Algorithmic fix (v17): make the newest landmark the reference edge for
    // the next distance-update round. Missing in the pre-v17 code.
    last_landmark_u = edges_vec[best_idx].first;
    last_landmark_v = edges_vec[best_idx].second;
    materialize_rows_for_pair(last_landmark_u, last_landmark_v);

    if (landmarks->size() == target_landmarks || landmarks->size() % 8 == 0) {
      cout << "[EdgeLandMark::H2] selected " << landmarks->size()
           << "/" << target_landmarks << " landmarks" << endl;
    }
  }
}

void EdgeLandMark::get_exact_landmarks() {
  // Exact greedy scoring needs candidate-root SP/LE rows for all vertices.
  // Materialize lazily here so the constructor no longer pays this cost
  // unless exact GEL is actually requested.
  for (unsigned int i = 0; i < no_nodes; ++i) {
    ensure_sp(i);
  }

  vector<double> estimated_nodes;
  for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
    estimated_nodes.push_back(0.);
  }

  while (landmarks->size() < no_landmarks) {
    if (landmarks->size() > 0) {
      cout << "Selected " << landmarks->size() << " number of landmarks"
           << endl;
    }
    double value_out = 0.0;
    unsigned int out_x = 0;
    unsigned int out_y = 0;
    double out_edge_len = 0.0;
    map<pair<unsigned int, unsigned int>, double> contributions;

    for (unsigned int i = 0; i < no_nodes; ++i) {
      for (unsigned int j = i + 1; j < no_nodes; ++j) {
        unsigned int u = min(i, j), v = max(i, j);
        if (edges_map->find(make_pair(u, v)) != edges_map->end()) {
          continue;
        }
        double lb_landmarks = this->lookup(u, v);
        for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
          estimated_nodes[index_i] =
              max(0., 2 * max(le_vector->at(index_i)->at(u),
                              le_vector->at(index_i)->at(v)) -
                          sp_vector->at(index_i)->at(u) -
                          sp_vector->at(index_i)->at(v));
        }

        for (auto itr = edges_map->begin(); itr != edges_map->end(); ++itr) {
          unsigned int node_x = itr->first.first, node_y = itr->first.second;
          pair<unsigned int, unsigned int> pair_xy = make_pair(node_x, node_y);
          if (landmarks->find(pair_xy) != landmarks->end()) {
            continue;
          }

          double edge_len_xy = itr->second;
          double ux = sp_vector->at(u)->at(node_x),
                 uy = sp_vector->at(u)->at(node_y);
          double vx = sp_vector->at(v)->at(node_x),
                 vy = sp_vector->at(v)->at(node_y);
          double lb_edge = 0.0;
          lb_edge = max({lb_edge, edge_len_xy - ux - vy, edge_len_xy - uy - vx,
                         estimated_nodes[node_x], estimated_nodes[node_y]});
          if (lb_edge > lb_landmarks) {
            if (contributions.find(pair_xy) == contributions.end()) {
              contributions.insert({pair_xy, 0.0});
            }
            contributions[pair_xy] += lb_edge - lb_landmarks;
            if (value_out < contributions[pair_xy]) {
              value_out = contributions[pair_xy];
              out_x = node_x;
              out_y = node_y;
              out_edge_len = edge_len_xy;
            }
          }
        }
      }
    }
    if (value_out == 0.0) {
      cout << "Something wrong here!" << endl;
    } else {
      pair<unsigned int, unsigned int> pair_out = make_pair(out_x, out_y);
      landmarks->insert({pair_out, out_edge_len});
    }
  }
  clean_unwanted_shortest_paths();
}

void EdgeLandMark::clean_unwanted_shortest_paths() {
  for (unsigned int index_i = 0; index_i < no_nodes; ++index_i) {
    bool present = false;
    for (auto it = landmarks->begin(); it != landmarks->end(); ++it) {
      if (it->first.first == index_i || it->first.second == index_i) {
        present = true;
        break;
      }
    }
    if (!present) {
      delete sp_vector->at(index_i);
      delete le_vector->at(index_i);
      sp_vector->at(index_i) = NULL;
      le_vector->at(index_i) = NULL;
    }
  }
}

void EdgeLandMark::get_landmarks() {
  vector<double> estimated_nodes;
  for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
    estimated_nodes.push_back(0.);
  }

  while (landmarks->size() < no_landmarks) {
    if (landmarks->size() > 0) {
      cout << "Selected " << landmarks->size() << " number of landmarks"
           << endl;
    }
    double value_out = 0.0;
    unsigned int out_x = 0;
    unsigned int out_y = 0;
    double out_edge_len = 0.0;
    map<pair<unsigned int, unsigned int>, double> contributions;
    for (unsigned int i = 0; i < no_samples; ++i) {
      unsigned int u = rand() % this->no_nodes;
      unsigned int v = rand() % this->no_nodes;

      while (edges_map->find(make_pair(min(u, v), max(u, v))) !=
             edges_map->end()) {
        u = rand() % this->no_nodes;
        v = rand() % this->no_nodes;
      }

      // Lazy SPT materialization for the sampled query endpoints.
      // This supports the binary-lifting/SPT ablation measurement without
      // forcing all-source SP/LE construction in the constructor.
      ensure_sp(u);
      ensure_sp(v);
      for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
        estimated_nodes[index_i] =
            max(0., 2 * max(le_vector->at(u)->at(index_i),
                            le_vector->at(v)->at(index_i)) -
                        sp_vector->at(u)->at(index_i) -
                        sp_vector->at(v)->at(index_i));
      }

      double lb_landmarks = this->lookup(u, v);
      // for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
      //   unsigned int node_x = index_i;
      //   for (list<pair<unsigned int, double>>::iterator it =
      //            adj_list->at(index_i)->begin();
      //        it != adj_list->at(index_i)->end(); ++it) {
      //     if (it->first <= index_i) {
      //       continue;
      //     }
      //     unsigned int node_y = it->first;
      //     pair<unsigned int, unsigned int> pair_xy = make_pair(node_x, node_y);
      //     if (landmarks->find(pair_xy) != landmarks->end()) {
      //       continue;
      //     }

      //     double edge_len_xy = it->second;
      //     double ux = sp_vector->at(u)->at(node_x),
      //            uy = sp_vector->at(u)->at(node_y);
      //     double vx = sp_vector->at(v)->at(node_x),
      //            vy = sp_vector->at(v)->at(node_y);
      //     double lb_edge = 0.0;
      //     lb_edge = max({lb_edge, edge_len_xy - ux - vy, edge_len_xy - uy - vx,
      //                    estimated_nodes[node_x], estimated_nodes[node_y]});
      //     if (lb_edge > lb_landmarks) {
      //       if (contributions.find(pair_xy) == contributions.end()) {
      //         contributions.emplace(pair_xy, 0.0);
      //       }
      //       contributions[pair_xy] += lb_edge - lb_landmarks;
      //       if (value_out < contributions[pair_xy]) {
      //         value_out = contributions[pair_xy];
      //         out_x = node_x;
      //         out_y = node_y;
      //         out_edge_len = edge_len_xy;
      //       }
      //     }
      //   }
      // }
      for (auto itr = edges_map->begin(); itr != edges_map->end(); ++itr) {
        unsigned int node_x = itr->first.first, node_y = itr->first.second;
        pair<unsigned int, unsigned int> pair_xy = make_pair(node_x, node_y);
        if (landmarks->find(pair_xy) != landmarks->end()) {
          continue;
        }

        double edge_len_xy = itr->second;
        double ux = sp_vector->at(u)->at(node_x),
                uy = sp_vector->at(u)->at(node_y);
        double vx = sp_vector->at(v)->at(node_x),
                vy = sp_vector->at(v)->at(node_y);
        double lb_edge = 0.0;
        lb_edge = max({lb_edge, edge_len_xy - ux - vy, edge_len_xy - uy - vx,
                        estimated_nodes[node_x], estimated_nodes[node_y]});
        if (lb_edge > lb_landmarks) {
          if (contributions.find(pair_xy) == contributions.end()) {
            contributions.insert({pair_xy, 0.0});
          }
          contributions[pair_xy] += lb_edge - lb_landmarks;
          if (value_out < contributions[pair_xy]) {
            value_out = contributions[pair_xy];
            out_x = node_x;
            out_y = node_y;
            out_edge_len = edge_len_xy;
          }
        }
      }
    }
    if (value_out == 0.0) {
      cout << "Something wrong here!" << endl;
    } else {
      pair<unsigned int, unsigned int> pair_out = make_pair(out_x, out_y);
      // landmarks->emplace(pair_out, out_edge_len);
      landmarks->insert({pair_out, out_edge_len});
      ensure_sp(out_x);
      ensure_sp(out_y);
      cout << "Adding " << out_x << " and " << out_y << " into landmarks."
           << endl;
    }
  }
}

double EdgeLandMark::lookup(unsigned int u, unsigned int v) {
  if (edges_map->find(make_pair(min(u, v), max(u, v))) != edges_map->end()) {
    return edges_map->at(make_pair(min(u, v), max(u, v)));
  }
  double lb_landmarks = 0.0;
  for (auto it = landmarks->begin(); it != landmarks->end(); ++it) {
    unsigned int x = it->first.first, y = it->first.second;
    ensure_sp(x);
    ensure_sp(y);
    double edge_len_landmark = it->second;
    double ux = sp_vector->at(x)->at(u), uy = sp_vector->at(y)->at(u);
    double vx = sp_vector->at(x)->at(v), vy = sp_vector->at(y)->at(v);
    double max_along_x =
        2 * max(le_vector->at(x)->at(u), le_vector->at(x)->at(v)) -
        sp_vector->at(x)->at(u) - sp_vector->at(x)->at(v);
    double max_along_y =
        2 * max(le_vector->at(y)->at(u), le_vector->at(y)->at(v)) -
        sp_vector->at(y)->at(u) - sp_vector->at(y)->at(v);
    lb_landmarks = max({lb_landmarks, edge_len_landmark - ux - vy,
                        edge_len_landmark - uy - vx, max_along_x, max_along_y});
  }
  return lb_landmarks;
}

double EdgeLandMark::lookup_ub(unsigned int u, unsigned int v) {
  if (edges_map->find(make_pair(min(u, v), max(u, v))) != edges_map->end()) {
    return edges_map->at(make_pair(min(u, v), max(u, v)));
  }
  double ub_landmarks = 1.0;
  for (auto it = landmarks->begin(); it != landmarks->end(); ++it) {
    unsigned int x = it->first.first, y = it->first.second;
    ensure_sp(x);
    ensure_sp(y);
    double edge_len_landmark = it->second;
    double ux = sp_vector->at(x)->at(u), uy = sp_vector->at(y)->at(u);
    double vx = sp_vector->at(x)->at(v), vy = sp_vector->at(y)->at(v);
    ub_landmarks =
        min({ub_landmarks, ux + vx, uy + vy, ux + vy + edge_len_landmark,
             uy + vx + edge_len_landmark});
  }
  return ub_landmarks;
}

vector<double> *EdgeLandMark::lookup_multiple(unsigned int u, unsigned int v) {
  vector<double> *lb_vals = new vector<double>();
  double lb_landmarks = 0.0;
  for (auto it = landmarks->begin(); it != landmarks->end(); ++it) {
    unsigned int x = it->first.first, y = it->first.second;
    ensure_sp(x);
    ensure_sp(y);
    double edge_len_landmark = it->second;
    double ux = sp_vector->at(x)->at(u), uy = sp_vector->at(y)->at(u);
    double vx = sp_vector->at(x)->at(v), vy = sp_vector->at(y)->at(v);
    double max_along_x =
        2 * max(le_vector->at(x)->at(u), le_vector->at(x)->at(v)) -
        sp_vector->at(x)->at(u) - sp_vector->at(x)->at(v);
    double max_along_y =
        2 * max(le_vector->at(y)->at(u), le_vector->at(y)->at(v)) -
        sp_vector->at(y)->at(u) - sp_vector->at(y)->at(v);
    lb_landmarks = max({lb_landmarks, edge_len_landmark - ux - vy,
                        edge_len_landmark - uy - vx, max_along_x, max_along_y});
    lb_vals->push_back(lb_landmarks);
  }
  return lb_vals;
}

// ============================================================================
// New: dense-graph-tolerant sampler (paper §4.2 PDF-based alternative).
// ============================================================================
// Strategy: when missing edges are abundant (rejection probability of hitting
// a known edge is < ~50%), use plain rejection sampling. As the graph
// approaches dense, switch to a two-step sampler:
//   1. Pick first endpoint u with probability proportional to its number of
//      missing edges (n - 1 - deg(u)).
//   2. Pick second endpoint v uniformly from the missing neighbors of u.
// This avoids the 1/p blowup of rejection sampling as p -> 0.
pair<unsigned int, unsigned int> EdgeLandMark::sample_missing_edge_robust() {
  const uint64_t total_pairs = (uint64_t)no_nodes * (uint64_t)(no_nodes - 1) / 2ULL;
  const uint64_t known_pairs = (uint64_t)edges_map->size();
  const uint64_t missing_pairs = total_pairs - known_pairs;

  if (missing_pairs == 0) {
    // No missing edges: caller should have checked. Return (0,0) as a
    // sentinel; behavior is undefined but we don't infinite-loop.
    return {0u, 0u};
  }

  // Switch threshold: use rejection sampling when at least half the pairs are
  // missing (so expected attempts <= 2). Otherwise build the PDF.
  if (missing_pairs * 2 >= total_pairs) {
    std::uniform_int_distribution<unsigned int> uni(0, no_nodes - 1);
    while (true) {
      unsigned int u = uni(rng_);
      unsigned int v = uni(rng_);
      if (u == v) continue;
      auto key = make_pair(min(u, v), max(u, v));
      if (edges_map->find(key) == edges_map->end()) return {u, v};
    }
  }

  // Dense-graph path: construct per-vertex missing-degree weights once per
  // call (cheap relative to the alternative of long rejection chains).
  vector<double> weights(no_nodes, 0.0);
  for (unsigned int i = 0; i < no_nodes; ++i) {
    weights[i] = double(no_nodes - 1) - double(adj_list->at(i)->size());
    if (weights[i] < 0.0) weights[i] = 0.0;
  }
  std::discrete_distribution<unsigned int> pick_u(weights.begin(), weights.end());

  for (int attempt = 0; attempt < 64; ++attempt) {
    unsigned int u = pick_u(rng_);
    // Build the set of known neighbors of u (small for sparse rows).
    vector<bool> is_known(no_nodes, false);
    is_known[u] = true;
    for (const auto &pr : *adj_list->at(u)) is_known[pr.first] = true;
    // Reservoir-sample one missing neighbor.
    unsigned int chosen = no_nodes;
    unsigned int seen = 0;
    std::uniform_int_distribution<unsigned int> uni01(0, 0);
    for (unsigned int v = 0; v < no_nodes; ++v) {
      if (is_known[v]) continue;
      ++seen;
      std::uniform_int_distribution<unsigned int> pick(0, seen - 1);
      if (pick(rng_) == 0) chosen = v;
    }
    if (chosen != no_nodes) return {u, chosen};
  }
  // Extreme fallback: shouldn't happen unless the discrete distribution
  // picks vertices that are themselves fully-known (degree n-1).
  return {0u, no_nodes > 1 ? 1u : 0u};
}

void EdgeLandMark::get_landmarks_seeded(unsigned int seed) {
  rng_.seed(seed);

  vector<double> estimated_nodes(no_nodes, 0.0);

  while (landmarks->size() < no_landmarks) {
    if (landmarks->size() > 0) {
      cout << "Selected " << landmarks->size() << " number of landmarks (lazy seeded)" << endl;
    }
    double value_out = 0.0;
    unsigned int out_x = 0, out_y = 0;
    double out_edge_len = 0.0;
    map<pair<unsigned int, unsigned int>, double> contributions;

    for (unsigned int i = 0; i < no_samples; ++i) {
      auto sample = sample_missing_edge_robust();
      unsigned int u = sample.first, v = sample.second;
      if (u == v) continue;

      ensure_sp(u);
      ensure_sp(v);
      for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
        estimated_nodes[index_i] = max(0.,
            2 * max(le_vector->at(u)->at(index_i),
                    le_vector->at(v)->at(index_i))
            - sp_vector->at(u)->at(index_i)
            - sp_vector->at(v)->at(index_i));
      }

      double lb_landmarks = this->lookup(u, v);
      for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
        unsigned int node_x = index_i;
        for (auto it = adj_list->at(index_i)->begin();
             it != adj_list->at(index_i)->end(); ++it) {
          if (it->first <= index_i) continue;
          unsigned int node_y = it->first;
          auto pair_xy = make_pair(node_x, node_y);
          if (landmarks->find(pair_xy) != landmarks->end()) continue;

          double edge_len_xy = it->second;
          double ux = sp_vector->at(u)->at(node_x),
                 uy = sp_vector->at(u)->at(node_y);
          double vx = sp_vector->at(v)->at(node_x),
                 vy = sp_vector->at(v)->at(node_y);
          double lb_edge = max({0.0,
                                edge_len_xy - ux - vy,
                                edge_len_xy - uy - vx,
                                estimated_nodes[node_x],
                                estimated_nodes[node_y]});
          if (lb_edge > lb_landmarks) {
            auto it2 = contributions.find(pair_xy);
            if (it2 == contributions.end()) {
              it2 = contributions.emplace(pair_xy, 0.0).first;
            }
            it2->second += lb_edge - lb_landmarks;
            if (value_out < it2->second) {
              value_out = it2->second;
              out_x = node_x;
              out_y = node_y;
              out_edge_len = edge_len_xy;
            }
          }
        }
      }
    }

    if (value_out == 0.0) {
      cout << "[get_landmarks_seeded] No positive-gain candidate; stopping early." << endl;
      break;
    }
    landmarks->emplace(make_pair(out_x, out_y), out_edge_len);
    ensure_sp(out_x);
    ensure_sp(out_y);
  }
  clean_unwanted_shortest_paths();
}

// ============================================================================
// New: lazy / accelerated greedy (Minoux 1978) for GEL.
// ============================================================================
// Submodular monotone maximization admits the "lazy evaluation" speedup: the
// marginal gain of any candidate at iteration t is upper-bounded by its gain
// at any earlier iteration. So we keep a max-heap of (cached_gain,
// iter_when_cached, edge); each round we pop the top, refresh its gain only
// if its iter_when_cached < current_iter, and push back. The first popped
// candidate that's already up-to-date is the true argmax.
//
// Same (1 - 1/e) guarantee as get_exact_landmarks(); typically much faster.
namespace {
struct LazyGainEntry {
  double gain;
  unsigned int last_iter;
  pair<unsigned int, unsigned int> edge;
  double edge_len;
  bool operator<(const LazyGainEntry &o) const {
    // max-heap: top has the largest gain. Tie-break on edge for determinism.
    if (gain != o.gain) return gain < o.gain;
    if (edge.first != o.edge.first) return edge.first > o.edge.first;
    return edge.second > o.edge.second;
  }
};
} // namespace

void EdgeLandMark::get_exact_landmarks_lazy() {
  // Lazy exact GEL still needs candidate-root SP/LE rows for all candidate
  // endpoints. Materialize them here, not in the constructor.
  for (unsigned int i = 0; i < no_nodes; ++i) {
    ensure_sp(i);
  }

  // Cache current best LB per (i, j) missing pair so marginal gain
  // computation against an additional candidate is straightforward.
  vector<pair<unsigned int, unsigned int>> queries;
  queries.reserve(((size_t)no_nodes * (size_t)(no_nodes - 1)) / 2);
  for (unsigned int i = 0; i < no_nodes; ++i) {
    for (unsigned int j = i + 1; j < no_nodes; ++j) {
      auto key = make_pair(i, j);
      if (edges_map->find(key) == edges_map->end()) queries.push_back(key);
    }
  }
  if (queries.empty()) return;

  vector<double> current_best(queries.size(), 0.0);

  // Helper: marginal gain of edge (x, y, w) over current_best.
  auto marginal_gain = [&](unsigned int x, unsigned int y, double w) -> double {
    double gain = 0.0;
    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first, v = queries[qi].second;
      double ux = sp_vector->at(x)->at(u), uy = sp_vector->at(y)->at(u);
      double vx = sp_vector->at(x)->at(v), vy = sp_vector->at(y)->at(v);
      double max_along_x =
          2 * max(le_vector->at(x)->at(u), le_vector->at(x)->at(v))
          - sp_vector->at(x)->at(u) - sp_vector->at(x)->at(v);
      double max_along_y =
          2 * max(le_vector->at(y)->at(u), le_vector->at(y)->at(v))
          - sp_vector->at(y)->at(u) - sp_vector->at(y)->at(v);
      double cand = max({0.0, w - ux - vy, w - uy - vx, max_along_x, max_along_y});
      if (cand > current_best[qi]) gain += cand - current_best[qi];
    }
    return gain;
  };

  // Initial pass: compute every candidate's gain at iter 0 and push.
  std::priority_queue<LazyGainEntry> heap;
  for (auto &kv : *edges_map) {
    if (landmarks->find(kv.first) != landmarks->end()) continue;
    double g = marginal_gain(kv.first.first, kv.first.second, kv.second);
    heap.push({g, 0u, kv.first, kv.second});
  }

  for (unsigned int iter = 0; iter < no_landmarks; ++iter) {
    if (heap.empty()) break;
    LazyGainEntry top = heap.top();
    while (top.last_iter < iter) {
      // Stale: refresh and re-insert.
      heap.pop();
      double fresh = marginal_gain(top.edge.first, top.edge.second, top.edge_len);
      heap.push({fresh, iter, top.edge, top.edge_len});
      top = heap.top();
    }
    // Top is fresh and is the argmax by submodularity.
    if (top.gain <= 0.0) {
      cout << "[get_exact_landmarks_lazy] No positive gain; stopping early at iter "
           << iter << endl;
      break;
    }
    landmarks->insert({top.edge, top.edge_len});
    heap.pop();

    // Apply commitment: refresh current_best for this landmark.
    unsigned int x = top.edge.first, y = top.edge.second;
    double w = top.edge_len;
    for (size_t qi = 0; qi < queries.size(); ++qi) {
      unsigned int u = queries[qi].first, v = queries[qi].second;
      double ux = sp_vector->at(x)->at(u), uy = sp_vector->at(y)->at(u);
      double vx = sp_vector->at(x)->at(v), vy = sp_vector->at(y)->at(v);
      double max_along_x =
          2 * max(le_vector->at(x)->at(u), le_vector->at(x)->at(v))
          - sp_vector->at(x)->at(u) - sp_vector->at(x)->at(v);
      double max_along_y =
          2 * max(le_vector->at(y)->at(u), le_vector->at(y)->at(v))
          - sp_vector->at(y)->at(u) - sp_vector->at(y)->at(v);
      double cand = max({0.0, w - ux - vy, w - uy - vx, max_along_x, max_along_y});
      if (cand > current_best[qi]) current_best[qi] = cand;
    }
  }
}

// ============================================================================
// New: batched per-source LB fill (kNN-friendly).
// ============================================================================
void EdgeLandMark::lookup_from_source(unsigned int u, vector<double> &LB) {
  if (LB.size() != (size_t)no_nodes) LB.assign(no_nodes, 0.0);
  else std::fill(LB.begin(), LB.end(), 0.0);

  // Walk landmarks once; for each, update every destination v in a single pass.
  // This is dramatically more cache-friendly than calling lookup(u, v) n times.
  for (auto &lm : *landmarks) {
    unsigned int x = lm.first.first;
    unsigned int y = lm.first.second;
    ensure_sp(x);
    ensure_sp(y);
    double w = lm.second;
    const vector<double> *spx = sp_vector->at(x);
    const vector<double> *spy = sp_vector->at(y);
    const vector<double> *lex = le_vector->at(x);
    const vector<double> *ley = le_vector->at(y);
    double ux = spx->at(u), uy = spy->at(u);
    double leu_x = lex->at(u), leu_y = ley->at(u);
    for (unsigned int v = 0; v < no_nodes; ++v) {
      if (v == u) continue;
      double vx = spx->at(v), vy = spy->at(v);
      double lev_x = lex->at(v), lev_y = ley->at(v);
      double max_along_x = 2 * max(leu_x, lev_x) - ux - vx;
      double max_along_y = 2 * max(leu_y, lev_y) - uy - vy;
      double cand = max({0.0, w - ux - vy, w - uy - vx, max_along_x, max_along_y});
      if (cand > LB[v]) LB[v] = cand;
    }
  }
}

// ============================================================================
// LCA-based lookup (proper §4.1 bound, Eq. 7).
// ============================================================================
void EdgeLandMark::clear_lca_cache() {
  for (auto &kv : lca_spt_cache_) {
    if (!kv.second) continue;
    for (auto *node : *kv.second) delete node;
    delete kv.second;
  }
  lca_spt_cache_.clear();
}

// Internal helper: get (or build) a fully-preprocessed SPT rooted at `r`.
namespace {
inline vector<shortest_path_tree*>* edgelandmark_ensure_spt(
    std::map<unsigned int, std::vector<shortest_path_tree*>*> &cache,
    vector<list<pair<unsigned int, double>>*> *adj_list,
    unsigned int r) {
  auto it = cache.find(r);
  if (it != cache.end()) return it->second;
  // Build with full preprocessing (lifting + RMQ-LCA + Cartesian) so both
  // find_LCA (used here) and find_HCA (potential future user) work.
  auto *t = Dijkstra(adj_list, r, true);
  cache[r] = t;
  return t;
}

// Per-vertex LCA-based folding term using a tree rooted at `root_v`:
//   2 * max_edge(VV(u, v) in tau_root_v)  -  pathlen(u, v in tau_root_v)
inline double lca_vertex_term(vector<shortest_path_tree*> *tree,
                              unsigned int u, unsigned int v) {
  auto info = find_LCA(tree, u, v);
  shortest_path_tree *lca_node = info.first;
  if (!lca_node) return 0.0;
  double max_edge_on_path = info.second;
  double dist_root_u   = tree->at(u)->path_length;
  double dist_root_v   = tree->at(v)->path_length;
  double dist_root_lca = lca_node->path_length;
  double tree_dist_uv = dist_root_u + dist_root_v - 2.0 * dist_root_lca;
  return 2.0 * max_edge_on_path - tree_dist_uv;
}
} // namespace

double EdgeLandMark::lookup_lca(unsigned int u, unsigned int v) {
  // Known-edge fast path.
  if (edges_map->find(make_pair(min(u, v), max(u, v))) != edges_map->end()) {
    return edges_map->at(make_pair(min(u, v), max(u, v)));
  }
  double best = 0.0;
  for (auto &lm : *landmarks) {
    unsigned int a = lm.first.first, b = lm.first.second;
    double w = lm.second;
    auto *treeA = edgelandmark_ensure_spt(lca_spt_cache_, adj_list, a);
    auto *treeB = edgelandmark_ensure_spt(lca_spt_cache_, adj_list, b);

    double au = treeA->at(u)->path_length;
    double av = treeA->at(v)->path_length;
    double bu = treeB->at(u)->path_length;
    double bv = treeB->at(v)->path_length;

    best = std::max(best, w - au - bv);
    best = std::max(best, w - av - bu);
    best = std::max(best, lca_vertex_term(treeA, u, v));
    best = std::max(best, lca_vertex_term(treeB, u, v));
  }
  if (best < 0.0) best = 0.0;
  return best;
}

vector<double> *EdgeLandMark::lookup_multiple_lca(unsigned int u, unsigned int v) {
  vector<double> *out = new vector<double>();
  out->reserve(landmarks->size());
  if (edges_map->find(make_pair(min(u, v), max(u, v))) != edges_map->end()) {
    double known = edges_map->at(make_pair(min(u, v), max(u, v)));
    for (size_t i = 0; i < landmarks->size(); ++i) out->push_back(known);
    return out;
  }
  double running = 0.0;
  for (auto &lm : *landmarks) {
    unsigned int a = lm.first.first, b = lm.first.second;
    double w = lm.second;
    auto *treeA = edgelandmark_ensure_spt(lca_spt_cache_, adj_list, a);
    auto *treeB = edgelandmark_ensure_spt(lca_spt_cache_, adj_list, b);
    double au = treeA->at(u)->path_length;
    double av = treeA->at(v)->path_length;
    double bu = treeB->at(u)->path_length;
    double bv = treeB->at(v)->path_length;
    double cand = std::max({0.0, w - au - bv, w - av - bu,
                            lca_vertex_term(treeA, u, v),
                            lca_vertex_term(treeB, u, v)});
    if (cand > running) running = cand;
    out->push_back(running);
  }
  return out;
}



