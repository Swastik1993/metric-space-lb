#include "EdgeLandMark.h"

EdgeLandMark::EdgeLandMark(vector<list<pair<unsigned int, double>> *> *adj_list,
                           unsigned int n_nodes, unsigned int k,
                           unsigned int sampling_size) {
  cout << "No nodes : " << n_nodes << endl
       << "No landmarks : " << k << endl
       << "Samples : " << sampling_size << endl;
  no_nodes = n_nodes;
  no_landmarks = k;
  no_samples = sampling_size;
  this->adj_list = adj_list;
  store_map();
  sp_vector = new vector<vector<double> *>();
  le_vector = new vector<vector<double> *>();
  for (unsigned int i = 0; i < no_nodes; ++i) {
    pair<vector<double> *, vector<double> *> shortest_path_info =
        DijkstraELM(adj_list, i);
    sp_vector->push_back(shortest_path_info.first);
    le_vector->push_back(shortest_path_info.second);
  }
  landmarks = new map<pair<unsigned int, unsigned int>, double>();
}

long EdgeLandMark::_sizeof() {

  long size_val = sizeof(no_nodes) + sizeof(no_landmarks) + sizeof(no_samples) +
                  sizeof(*sp_vector) + sizeof(*le_vector) + sizeof(*landmarks) +
                  sizeof(*edges_map);
  for (unsigned int i = 0; i < sp_vector->size(); ++i) {
    size_val += sizeof(sp_vector->at(i));
  }
  for (unsigned int i = 0; i < le_vector->size(); ++i) {
    size_val += sizeof(le_vector->at(i));
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
  map<pair<unsigned int, unsigned int>, double> distance_map;
  unsigned int last_landmark_u, last_landmark_v;
  while (landmarks->size() != no_landmarks) {
    if (landmarks->size() == 0) {
      for (auto itr = edges_map->begin(); itr != edges_map->end(); ++itr) {
        if (itr == edges_map->begin()) {
          landmarks->insert({itr->first, itr->second});
          last_landmark_u = itr->first.first;
          last_landmark_v = itr->first.second;
        } else {
          distance_map.insert({itr->first, 0.0});
        }
      }
    } else {
      unsigned int min_u, min_v;
      double min_dist = no_nodes;
      for (auto itr = edges_map->begin(); itr != edges_map->end(); ++itr) {
        if (landmarks->find(itr->first) != landmarks->end()) {
          continue;
        }
        unsigned int x = itr->first.first, y = itr->first.second;
        distance_map.find(itr->first)->second +=
            min(sp_vector->at(last_landmark_u)->at(x) +
                    sp_vector->at(last_landmark_v)->at(y),
                sp_vector->at(last_landmark_v)->at(x) +
                    sp_vector->at(last_landmark_u)->at(y));
        if (distance_map.find(itr->first)->second < min_dist) {
          min_dist = distance_map.find(itr->first)->second;
          min_u = x;
          min_v = y;
        }
      }
      landmarks->insert(
          {make_pair(min_u, min_v), edges_map->at(make_pair(min_u, min_v))});
      distance_map.erase(make_pair(min_u, min_v));
    }
  }
}

void EdgeLandMark::get_exact_landmarks() {
  vector<double> estimated_nodes;
  for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
    estimated_nodes.push_back(0.);
  }

  while (landmarks->size() != no_landmarks) {
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

  while (landmarks->size() != no_landmarks) {
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

      for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
        estimated_nodes[index_i] =
            max(0., 2 * max(le_vector->at(index_i)->at(u),
                            le_vector->at(index_i)->at(v)) -
                        sp_vector->at(index_i)->at(u) -
                        sp_vector->at(index_i)->at(v));
      }

      double lb_landmarks = this->lookup(u, v);
      for (unsigned index_i = 0; index_i < no_nodes; ++index_i) {
        unsigned int node_x = index_i;
        for (list<pair<unsigned int, double>>::iterator it =
                 adj_list->at(index_i)->begin();
             it != adj_list->at(index_i)->end(); ++it) {
          if (it->first <= index_i) {
            continue;
          }
          unsigned int node_y = it->first;
          pair<unsigned int, unsigned int> pair_xy = make_pair(node_x, node_y);
          if (landmarks->find(pair_xy) != landmarks->end()) {
            continue;
          }

          double edge_len_xy = it->second;
          double ux = sp_vector->at(u)->at(node_x),
                 uy = sp_vector->at(u)->at(node_y);
          double vx = sp_vector->at(v)->at(node_x),
                 vy = sp_vector->at(v)->at(node_y);
          double lb_edge = 0.0;
          lb_edge = max({lb_edge, edge_len_xy - ux - vy, edge_len_xy - uy - vx,
                         estimated_nodes[node_x], estimated_nodes[node_y]});
          if (lb_edge > lb_landmarks) {
            if (contributions.find(pair_xy) == contributions.end()) {
              contributions.emplace(pair_xy, 0.0);
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
      landmarks->emplace(pair_out, out_edge_len);
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
