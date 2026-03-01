#include "GraphUtils.h"
#include <numbers>
#include <algorithm> 

int check_connected(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id) {
  int count = 0;
  unsigned int nodes = boost::num_vertices(g);
  if (nodes == 0) {
    return count;
  }
  list<unsigned int> queue;
  unsigned int total_visited = 0;
  vector<bool> visited(nodes, false);
  while (total_visited != nodes) {
    ++count;
    queue.clear();
    for (unsigned int i = 0; i < nodes; ++i) {
      if (!visited.at(i)) {
        ++total_visited;
        visited.at(i) = true;
        queue.push_back(i);
        break;
      }
    }
    while (!queue.empty()) {
      unsigned int node = queue.back();
      queue.pop_back();
      boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
      boost::graph_traits<Graph>::vertex_iterator i, end;
      for (tie(i, end) = vertices(g); i != end; ++i) {
        if (id[*i] == node) {
          break;
        }
      }
      for (tie(ei, edge_end) = out_edges(*i, g); ei != edge_end; ++ei) {
        unsigned int neighbour = id[target(*ei, g)];
        if (!visited.at(neighbour)) {
          visited.at(neighbour) = true;
          queue.push_back(neighbour);
          ++total_visited;
        }
      }
    }
  }
  return count;
}

unsigned long count_better(vector<double> *lower_bounds,
                           vector<double> *upper_bounds) {
  unsigned long total = 0;
  unsigned int length = lower_bounds->size();
  sort(lower_bounds->begin(), lower_bounds->end());
  sort(upper_bounds->begin(), upper_bounds->end());
  cout << lower_bounds->at(0) << lower_bounds->at(1) << endl;
  cout << upper_bounds->at(0) << upper_bounds->at(1) << endl;
  unsigned int lb_i = 0, ub_j = 0;
  while (ub_j != length && lb_i != length) {
    if (upper_bounds->at(ub_j) > lower_bounds->at(lb_i)) {
      ++lb_i;
    } else {
      total += length - lb_i;
      ++ub_j;
    }
  }
  return total;
}

unsigned long long count_better_sampled(const vector<double> &sampled_lower_bounds,
                                        const vector<double> &sorted_upper_bounds,
                                        bool strict) {
  unsigned long long total = 0ULL;
  if (sorted_upper_bounds.empty() || sampled_lower_bounds.empty()) return 0ULL;

  for (double lb : sampled_lower_bounds) {
    if (strict) {
      // counts ub < lb  (matches main.cpp's current saved logic using lower_bound)
      total += (unsigned long long)std::distance(
          sorted_upper_bounds.begin(),
          std::lower_bound(sorted_upper_bounds.begin(), sorted_upper_bounds.end(), lb));
    } else {
      // counts ub <= lb (matches count_better()'s >= behavior more closely)
      total += (unsigned long long)std::distance(
          sorted_upper_bounds.begin(),
          std::upper_bound(sorted_upper_bounds.begin(), sorted_upper_bounds.end(), lb));
    }
  }
  return total;
}

vector<list<pair<unsigned int, double>> *> *get_adjacency_list(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id,
    vector<vector<double> *> *dist) {

  vector<list<pair<unsigned int, double>> *> *adj_list =
      new vector<list<pair<unsigned int, double>> *>();
  boost::graph_traits<Graph>::vertex_iterator i, end;
  boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;

  for (tie(i, end) = vertices(g); i != end; ++i) {
    unsigned int end_point1 = id[*i];
    adj_list->push_back(new list<pair<unsigned int, double>>());
    for (tie(ei, edge_end) = out_edges(*i, g); ei != edge_end; ++ei) {
      unsigned int end_point2 = id[target(*ei, g)];
      if (dist == nullptr) {
        cout << "dist is null" << endl;
        exit(1);
      }
      if (end_point1 >= dist->size()) {
        cout << "end_point1 out of bounds: " << end_point1
             << " size: " << dist->size() << endl;
        exit(1);
      }
      if (dist->at(end_point1) == nullptr) {
        cout << "dist->at(" << end_point1 << ") is null" << endl;
        exit(1);
      }
      if (end_point2 >= dist->at(end_point1)->size()) {
        cout << "end_point2 out of bounds: " << end_point2
             << " size: " << dist->at(end_point1)->size() << endl;
        exit(1);
      }
      double distance = dist->at(end_point1)->at(end_point2);
      std::pair<unsigned int, double> pr = make_pair(end_point2, distance);
      adj_list->at(end_point1)->push_back(pr);
    }
  }
  return adj_list;
}

vector<vector<double> *> *get_adjacency_matrix(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id,
    vector<vector<double> *> *dist, double default_missing) {

  vector<vector<double> *> *adj_matrix = new vector<vector<double> *>();
  boost::graph_traits<Graph>::vertex_iterator i, end;
  boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
  unsigned int nodes = num_vertices(g);
  for (unsigned int i = 0; i < nodes; ++i) {
    adj_matrix->push_back(new vector<double>());
    for (unsigned int j = 0; j < nodes; ++j) {
      adj_matrix->at(i)->push_back(default_missing);
    }
    adj_matrix->at(i)->at(i) = 0.;
  }
  for (tie(i, end) = vertices(g); i != end; ++i) {
    unsigned int end_point1 = id[*i];
    for (tie(ei, edge_end) = out_edges(*i, g); ei != edge_end; ++ei) {
      unsigned int end_point2 = id[target(*ei, g)];
      double distance = dist->at(end_point1)->at(end_point2);
      adj_matrix->at(end_point1)->at(end_point2) = distance;
    }
  }
  return adj_matrix;
}

vector<vector<double> *> *distance_matrix(unsigned int nodes, unsigned int dims,
                                          unsigned p) {
  uniform_real_distribution<double> runif(-1., 1.);
  uniform_real_distribution<double> aunif(
      0., 2. * boost::math::constants::pi<double>());
  default_random_engine re;
  vector<vector<double> *> *dist = new vector<vector<double> *>();
  vector<vector<double> *> *points = new vector<vector<double> *>();
  double max_val = 0;
  for (unsigned int i = 0; i < nodes; ++i) {
    double r = runif(re);
    double entity = r;
    double angle = aunif(re);
    points->push_back(new vector<double>());
    for (unsigned int j = 1; j < dims; ++j) {
      points->at(i)->push_back(entity * cos(angle));
      entity *= sin(angle);
      angle = aunif(re);
    }
    points->at(i)->push_back(entity);
  }

  cout << "Allocated points" << endl;

  for (unsigned int i = 0; i < nodes; ++i) {
    dist->push_back(new vector<double>());
    for (unsigned int j = 0; j < nodes; ++j) {
      if (i == j) {
        dist->at(i)->push_back((double)0.);
      } else {
        double total = 0;
        for (unsigned int k = 0; k < dims; ++k) {
          double val = points->at(i)->at(k) - points->at(j)->at(k);
          if (val < 0)
            val *= -1.;
          total += pow(val, (double)p);
        }
        total = pow(total, (double)1. / p) / 2;
        max_val = max(max_val, total);
        dist->at(i)->push_back(total);
      }
    }
  }

  for (unsigned int i = 0; i < nodes; ++i) {
    delete points->at(i);
  }
  delete points;
  cout << "Maximum value encountered is " << max_val << endl;
  return dist;
}

void clean_up_adj_matrix(vector<vector<double> *> *adj_matrix) {
  unsigned int adj_matrix_size = adj_matrix->size();
  while (adj_matrix_size-- > 0) {
    delete adj_matrix->at(adj_matrix_size);
  }
  delete adj_matrix;
}

void clean_up_adj_list(vector<list<pair<unsigned int, double>> *> *adj_lst) {
  unsigned int adj_lst_size = adj_lst->size();
  while (adj_lst_size-- > 0) {
    delete adj_lst->at(adj_lst_size);
  }
  delete adj_lst;
}

map<pair<unsigned int, unsigned int>, double> *convert_adjList_to_knownEdges(
    vector<list<pair<unsigned int, double>> *> *adj_lst) {
  map<pair<unsigned int, unsigned int>, double> *known_edges =
      new map<pair<unsigned int, unsigned int>, double>();
  ofstream myfile;
  std::string file_name = "graph_" + to_string(adj_lst->size()) + ".txt";
  myfile.open(file_name);
  unsigned int ctr = 0;
  for (unsigned int u = 0; u < adj_lst->size(); ++u) {
    list<pair<unsigned int, double>> *lst = adj_lst->at(u);
    for (auto it = lst->begin(); it != lst->end(); ++it) {
      pair<unsigned int, double> p = *it;
      unsigned int v = p.first, u1 = u;
      if (v < u1) {
        v = u1;
        u1 = p.first;
      }
      std::stringstream u_;
      std::stringstream v_;
      std::stringstream p_;
      u_ << u1;
      v_ << v;
      p_ << p.second;
      std::string out_u = u_.str();
      std::string out_v = v_.str();
      std::string out_p = p_.str();
      std::string outline = out_u + " " + out_v + " " + out_p + "\n";
      myfile << outline;
      known_edges->insert(make_pair(make_pair(u1, v), p.second));
      ++ctr;
    }
  }
  myfile.close();
  cout << "Input (Graph)File is written; You can start Python code" << endl;

  return known_edges;
}

vector<list<pair<unsigned int, double>> *> *get_adj_list_file(char *filename) {
  vector<list<pair<unsigned int, double>> *> *adj_list =
      new vector<list<pair<unsigned int, double>> *>();
  ifstream ifs;
  ifs.open(filename, ifstream::in);
  unsigned int nodes;
  ifs >> nodes;
  for (int i = 0; i < nodes; i++) {
    adj_list->push_back(new list<pair<unsigned int, double>>());
  }
  unsigned int u, v, edge_numbers = 0;
  double dist;
  while (ifs >> u >> v >> dist) {
    adj_list->at(u)->push_back(make_pair(v, dist));
    adj_list->at(v)->push_back(make_pair(u, dist));
    ++edge_numbers;
    unsigned int u1 = min(u, v);
    unsigned int v1 = max(u, v);
    cout << u1 << " " << v1 << endl;
  }
  cout << "From read file" << edge_numbers << endl;
  return adj_list;
}

vector<vector<double> *> *get_adj_matrix_file(char *filename,
                                              double default_val) {
  vector<vector<double> *> *adj_mat = new vector<vector<double> *>();
  ifstream ifs;
  ifs.open(filename, ifstream::in);
  unsigned int nodes;
  ifs >> nodes;
  for (int i = 0; i < nodes; i++) {
    adj_mat->push_back(new vector<double>());
    for (int j = 0; j < nodes; j++) {
      adj_mat->at(i)->push_back(default_val);
    }
    adj_mat->at(i)->at(i) = 0.;
  }
  unsigned int u, v;
  double dist;
  while (ifs >> u >> v >> dist) {
    adj_mat->at(u)->at(v) = adj_mat->at(v)->at(u) = dist;
  }
  return adj_mat;
}

Graph waxman_graph_distances(unsigned int n, vector<vector<double> *> *distance,
                             double beta, double alpha) {
  cout << "Using WAXMAN based model " << endl;
  Graph g;
  vector<Traits::vertex_descriptor> node_desc;
  for (unsigned int i = 0; i < n; i++) {
    node_desc.push_back(boost::add_vertex(g));
  }

  double L = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      L = max(L, distance->at(i)->at(j));
    }
  }

  long total = 0;
  for (unsigned int i = 0; i < n - 1; i++) {
    for (unsigned int j = i + 1; j < n; j++) {
      // Use exp() for portability (std::numbers::e is C++20).
      double prob = beta * exp(-distance->at(i)->at(j) / (L * alpha));
      double sample_p = (double)rand() / RAND_MAX;
      if (sample_p <= prob) {
        cout << "Adding edge " << i << " " << j << " with distance "
             << distance->at(i)->at(j) << endl;
        boost::add_edge(node_desc.at(i), node_desc.at(j), g);
        boost::add_edge(node_desc.at(j), node_desc.at(i), g);
        ++total;
      }
    }
  }
  cout << "Total edges = " << total << endl;
  return g;
}

pair<Graph, vector<vector<double> *> *> waxman_graph(unsigned int n,
                                                     double beta, double alpha,
                                                     unsigned int p,
                                                     unsigned int dims) {

  cout << "Dimensions WAXMAN " << dims << endl;
  Graph g;
  vector<Traits::vertex_descriptor> node_desc;
  for (unsigned int i = 0; i < n; i++) {
    node_desc.push_back(boost::add_vertex(g));
  }

  vector<vector<double>> locations;
  for (unsigned int i = 0; i < n; i++) {
    vector<double> loc_i;
    for (unsigned int j = 0; j < dims; ++j) {
      double rn_gen = (double)rand() / RAND_MAX;
      loc_i.push_back(rn_gen);
    }
    locations.push_back(loc_i);
  }

  vector<vector<double> *> *distance = new vector<vector<double> *>();
  double L = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    vector<double> *new_distance = new vector<double>();
    for (unsigned int j = 0; j < n; j++) {
      if (i > j) {
        new_distance->push_back(distance->at(j)->at(i));
      } else if (i == j) {
        new_distance->push_back(0.);
      } else {
        double distance_ij = 0.0;
        for (unsigned int loc = 0; loc < dims; ++loc) {
          distance_ij +=
              abs(pow(locations.at(i).at(loc) - locations.at(j).at(loc), p));
        }
        distance_ij = pow(distance_ij, 1. / p);
        L = max(L, distance_ij);
        new_distance->push_back(distance_ij);
      }
    }
    distance->push_back(new_distance);
  }

  long total = 0;
  for (unsigned int i = 0; i < n - 1; i++) {
    for (unsigned int j = i + 1; j < n; j++) {
      // Use exp() for portability (std::numbers::e is C++20).
      double prob = beta * exp(-distance->at(i)->at(j) / (L * alpha));
      double sample_p = (double)rand() / RAND_MAX;
      if (sample_p <= prob) {
        boost::add_edge(node_desc.at(i), node_desc.at(j), g);
        boost::add_edge(node_desc.at(j), node_desc.at(i), g);
        ++total;
      }
    }
  }
  cout << "Total edges = " << total << endl;
  return make_pair(g, distance);
}
