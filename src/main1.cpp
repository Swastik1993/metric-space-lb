#include <random>
#include <sys/resource.h>
#include "Dijkstra.h"
#include "ELM_SPTree.h"
#include "EdgeLandMark.h"
#include "GraphDefinitions.h"
#include "GraphUtils.h"
#include "LowerBound.h"
#include "ShortestPathTree.h"
#include "TriSearch.h"
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>  // sort, lower_bound
#include <map>        // map
#include <set>        // set
#include <cmath>      // sqrt

using namespace std;

const int DEBUG = 0;
const bool EXACT = true;
const bool RENYI_ERDOS = true;
const bool WAXMAN = !RENYI_ERDOS;

static size_t get_peak_rss_kb() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_maxrss; // KB on Linux
}

inline pair<unsigned int, unsigned int>
sample_missing_edge(unsigned int nodes,
                    map<pair<unsigned int, unsigned int>, double> *edges_map) {
  unsigned int u = rand() % nodes;
  unsigned int v = rand() % nodes;

  while (u == v ||
         edges_map->find(make_pair(min(u, v), max(u, v))) != edges_map->end()) {
    u = rand() % nodes;
    v = rand() % nodes;
  }
  return {u, v};
}

static double splub(vector<list<pair<unsigned int, double>> *> *adj_lst, unsigned int a, unsigned int b) {
  double lb = 0;
  pair<vector<double>*, vector<double>*> a_shortest_path = DijkstraELM(adj_lst, a);
  pair<vector<double>*, vector<double>*> b_shortest_path = DijkstraELM(adj_lst, b); 

  for(int i =0; i < adj_lst->size(); ++i) {
    for (pair<unsigned int, double>& item : *adj_lst->at(i)) {
      int  j= item.first;
      double length = item.second;
      double sp_ai = a_shortest_path.first->at(i);
      double sp_aj = a_shortest_path.first->at(j);
      double sp_bi = b_shortest_path.first->at(i);
      double sp_bj = b_shortest_path.first->at(j);

      lb = max(lb, length - sp_ai, sp_bj);
      lb = max(lb, length - sp_aj, sp_bi);
    }
  }
  return lb; 
}

int main(int argc, char **argv) {
  unsigned int init = time(NULL);
  unsigned int nodes = 1000, k = 10, sampling_size = 2000;
  double prob = 0.05;

  cout << "This program can be run with following options" << endl;
  cout << "./edge_landmarks <random seed> <number of nodes> <edge density> "
          "<size of landmarks> <number of samples> <filename>"
       << endl;
  cout << "random seed - positive integer" << endl;
  cout << "number of nodes - positive integer" << endl;
  cout << "edge density - floating point number. Indicates the density of the "
          "graph."
       << endl;
  cout << "size of landmarks - positive integer" << endl;
  cout << "number of samples - positive integer" << endl;
  cout << "number of samples - positive integer" << endl;
  cout << "filename - filename to read the distance " << endl;
  if (argc > 1) {
    init = stoi(argv[1]);
  }
  if (argc > 2) {
    nodes = stoi(argv[2]);
  }
  if (argc > 3) {
    prob = stof(argv[3]);
  }
  if (argc > 4) {
    k = stoi(argv[4]);
  }
  if (argc > 5) {
    sampling_size = stoi(argv[5]);
  }

  const bool LARGE_GRAPH = (nodes >= 2000);
  const unsigned int RMSE_SAMPLES = sampling_size; // with repetition
  const unsigned int ERC_SAMPLES  = sampling_size; // without repetition

  // Create graph with 100 nodes and edges with probability 0.05
  srand(init);
  vector<vector<double> *> *distance;
  Graph g;
  if (argc > 6) {
    // Read distance matrix
    distance = get_adj_matrix_file(argv[6], 0.);
    if (RENYI_ERDOS) {

      boost::minstd_rand gen;
      g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);
    } else {
      // Generate graph based on the distances
      g = waxman_graph_distances(nodes, distance, 0.4, 0.2);
    }
  } else {
    if (RENYI_ERDOS) {
      distance = distance_matrix(nodes, 6);

      boost::minstd_rand gen;
      g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);
    } else if (WAXMAN) {
      pair<Graph, vector<vector<double> *> *> graph_info = waxman_graph(nodes);
      g = graph_info.first;
      distance = graph_info.second;
    }
  }

  boost::property_map<Graph, unsigned int VertexProperties::*>::type id =
      get(&VertexProperties::index, g);
  boost::graph_traits<Graph>::vertex_iterator vi, viend;

  int vert_num = 0;
  for (tie(vi, viend) = vertices(g); vi != viend; ++vi) {
    id[*vi] = vert_num++;
  }

  boost::graph_traits<Graph>::vertex_iterator i, end;
  boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;

  int total = 0;
  for (tie(i, end) = vertices(g); i != end; ++i) {
    if (DEBUG) {
      cout << id[*i] << " ";
    }
    int count = 0;
    for (tie(ei, edge_end) = out_edges(*i, g); ei != edge_end; ++ei) {
      if (DEBUG) {
        cout << id[target(*ei, g)] << "  ";
      }
      ++count;
    }
    total += count;
    if (DEBUG) {
      cout << count << endl;
    }
  }
  cout << "Total nodes : " << nodes << endl;
  cout << "Total edges : " << total / 2. << endl;
  cout << "Components " << check_connected(g, id) << endl;
  cout << "Adjacency matrix" << endl;

  vector<list<pair<unsigned int, double>> *> *adj_lst =
      get_adjacency_list(g, id, distance);
  vector<vector<double> *> *adj_mat =
      get_adjacency_matrix(g, id, distance, (double)-1.);
  vector<vector<double> *> *lb =
      get_adjacency_matrix(g, id, distance, (double)0.);
  vector<vector<double> *> *ub =
      get_adjacency_matrix(g, id, distance, (double)1.);
  vector<vector<double> *> *lb_elm =
      get_adjacency_matrix(g, id, distance, (double)0.);

  map<pair<unsigned int, unsigned int>, double> *known_edges =
      convert_adjList_to_knownEdges(adj_lst);

  auto start_lb_elm = chrono::high_resolution_clock::now();

  size_t mem0 = get_peak_rss_kb();
  EdgeLandMark *elm_sampling = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
  size_t mem1 = get_peak_rss_kb();
  // cout << "Sampling edges" << endl;
  elm_sampling->get_landmarks();
  auto stop_lb_elm_sampling = chrono::high_resolution_clock::now();

  EdgeLandMark *elm_exact = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
  size_t mem2 = get_peak_rss_kb();
  // cout << "Sampling edges" << endl;
  if (EXACT) {
    elm_exact->get_exact_landmarks();
  }
  auto stop_lb_elm_exact = chrono::high_resolution_clock::now();

  EdgeLandMark *elm_heuristic1 = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
  size_t mem3 = get_peak_rss_kb();
  // cout << "Sampling edges" << endl;
  elm_heuristic1->large_edge_heuristic();
  auto stop_lb_elm_h1 = chrono::high_resolution_clock::now();

  EdgeLandMark *elm_heuristic2 = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
  size_t mem4 = get_peak_rss_kb();
  // cout << "Sampling edges" << endl;
  elm_heuristic2->far_away_heuristic();
  auto stop_lb_elm_h2 = chrono::high_resolution_clock::now();

  TriSearch tri = TriSearch(adj_lst, nodes);
  size_t mem5 = get_peak_rss_kb();
  auto stop_lb_elm_ts_end = chrono::high_resolution_clock::now();
  cout << "ELM Sampling Time : "
       << chrono::duration_cast<chrono::microseconds>(stop_lb_elm_sampling -
                                                      start_lb_elm)
              .count()
       << endl;
  cout << "ELM Exact Time : "
       << chrono::duration_cast<chrono::microseconds>(stop_lb_elm_exact -
                                                      stop_lb_elm_sampling)
              .count()
       << endl;
  cout << "ELM H1 Time : "
       << chrono::duration_cast<chrono::microseconds>(stop_lb_elm_h1 -
                                                      stop_lb_elm_exact)
              .count()
       << endl;
  cout << "ELM H2 Time : "
       << chrono::duration_cast<chrono::microseconds>(stop_lb_elm_h2 -
                                                      stop_lb_elm_h1)
              .count()
       << endl;
  
  cout << "PeakRSS(KB) after elm_sampling: " << mem1 << "\n";
  cout << "PeakRSS(KB) after elm_exact:    " << mem2 << "\n";
  cout << "PeakRSS(KB) after elm_h1:       " << mem3 << "\n";
  cout << "PeakRSS(KB) after elm_h2:       " << mem4 << "\n";
  cout << "PeakRSS(KB) after TriSearch:    " << mem5 << "\n";

  double total_lb_sw = 0.;
  double graph_weight_orig = 0.;
  auto start_lb_sw = chrono::high_resolution_clock::now();
  SashaWang(lb, ub);
  compute_lb(adj_lst, adj_mat, lb);
  auto stop_lb_sw = chrono::high_resolution_clock::now();
  auto duration_lb_sw = chrono::duration_cast<chrono::microseconds>(stop_lb_sw - start_lb_sw);

  // REQUIRED: in LARGE_GRAPH branch you no longer hit the "else { graph_weight_orig += ... }"
  // so compute graph_weight_orig once here.
  for (unsigned int a = 0; a < nodes; ++a) {
    for (unsigned int b = a + 1; b < nodes; ++b) {
      if (adj_mat->at(a)->at(b) >= -0.1) {
        graph_weight_orig += lb->at(a)->at(b);
      }
    }
  }

  double relative = 0.;
  unsigned int relative_count = 0;
  double rmse_tri = 0.0;
  vector<double> rmse_exact, rmse_h1, rmse_h2, rmse_sampling;
  unsigned long long sw_saved = 0, tri_saved = 0;
  vector<unsigned long long> elm_exact_saved, elm_sampling_saved, elm_h1_saved,
      elm_h2_saved;
  for (unsigned int i = 0; i < k; ++i) {
    rmse_exact.push_back(0.0);
    rmse_h1.push_back(0.0);
    rmse_h2.push_back(0.0);
    rmse_sampling.push_back(0.0);
    elm_exact_saved.push_back(0);
    elm_sampling_saved.push_back(0);
    elm_h1_saved.push_back(0);
    elm_h2_saved.push_back(0);
  }
  vector<double> sorted_ub;
  for (unsigned int i = 0; i < nodes; ++i) {
    for (unsigned int j = i + 1; j < nodes; ++j) {
      sorted_ub.push_back(ub->at(i)->at(j));
    }
  }
  sort(sorted_ub.begin(), sorted_ub.end());
  std::chrono::seconds dummy_s(0);
  auto duration_tri_rt = chrono::duration_cast<chrono::microseconds>(dummy_s);
  auto duration_elm_rt = chrono::duration_cast<chrono::microseconds>(dummy_s);

  int uv_count = 0;

  if (!LARGE_GRAPH) {
    for (unsigned int i = 0; i < nodes; ++i) {
      for (unsigned int j = i + 1; j < nodes; ++j) {
        if (adj_mat->at(i)->at(j) < -0.1) {
          if (lb->at(i)->at(j) > 0) {
            vector<double> *temp_exact;
            if (EXACT) {
              temp_exact = elm_exact->lookup_multiple(i, j);
            }
            vector<double> *temp_sampling = elm_sampling->lookup_multiple(i, j);
            vector<double> *temp_h1 = elm_heuristic1->lookup_multiple(i, j);
            vector<double> *temp_h2 = elm_heuristic2->lookup_multiple(i, j);
            for (unsigned index_i = 0; index_i < k; ++index_i) {
              if (EXACT) {
                rmse_exact.at(index_i) +=
                    ((temp_exact->at(index_i) - lb->at(i)->at(j)) *
                    (temp_exact->at(index_i) - lb->at(i)->at(j)));
                elm_exact_saved.at(index_i) +=
                    std::distance(sorted_ub.begin(),
                                  lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                              temp_exact->at(index_i)));
              }
              rmse_sampling.at(index_i) +=
                  (temp_sampling->at(index_i) - lb->at(i)->at(j)) *
                  (temp_sampling->at(index_i) - lb->at(i)->at(j));
              elm_sampling_saved.at(index_i) +=
                  std::distance(sorted_ub.begin(),
                                lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                            temp_sampling->at(index_i)));
              rmse_h1.at(index_i) += ((temp_h1->at(index_i) - lb->at(i)->at(j)) *
                                      (temp_h1->at(index_i) - lb->at(i)->at(j)));
              elm_h1_saved.at(index_i) +=
                  std::distance(sorted_ub.begin(),
                                lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                            temp_h1->at(index_i)));
              rmse_h2.at(index_i) += ((temp_h2->at(index_i) - lb->at(i)->at(j)) *
                                      (temp_h2->at(index_i) - lb->at(i)->at(j)));
              elm_h2_saved.at(index_i) +=
                  std::distance(sorted_ub.begin(),
                                lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                            temp_h2->at(index_i)));
            }
            if (EXACT) {
              delete temp_exact;
            }
            delete temp_h1;
            delete temp_h2;
            delete temp_sampling;
            double dummy = 0.;
            auto start_tri_rt = chrono::high_resolution_clock::now();
            dummy = tri.lookup(i, j);
            auto stop_tri_rt = chrono::high_resolution_clock::now();
            duration_tri_rt += chrono::duration_cast<chrono::microseconds>(
                stop_tri_rt - start_tri_rt);
            auto start_elm_rt = chrono::high_resolution_clock::now();
            dummy = elm_sampling->lookup(i, j);
            auto stop_elm_rt = chrono::high_resolution_clock::now();
            duration_elm_rt += chrono::duration_cast<chrono::microseconds>(
                stop_elm_rt - start_elm_rt);
            rmse_tri += ((tri.lookup(i, j) - lb->at(i)->at(j)) *
                        (tri.lookup(i, j) - lb->at(i)->at(j)));
            tri_saved += std::distance(
                sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                              tri.lookup(i, j)));
            sw_saved += std::distance(
                sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                              lb->at(i)->at(j)));
            relative += 1 - (lb_elm->at(i)->at(j)) / (lb->at(i)->at(j));
            relative_count++;
            total_lb_sw += lb->at(i)->at(j);
          }
        } else {
          graph_weight_orig += lb->at(i)->at(j);
        }
      }
    }
  } else {
    // =========================
    // LARGE GRAPH:
    // RMSE/Relative -> sampling WITH repetition using sample_missing_edge
    // Edge-resolution -> sampling WITHOUT repetition using set + sample_missing_edge
    // =========================

    // ---- RMSE / Relative WITH repetition ----
    for (unsigned int s = 0; s < RMSE_SAMPLES; ++s) {
      auto [a, b] = sample_missing_edge(nodes, known_edges);
      double lb_val = splub(adj_lst, a, b);

      // Keep sample budget consistent: if lb invalid, resample
      if (lb_val <= 0) { --s; continue; }

      vector<double> *temp_exact;
      if (EXACT) temp_exact = elm_exact->lookup_multiple(a, b);
      vector<double> *temp_sampling = elm_sampling->lookup_multiple(a, b);
      vector<double> *temp_h1 = elm_heuristic1->lookup_multiple(a, b);
      vector<double> *temp_h2 = elm_heuristic2->lookup_multiple(a, b);

      for (unsigned index_i = 0; index_i < k; ++index_i) {
        if (EXACT) {
          rmse_exact.at(index_i) +=
              (temp_exact->at(index_i) - lb_val) *
              (temp_exact->at(index_i) - lb_val);
        }
        rmse_sampling.at(index_i) +=
            (temp_sampling->at(index_i) - lb_val) *
            (temp_sampling->at(index_i) - lb_val);
        rmse_h1.at(index_i) +=
            (temp_h1->at(index_i) - lb_val) *
            (temp_h1->at(index_i) - lb_val);
        rmse_h2.at(index_i) +=
            (temp_h2->at(index_i) - lb_val) *
            (temp_h2->at(index_i) - lb_val);
      }

      if (EXACT) delete temp_exact;
      delete temp_h1;
      delete temp_h2;
      delete temp_sampling;

      double tri_lb = tri.lookup(a, b);
      rmse_tri += (tri_lb - lb_val) * (tri_lb - lb_val);

      relative += 1.0 - (lb_elm->at(a)->at(b)) / (lb_val);
      relative_count++;
      total_lb_sw += lb_val;

      // query time tri
      auto start_tri_rt = chrono::high_resolution_clock::now();
      double tri_lb = tri.lookup(a, b);
      auto stop_tri_rt = chrono::high_resolution_clock::now();
      duration_tri_rt += chrono::duration_cast<chrono::microseconds>(stop_tri_rt - start_tri_rt);

      // query time elm
      auto start_elm_rt = chrono::high_resolution_clock::now();
      double dummy = elm_sampling->lookup(a, b);
      auto stop_elm_rt = chrono::high_resolution_clock::now();
      duration_elm_rt += chrono::duration_cast<chrono::microseconds>(stop_elm_rt - start_elm_rt);

    }

    // ---- Edge resolution WITHOUT repetition ----
    // set<pair<unsigned int, unsigned int>> used;
    for (int res_i = 0; res_i < ERC_SAMPLES; res_i++) {
      // auto [a, b] = sample_missing_edge(nodes, known_edges);
      
      unsigned int u1 = rand() % nodes;
      unsigned int v1 = rand() % nodes;
      while(u1 == v1) {
        unsigned int v1 = rand() % nodes;
      }

      unsigned int u2 = rand() % nodes;
      unsigned int v2 = rand() % nodes;
      while(u2 == v2) {
        unsigned int v2 = rand() % nodes;
      }

      while(u1==u2 & v1==v2) {
        u2 = rand() % nodes;
        v2 = rand() % nodes;
          while(u2 == v2) {
            v2 = rand() % nodes;
        }
      }

      auto key_u1v1 = make_pair(min(u1, v1), max(u1, v1));
      auto key_u2v2 = make_pair(min(u2, v2), max(u2, v2));
      // if (!used.insert(key).second) continue;
      // if (lb->at(a)->at(b) <= 0) continue;

      u1 = key_u1v1.first;
      u2 = key_u2v2.first;
      v1 = key_u1v1.second;
      v2 = key_u2v2.second;

      vector<double> *temp_exact_u1v1;
      if (EXACT) temp_exact_u1v1 = elm_exact->lookup_multiple(u1, v1);
      vector<double> *temp_sampling_u1v1 = elm_sampling->lookup_multiple(u1, v1);
      vector<double> *temp_h1_u1v1 = elm_heuristic1->lookup_multiple(u1, v1);
      vector<double> *temp_h2_u1v1 = elm_heuristic2->lookup_multiple(u1, v1);

      vector<double> *temp_exact_u2v2;
      if (EXACT) temp_exact_u2v2 = elm_exact->lookup_multiple(u2, v2);
      vector<double> *temp_sampling_u2v2 = elm_sampling->lookup_multiple(u2, v2);
      vector<double> *temp_h1_u2v2 = elm_heuristic1->lookup_multiple(u2, v2);
      vector<double> *temp_h2_u2v2 = elm_heuristic2->lookup_multiple(u2, v2);

      pair<vector<double>*, vector<double>*> u1v1_shortest_path = DijkstraELM(adj_lst, u1);
      pair<vector<double>*, vector<double>*> u2v2_shortest_path = DijkstraELM(adj_lst, u2); 

      double ub_u1v1 = u1v1_shortest_path.first->at(v1);
      double ub_u2v2 = u2v2_shortest_path.first->at(v2);

      for (unsigned index_i = 0; index_i < k; ++index_i) {
        if (EXACT) {
          if (ub_u1v1 < temp_exact_u2v2->at(index_i) || ub_u2v2 < temp_exact_u1v1->at(index_i)) {
            elm_exact_saved.at(index_i) += 1;
          }
        }
        if (ub_u1v1 < temp_sampling_u2v2->at(index_i) || ub_u2v2 < temp_sampling_u1v1->at(index_i)) {
          elm_sampling_saved.at(index_i) += 1;
        }
        if (ub_u1v1 < temp_h1_u2v2->at(index_i) || ub_u2v2 < temp_h1_u1v1->at(index_i)) {
          elm_h1_saved.at(index_i) += 1;
        }
        if (ub_u1v1 < temp_h2_u2v2->at(index_i) || ub_u2v2 < temp_h2_u1v1->at(index_i)) {
          elm_h2_saved.at(index_i) += 1;
        }
      }

      if (EXACT) delete temp_exact_u1v1;
      if (EXACT) delete temp_exact_u2v2;
      delete temp_h1_u1v1;
      delete temp_h1_u2v2;
      delete temp_h2_u1v1;
      delete temp_h2_u2v2;
      delete temp_sampling_u1v1;
      delete temp_sampling_u2v2;

      tri_saved += std::distance(
          sorted_ub.begin(),
          lower_bound(sorted_ub.begin(), sorted_ub.end(), tri_lb));
    }
  }

  cout << "Query time elb : " << duration_elm_rt.count() << endl;
  cout << "Query time tri : " << duration_tri_rt.count() << endl;

  std::cout << "Duration SW LB:" << duration_lb_sw.count() / 1000000.0 << endl;
  cout << "Total Original Graph weight : " << graph_weight_orig << endl;
  cout << "Total SW LB : " << total_lb_sw << endl;
  cout << "Sum Relative Error on edge : " << relative / relative_count << endl;
  for (unsigned index_i = 0; index_i < k; ++index_i) {
    cout << "root mean square error sampling for " << index_i << ": "
         << sqrt(rmse_sampling.at(index_i) / relative_count) << endl;
    cout << "Saved Sampling " << index_i << ": "
         << elm_sampling_saved.at(index_i) << endl;
    if (EXACT) {
      cout << "root mean square error exact " << index_i << ": "
           << sqrt(rmse_exact.at(index_i) / relative_count) << endl;
      cout << "Saved Exact " << index_i << ": " << elm_exact_saved.at(index_i)
           << endl;
    }
    cout << "root mean square error H1 " << index_i << ": "
         << sqrt(rmse_h1.at(index_i) / relative_count) << endl;
    cout << "Saved H1 " << index_i << ": " << elm_h1_saved.at(index_i) << endl;
    cout << "root mean square error H2 " << index_i << ": "
         << sqrt(rmse_h2.at(index_i) / relative_count) << endl;
    cout << "Saved H2 " << index_i << ": " << elm_h2_saved.at(index_i) << endl;
  }
  cout << "root mean square error tri: " << sqrt(rmse_tri / relative_count)
       << endl;
  cout << "Saved TriSearch " << tri_saved << endl;
  cout << "Saved SW " << sw_saved << endl;
  cout << "Size of TriSearch " << tri._sizeof() << std::endl;
  cout << "Size of ELM " << elm_sampling->_sizeof() << std::endl;

  delete elm_exact;
  delete elm_heuristic1;
  delete elm_heuristic2;
  delete elm_sampling;
  delete known_edges;
  clean_up_adj_list(adj_lst);
  clean_up_adj_matrix(adj_mat);
  clean_up_adj_matrix(lb);
  clean_up_adj_matrix(ub);
  clean_up_adj_matrix(lb_elm);
  for (unsigned int i = 0; i < nodes; ++i) {
    if (DEBUG) {
      for (unsigned int j = 0; j < nodes; ++j) {
        cout << " " << distance->at(i)->at(j);
      }
      cout << endl;
    }
    delete distance->at(i);
  }
  delete distance;
  return 0;
}
