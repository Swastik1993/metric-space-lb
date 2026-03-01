#include <time.h>
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
#include "KnnMetrics.h"

using namespace std;

bool EXACT;
const int DEBUG = 0;
const bool RENYI_ERDOS = true;
const bool WAXMAN = !RENYI_ERDOS;
const unsigned int RMSE_SAMPLES = 3000; // with repetition
const unsigned int ERC_SAMPLES  = 3000; // without repetition


static inline uint64_t cpu_time_us() {
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return uint64_t(ts.tv_sec) * 1000000ULL + ts.tv_nsec / 1000ULL;
}


static size_t get_peak_rss_kb() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_maxrss; // KB on Linux
}


inline pair<unsigned int, unsigned int>
sample_missing_edge(unsigned int nodes, map<pair<unsigned int, unsigned int>, double> *edges_map) {
  unsigned int u = rand() % nodes;
  unsigned int v = rand() % nodes;

  while (u == v ||
         edges_map->find(make_pair(min(u, v), max(u, v))) != edges_map->end()) {
    u = rand() % nodes;
    v = rand() % nodes;
  }
  return {u, v};
}


static double splub(vector<list<pair<unsigned int, double>> *> *adj_lst,
                    unsigned int a, unsigned int b) {
  double lb = 0.0;

  auto a_shortest_path = DijkstraELM(adj_lst, a);
  auto b_shortest_path = DijkstraELM(adj_lst, b);

  for (size_t i = 0; i < adj_lst->size(); ++i) {
    for (const auto& item : *adj_lst->at(i)) {
      unsigned j = item.first;
      double length = item.second;

      double sp_ai = a_shortest_path.first->at(i);
      double sp_aj = a_shortest_path.first->at(j);
      double sp_bi = b_shortest_path.first->at(i);
      double sp_bj = b_shortest_path.first->at(j);

      lb = std::max(lb, length - sp_ai - sp_bj);
      lb = std::max(lb, length - sp_aj - sp_bi);
    }
  }

  // free Dijkstra outputs
  delete a_shortest_path.first;
  if (a_shortest_path.second) delete a_shortest_path.second;
  delete b_shortest_path.first;
  if (b_shortest_path.second) delete b_shortest_path.second;

  return lb;
}



static void multiple_sample_exec(
  unsigned int init, 
  bool LARGE_GRAPH,
  vector<list<pair<unsigned int, double>> *> *adj_lst, 
  vector<vector<double> *> *adj_mat, 
  vector<vector<double> *> *lb, 
  vector<vector<double> *> *ub, 
  vector<vector<double> *> *lb_elm, 
  TriSearch &tri, 
  map<pair<unsigned int, unsigned int>, double> *known_edges, 
  vector<vector<double> *> *distance, 
  unsigned int nodes, 
  unsigned int k, 
  unsigned int sampling_size, 
  unsigned int knn_k, 
  unsigned int knn_queries, 
  unsigned int knn_pair_samples) {
    uint64_t start_lb_time = cpu_time_us();
    size_t mem0 = get_peak_rss_kb();
    ELM_SPTree *elm_sampling = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
    size_t sampling_mem = get_peak_rss_kb();
    cout << "Sampling edges" << endl;
    elm_sampling->get_sampling_landmarks();
    uint64_t stop_sampling_end = cpu_time_us();
    cout << "Preprocessing Sampled edges" << endl;
    elm_sampling->preprocess_landmark_roots();

    size_t exact_mem = 0;
    cout << "Exact edges" << endl;
    ELM_SPTree *elm_exact = nullptr;
    if (EXACT) {
      elm_exact = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
      exact_mem = get_peak_rss_kb();
      elm_exact->get_exact_landmarks();
    }
    uint64_t stop_exact_end = cpu_time_us();
    cout << "Preprocessing Exact edges" << endl;
    if (elm_exact) elm_exact->preprocess_landmark_roots();

    EdgeLandMark *elm_heuristic1 = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
    size_t h1_mem = get_peak_rss_kb();
    cout << "H1 edges" << endl;
    elm_heuristic1->large_edge_heuristic();
    uint64_t stop_h1_end = cpu_time_us();

    EdgeLandMark *elm_heuristic2 = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
    size_t h2_mem = get_peak_rss_kb();
    cout << "H2 edges" << endl;
    elm_heuristic2->far_away_heuristic();
    uint64_t stop_h2_end = cpu_time_us();

    cout << "ELM Sampling Time : " << (stop_sampling_end - start_lb_time) << endl;
    cout << "ELM Exact Time : " << (stop_exact_end - stop_sampling_end) << endl;
    cout << "ELM H1 Time : " << (stop_h1_end - stop_exact_end) << endl;
    cout << "ELM H2 Time : " << (stop_h2_end - stop_h1_end) << endl;
    
    cout << "PeakRSS(KB) after elm_sampling: " << sampling_mem << "\n";
    cout << "PeakRSS(KB) after elm_exact:    " << exact_mem << "\n";
    cout << "PeakRSS(KB) after elm_h1:       " << h1_mem << "\n";
    cout << "PeakRSS(KB) after elm_h2:       " << h2_mem << "\n";

    if (knn_queries > 0) {
      cout << "\n====================\n";
      cout << "[RUNNING KNN METRICS] (large_main)\n";
      cout << "====================\n";
      run_knn_metrics(
        adj_lst,
        tri,
        elm_exact,    
        elm_sampling, 
        elm_heuristic1,
        elm_heuristic2,
        nodes,
        knn_k,
        knn_queries,
        knn_pair_samples,
        init);
    }

    double total_lb_sw = 0.;
    double graph_weight_orig = 0.;
    uint64_t duration_lb_sw = 0;
    
    if (!LARGE_GRAPH) {
      uint64_t start_lb_sw = cpu_time_us();
      SashaWang(lb, ub);
      compute_lb(adj_lst, adj_mat, lb);
      uint64_t stop_lb_sw = cpu_time_us();
      duration_lb_sw = (stop_lb_sw - start_lb_sw);
    } else {
      cout << "Skipping SashaWang/compute_lb for LARGE_GRAPH." << endl;
    }

    if (!LARGE_GRAPH) {
    // REQUIRED: in LARGE_GRAPH branch you no longer hit the "else { graph_weight_orig += ... }"
    // so compute graph_weight_orig once here.
      for (unsigned int a = 0; a < nodes; ++a) {
        for (unsigned int b = a + 1; b < nodes; ++b) {
          if (adj_mat->at(a)->at(b) >= -0.1) {
            graph_weight_orig += lb->at(a)->at(b);
          }
        }
      }
    }

    cout << " =====>>> LB begins:   " << "\n";

    double relative = 0.;
    unsigned int relative_count = 0;
    double rmse_tri = 0.0;
    vector<double> rmse_exact, rmse_h1, rmse_h2, rmse_sampling;
    unsigned long long sw_saved = 0, tri_saved = 0;
    vector<unsigned long long> elm_exact_saved, elm_sampling_saved, elm_h1_saved, elm_h2_saved;
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

    if (!LARGE_GRAPH) {
      for (unsigned int i = 0; i < nodes; ++i) {
        for (unsigned int j = i + 1; j < nodes; ++j) {
          sorted_ub.push_back(ub->at(i)->at(j));
        }
      }
      sort(sorted_ub.begin(), sorted_ub.end());
    }
    
    uint64_t duration_tri_rt = 0;
    uint64_t duration_exact_rt = 0;
    uint64_t duration_elm_rt = 0;
    uint64_t duration_h1_rt = 0;
    uint64_t duration_h2_rt = 0;

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
              
              uint64_t start_tri_rt = cpu_time_us();
              double tri_val_ij = tri.lookup(i, j);
              uint64_t stop_tri_rt = cpu_time_us();
              duration_tri_rt += (stop_tri_rt - start_tri_rt);
              
              uint64_t start_exact_rt = cpu_time_us();
              dummy = elm_exact->lookup(i, j);
              uint64_t stop_exact_rt = cpu_time_us();
              duration_exact_rt += (stop_exact_rt - start_exact_rt);
              
              uint64_t start_elm_rt = cpu_time_us();
              dummy = elm_sampling->lookup(i, j);
              uint64_t stop_elm_rt = cpu_time_us();
              duration_elm_rt += (stop_elm_rt - start_elm_rt);
              
              // query time h1
              uint64_t start_h1_rt = cpu_time_us();
              dummy = elm_heuristic1->lookup(i, j);
              uint64_t stop_h1_rt = cpu_time_us();
              duration_h1_rt += (stop_h1_rt - start_h1_rt);

              // query time h2
              uint64_t start_h2_rt = cpu_time_us();
              dummy = elm_heuristic2->lookup(i, j);
              uint64_t stop_h2_rt = cpu_time_us();
              duration_h2_rt += (stop_h2_rt - start_h2_rt);


              rmse_tri += ((tri_val_ij - lb->at(i)->at(j)) *
                          (tri_val_ij - lb->at(i)->at(j)));
              tri_saved += std::distance(
                  sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                                tri_val_ij));
              sw_saved += std::distance(
                  sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                                lb->at(i)->at(j)));
              relative += 1 - (lb_elm->at(i)->at(j)) / (lb->at(i)->at(j));
              relative_count++;
              total_lb_sw += lb->at(i)->at(j);
            }
          // } else {
          //   graph_weight_orig += lb->at(i)->at(j);
          }
        }
      }
    } else {
      // =========================
      // LARGE GRAPH:
      // RMSE/Relative -> sampling WITH repetition 
      // Edge-resolution -> sampling WITHOUT repetition 
      // =========================

      double tri_lb = 0.0;

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
            rmse_exact.at(index_i) += (temp_exact->at(index_i) - lb_val) * (temp_exact->at(index_i) - lb_val);
          }
          rmse_sampling.at(index_i) += (temp_sampling->at(index_i) - lb_val) * (temp_sampling->at(index_i) - lb_val);
          rmse_h1.at(index_i) += (temp_h1->at(index_i) - lb_val) * (temp_h1->at(index_i) - lb_val);
          rmse_h2.at(index_i) += (temp_h2->at(index_i) - lb_val) * (temp_h2->at(index_i) - lb_val);
          // cout << "=====> RMSE Sampling, H1 and H2: " << rmse_sampling.at(index_i) << ", " << rmse_h1.at(index_i) << ", " << rmse_h2.at(index_i) << ", ";
        }

        if (EXACT) delete temp_exact;
        delete temp_h1;
        delete temp_h2;
        delete temp_sampling;

        // query time tri
        uint64_t start_tri_rt = cpu_time_us();
        tri_lb = tri.lookup(a, b);
        uint64_t stop_tri_rt = cpu_time_us();
        duration_tri_rt += (stop_tri_rt - start_tri_rt);
        rmse_tri += (tri_lb - lb_val) * (tri_lb - lb_val);

        relative += 1.0 - (lb_elm->at(a)->at(b)) / (lb_val);
        relative_count++;
        total_lb_sw += lb_val;

        // query time elm
        uint64_t start_elm_rt = cpu_time_us();
        double dummy = elm_sampling->lookup(a, b);
        uint64_t stop_elm_rt = cpu_time_us();
        duration_elm_rt += (stop_elm_rt - start_elm_rt);

        // query time h1
        uint64_t start_h1_rt = cpu_time_us();
        dummy = elm_heuristic1->lookup(a, b);
        uint64_t stop_h1_rt = cpu_time_us();
        duration_h1_rt += (stop_h1_rt - start_h1_rt);

        // query time h2
        uint64_t start_h2_rt = cpu_time_us();
        dummy = elm_heuristic2->lookup(a, b);
        uint64_t stop_h2_rt = cpu_time_us();
        duration_h2_rt += (stop_h2_rt - start_h2_rt);

      }

      // ---- Edge resolution WITHOUT repetition ----
      // set<pair<unsigned int, unsigned int>> used;
      for (int res_i = 0; res_i < ERC_SAMPLES; res_i++) {
        // auto [a, b] = sample_missing_edge(nodes, known_edges);
        
        unsigned int u1 = rand() % nodes;
        unsigned int v1 = rand() % nodes;
        while(u1 == v1) {
          v1 = rand() % nodes;
        }

        unsigned int u2 = rand() % nodes;
        unsigned int v2 = rand() % nodes;
        while(u2 == v2) {
          v2 = rand() % nodes;
        }

        while(u1==u2 && v1==v2) {
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

        // Fix Memory leak:
        delete u1v1_shortest_path.first;
        if (u1v1_shortest_path.second) delete u1v1_shortest_path.second;
        delete u2v2_shortest_path.first;
        if (u2v2_shortest_path.second) delete u2v2_shortest_path.second;

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

          // cout << "=====> ELM Sampling, H1 and H2: " << elm_sampling_saved.at(index_i) << ", " << elm_h1_saved.at(index_i) << ", " << elm_h2_saved.at(index_i) << ", ";
        }

        // query time tri
        uint64_t t0 = cpu_time_us();
        double tri_lb_u1v1 = tri.lookup(u1, v1);
        double tri_lb_u2v2 = tri.lookup(u2, v2);
        uint64_t t1 = cpu_time_us();
        duration_tri_rt += (t1 - t0);

        if (EXACT) delete temp_exact_u1v1;
        if (EXACT) delete temp_exact_u2v2;
        delete temp_h1_u1v1;
        delete temp_h1_u2v2;
        delete temp_h2_u1v1;
        delete temp_h2_u2v2;
        delete temp_sampling_u1v1;
        delete temp_sampling_u2v2;

        // tri_saved += std::distance(sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(), tri_lb));
        if (ub_u1v1 < tri_lb_u2v2 || ub_u2v2 < tri_lb_u1v1) {
          tri_saved += 1;
        }
      }
    }

    cout << "Query time tri : " << duration_tri_rt << endl;
    cout << "Query time exact : " << duration_exact_rt << endl;
    cout << "Query time sampling : " << duration_elm_rt << endl;
    cout << "Query time heuristic1: " << duration_h1_rt << endl;
    cout << "Query time heuristic2: " << duration_h2_rt << endl;

    cout << "Duration SW LB: "  << duration_lb_sw << endl; 
    cout << "Total Original Graph weight : " << graph_weight_orig << endl;
    cout << "Total SW LB : " << total_lb_sw << endl;
    cout << "Sum Relative Error on edge : " << relative / relative_count << endl;
    for (unsigned index_i = 0; index_i < k; ++index_i) {
      cout << "root mean square error sampling for " << index_i << ": " << sqrt(rmse_sampling.at(index_i) / relative_count) << endl;
      cout << "Saved Sampling " << index_i << ": " << elm_sampling_saved.at(index_i) << endl;
      if (EXACT) {
        cout << "root mean square error exact " << index_i << ": " << sqrt(rmse_exact.at(index_i) / relative_count) << endl;
        cout << "Saved Exact " << index_i << ": " << elm_exact_saved.at(index_i) << endl;
      }
      cout << "root mean square error H1 " << index_i << ": " << sqrt(rmse_h1.at(index_i) / relative_count) << endl;
      cout << "Saved H1 " << index_i << ": " << elm_h1_saved.at(index_i) << endl;
      cout << "root mean square error H2 " << index_i << ": " << sqrt(rmse_h2.at(index_i) / relative_count) << endl;
      cout << "Saved H2 " << index_i << ": " << elm_h2_saved.at(index_i) << endl;
    }
    cout << "root mean square error tri: " << sqrt(rmse_tri / relative_count) << endl;
    cout << "Saved TriSearch " << tri_saved << endl;
    cout << "Saved SW " << sw_saved << endl;
    cout << "Size of TriSearch " << tri._sizeof() << endl;
    cout << "Size of ELM " << elm_sampling->_sizeof() << endl;

    if (elm_exact) delete elm_exact;
    delete elm_heuristic1;
    delete elm_heuristic2;
    delete elm_sampling;
}






int main(int argc, char **argv) {
  unsigned int init = time(NULL);
  unsigned int nodes = 1000, k = 10, sampling_size = 2000;
  double prob = 0.05;

  cout << "This program can be run with following options" << endl;
  cout << "./edge_landmarks <random seed> <number of nodes> <edge density> <size of landmarks> <number of samples> <filename>" << endl;
  cout << "random seed - positive integer" << endl;
  cout << "number of nodes - positive integer" << endl;
  cout << "edge density - floating point number. Indicates the density of the graph." << endl;
  cout << "size of landmarks - positive integer" << endl;
  cout << "number of samples - positive integer" << endl;
  cout << "number of samples - positive integer" << endl;
  cout << "filename - filename to read the distance " << endl;
  
  if (argc > 1) init = stoi(argv[1]);
  if (argc > 2) nodes = stoi(argv[2]);
  if (argc > 3) prob = stof(argv[3]);
  if (argc > 4) k = stoi(argv[4]);
  if (argc > 5) sampling_size = stoi(argv[5]);

  const bool LARGE_GRAPH = (nodes >= 2000);
  EXACT = !(LARGE_GRAPH);

  unsigned int knn_queries = 50;          // 0 = disabled
  unsigned int knn_pair_samples = 2000;   // default
  unsigned int knn_k = 50;                // 0 => "use landmark_k" ONLY if you want fallback; otherwise set a real default like 50

  for (int ai = 6; ai < argc; ++ai) {
    if (string(argv[ai]) == "--knn") {
      // usage:
      //   --knn <num_queries> <pair_samples> [topk]
      if (ai + 2 >= argc) { cerr << "Bad --knn usage\n"; exit(2); }

      knn_queries = (unsigned)stoul(argv[ai + 1]);
      knn_pair_samples = (unsigned)stoul(argv[ai + 2]);

      // Optional 3rd arg: knn_k (Top-K)
      if (ai + 3 < argc) {
        // Only treat argv[ai+3] as knn_k if it does NOT look like a flag.
        string maybe = argv[ai + 3];
        if (!maybe.empty() && maybe[0] != '-') {
          knn_k = (unsigned)stoul(maybe);
          ai += 3;
        } else {
          ai += 2;
        }
      } else {
        ai += 2;
      }
    }
  }

  // Keep your original print, but extend it so you can SEE topk:
  cout << "knn_queries = " << knn_queries
      << ", knn_pair_samples = " << knn_pair_samples
      << ", knn_k = " << knn_k
      << endl;
  
  // Create graph with 100 nodes and edges with probability 0.05
  srand(init);
  
  vector<vector<double> *> *distance;
  Graph g;
  if (argc > 6) {
    // Read distance matrix
    distance = get_adj_matrix_file(argv[6], 0.);
    if (RENYI_ERDOS) {

      boost::minstd_rand gen(init);
      g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);
    } else {
      // Generate graph based on the distances
      g = waxman_graph_distances(nodes, distance, 0.4, 0.2);
    }
  } else {
    if (RENYI_ERDOS) {
      distance = distance_matrix(nodes, 6);

      boost::minstd_rand gen(init);
      g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);
    } else if (WAXMAN) {
      pair<Graph, vector<vector<double> *> *> graph_info = waxman_graph(nodes);
      g = graph_info.first;
      distance = graph_info.second;
    }
  }

  boost::property_map<Graph, unsigned int VertexProperties::*>::type id = get(&VertexProperties::index, g);
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

  vector<list<pair<unsigned int, double>> *> *adj_lst = get_adjacency_list(g, id, distance);
  vector<vector<double> *> *adj_mat; 
  vector<vector<double> *> *lb;
  vector<vector<double> *> *ub;
  vector<vector<double> *> *lb_elm = get_adjacency_matrix(g, id, distance, (double)0.);;

  if (!LARGE_GRAPH) {
    adj_mat = get_adjacency_matrix(g, id, distance, (double)-1.);
    lb = get_adjacency_matrix(g, id, distance, (double)0.);
    ub = get_adjacency_matrix(g, id, distance, (double)1.);
  }

  map<pair<unsigned int, unsigned int>, double> *known_edges =
      convert_adjList_to_knownEdges(adj_lst);

    cout << " =====>> Proceeding for Sampling " << endl;

  uint64_t start_tri_time = cpu_time_us();

  TriSearch tri = TriSearch(adj_lst, nodes);
  size_t tri_mem = get_peak_rss_kb();
  uint64_t stop_tri_end = cpu_time_us();

  cout << "TriSearch Time : " << (stop_tri_end - start_tri_time) << endl;
  cout << "PeakRSS(KB) after TriSearch:    " << tri_mem << "\n";

  // Declare and initialize an array of 4 samples
  int numbers[4] = {50, 100, 200, 300};
  // Loop from index 0 to 4 (5 elements)
  for (int i = 0; i < 4; i++) {
      sampling_size = numbers[i];
      srand(init);
      multiple_sample_exec(init, LARGE_GRAPH, 
        adj_lst, adj_mat, lb, ub, lb_elm, tri, known_edges, distance, 
        nodes, k, sampling_size, 
        knn_k, knn_queries, knn_pair_samples);
  }

  // multiple_sample_exec(init, LARGE_GRAPH, 
  //   adj_lst, adj_mat, lb, ub, lb_elm, tri, known_edges, distance, 
  //   nodes, k, sampling_size, 
  //   knn_k, knn_queries, knn_pair_samples);
  
  delete known_edges;
  clean_up_adj_list(adj_lst);
  clean_up_adj_matrix(lb_elm);
  if (!LARGE_GRAPH) {
    clean_up_adj_matrix(adj_mat);
    clean_up_adj_matrix(lb);
    clean_up_adj_matrix(ub);
  }

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



// g++ -std=c++23 -O3 -DNDEBUG -march=native -pthread *.cpp -o large_edge_landmarks_sampling


