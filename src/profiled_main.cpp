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
#include <iomanip>    // fixed, setprecision
#include <sstream>

using namespace std;

const int DEBUG = 0;
const bool EXACT = false;
const bool RENYI_ERDOS = true;
const bool WAXMAN = !RENYI_ERDOS;

static size_t get_peak_rss_kb() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_maxrss; // KB on Linux
}

// =========================
// Profiling helpers
// =========================
using Clock = std::chrono::steady_clock;

static inline double seconds_since(const Clock::time_point& t0) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - t0).count();
}

struct ScopeTimer {
  const char* label;
  Clock::time_point t0;
  ScopeTimer(const char* lbl) : label(lbl), t0(Clock::now()) {}
  ~ScopeTimer() {
    double s = seconds_since(t0);
    std::cout << "[T] " << label << " took " << std::fixed << std::setprecision(6) << s << " s\n";
  }
};

#define SCOPE_TIMER(lbl) ScopeTimer _scope_timer_##__LINE__(lbl)

static inline void checkpoint(const char* lbl, const Clock::time_point& t0_global) {
  std::cout << "[T] checkpoint: " << lbl << " @ +" << std::fixed << std::setprecision(6)
            << seconds_since(t0_global) << " s\n";
}

// Optional: accumulate microseconds into a counter
static inline void add_us(std::chrono::microseconds& acc, const Clock::time_point& t0) {
  acc += std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - t0);
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
  // NOTE: splub is expensive (2x Dijkstra + scan all edges). We time it, but avoid spamming output.
  static int splub_prints = 0;
  const bool do_print = (splub_prints++ < 5);

  auto t0 = Clock::now();

  double lb = 0.0;

  auto t_dij = Clock::now();
  pair<vector<double>*, vector<double>*> a_shortest_path = DijkstraELM(adj_lst, a);
  pair<vector<double>*, vector<double>*> b_shortest_path = DijkstraELM(adj_lst, b);
  if (do_print) {
    std::cout << "[T] splub: DijkstraELM pair took " << std::fixed << std::setprecision(6)
              << seconds_since(t_dij) << " s\n";
  }

  auto t_scan = Clock::now();
  for (size_t i = 0; i < adj_lst->size(); ++i) {
    for (const auto& item : *adj_lst->at(i)) {
      int  j = item.first;
      double length = item.second;
      double sp_ai = a_shortest_path.first->at(i);
      double sp_aj = a_shortest_path.first->at(j);
      double sp_bi = b_shortest_path.first->at(i);
      double sp_bj = b_shortest_path.first->at(j);

      lb = std::max(lb, length - sp_ai - sp_bj);
      lb = std::max(lb, length - sp_aj - sp_bi);
    }
  }
  if (do_print) {
    std::cout << "[T] splub: scan all edges took " << std::fixed << std::setprecision(6)
              << seconds_since(t_scan) << " s\n";
    std::cout << "[T] splub total took " << std::fixed << std::setprecision(6)
              << seconds_since(t0) << " s\n";
  }

  // NOTE: DijkstraELM likely allocates vectors; original code didn't free them.
  // We do NOT change logic/ownership assumptions here.

  return lb;
}

int main(int argc, char **argv) {
  auto T_GLOBAL = Clock::now();
  SCOPE_TIMER("TOTAL main()");

  unsigned int init = time(NULL);
  unsigned int nodes = 1000, k = 10, sampling_size = 2000;
  double prob = 0.05;

  {
    SCOPE_TIMER("parse args + usage prints");
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

    if (argc > 1) init = stoi(argv[1]);
    if (argc > 2) nodes = stoi(argv[2]);
    if (argc > 3) prob = stof(argv[3]);
    if (argc > 4) k = stoi(argv[4]);
    if (argc > 5) sampling_size = stoi(argv[5]);
  }
  checkpoint("after args", T_GLOBAL);

  const bool LARGE_GRAPH = (nodes >= 2000);
  const unsigned int RMSE_SAMPLES = 10000; sampling_size; // with repetition
  const unsigned int ERC_SAMPLES  = 10000; sampling_size; // without repetition

  // Create graph with nodes and edges with probability prob
  srand(init);

  vector<vector<double> *> *distance = nullptr;
  Graph g;

  {
    SCOPE_TIMER("graph + distance build");
    if (argc > 6) {
      {
        SCOPE_TIMER("read distance matrix from file");
        distance = get_adj_matrix_file(argv[6], 0.);
      }
      if (RENYI_ERDOS) {
        SCOPE_TIMER("generate Erdos-Renyi graph (file distance present)");
        boost::minstd_rand gen;
        g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);
      } else {
        SCOPE_TIMER("generate Waxman graph from distances");
        g = waxman_graph_distances(nodes, distance, 0.4, 0.2);
      }
    } else {
      if (RENYI_ERDOS) {
        {
          SCOPE_TIMER("compute distance matrix (synthetic)");
          distance = distance_matrix(nodes, 6);
        }
        {
          SCOPE_TIMER("generate Erdos-Renyi graph (synthetic distance)");
          boost::minstd_rand gen;
          g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);
        }
      } else if (WAXMAN) {
        SCOPE_TIMER("generate Waxman graph (with generated distances)");
        pair<Graph, vector<vector<double> *> *> graph_info = waxman_graph(nodes);
        g = graph_info.first;
        distance = graph_info.second;
      }
    }
  }
  checkpoint("after graph build", T_GLOBAL);

  boost::property_map<Graph, unsigned int VertexProperties::*>::type id =
      get(&VertexProperties::index, g);

  {
    SCOPE_TIMER("assign vertex ids");
    boost::graph_traits<Graph>::vertex_iterator vi, viend;
    int vert_num = 0;
    for (tie(vi, viend) = vertices(g); vi != viend; ++vi) {
      id[*vi] = vert_num++;
    }
  }
  checkpoint("after vertex ids", T_GLOBAL);

  {
    SCOPE_TIMER("count edges");
    boost::graph_traits<Graph>::vertex_iterator i, end;
    boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;

    int total = 0;
    for (tie(i, end) = vertices(g); i != end; ++i) {
      if (DEBUG) cout << id[*i] << " ";
      int count = 0;
      for (tie(ei, edge_end) = out_edges(*i, g); ei != edge_end; ++ei) {
        if (DEBUG) cout << id[target(*ei, g)] << "  ";
        ++count;
      }
      total += count;
      if (DEBUG) cout << count << endl;
    }

    cout << "Total nodes : " << nodes << endl;
    cout << "Total edges : " << total / 2. << endl;

    {
      SCOPE_TIMER("check_connected()");
      cout << "Components " << check_connected(g, id) << endl;
    }
  }
  checkpoint("after graph stats", T_GLOBAL);

  cout << "Adjacency matrix" << endl;

  vector<list<pair<unsigned int, double>> *> *adj_lst = nullptr;
  vector<vector<double> *> *adj_mat = nullptr;
  vector<vector<double> *> *lb = nullptr;
  vector<vector<double> *> *ub = nullptr;
  vector<vector<double> *> *lb_elm = nullptr;

  {
    SCOPE_TIMER("build adjacency list/matrices");
    {
      SCOPE_TIMER("get_adjacency_list()");
      adj_lst = get_adjacency_list(g, id, distance);
    }
    {
      SCOPE_TIMER("get_adjacency_matrix() adj_mat");
      adj_mat = get_adjacency_matrix(g, id, distance, (double)-1.);
    }
    {
      SCOPE_TIMER("get_adjacency_matrix() lb");
      lb = get_adjacency_matrix(g, id, distance, (double)0.);
    }
    {
      SCOPE_TIMER("get_adjacency_matrix() ub");
      ub = get_adjacency_matrix(g, id, distance, (double)1.);
    }
    {
      SCOPE_TIMER("get_adjacency_matrix() lb_elm");
      lb_elm = get_adjacency_matrix(g, id, distance, (double)0.);
    }
  }
  checkpoint("after adjacency structures", T_GLOBAL);

  map<pair<unsigned int, unsigned int>, double> *known_edges = nullptr;
  {
    SCOPE_TIMER("convert_adjList_to_knownEdges()");
    known_edges = convert_adjList_to_knownEdges(adj_lst);
  }
  checkpoint("after known_edges map", T_GLOBAL);

  cout << " =====>> Proceeding for Sampling " << endl;

  // Memory checkpoints (peak RSS)
  size_t mem0 = get_peak_rss_kb();

  ELM_SPTree *elm_sampling = nullptr;
  ELM_SPTree *elm_exact = nullptr;
  EdgeLandMark *elm_heuristic1 = nullptr;
  EdgeLandMark *elm_heuristic2 = nullptr;

  {
    SCOPE_TIMER("ELM_SPTree ctor (sampling)");
    elm_sampling = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
  }
  size_t mem1 = get_peak_rss_kb();
  {
    SCOPE_TIMER("elm_sampling->get_sampling_landmarks()");
    elm_sampling->get_sampling_landmarks();
  }
  checkpoint("after elm_sampling", T_GLOBAL);

  {
    SCOPE_TIMER("ELM_SPTree ctor (exact)");
    elm_exact = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
  }
  size_t mem2 = get_peak_rss_kb();
  if (EXACT) {
    SCOPE_TIMER("elm_exact->get_exact_landmarks()");
    elm_exact->get_exact_landmarks();
  }
  checkpoint("after elm_exact", T_GLOBAL);

  {
    SCOPE_TIMER("EdgeLandMark ctor (H1)");
    elm_heuristic1 = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
  }
  size_t mem3 = get_peak_rss_kb();
  {
    SCOPE_TIMER("elm_heuristic1->large_edge_heuristic()");
    elm_heuristic1->large_edge_heuristic();
  }
  checkpoint("after elm_h1", T_GLOBAL);

  {
    SCOPE_TIMER("EdgeLandMark ctor (H2)");
    elm_heuristic2 = new EdgeLandMark(adj_lst, nodes, k, sampling_size);
  }
  size_t mem4 = get_peak_rss_kb();
  {
    SCOPE_TIMER("elm_heuristic2->far_away_heuristic()");
    elm_heuristic2->far_away_heuristic();
  }
  checkpoint("after elm_h2", T_GLOBAL);

  TriSearch tri(adj_lst, nodes);
  size_t mem5 = get_peak_rss_kb();
  checkpoint("after TriSearch ctor", T_GLOBAL);

  cout << "PeakRSS(KB) after baseline (mem0): " << mem0 << "\n";
  cout << "PeakRSS(KB) after elm_sampling:    " << mem1 << "\n";
  cout << "PeakRSS(KB) after elm_exact:       " << mem2 << "\n";
  cout << "PeakRSS(KB) after elm_h1:          " << mem3 << "\n";
  cout << "PeakRSS(KB) after elm_h2:          " << mem4 << "\n";
  cout << "PeakRSS(KB) after TriSearch:       " << mem5 << "\n";

  double total_lb_sw = 0.;
  double graph_weight_orig = 0.;

  std::chrono::microseconds duration_lb_sw_us(0);

  {
    SCOPE_TIMER("SashaWang/compute_lb block");
    auto t0 = Clock::now();
    if (!LARGE_GRAPH) {
      {
        SCOPE_TIMER("SashaWang(lb, ub)");
        SashaWang(lb, ub);
      }
      {
        SCOPE_TIMER("compute_lb(adj_lst, adj_mat, lb)");
        compute_lb(adj_lst, adj_mat, lb);
      }
    } else {
      cout << "Skipping SashaWang/compute_lb for LARGE_GRAPH." << endl;
    }
    duration_lb_sw_us = std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - t0);
  }
  checkpoint("after SW/compute_lb", T_GLOBAL);

  cout << " =====>>> LB begins:   " << "\n";

  double relative = 0.;
  unsigned int relative_count = 0;
  double rmse_tri = 0.0;

  vector<double> rmse_exact, rmse_h1, rmse_h2, rmse_sampling;
  unsigned long long sw_saved = 0, tri_saved = 0;
  vector<unsigned long long> elm_exact_saved, elm_sampling_saved, elm_h1_saved, elm_h2_saved;

  {
    SCOPE_TIMER("init metric vectors");
    rmse_exact.assign(k, 0.0);
    rmse_h1.assign(k, 0.0);
    rmse_h2.assign(k, 0.0);
    rmse_sampling.assign(k, 0.0);

    elm_exact_saved.assign(k, 0);
    elm_sampling_saved.assign(k, 0);
    elm_h1_saved.assign(k, 0);
    elm_h2_saved.assign(k, 0);
  }

  vector<double> sorted_ub;
  {
    SCOPE_TIMER("build sorted_ub (upper triangle)");
    sorted_ub.reserve((size_t)nodes * (nodes - 1) / 2);
    for (unsigned int i = 0; i < nodes; ++i) {
      for (unsigned int j = i + 1; j < nodes; ++j) {
        sorted_ub.push_back(ub->at(i)->at(j));
      }
    }
  }
  {
    SCOPE_TIMER("sort(sorted_ub)");
    sort(sorted_ub.begin(), sorted_ub.end());
  }

  // Query-time accumulators
  std::chrono::microseconds duration_tri_rt(0);
  std::chrono::microseconds duration_elm_rt(0);

  // Additional accumulators for lookup_multiple calls (optional)
  std::chrono::microseconds dur_lookup_sampling(0), dur_lookup_h1(0), dur_lookup_h2(0), dur_lookup_exact(0);
  std::chrono::microseconds dur_tri_lookup_only(0);

  if (!LARGE_GRAPH) {
    SCOPE_TIMER("SMALL_GRAPH evaluation (full pair loops)");
    auto t_pairs = Clock::now();
    unsigned long long pair_counter = 0;

    for (unsigned int i = 0; i < nodes; ++i) {
      for (unsigned int j = i + 1; j < nodes; ++j) {
        ++pair_counter;

        // periodic progress (avoid flooding)
        if ((pair_counter % 50000ULL) == 0ULL) {
          cout << "[T] small-graph pair loop progress pairs=" << pair_counter
               << " elapsed=" << std::fixed << std::setprecision(6)
               << seconds_since(t_pairs) << " s\n";
        }

        if (adj_mat->at(i)->at(j) < -0.1) {
          if (lb->at(i)->at(j) > 0) {
            vector<double> *temp_exact = nullptr;

            if (EXACT) {
              auto t0 = Clock::now();
              temp_exact = elm_exact->lookup_multiple(i, j);
              add_us(dur_lookup_exact, t0);
            }

            vector<double> *temp_sampling = nullptr;
            {
              auto t0 = Clock::now();
              temp_sampling = elm_sampling->lookup_multiple(i, j);
              add_us(dur_lookup_sampling, t0);
            }

            vector<double> *temp_h1 = nullptr;
            {
              auto t0 = Clock::now();
              temp_h1 = elm_heuristic1->lookup_multiple(i, j);
              add_us(dur_lookup_h1, t0);
            }

            vector<double> *temp_h2 = nullptr;
            {
              auto t0 = Clock::now();
              temp_h2 = elm_heuristic2->lookup_multiple(i, j);
              add_us(dur_lookup_h2, t0);
            }

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
              rmse_h1.at(index_i) +=
                  ((temp_h1->at(index_i) - lb->at(i)->at(j)) *
                   (temp_h1->at(index_i) - lb->at(i)->at(j)));
              elm_h1_saved.at(index_i) +=
                  std::distance(sorted_ub.begin(),
                                lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                            temp_h1->at(index_i)));
              rmse_h2.at(index_i) +=
                  ((temp_h2->at(index_i) - lb->at(i)->at(j)) *
                   (temp_h2->at(index_i) - lb->at(i)->at(j)));
              elm_h2_saved.at(index_i) +=
                  std::distance(sorted_ub.begin(),
                                lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                            temp_h2->at(index_i)));
            }

            if (EXACT) delete temp_exact;
            delete temp_h1;
            delete temp_h2;
            delete temp_sampling;

            // tri query time
            double dummy = 0.;
            auto start_tri_rt = Clock::now();
            dummy = tri.lookup(i, j);
            auto stop_tri_rt = Clock::now();
            duration_tri_rt += std::chrono::duration_cast<std::chrono::microseconds>(stop_tri_rt - start_tri_rt);

            // elm query time
            auto start_elm_rt = Clock::now();
            dummy = elm_sampling->lookup(i, j);
            auto stop_elm_rt = Clock::now();
            duration_elm_rt += std::chrono::duration_cast<std::chrono::microseconds>(stop_elm_rt - start_elm_rt);

            // rmse tri / saved
            auto t_tri_only = Clock::now();
            double tri_val = tri.lookup(i, j);
            add_us(dur_tri_lookup_only, t_tri_only);

            rmse_tri += ((tri_val - lb->at(i)->at(j)) *
                         (tri_val - lb->at(i)->at(j)));
            tri_saved += std::distance(
                sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(),
                                               tri_val));
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

    cout << "[T] small-graph pair loop total elapsed=" << std::fixed << std::setprecision(6)
         << seconds_since(t_pairs) << " s\n";

  } else {
    SCOPE_TIMER("LARGE_GRAPH evaluation (sampling loops)");

    double tri_lb = 0.0;

    // ---- RMSE / Relative WITH repetition ----
    auto t_rmse_loop = Clock::now();
    for (unsigned int s = 0; s < RMSE_SAMPLES; ++s) {
      if ((s % 50) == 0 && s > 0) {
        cout << "[T] RMSE loop progress s=" << s
             << " elapsed=" << std::fixed << std::setprecision(6)
             << seconds_since(t_rmse_loop) << " s\n";
      }

      auto ab = sample_missing_edge(nodes, known_edges);
      unsigned int a = ab.first, b = ab.second;

      auto t_splub = Clock::now();
      double lb_val = splub(adj_lst, a, b);
      if ((s % 50) == 0 && s > 0) {
        cout << "[T] splub sample s=" << s << " took "
             << std::fixed << std::setprecision(6) << seconds_since(t_splub) << " s\n";
      }

      // Keep sample budget consistent: if lb invalid, resample
      if (lb_val <= 0) { --s; continue; }

      vector<double> *temp_exact = nullptr;
      if (EXACT) {
        auto t0 = Clock::now();
        temp_exact = elm_exact->lookup_multiple(a, b);
        add_us(dur_lookup_exact, t0);
      }

      vector<double> *temp_sampling = nullptr;
      {
        auto t0 = Clock::now();
        temp_sampling = elm_sampling->lookup_multiple(a, b);
        add_us(dur_lookup_sampling, t0);
      }

      vector<double> *temp_h1 = nullptr;
      {
        auto t0 = Clock::now();
        temp_h1 = elm_heuristic1->lookup_multiple(a, b);
        add_us(dur_lookup_h1, t0);
      }

      vector<double> *temp_h2 = nullptr;
      {
        auto t0 = Clock::now();
        temp_h2 = elm_heuristic2->lookup_multiple(a, b);
        add_us(dur_lookup_h2, t0);
      }

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

      tri_lb = tri.lookup(a, b);
      rmse_tri += (tri_lb - lb_val) * (tri_lb - lb_val);

      relative += 1.0 - (lb_elm->at(a)->at(b)) / (lb_val);
      relative_count++;
      total_lb_sw += lb_val;

      // query time tri
      auto start_tri_rt = Clock::now();
      tri_lb = tri.lookup(a, b);
      auto stop_tri_rt = Clock::now();
      duration_tri_rt += std::chrono::duration_cast<std::chrono::microseconds>(stop_tri_rt - start_tri_rt);

      // query time elm
      auto start_elm_rt = Clock::now();
      double dummy = elm_sampling->lookup(a, b);
      auto stop_elm_rt = Clock::now();
      duration_elm_rt += std::chrono::duration_cast<std::chrono::microseconds>(stop_elm_rt - start_elm_rt);
    }
    cout << "[T] RMSE loop total elapsed=" << std::fixed << std::setprecision(6)
         << seconds_since(t_rmse_loop) << " s\n";

    // ---- Edge resolution WITHOUT repetition ----
    auto t_erc_loop = Clock::now();
    for (int res_i = 0; res_i < (int)ERC_SAMPLES; res_i++) {
      if ((res_i % 50) == 0 && res_i > 0) {
        cout << "[T] ERC loop progress res_i=" << res_i
             << " elapsed=" << std::fixed << std::setprecision(6)
             << seconds_since(t_erc_loop) << " s\n";
      }

      unsigned int u1 = rand() % nodes;
      unsigned int v1 = rand() % nodes;
      while (u1 == v1) {
        v1 = rand() % nodes;
      }

      unsigned int u2 = rand() % nodes;
      unsigned int v2 = rand() % nodes;
      while (u2 == v2) {
        v2 = rand() % nodes;
      }

      while ((u1 == u2) && (v1 == v2)) {
        u2 = rand() % nodes;
        v2 = rand() % nodes;
        while (u2 == v2) {
          v2 = rand() % nodes;
        }
      }

      auto key_u1v1 = make_pair(min(u1, v1), max(u1, v1));
      auto key_u2v2 = make_pair(min(u2, v2), max(u2, v2));

      u1 = key_u1v1.first;
      u2 = key_u2v2.first;
      v1 = key_u1v1.second;
      v2 = key_u2v2.second;

      vector<double> *temp_exact_u1v1 = nullptr;
      if (EXACT) {
        auto t0 = Clock::now();
        temp_exact_u1v1 = elm_exact->lookup_multiple(u1, v1);
        add_us(dur_lookup_exact, t0);
      }
      vector<double> *temp_sampling_u1v1 = nullptr;
      {
        auto t0 = Clock::now();
        temp_sampling_u1v1 = elm_sampling->lookup_multiple(u1, v1);
        add_us(dur_lookup_sampling, t0);
      }
      vector<double> *temp_h1_u1v1 = nullptr;
      {
        auto t0 = Clock::now();
        temp_h1_u1v1 = elm_heuristic1->lookup_multiple(u1, v1);
        add_us(dur_lookup_h1, t0);
      }
      vector<double> *temp_h2_u1v1 = nullptr;
      {
        auto t0 = Clock::now();
        temp_h2_u1v1 = elm_heuristic2->lookup_multiple(u1, v1);
        add_us(dur_lookup_h2, t0);
      }

      vector<double> *temp_exact_u2v2 = nullptr;
      if (EXACT) {
        auto t0 = Clock::now();
        temp_exact_u2v2 = elm_exact->lookup_multiple(u2, v2);
        add_us(dur_lookup_exact, t0);
      }
      vector<double> *temp_sampling_u2v2 = nullptr;
      {
        auto t0 = Clock::now();
        temp_sampling_u2v2 = elm_sampling->lookup_multiple(u2, v2);
        add_us(dur_lookup_sampling, t0);
      }
      vector<double> *temp_h1_u2v2 = nullptr;
      {
        auto t0 = Clock::now();
        temp_h1_u2v2 = elm_heuristic1->lookup_multiple(u2, v2);
        add_us(dur_lookup_h1, t0);
      }
      vector<double> *temp_h2_u2v2 = nullptr;
      {
        auto t0 = Clock::now();
        temp_h2_u2v2 = elm_heuristic2->lookup_multiple(u2, v2);
        add_us(dur_lookup_h2, t0);
      }

      // Dijkstra for UB comparisons
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

      tri_saved += std::distance(sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(), tri_lb));
    }
    cout << "[T] ERC loop total elapsed=" << std::fixed << std::setprecision(6)
         << seconds_since(t_erc_loop) << " s\n";
  }

  cout << "Query time elb : " << duration_elm_rt.count() << " us" << endl;
  cout << "Query time tri : " << duration_tri_rt.count() << " us" << endl;

  cout << "[T] lookup_multiple totals (us): "
       << "sampling=" << dur_lookup_sampling.count()
       << " h1=" << dur_lookup_h1.count()
       << " h2=" << dur_lookup_h2.count();
  if (EXACT) cout << " exact=" << dur_lookup_exact.count();
  cout << endl;

  cout << "Duration SW LB: " << duration_lb_sw_us.count() / 1000000.0 << " s" << endl;
  cout << "Total Original Graph weight : " << graph_weight_orig << endl;
  cout << "Total SW LB : " << total_lb_sw << endl;
  cout << "Sum Relative Error on edge : " << (relative_count ? (relative / relative_count) : 0.0) << endl;

  for (unsigned index_i = 0; index_i < k; ++index_i) {
    cout << "root mean square error sampling for " << index_i << ": "
         << sqrt(rmse_sampling.at(index_i) / (relative_count ? relative_count : 1)) << endl;
    cout << "Saved Sampling " << index_i << ": "
         << elm_sampling_saved.at(index_i) << endl;

    if (EXACT) {
      cout << "root mean square error exact " << index_i << ": "
           << sqrt(rmse_exact.at(index_i) / (relative_count ? relative_count : 1)) << endl;
      cout << "Saved Exact " << index_i << ": " << elm_exact_saved.at(index_i) << endl;
    }

    cout << "root mean square error H1 " << index_i << ": "
         << sqrt(rmse_h1.at(index_i) / (relative_count ? relative_count : 1)) << endl;
    cout << "Saved H1 " << index_i << ": " << elm_h1_saved.at(index_i) << endl;

    cout << "root mean square error H2 " << index_i << ": "
         << sqrt(rmse_h2.at(index_i) / (relative_count ? relative_count : 1)) << endl;
    cout << "Saved H2 " << index_i << ": " << elm_h2_saved.at(index_i) << endl;
  }

  cout << "root mean square error tri: " << sqrt(rmse_tri / (relative_count ? relative_count : 1))
       << endl;
  cout << "Saved TriSearch " << tri_saved << endl;
  cout << "Saved SW " << sw_saved << endl;
  cout << "Size of TriSearch " << tri._sizeof() << std::endl;
  cout << "Size of ELM " << elm_sampling->_sizeof() << std::endl;

  {
    SCOPE_TIMER("cleanup + deletes");
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

    {
      SCOPE_TIMER("delete distance matrix rows");
      for (unsigned int i = 0; i < nodes; ++i) {
        if (DEBUG) {
          for (unsigned int j = 0; j < nodes; ++j) {
            cout << " " << distance->at(i)->at(j);
          }
          cout << endl;
        }
        delete distance->at(i);
      }
    }
    delete distance;
  }

  checkpoint("end of main", T_GLOBAL);
  return 0;
}




