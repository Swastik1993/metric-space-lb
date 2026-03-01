// knn_eval_main.cpp
// ------------------------------------------------------------
// Bound-based k-NN evaluation with PRINTED metrics (no plots).
//
// Per query source node X:
//   UB(X,v) = exact shortest path distance via DijkstraELM
//   LB_tri(X,v) = TriSearch.lookup
//   LB_elm(X,v) = ELM_SPTree.lookup (optional; can be disabled)
//   LB_comb(X,v) = max(LB_tri, LB_elm)  (tighter lower bound)
//
// Metrics (printed):
//   (A) "Hasse-resolved %" (sampled): fraction of node pairs (a,b) where
//          UB(a) < LB(b)  OR  UB(b) < LB(a)
//       This measures how many pairwise order relations are *certified*.
//   (B) Top-K certification counts:
//          definite_in  : UB(v) < Kth-smallest LB among all nodes (excluding X)
//          definite_out : LB(v) > Kth-smallest UB among all nodes (excluding X)
//          unresolved   : everything else
//       plus precision of definite_in vs ground-truth top-K.
//
// Build (from rmq_cartesian/ folder):
//   g++ -O3 -std=c++20 -I. knn_eval_main.cpp \
//       Dijkstra.cpp TriSearch.cpp ELM_SPTree.cpp EdgeLandMark.cpp \
//       ShortestPathTree.cpp spt.cpp LowerBound.cpp GraphUtils.cpp \
//       -o knn_eval
//
// Example:
//   ./knn_eval 1 500 0.03 10 2000 25 --elm 1 --elm-mode sampling --elm-landmarks 10 --elm-samples 2000
//
// Positional args:
//   <seed> <nodes> <edge_density> <K> <pair_samples> <num_queries>
//
// Flags:
//   --elm 0|1                 (default 0)
//   --elm-mode sampling|exact  (default sampling)
//   --elm-landmarks <int>      (default 10)
//   --elm-samples <int>        (default 2000)
//   --distance-matrix <file>   (optional; uses get_adj_matrix_file)
//
// Notes:
//   - If you hit crashes with ELM, try higher edge_density or run with --elm 0 first.
//
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Dijkstra.h"
#include "ELM_SPTree.h"
#include "GraphDefinitions.h"
#include "GraphUtils.h"
#include "TriSearch.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::vector;

static inline bool is_finite(double x) { return std::isfinite(x); }

static inline unsigned int clampK(unsigned int K, unsigned int n_minus_one) {
  return std::min(K, n_minus_one);
}

static vector<unsigned int> topk_truth(const vector<double>& dist, unsigned int src, unsigned int K) {
  vector<unsigned int> ids;
  ids.reserve(dist.size() - 1);
  for (unsigned int v = 0; v < dist.size(); ++v) if (v != src) ids.push_back(v);

  auto cmp = [&](unsigned int a, unsigned int b) {
    if (dist[a] != dist[b]) return dist[a] < dist[b];
    return a < b;
  };

  if (K > ids.size()) K = (unsigned int)ids.size();
  std::nth_element(ids.begin(), ids.begin() + K, ids.end(), cmp);
  ids.resize(K);
  std::sort(ids.begin(), ids.end(), cmp);
  return ids;
}

static double kth_smallest_excluding_self(const vector<double>& arr, unsigned int self, unsigned int K) {
  vector<double> tmp;
  tmp.reserve(arr.size() - 1);
  for (unsigned int i = 0; i < arr.size(); ++i) {
    if (i == self) continue;
    double x = arr[i];
    if (!is_finite(x)) x = std::numeric_limits<double>::infinity();
    tmp.push_back(x);
  }
  if (tmp.empty()) return 0.0;
  K = clampK(K, (unsigned int)tmp.size());
  std::nth_element(tmp.begin(), tmp.begin() + (K - 1), tmp.end());
  return tmp[K - 1];
}

static inline void sample_resolved_pairs(
    const vector<double>& LB, const vector<double>& UB,
    unsigned int src, unsigned int samples, std::mt19937_64& rng,
    unsigned long long& resolved_out, unsigned long long& tried_out) {

  std::uniform_int_distribution<unsigned int> uni(0, (unsigned int)LB.size() - 1);

  for (unsigned int t = 0; t < samples; ++t) {
    unsigned int a = uni(rng);
    unsigned int b = uni(rng);
    if (a == src || b == src || a == b) continue;

    double lba = LB[a], lbb = LB[b];
    double uba = UB[a], ubb = UB[b];

    if (!is_finite(lba) || lba < 0) lba = 0.0;
    if (!is_finite(lbb) || lbb < 0) lbb = 0.0;
    if (!is_finite(uba)) uba = std::numeric_limits<double>::infinity();
    if (!is_finite(ubb)) ubb = std::numeric_limits<double>::infinity();

    ++tried_out;
    if (uba < lbb || ubb < lba) ++resolved_out;
  }
}

struct TopKStats {
  unsigned long long definite_in = 0;
  unsigned long long definite_out = 0;
  unsigned long long unresolved = 0;

  unsigned long long definite_in_correct = 0; // among definite_in
  unsigned long long definite_in_total = 0;

  unsigned long long queries = 0;
};

static TopKStats eval_topk(
    const vector<double>& LB, const vector<double>& UB,
    const vector<double>& exact, unsigned int src, unsigned int K) {

  TopKStats st;
  st.queries = 1;

  unsigned int n = (unsigned int)LB.size();
  if (n <= 1) return st;
  K = clampK(K, n - 1);

  vector<unsigned int> truth = topk_truth(exact, src, K);
  std::unordered_set<unsigned int> truth_set(truth.begin(), truth.end());

  const double kthLB = kth_smallest_excluding_self(LB, src, K);
  const double kthUB = kth_smallest_excluding_self(UB, src, K);

  for (unsigned int v = 0; v < n; ++v) {
    if (v == src) continue;

    double lb = LB[v], ub = UB[v];
    if (!is_finite(lb) || lb < 0) lb = 0.0;
    if (!is_finite(ub)) ub = std::numeric_limits<double>::infinity();

    bool in_cert = (ub < kthLB);
    bool out_cert = (lb > kthUB);

    if (in_cert) {
      st.definite_in++;
      st.definite_in_total++;
      if (truth_set.find(v) != truth_set.end()) st.definite_in_correct++;
    } else if (out_cert) {
      st.definite_out++;
    } else {
      st.unresolved++;
    }
  }

  return st;
}

static void usage() {
  cout << "Usage:\n"
       << "  ./knn_eval <seed> <nodes> <edge_density> <K> <pair_samples> <num_queries>\n"
       << "           [--elm 0|1] [--elm-mode sampling|exact] [--elm-landmarks L] [--elm-samples S]\n"
       << "           [--distance-matrix <file>]\n";
}

static bool pop_kv(int& i, int argc, char** argv, const string& key, string& out) {
  if (string(argv[i]) != key) return false;
  if (i + 1 >= argc) {
    cerr << "[ERROR] Missing value for " << key << endl;
    std::exit(2);
  }
  out = argv[i + 1];
  i += 1;
  return true;
}

int main(int argc, char** argv) {
  if (argc < 7) {
    usage();
    cerr << "\n[ERROR] Need at least 6 positional args.\n";
    return 2;
  }

  unsigned int seed = (unsigned int)std::stoul(argv[1]);
  unsigned int nodes = (unsigned int)std::stoul(argv[2]);
  double prob = std::stod(argv[3]);
  unsigned int K = (unsigned int)std::stoul(argv[4]);
  unsigned int pair_samples = (unsigned int)std::stoul(argv[5]);
  unsigned int num_queries = (unsigned int)std::stoul(argv[6]);

  bool use_elm = false;
  string elm_mode = "sampling";
  unsigned int elm_landmarks = 10;
  unsigned int elm_samples = 2000;
  string dist_file = "";

  for (int i = 7; i < argc; ++i) {
    string v;
    if (pop_kv(i, argc, argv, "--elm", v)) { use_elm = (v != "0"); continue; }
    if (pop_kv(i, argc, argv, "--elm-mode", v)) { elm_mode = v; continue; }
    if (pop_kv(i, argc, argv, "--elm-landmarks", v)) { elm_landmarks = (unsigned int)std::stoul(v); continue; }
    if (pop_kv(i, argc, argv, "--elm-samples", v)) { elm_samples = (unsigned int)std::stoul(v); continue; }
    if (pop_kv(i, argc, argv, "--distance-matrix", v)) { dist_file = v; continue; }
    cerr << "[WARN] Unrecognized arg: " << argv[i] << endl;
  }

  srand(seed);
  std::mt19937_64 rng(seed);

  // Distance matrix
  vector<vector<double>*>* distance = nullptr;
  if (!dist_file.empty()) distance = get_adj_matrix_file(const_cast<char*>(dist_file.c_str()), 0.);
  else distance = distance_matrix(nodes, 6);

  // Build ER graph
  Graph g;
  boost::minstd_rand gen;
  g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);

  // Initialize vertex indices 0..n-1 (GraphUtils assumes this)
  auto id = get(&VertexProperties::index, g);
  boost::graph_traits<Graph>::vertex_iterator vi, viend;
  unsigned int vert_num = 0;
  for (tie(vi, viend) = vertices(g); vi != viend; ++vi) id[*vi] = vert_num++;

  // adjacency list with edge weights from distance matrix
  vector<std::list<pair<unsigned int, double>>*>* adj_list = get_adjacency_list(g, id, distance);
  cout << "[INFO] Built adjacency list. nodes=" << adj_list->size() << " prob=" << prob << endl;

  // Lower bound structures
  TriSearch tri(adj_list, (unsigned int)adj_list->size());

  ELM_SPTree* elm = nullptr;
  if (use_elm) {
    elm = new ELM_SPTree(adj_list, (unsigned int)adj_list->size(), elm_landmarks, elm_samples);
    if (elm_mode == "exact") elm->get_exact_landmarks();
    else elm->get_sampling_landmarks();
    cout << "[INFO] ELM enabled. mode=" << elm_mode
         << " landmarks=" << elm_landmarks
         << " samples=" << elm_samples << endl;
  } else {
    cout << "[INFO] ELM disabled (use TriSearch only). Enable with --elm 1" << endl;
  }

  std::uniform_int_distribution<unsigned int> qdist(0, (unsigned int)adj_list->size() - 1);

  unsigned long long tri_resolved = 0, tri_tried = 0;
  unsigned long long elm_resolved = 0, elm_tried = 0;
  unsigned long long comb_resolved = 0, comb_tried = 0;

  TopKStats tri_topk, elm_topk, comb_topk;

  for (unsigned int qi = 0; qi < num_queries; ++qi) {
    unsigned int src = qdist(rng);

    // exact UB
    auto exact_pair = DijkstraELM(adj_list, src);
    vector<double>* exact_dist = exact_pair.first;

    vector<double> UB(adj_list->size(), std::numeric_limits<double>::infinity());
    vector<double> LB_tri(adj_list->size(), 0.0);
    vector<double> LB_elm(adj_list->size(), 0.0);
    vector<double> LB_comb(adj_list->size(), 0.0);

    for (unsigned int v = 0; v < adj_list->size(); ++v) {
      UB[v] = exact_dist->at(v);
      if (v == src) continue;

      double lbt = tri.lookup(src, v);
      if (!is_finite(lbt) || lbt < 0) lbt = 0.0;
      LB_tri[v] = lbt;

      double lbe = 0.0;
      if (elm) {
        lbe = elm->lookup(src, v);
        if (!is_finite(lbe) || lbe < 0) lbe = 0.0;
      }
      LB_elm[v] = lbe;

      LB_comb[v] = std::max(LB_tri[v], LB_elm[v]);
    }

    // (A) Hasse-style resolved comparisons (sampled)
    unsigned long long r = 0, t = 0;
    sample_resolved_pairs(LB_tri, UB, src, pair_samples, rng, r, t);
    tri_resolved += r; tri_tried += t;

    r = 0; t = 0;
    sample_resolved_pairs(LB_elm, UB, src, pair_samples, rng, r, t);
    elm_resolved += r; elm_tried += t;

    r = 0; t = 0;
    sample_resolved_pairs(LB_comb, UB, src, pair_samples, rng, r, t);
    comb_resolved += r; comb_tried += t;

    // (B) Top-K certification
    TopKStats t1 = eval_topk(LB_tri, UB, *exact_dist, src, K);
    tri_topk.definite_in += t1.definite_in;
    tri_topk.definite_out += t1.definite_out;
    tri_topk.unresolved += t1.unresolved;
    tri_topk.definite_in_correct += t1.definite_in_correct;
    tri_topk.definite_in_total += t1.definite_in_total;
    tri_topk.queries++;

    TopKStats t2 = eval_topk(LB_elm, UB, *exact_dist, src, K);
    elm_topk.definite_in += t2.definite_in;
    elm_topk.definite_out += t2.definite_out;
    elm_topk.unresolved += t2.unresolved;
    elm_topk.definite_in_correct += t2.definite_in_correct;
    elm_topk.definite_in_total += t2.definite_in_total;
    elm_topk.queries++;

    TopKStats t3 = eval_topk(LB_comb, UB, *exact_dist, src, K);
    comb_topk.definite_in += t3.definite_in;
    comb_topk.definite_out += t3.definite_out;
    comb_topk.unresolved += t3.unresolved;
    comb_topk.definite_in_correct += t3.definite_in_correct;
    comb_topk.definite_in_total += t3.definite_in_total;
    comb_topk.queries++;

    auto prec = [](unsigned long long ok, unsigned long long tot) -> double {
      return tot ? (double)ok / (double)tot : 0.0;
    };

    cout << "\n[QUERY " << (qi + 1) << "/" << num_queries << "] src=" << src << endl;
    cout << "  TopK(Tri):  in=" << t1.definite_in << " out=" << t1.definite_out
         << " unk=" << t1.unresolved
         << " prec=" << prec(t1.definite_in_correct, t1.definite_in_total) << endl;
    cout << "  TopK(ELM):  in=" << t2.definite_in << " out=" << t2.definite_out
         << " unk=" << t2.unresolved
         << " prec=" << prec(t2.definite_in_correct, t2.definite_in_total) << endl;
    cout << "  TopK(Comb): in=" << t3.definite_in << " out=" << t3.definite_out
         << " unk=" << t3.unresolved
         << " prec=" << prec(t3.definite_in_correct, t3.definite_in_total) << endl;

    delete exact_pair.first;
    delete exact_pair.second;
  }

  auto pct = [](unsigned long long a, unsigned long long b) -> double {
    return b ? (100.0 * (double)a / (double)b) : 0.0;
  };

  auto print_topk_agg = [&](const char* name, const TopKStats& st) {
    double prec = st.definite_in_total ? (double)st.definite_in_correct / (double)st.definite_in_total : 0.0;
    double q = st.queries ? (double)st.queries : 1.0;
    cout << "TopK[" << name << "] avg_in=" << (double)st.definite_in / q
         << " avg_out=" << (double)st.definite_out / q
         << " avg_unk=" << (double)st.unresolved / q
         << " definite_in_precision=" << prec
         << " (correct=" << st.definite_in_correct << "/" << st.definite_in_total << ")"
         << endl;
  };

  cout << "\n====================\n";
  cout << "AGGREGATE (" << num_queries << " queries, K=" << K << ")\n";
  cout << "====================\n";
  cout << "Hasse-sampled resolved%:\n";
  cout << "  Tri:  " << pct(tri_resolved, tri_tried) << "% (resolved=" << tri_resolved << " tried=" << tri_tried << ")\n";
  cout << "  ELM:  " << pct(elm_resolved, elm_tried) << "% (resolved=" << elm_resolved << " tried=" << elm_tried << ")\n";
  cout << "  Comb: " << pct(comb_resolved, comb_tried) << "% (resolved=" << comb_resolved << " tried=" << comb_tried << ")\n";

  print_topk_agg("Tri", tri_topk);
  print_topk_agg("ELM", elm_topk);
  print_topk_agg("Comb", comb_topk);

  if (elm) delete elm;
  clean_up_adj_list(adj_list);
  clean_up_adj_matrix(distance);
  return 0;
}



