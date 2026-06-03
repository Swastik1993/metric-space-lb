#include <time.h>
#include <random>
#include <sys/resource.h>
#include "Dijkstra.h"
#include "ELM_SPTree.h"
#include "spt.h"
#include "EdgeLandMark.h"
#include "GraphDefinitions.h"
#include "GraphUtils.h"
#include "LowerBound.h"
#include "ShortestPathTree.h"
#include "TriSearch.h"
#include "Downstream.h"
#include "Subsample.h"
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <atomic>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <map>
#include <set>
#include <cmath>
#include <tuple>
#include <unordered_set>
#include <iomanip>
#include <limits>
#include <cassert>
#include <fstream>
#include <sstream>
#include <filesystem>

using namespace std;

bool EXACT;
const int DEBUG = 0;
const bool RENYI_ERDOS = true;
const bool WAXMAN = !RENYI_ERDOS;
// Promoted from `const` to runtime-tunable so the CLI can override them.
// Defaults preserved exactly as before.
unsigned int RMSE_SAMPLES = 10000;
unsigned int ERC_SAMPLES  = 10000;
double STOCH_EPS = 0.10;
// Default candidate cap for the stochastic-greedy variants (ELM Stochastic
// and SPT Stochastic). 1000 is large enough that on n=50000 graphs the
// Mirzasoleiman (1 - 1/e - eps) bound is rarely the binding constraint at
// realistic k (>= 64), while keeping each iteration to ~1000 SPT builds.
// Override at the CLI with --stoch-cand-cap N.
unsigned int STOCH_CAND_CAP = 1000;

// Ablation flags (Reviewer R3.D7). Set by --ablation-no-hca,
// --ablation-no-rmq, --ablation-no-spt on the CLI. Propagated to
// every ELM_SPTree built in multiple_sample_exec. Default: all false (=
// paper behavior).
bool ABLATION_NO_HCA = false;
bool ABLATION_NO_RMQ = false;
bool ABLATION_NO_SPT = false;

// Single-run ablation timing+quality report (--ablation-report). When set,
// multiple_sample_exec runs a self-contained component breakdown instead of
// the normal benchmark pipeline: it measures per-component build and query
// times (SPT, RMQ, HCA) AND evaluates the same fixed query set with each
// component optimized vs unoptimized, all in ONE process. The number of
// fixed evaluation queries is configurable.
bool ABLATION_REPORT = false;
unsigned int ABLATION_REPORT_QUERIES = 5000;
// SGEL objective for the HCA/SGEL sampling model.
// raw      = original L1 objective: sum of marginal lower-bound increases.
// relative = Zoom-discussed ARLB-aligned objective: sum of marginal
//            lower-bound increases divided by exact/SPLUB OPT per query.
// both     = build both models for the requested ablation.
std::string SGEL_OBJECTIVE_MODE = "relative";

// Downstream tasks (Reviewer R3.D2). Enabled by --downstream. Runs an MST
// (Prim's) + PAM K-Medoid per method using the materialized LB matrix.
// Compares to true-distance reference; reports edge overlap + Jaccard for
// MST and intra-cluster cost for PAM. Requires legacy matrix mode
// (distance != nullptr) -- skipped automatically in --sparse-graph runs.
bool RUN_DOWNSTREAM = false;
// Persistent landmark cache (saves the expensive selection step across runs).
// --cache-dir <path> enables it. Cache file name is keyed on dataset, seed, n,
// density, K, sampling_size; the file's metadata header additionally encodes
// the objective mode + ablation flags that affect landmark selection so
// landmarks selected under a different objective/variant are never silently reused.
std::string CACHE_DIR = "";
// Persistent KNN bound-row cache. Separate from the landmark cache. When
// enabled, KNN evaluation loads/saves UB rows and per-method LB rows by query
// source, so a later run with more query nodes reuses rows computed earlier.
std::string KNN_BOUND_CACHE_DIR = "";
// SIFT/high-D diagnostic: select query nodes by LB/UB separability instead of
// a random query subset. Enabled by --knn-cert-select X K. The program scans
// every node, computes UB(q,*) and Combined(max) LB(q,*), scores each query by
// how many candidates satisfy LB(q,v) > kth-smallest UB(q,*), then evaluates
// KNN only on the top-X most certifiable query nodes.
bool KNN_CERT_SELECT = false;
std::string CACHE_DATASET = "";   // basename of the dataset file
double CACHE_DENSITY = 0.0;       // density (argv[3])
// Non-paper quality-degradation metrics (LB used AS the distance, producing a
// different/wrong MST and PAM cost per method). OFF by default: the default
// downstream output is the paper-faithful Augustine protocol only (exact
// output + oracle-call savings). Enable with --downstream-lb-as-distance.
bool DOWNSTREAM_LB_AS_DISTANCE = false;
std::vector<unsigned int> CLUSTER_KS = {3, 5, 10};
// v25: include SPLUB as an explicit downstream baseline only for small
// matrix-mode datasets by default. Materializing SPLUB LB for all n^2 pairs
// runs Dijkstra for many sources and is too expensive for Porto/SIFT unless
// the user explicitly raises this cap.
unsigned int DOWNSTREAM_SPLUB_MAX_N = 500;

// SPLUB as an explicit baseline INSIDE the KNN / ERC (PairResolved) metrics
// (Reviewer R1.O2 / R3.D6). SPLUB is the exact lower bound, so it is the
// upper envelope any approximate method can reach for both edge resolution
// and kNN certification. Producing a SPLUB row requires all-pairs shortest
// paths (one Dijkstra per source) so the per-(q,v) exact LB can be folded;
// that is only affordable on small graphs. KNN_SPLUB_ENABLE turns it on,
// gated by KNN_SPLUB_MAX_N (graphs with more nodes skip SPLUB, matching the
// paper, which does not report SPLUB on Porto/SIFT). Override with
// --knn-splub / --no-knn-splub and --knn-splub-max-n N.
bool KNN_SPLUB_ENABLE = true;
unsigned int KNN_SPLUB_MAX_N = 2000;

// Node subsampling (v16, Reviewer R1.W3 / R3.D5 + variance-across-seeds).
// When --subsample-nodes N is passed in --sparse-graph mode, the loaded
// full graph is replaced with a random N-vertex induced subgraph BEFORE
// any preprocessing. --subsample-seed sets the RNG seed (independent of
// the binary's main `init` seed); on connectivity failure the seed is
// incremented by 1 and the sample is redrawn, up to
// --subsample-max-retries times before a hard error.
unsigned int SUBSAMPLE_NODES = 0;       // 0 = disabled
unsigned int SUBSAMPLE_SEED = 1;
unsigned int SUBSAMPLE_MAX_RETRIES = 50;
// v16.1: when false (default), the binary uses an LCC fallback when the
// induced subgraph has > 1 component, dropping fringe nodes outside the
// giant component. When true (opt-in via --subsample-strict-connectivity),
// the binary instead retries with seed+1 up to SUBSAMPLE_MAX_RETRIES times
// and hard-errors if no single-component sample is found. The default was
// flipped in v16.1 after the v16 strict policy failed on every seed for the
// production SIFT n=100k, M=24 graph (consistently 10-25 small fringe
// components per random 50k subsample; the giant component still covers
// 99.97% of nodes, so LCC fallback is the right default for this data).
bool SUBSAMPLE_STRICT_CONNECTIVITY = false;
// Sampling-size sweep is now configurable via --sweep "50,100,200,300".
// std::vector<unsigned int> SWEEP_SAMPLE_SIZES = {50, 100, 200, 300};
std::vector<unsigned int> SWEEP_SAMPLE_SIZES = {300};

static inline uint64_t cpu_time_us() {
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return uint64_t(ts.tv_sec) * 1000000ULL + ts.tv_nsec / 1000ULL;
}

static size_t get_peak_rss_kb() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_maxrss;
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

  delete a_shortest_path.first;
  if (a_shortest_path.second) delete a_shortest_path.second;
  delete b_shortest_path.first;
  if (b_shortest_path.second) delete b_shortest_path.second;

  return lb;
}

// Cached SPLUB: caches Dijkstra-from-source results across calls. When the
// same endpoint reappears in many sampled queries (the common case for
// thousands of RMSE / ERC samples), this cuts ground-truth computation from
// 2K Dijkstras to one per unique endpoint. Equivalent semantics to splub();
// shares its arithmetic exactly.
static double splub_cached(SplubCache &cache, unsigned int a, unsigned int b) {
  return cache.splub(a, b);
}

// ---- SPLUB support for the KNN / ERC metrics --------------------------------
// Build all-pairs shortest-path distance rows (one Dijkstra per source). Only
// called on small graphs (gated by KNN_SPLUB_MAX_N) because it is O(n) Dijkstras
// and stores n*n doubles. Parallelized over sources; each thread allocates its
// own Dijkstra scratch, so the writes into distinct SP[src] rows are race-free.
static void build_all_pairs_sp(
    vector<list<pair<unsigned int, double>> *> *adj_lst,
    unsigned int nodes,
    vector<vector<double>> &SP) {
  SP.assign(nodes, vector<double>());
  #pragma omp parallel for schedule(dynamic, 16)
  for (long long s = 0; s < (long long)nodes; ++s) {
    auto r = DijkstraELM(adj_lst, (unsigned int)s);
    SP[(size_t)s] = *r.first;        // copy the distance row
    delete r.first;
    if (r.second) delete r.second;
  }
}

// Exact SPLUB lower-bound row for a fixed source `src` against every node v,
// folded from the precomputed all-pairs SP rows. This reproduces splub()'s
// arithmetic exactly: LB(src,v) = max over known directed edges (i,j) with
// length L of { L - sp(src,i) - sp(v,j), L - sp(src,j) - sp(v,i) }, floored at 0.
// O(n * m) for the whole row; only used on small graphs.
static void splub_lb_row(
    vector<list<pair<unsigned int, double>> *> *adj_lst,
    unsigned int nodes,
    const vector<vector<double>> &SP,
    unsigned int src,
    vector<double> &LB_out) {
  const vector<double> &sa = SP[src];
  for (unsigned int v = 0; v < nodes; ++v) {
    if (v == src) { LB_out[v] = 0.0; continue; }
    const vector<double> &sb = SP[v];
    double best = 0.0;
    for (unsigned int i = 0; i < nodes; ++i) {
      const double sa_i = sa[i], sb_i = sb[i];
      for (const auto &item : *adj_lst->at(i)) {
        const unsigned int j = item.first;
        const double L = item.second;
        const double t1 = L - sa_i - sb[j];
        if (t1 > best) best = t1;
        const double t2 = L - sa[j] - sb_i;
        if (t2 > best) best = t2;
      }
    }
    LB_out[v] = best;
  }
}

static double matrix_max_offdiag(vector<vector<double> *> *mat) {
  if (mat == nullptr) return 1.0;
  double mx = 0.0;
  for (unsigned int i = 0; i < mat->size(); ++i) {
    if (mat->at(i) == nullptr) continue;
    for (unsigned int j = 0; j < mat->at(i)->size(); ++j) {
      if (i == j) continue;
      double v = mat->at(i)->at(j);
      if (std::isfinite(v)) mx = std::max(mx, v);
    }
  }
  return (mx > 0.0) ? mx : 1.0;
}

static SimpleGraph build_simple_graph(
    vector<list<pair<unsigned int, double>> *> *adj_lst,
    unsigned int nodes) {
  SimpleGraph g((int)nodes);
  for (unsigned int u = 0; u < nodes; ++u) {
    for (const auto &pr : *adj_lst->at(u)) {
      unsigned int v = pr.first;
      double w = pr.second;
      if (u < v) g.addEdge((int)u, (int)v, w);
    }
  }
  return g;
}

struct CurveModel {
  string name;
  function<double(unsigned int, unsigned int)> lookup;
  function<vector<double>*(unsigned int, unsigned int)> lookup_multiple;
  function<size_t()> size_of;
  uint64_t build_time_us = 0;
  size_t peak_rss_kb = 0;
};

static inline bool km_isfinite(double x) { return std::isfinite(x); }

struct KmAgg {
  unsigned long long definite_in = 0;
  unsigned long long definite_out = 0;
  unsigned long long unresolved = 0;
  unsigned long long definite_in_total = 0;
  unsigned long long definite_in_correct = 0;
  unsigned long long pair_resolved = 0;
  unsigned long long pair_tried = 0;
};

struct KnnQuerySelectScore {
  unsigned src = 0;
  unsigned long long certified_out_like = 0; // count LB(q,v) > kth-smallest UB(q,*)
  unsigned long long denom = 0;              // candidates_per_query - min(K,candidates)
  double score_pct = 0.0;
  double tau_k = 0.0;                        // kth-smallest UB threshold
  double lb_boundary = 0.0;                  // (C-K)th largest LB, if defined
  double gap = 0.0;                          // lb_boundary - tau_k
};

static inline vector<unsigned> km_truth_topk(const vector<double>& exact, unsigned src, unsigned K) {
  const unsigned n = (unsigned)exact.size();
  vector<unsigned> idx;
  idx.reserve(n > 0 ? n - 1 : 0);
  for (unsigned v = 0; v < n; ++v) if (v != src) idx.push_back(v);
  auto cmp = [&](unsigned a, unsigned b) {
    double da = exact[a], db = exact[b];
    if (!km_isfinite(da)) da = numeric_limits<double>::infinity();
    if (!km_isfinite(db)) db = numeric_limits<double>::infinity();
    if (da != db) return da < db;
    return a < b;
  };
  if (K > idx.size()) K = (unsigned)idx.size();
  if (K == 0) return {};
  nth_element(idx.begin(), idx.begin() + (K - 1), idx.end(), cmp);
  idx.resize(K);
  sort(idx.begin(), idx.end(), cmp);
  return idx;
}


// Deterministic candidate-pool construction for candidate-limited KNN.
// If candidate_limit==0 or candidate_limit>=n-1, this returns all non-src
// candidates. Otherwise it returns a seeded, reproducible subset of exactly
// candidate_limit non-src nodes.
static inline vector<unsigned> km_candidate_pool(unsigned n,
                                                 unsigned src,
                                                 unsigned candidate_limit,
                                                 unsigned long long seed) {
  vector<unsigned> cand;
  if (n < 2) return cand;
  cand.reserve(n - 1);
  for (unsigned v = 0; v < n; ++v) if (v != src) cand.push_back(v);
  if (candidate_limit == 0 || candidate_limit >= cand.size()) return cand;
  uint64_t x = (uint64_t)seed ^ (0x9e3779b97f4a7c15ULL + ((uint64_t)src << 6) + ((uint64_t)src >> 2));
  x ^= (x >> 30); x *= 0xbf58476d1ce4e5b9ULL;
  x ^= (x >> 27); x *= 0x94d049bb133111ebULL;
  x ^= (x >> 31);
  std::mt19937 rng((uint32_t)(x ^ (x >> 32)));
  std::shuffle(cand.begin(), cand.end(), rng);
  cand.resize(candidate_limit);
  std::sort(cand.begin(), cand.end());
  return cand;
}

static inline vector<unsigned> km_truth_topk_candidates(const vector<double>& exact,
                                                        const vector<unsigned>& candidates,
                                                        unsigned K) {
  vector<unsigned> idx = candidates;
  auto cmp = [&](unsigned a, unsigned b) {
    double da = exact[a], db = exact[b];
    if (!km_isfinite(da)) da = numeric_limits<double>::infinity();
    if (!km_isfinite(db)) db = numeric_limits<double>::infinity();
    if (da != db) return da < db;
    return a < b;
  };
  if (K > idx.size()) K = (unsigned)idx.size();
  if (K == 0) return {};
  nth_element(idx.begin(), idx.begin() + (K - 1), idx.end(), cmp);
  idx.resize(K);
  sort(idx.begin(), idx.end(), cmp);
  return idx;
}

static inline void km_eval_topk(const vector<double>& LB,
                                const vector<double>& UB,
                                const vector<double>& exact,
                                unsigned src,
                                unsigned K,
                                KmAgg& agg) {
  if (LB.size() < 2) return;
  const unsigned n = (unsigned)LB.size();
  K = std::min<unsigned>(K, n - 1);
  if (K == 0) return;

  vector<unsigned> truth = km_truth_topk(exact, src, K);
  unordered_set<unsigned> truth_set(truth.begin(), truth.end());

  for (unsigned v = 0; v < n; ++v) {
    if (v == src) continue;
    double lbv = LB[v], ubv = UB[v];
    if (!km_isfinite(lbv) || lbv < 0) lbv = 0.0;
    if (!km_isfinite(ubv)) ubv = numeric_limits<double>::infinity();

    unsigned definitely_closer = 0;
    for (unsigned u = 0; u < n; ++u) {
      if (u == src || u == v) continue;
      double ubu = UB[u];
      if (!km_isfinite(ubu)) ubu = numeric_limits<double>::infinity();
      if (ubu < lbv) {
        definitely_closer++;
        if (definitely_closer >= K) break;
      }
    }
    bool out_cert = (definitely_closer >= K);

    unsigned possibly_closer = 0;
    for (unsigned u = 0; u < n; ++u) {
      if (u == src || u == v) continue;
      double lbu = LB[u];
      if (!km_isfinite(lbu) || lbu < 0) lbu = 0.0;
      if (lbu < ubv) {
        possibly_closer++;
        if (possibly_closer >= K) break;
      }
    }
    bool in_cert = (possibly_closer <= K - 1);

    if (in_cert) {
      agg.definite_in++;
      agg.definite_in_total++;
      if (truth_set.count(v)) agg.definite_in_correct++;
    } else if (out_cert) {
      agg.definite_out++;
    } else {
      agg.unresolved++;
    }
  }
}

// Fast O(n log n) equivalent of km_eval_topk. Produces IDENTICAL certification
// counts (definite_in / definite_out / unresolved / precision) but replaces
// the per-candidate O(n) inner scans with two sorted-array binary searches:
//
//   definitely_closer(v) = #{u != src : UB[u] <  LB[v]}      (out-cert if >= K)
//   possibly_closer(v)   = #{u != src, u != v : LB[u] < UB[v]} (in-cert if <= K-1)
//
// Build sorted copies of UB and LB over all u != src once per query (O(n log n)),
// then each candidate v is O(log n). This is what makes all-N-queries feasible
// at n=50k; the O(n^2)-per-query original would take days over 50k queries.
//
// Correctness vs the original:
//  - definitely_closer: u==v never counts because UB[v] >= LB[v] (LB is a valid
//    lower bound, UB the exact SP), so UB[v] < LB[v] is false. No self-exclusion
//    needed. src is excluded by construction (omitted from the sorted array).
//  - possibly_closer: u==v WOULD be counted when LB[v] < UB[v], so we subtract 1
//    in that case. src excluded by construction.
static inline void km_eval_topk_fast(const vector<double>& LB,
                                     const vector<double>& UB,
                                     const vector<double>& exact,
                                     unsigned src,
                                     unsigned K,
                                     KmAgg& agg) {
  const unsigned n = (unsigned)LB.size();
  if (n < 2) return;
  K = std::min<unsigned>(K, n - 1);
  if (K == 0) return;

  // Sorted UB[] and LB[] over all u != src (sanitized like the original).
  vector<double> sUB, sLB;
  sUB.reserve(n - 1);
  sLB.reserve(n - 1);
  for (unsigned u = 0; u < n; ++u) {
    if (u == src) continue;
    double ubu = UB[u];
    if (!km_isfinite(ubu)) ubu = numeric_limits<double>::infinity();
    double lbu = LB[u];
    if (!km_isfinite(lbu) || lbu < 0) lbu = 0.0;
    sUB.push_back(ubu);
    sLB.push_back(lbu);
  }
  sort(sUB.begin(), sUB.end());
  sort(sLB.begin(), sLB.end());

  vector<unsigned> truth = km_truth_topk(exact, src, K);
  unordered_set<unsigned> truth_set(truth.begin(), truth.end());

  for (unsigned v = 0; v < n; ++v) {
    if (v == src) continue;
    double lbv = LB[v], ubv = UB[v];
    if (!km_isfinite(lbv) || lbv < 0) lbv = 0.0;
    if (!km_isfinite(ubv)) ubv = numeric_limits<double>::infinity();

    // definitely_closer = #{u != src : UB[u] < lbv}
    unsigned long long definitely_closer =
        (unsigned long long)(lower_bound(sUB.begin(), sUB.end(), lbv) - sUB.begin());
    bool out_cert = (definitely_closer >= K);

    // possibly_closer = #{u != src, u != v : LB[u] < ubv}
    unsigned long long possibly_closer =
        (unsigned long long)(lower_bound(sLB.begin(), sLB.end(), ubv) - sLB.begin());
    if (lbv < ubv) {            // u==v was counted above; remove it
      if (possibly_closer > 0) possibly_closer--;
    }
    bool in_cert = (possibly_closer <= (unsigned long long)(K - 1));

    if (in_cert) {
      agg.definite_in++;
      agg.definite_in_total++;
      if (truth_set.count(v)) agg.definite_in_correct++;
    } else if (out_cert) {
      agg.definite_out++;
    } else {
      agg.unresolved++;
    }
  }
}


// Candidate-limited version of km_eval_topk_fast. The certification universe is
// the supplied candidate pool rather than all nodes.
static inline void km_eval_topk_fast_candidates(const vector<double>& LB,
                                                const vector<double>& UB,
                                                const vector<double>& exact,
                                                unsigned src,
                                                unsigned K,
                                                const vector<unsigned>& candidates,
                                                KmAgg& agg) {
  if (candidates.empty()) return;
  K = std::min<unsigned>(K, (unsigned)candidates.size());
  if (K == 0) return;

  vector<double> sUB, sLB;
  sUB.reserve(candidates.size());
  sLB.reserve(candidates.size());
  for (unsigned u : candidates) {
    double ubu = UB[u];
    if (!km_isfinite(ubu)) ubu = numeric_limits<double>::infinity();
    double lbu = LB[u];
    if (!km_isfinite(lbu) || lbu < 0) lbu = 0.0;
    sUB.push_back(ubu);
    sLB.push_back(lbu);
  }
  sort(sUB.begin(), sUB.end());
  sort(sLB.begin(), sLB.end());

  vector<unsigned> truth = km_truth_topk_candidates(exact, candidates, K);
  unordered_set<unsigned> truth_set(truth.begin(), truth.end());

  for (unsigned v : candidates) {
    double lbv = LB[v], ubv = UB[v];
    if (!km_isfinite(lbv) || lbv < 0) lbv = 0.0;
    if (!km_isfinite(ubv)) ubv = numeric_limits<double>::infinity();

    unsigned long long definitely_closer =
        (unsigned long long)(lower_bound(sUB.begin(), sUB.end(), lbv) - sUB.begin());
    bool out_cert = (definitely_closer >= K);

    unsigned long long possibly_closer =
        (unsigned long long)(lower_bound(sLB.begin(), sLB.end(), ubv) - sLB.begin());
    if (lbv < ubv && possibly_closer > 0) possibly_closer--; // remove v itself
    bool in_cert = (possibly_closer <= (unsigned long long)(K - 1));

    if (in_cert) {
      agg.definite_in++;
      agg.definite_in_total++;
      if (truth_set.count(v)) agg.definite_in_correct++;
    } else if (out_cert) {
      agg.definite_out++;
    } else {
      agg.unresolved++;
    }
  }
}

static inline bool km_pair_resolved(const vector<double>& LB,
                                    const vector<double>& UB,
                                    unsigned a,
                                    unsigned b) {
  double lba = LB[a], lbb = LB[b];
  double uba = UB[a], ubb = UB[b];
  if (!km_isfinite(lba) || lba < 0) lba = 0.0;
  if (!km_isfinite(lbb) || lbb < 0) lbb = 0.0;
  if (!km_isfinite(uba)) uba = numeric_limits<double>::infinity();
  if (!km_isfinite(ubb)) ubb = numeric_limits<double>::infinity();
  return (uba < lbb) || (ubb < lba);
}

static inline void km_fill_lb_from_lookup(vector<double>& LB,
                                          const function<double(unsigned int, unsigned int)>& lookup_fn,
                                          unsigned src) {
  const unsigned n = (unsigned)LB.size();
  for (unsigned v = 0; v < n; ++v) {
    if (v == src) {
      LB[v] = 0.0;
      continue;
    }
    double x = lookup_fn(src, v);
    if (!km_isfinite(x) || x < 0) x = 0.0;
    LB[v] = x;
  }
}


// -----------------------------------------------------------------------------
// Persistent KNN bound-row cache.
//
// Cache granularity is row-based rather than all-or-nothing:
//   UB/source q: shortest-path upper-bound row UB(q, *)
//   LB/method/source q: lower-bound row LB_method(q, *)
// A later KNN run with more query sources reuses previously materialized rows
// and appends only missing rows. If --knn-all N is used with N >= |V|, this
// materializes complete all-pairs UB/LB rows for the configured methods.
// -----------------------------------------------------------------------------
static std::string cache_sanitize_token(std::string s) {
  for (char &c : s) {
    unsigned char uc = static_cast<unsigned char>(c);
    if (!std::isalnum(uc) && c != '-' && c != '_' && c != '.') c = '_';
  }
  if (s.empty()) s = "empty";
  return s;
}

struct KnnBoundCacheStats {
  std::atomic<unsigned long long> ub_hits{0}, ub_misses{0};
  std::atomic<unsigned long long> lb_hits{0}, lb_misses{0};
  std::atomic<unsigned long long> ub_saved{0}, lb_saved{0};
};

static bool knn_bound_cache_enabled(const std::string &dir, const std::string &key) {
  return !dir.empty() && !key.empty();
}

static std::string knn_bound_cache_method_dir(const std::string &root,
                                              const std::string &key,
                                              const std::string &method_tag) {
  namespace fs = std::filesystem;
  fs::path p(root);
  p /= cache_sanitize_token(key);
  p /= cache_sanitize_token(method_tag);
  return p.string();
}

static std::string knn_bound_cache_row_path(const std::string &root,
                                            const std::string &key,
                                            const std::string &method_tag,
                                            unsigned src) {
  namespace fs = std::filesystem;
  fs::path p(knn_bound_cache_method_dir(root, key, method_tag));
  p /= ("src_" + std::to_string(src) + ".bin");
  return p.string();
}

static void knn_bound_cache_prepare_dir(const std::string &root,
                                        const std::string &key,
                                        const std::string &method_tag) {
  if (!knn_bound_cache_enabled(root, key)) return;
  namespace fs = std::filesystem;
  try {
    fs::create_directories(knn_bound_cache_method_dir(root, key, method_tag));
  } catch (const std::exception &e) {
    std::cerr << "[KNN-CACHE] WARN could not create cache directory for "
              << method_tag << ": " << e.what() << std::endl;
  }
}

static bool knn_bound_cache_load_row(const std::string &path,
                                     unsigned expected_n,
                                     unsigned expected_src,
                                     std::vector<double> &row) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  char magic[8] = {0};
  uint32_t version = 0;
  uint64_t n64 = 0, src64 = 0, count64 = 0;
  in.read(magic, sizeof(magic));
  in.read(reinterpret_cast<char*>(&version), sizeof(version));
  in.read(reinterpret_cast<char*>(&n64), sizeof(n64));
  in.read(reinterpret_cast<char*>(&src64), sizeof(src64));
  in.read(reinterpret_cast<char*>(&count64), sizeof(count64));
  if (!in || std::string(magic, magic + 7) != "KNNROW1" || version != 1 ||
      n64 != expected_n || src64 != expected_src || count64 != expected_n) {
    return false;
  }
  row.assign(expected_n, 0.0);
  in.read(reinterpret_cast<char*>(row.data()), sizeof(double) * (size_t)expected_n);
  if (!in) {
    row.clear();
    return false;
  }
  return true;
}

static bool knn_bound_cache_save_row(const std::string &path,
                                     unsigned n,
                                     unsigned src,
                                     const std::vector<double> &row) {
  if (row.size() != n) return false;
  namespace fs = std::filesystem;
  try {
    fs::create_directories(fs::path(path).parent_path());
  } catch (...) {
    return false;
  }
  const std::string tmp = path + ".tmp." + std::to_string(src) + "." +
                          std::to_string((unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count());
  {
    std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    char magic[8] = {'K','N','N','R','O','W','1','\0'};
    uint32_t version = 1;
    uint64_t n64 = n, src64 = src, count64 = n;
    out.write(magic, sizeof(magic));
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&n64), sizeof(n64));
    out.write(reinterpret_cast<const char*>(&src64), sizeof(src64));
    out.write(reinterpret_cast<const char*>(&count64), sizeof(count64));
    out.write(reinterpret_cast<const char*>(row.data()), sizeof(double) * (size_t)n);
    if (!out) return false;
  }
  try {
    std::error_code ec;
    fs::rename(tmp, path, ec);
    if (ec) {
      fs::remove(path, ec);
      ec.clear();
      fs::rename(tmp, path, ec);
      if (ec) {
        fs::remove(tmp);
        return false;
      }
    }
  } catch (...) {
    return false;
  }
  return true;
}

static bool knn_cache_load_or_compute_ub_row(
    vector<list<pair<unsigned, double>>*>* adj_lst,
    unsigned nodes,
    unsigned src,
    std::vector<double> &exact,
    std::vector<double> &UB,
    const std::string &cache_dir,
    const std::string &cache_key,
    KnnBoundCacheStats *stats) {
  bool use_cache = knn_bound_cache_enabled(cache_dir, cache_key);
  std::string path;
  if (use_cache) {
    path = knn_bound_cache_row_path(cache_dir, cache_key, "UB", src);
    if (knn_bound_cache_load_row(path, nodes, src, UB)) {
      exact = UB;
      if (stats) stats->ub_hits++;
      return true;
    }
    if (stats) stats->ub_misses++;
  }

  auto sp = DijkstraELM(adj_lst, src);
  assert(sp.first && "DijkstraELM returned null dist vector");
  exact.assign(sp.first->begin(), sp.first->end());
  delete sp.first;
  if (sp.second) delete sp.second;
  UB = exact;
  for (unsigned v = 0; v < nodes; ++v) {
    if (!km_isfinite(UB[v])) UB[v] = numeric_limits<double>::infinity();
    if (UB[v] < 0) UB[v] = 0.0;
  }
  exact = UB;
  if (use_cache && knn_bound_cache_save_row(path, nodes, src, UB)) {
    if (stats) stats->ub_saved++;
  }
  return false;
}

static bool knn_cache_load_or_compute_lb_row(
    std::vector<double> &LB,
    const function<double(unsigned int, unsigned int)> &lookup_fn,
    unsigned nodes,
    unsigned src,
    const std::string &method_tag,
    const std::string &cache_dir,
    const std::string &cache_key,
    KnnBoundCacheStats *stats) {
  bool use_cache = knn_bound_cache_enabled(cache_dir, cache_key);
  std::string path;
  if (use_cache) {
    path = knn_bound_cache_row_path(cache_dir, cache_key, method_tag, src);
    if (knn_bound_cache_load_row(path, nodes, src, LB)) {
      if (stats) stats->lb_hits++;
      return true;
    }
    if (stats) stats->lb_misses++;
  }
  km_fill_lb_from_lookup(LB, lookup_fn, src);
  if (use_cache && knn_bound_cache_save_row(path, nodes, src, LB)) {
    if (stats) stats->lb_saved++;
  }
  return false;
}

static void run_knn_metrics_combined(
    vector<list<pair<unsigned, double>>*>* adj_lst,
    TriSearch& tri,
    const vector<CurveModel>& models,
    unsigned nodes,
    unsigned K,
    unsigned num_queries,
    unsigned pair_samples_per_query,
    unsigned seed,
    bool full_mode = false,
    unsigned candidate_limit = 0,
    bool cert_select_mode = false,
    const std::string &knn_bound_cache_dir = "",
    const std::string &knn_bound_cache_key = "") {
  if (num_queries == 0 || nodes < 2) return;

  KnnBoundCacheStats knn_cache_stats;
  const bool use_knn_cache = knn_bound_cache_enabled(knn_bound_cache_dir, knn_bound_cache_key);
  if (use_knn_cache) {
    knn_bound_cache_prepare_dir(knn_bound_cache_dir, knn_bound_cache_key, "UB");
    knn_bound_cache_prepare_dir(knn_bound_cache_dir, knn_bound_cache_key, "Tri");
    for (const auto &m : models)
      knn_bound_cache_prepare_dir(knn_bound_cache_dir, knn_bound_cache_key, m.name);
    cout << "[KNN-CACHE] ENABLED dir=" << knn_bound_cache_dir
         << " key=" << knn_bound_cache_key
         << " candidate_limit=" << candidate_limit
         << " (row-granular full UB/LB rows reused for full or candidate-limited KNN)" << endl;
  }

  // SPLUB exact-LB baseline (Reviewer R1.O2 / R3.D6): build all-pairs SP once,
  // gated to small graphs. When enabled, an extra "SPLUB" method row is
  // reported alongside Tri / models / Combined in both KNN certification and
  // PairResolved (ERC).
  const bool do_splub = KNN_SPLUB_ENABLE && (nodes <= KNN_SPLUB_MAX_N);
  vector<vector<double>> SP_all;
  if (do_splub) {
    cout << "[KNN-SPLUB] building all-pairs shortest paths for the exact SPLUB "
            "baseline (n=" << nodes << ")..." << endl;
    build_all_pairs_sp(adj_lst, nodes, SP_all);
    cout << "[KNN-SPLUB] ready: SPLUB reported as an extra method row "
            "(exact LB upper-envelope for ERC and kNN certification)." << endl;
  } else if (KNN_SPLUB_ENABLE) {
    cout << "[KNN-SPLUB] skipped: n=" << nodes << " > KNN_SPLUB_MAX_N="
         << KNN_SPLUB_MAX_N << " (exact SPLUB over all pairs is infeasible at "
            "this scale; raise --knn-splub-max-n to force)." << endl;
  } else {
    cout << "[KNN-SPLUB] disabled (--no-knn-splub)." << endl;
  }

  // -------------------------------------------------------------------------
  // FULL-CANDIDATE MODE: a USER-PROVIDED number of query nodes (num_queries),
  // each evaluated against ALL N candidate nodes (no pair sampling). The
  // exhaustive TopK certification runs per query via the O(n log n) fast path,
  // parallelized across queries. This surfaces non-zero in/out certification
  // COUNTS (vs the ~0 from the tiny sampled mode) while keeping the query
  // count controllable so the run is affordable. Pass num_queries >= N to use
  // every node. Query nodes are a seeded, distinct sample for reproducibility.
  // PairResolved is not computed here (exhaustive all-pairs is O(n^3)); the
  // exhaustive TopK over all candidates IS the unsampled evaluation.
  // -------------------------------------------------------------------------
  if (full_mode) {
    const unsigned Q = std::min(num_queries, nodes);
    const unsigned candidates_per_query =
        (candidate_limit == 0 || candidate_limit >= nodes - 1) ? (nodes - 1) : candidate_limit;
    const bool candidate_limited = (candidates_per_query < nodes - 1);
    const size_t M = models.size();

    // Select query nodes. Default: if Q>=N use every node, otherwise a seeded
    // random subset. With --knn-cert-select, score every possible query by the
    // separation between decreasing lower bounds and increasing upper bounds:
    //   tau_k(q) = kth-smallest UB(q,*)
    //   score(q) = #{v: LB_combined(q,v) > tau_k(q)} / (C-k)
    // and evaluate the top-Q most certifiable query nodes. This is intended
    // only for high-dimensional/SIFT diagnostics where random query nodes often
    // have zero useful lower-bound separation.
    vector<unsigned> query_nodes;
    query_nodes.reserve(Q);

    if (cert_select_mode) {
      cout << "[KNN-SELECT] mode=CERTIFIABLE_QUERY_SELECT scanning=" << nodes
           << " candidate_queries to select top " << Q
           << " using Combined(max) LB, candidates_per_query=" << candidates_per_query
           << " K=" << K << endl;

      const unsigned Kc_scan = std::min<unsigned>(K, candidates_per_query);
      if (Kc_scan == 0 || Kc_scan >= candidates_per_query) {
        cout << "[KNN-SELECT] WARNING: K >= candidates_per_query, query selection is degenerate; "
             << "all candidates are inside top-K and certified-out denominator is zero.\n";
      }

      vector<KnnQuerySelectScore> scores(nodes);
      std::atomic<unsigned long long> select_done{0};
      int select_threads = 1;
      #ifdef _OPENMP
      select_threads = omp_get_max_threads();
      #endif

      #pragma omp parallel
      {
        vector<double> exact(nodes), UB(nodes), LB_tri(nodes), LB_comb(nodes);
        vector<vector<double>> modelLB(M, vector<double>(nodes));

        #pragma omp for schedule(dynamic, 16)
        for (long long src_ll = 0; src_ll < (long long)nodes; ++src_ll) {
          unsigned src = (unsigned)src_ll;
          knn_cache_load_or_compute_ub_row(adj_lst, nodes, src, exact, UB,
                                           knn_bound_cache_dir, knn_bound_cache_key,
                                           use_knn_cache ? &knn_cache_stats : nullptr);
          knn_cache_load_or_compute_lb_row(LB_tri,
              [&](unsigned u, unsigned v){ return tri.lookup(u, v); },
              nodes, src, "Tri", knn_bound_cache_dir, knn_bound_cache_key,
              use_knn_cache ? &knn_cache_stats : nullptr);
          for (size_t mi = 0; mi < M; ++mi) {
            knn_cache_load_or_compute_lb_row(modelLB[mi], models[mi].lookup,
                nodes, src, models[mi].name, knn_bound_cache_dir, knn_bound_cache_key,
                use_knn_cache ? &knn_cache_stats : nullptr);
          }

          for (unsigned v = 0; v < nodes; ++v) {
            double x = LB_tri[v];
            for (size_t mi = 0; mi < M; ++mi) x = std::max(x, modelLB[mi][v]);
            if (!km_isfinite(x) || x < 0) x = 0.0;
            LB_comb[v] = x;
          }

          vector<unsigned> cand = candidate_limited
              ? km_candidate_pool(nodes, src, candidates_per_query, seed + 99173ULL)
              : km_candidate_pool(nodes, src, 0, seed + 99173ULL);
          const unsigned C = (unsigned)cand.size();
          const unsigned Kc = std::min<unsigned>(K, C);

          KnnQuerySelectScore sc;
          sc.src = src;
          if (C > 0 && Kc > 0) {
            vector<double> sUB;
            sUB.reserve(C);
            for (unsigned v : cand) {
              double u = UB[v];
              if (!km_isfinite(u)) u = numeric_limits<double>::infinity();
              sUB.push_back(u);
            }
            nth_element(sUB.begin(), sUB.begin() + (Kc - 1), sUB.end());
            sc.tau_k = sUB[Kc - 1];

            unsigned long long cnt = 0;
            vector<double> sLB;
            sLB.reserve(C);
            for (unsigned v : cand) {
              double l = LB_comb[v];
              if (!km_isfinite(l) || l < 0) l = 0.0;
              sLB.push_back(l);
              if (l > sc.tau_k) cnt++;
            }
            sc.certified_out_like = cnt;
            sc.denom = (C > Kc) ? (unsigned long long)(C - Kc) : 0ULL;
            sc.score_pct = (sc.denom == 0) ? 0.0 : 100.0 * (double)cnt / (double)sc.denom;
            if (sc.denom > 0 && sc.denom <= sLB.size()) {
              // Boundary of the decreasing-LB list at rank C-K. If this value
              // exceeds tau_k, then every non-topK candidate can be separated.
              nth_element(sLB.begin(), sLB.begin() + (size_t)(sc.denom - 1), sLB.end(), greater<double>());
              sc.lb_boundary = sLB[(size_t)(sc.denom - 1)];
              sc.gap = sc.lb_boundary - sc.tau_k;
            }
          }
          scores[src] = sc;

          unsigned long long dsel = ++select_done;
          if ((dsel % 1000ULL) == 0ULL || dsel == nodes) {
            #pragma omp critical
            { cout << "[KNN-SELECT progress] " << dsel << "/" << nodes << " query nodes scored\n"; cout.flush(); }
          }
        }
      }

      sort(scores.begin(), scores.end(), [](const KnnQuerySelectScore &a, const KnnQuerySelectScore &b) {
        if (a.certified_out_like != b.certified_out_like) return a.certified_out_like > b.certified_out_like;
        if (a.score_pct != b.score_pct) return a.score_pct > b.score_pct;
        if (a.gap != b.gap) return a.gap > b.gap;
        return a.src < b.src;
      });

      query_nodes.clear();
      for (unsigned i = 0; i < Q && i < scores.size(); ++i) query_nodes.push_back(scores[i].src);

      unsigned long long sum_cert = 0, sum_denom = 0;
      for (unsigned i = 0; i < Q && i < scores.size(); ++i) {
        sum_cert += scores[i].certified_out_like;
        sum_denom += scores[i].denom;
      }
      cout << "[KNN-SELECT SUMMARY] selected_queries=" << query_nodes.size()
           << " avg_selected_out_score=" << ((sum_denom == 0) ? 0.0 : 100.0 * (double)sum_cert / (double)sum_denom)
           << "% total_selected_certified_out_like=" << sum_cert
           << " / " << sum_denom << endl;
      cout << "[KNN-SELECT TOP]" << endl;
      for (unsigned i = 0; i < std::min<unsigned>(10, (unsigned)scores.size()) && i < Q; ++i) {
        cout << "  rank=" << (i + 1)
             << " src=" << scores[i].src
             << " out_like=" << scores[i].certified_out_like
             << "/" << scores[i].denom
             << " score=" << scores[i].score_pct
             << "% tau_k=" << scores[i].tau_k
             << " lb_boundary=" << scores[i].lb_boundary
             << " gap=" << scores[i].gap << endl;
      }
    } else if (Q >= nodes) {
      query_nodes.resize(nodes);
      std::iota(query_nodes.begin(), query_nodes.end(), 0u);
    } else {
      vector<unsigned> all(nodes);
      std::iota(all.begin(), all.end(), 0u);
      std::mt19937 qrng(seed);
      std::shuffle(all.begin(), all.end(), qrng);
      query_nodes.assign(all.begin(), all.begin() + Q);
    }

    cout << "[KNN-FULL] mode="
         << (cert_select_mode ? "CERT_SELECTED" : (candidate_limited ? "CANDIDATE_LIMITED" : "FULL_CANDIDATE"))
         << " queries=" << query_nodes.size()
         << (cert_select_mode ? " (top certifiable query nodes)" : (Q >= nodes ? " (ALL nodes)" : " (seeded random subset)"))
         << " candidates_per_query=" << candidates_per_query
         << " K=" << K
         << " (no pair sampling; exhaustive TopK over "
         << (candidate_limited ? "candidate pool" : "all candidates") << " per query)"
         << endl;

    // Per-thread accumulators to avoid contention; merged after the loop.
    int max_threads = 1;
    #ifdef _OPENMP
    max_threads = omp_get_max_threads();
    #endif
    vector<vector<KmAgg>> tl_modelAgg(max_threads, vector<KmAgg>(M));
    vector<KmAgg> tl_triAgg(max_threads), tl_combAgg(max_threads);
    vector<KmAgg> tl_splubAgg(max_threads);

    std::atomic<unsigned long long> done{0};

    #pragma omp parallel
    {
      int tid = 0;
      #ifdef _OPENMP
      tid = omp_get_thread_num();
      #endif
      // Thread-local scratch buffers (reused across this thread's queries).
      vector<double> exact(nodes), UB(nodes), LB_tri(nodes), LB_comb(nodes);
      vector<vector<double>> modelLB(M, vector<double>(nodes));
      vector<double> LB_splub(do_splub ? nodes : 0);

      #pragma omp for schedule(dynamic, 64)
      for (long long qi = 0; qi < (long long)Q; ++qi) {
        unsigned src = query_nodes[(size_t)qi];

        knn_cache_load_or_compute_ub_row(adj_lst, nodes, src, exact, UB,
                                         knn_bound_cache_dir, knn_bound_cache_key,
                                         use_knn_cache ? &knn_cache_stats : nullptr);

        knn_cache_load_or_compute_lb_row(LB_tri,
            [&](unsigned u, unsigned v){ return tri.lookup(u, v); },
            nodes, src, "Tri", knn_bound_cache_dir, knn_bound_cache_key,
            use_knn_cache ? &knn_cache_stats : nullptr);
        for (size_t mi = 0; mi < M; ++mi)
          knn_cache_load_or_compute_lb_row(modelLB[mi], models[mi].lookup,
              nodes, src, models[mi].name, knn_bound_cache_dir, knn_bound_cache_key,
              use_knn_cache ? &knn_cache_stats : nullptr);

        for (unsigned v = 0; v < nodes; ++v) {
          double x = LB_tri[v];
          for (size_t mi = 0; mi < M; ++mi) x = std::max(x, modelLB[mi][v]);
          LB_comb[v] = x;
        }

        if (candidate_limited) {
          vector<unsigned> cand = km_candidate_pool(nodes, src, candidates_per_query, seed + 99173ULL);
          km_eval_topk_fast_candidates(LB_tri, UB, exact, src, K, cand, tl_triAgg[tid]);
          for (size_t mi = 0; mi < M; ++mi)
            km_eval_topk_fast_candidates(modelLB[mi], UB, exact, src, K, cand, tl_modelAgg[tid][mi]);
          km_eval_topk_fast_candidates(LB_comb, UB, exact, src, K, cand, tl_combAgg[tid]);
          if (do_splub) {
            splub_lb_row(adj_lst, nodes, SP_all, src, LB_splub);
            km_eval_topk_fast_candidates(LB_splub, UB, exact, src, K, cand, tl_splubAgg[tid]);
          }
        } else {
          km_eval_topk_fast(LB_tri, UB, exact, src, K, tl_triAgg[tid]);
          for (size_t mi = 0; mi < M; ++mi)
            km_eval_topk_fast(modelLB[mi], UB, exact, src, K, tl_modelAgg[tid][mi]);
          km_eval_topk_fast(LB_comb, UB, exact, src, K, tl_combAgg[tid]);
          if (do_splub) {
            splub_lb_row(adj_lst, nodes, SP_all, src, LB_splub);
            km_eval_topk_fast(LB_splub, UB, exact, src, K, tl_splubAgg[tid]);
          }
        }

        unsigned long long d2 = ++done;
        if ((d2 % 5000ULL) == 0ULL) {
          #pragma omp critical
          { cout << "[KNN-FULL progress] " << d2 << "/" << Q << " queries\n"; cout.flush(); }
        }
      }
    }

    // Merge thread-local accumulators.
    auto merge = [](KmAgg& dst, const KmAgg& s) {
      dst.definite_in += s.definite_in;
      dst.definite_out += s.definite_out;
      dst.unresolved += s.unresolved;
      dst.definite_in_total += s.definite_in_total;
      dst.definite_in_correct += s.definite_in_correct;
      dst.pair_resolved += s.pair_resolved;
      dst.pair_tried += s.pair_tried;
    };
    vector<KmAgg> modelAgg(M);
    KmAgg triAgg, combAgg, splubAgg;
    for (int t = 0; t < max_threads; ++t) {
      merge(triAgg, tl_triAgg[t]);
      merge(combAgg, tl_combAgg[t]);
      if (do_splub) merge(splubAgg, tl_splubAgg[t]);
      for (size_t mi = 0; mi < M; ++mi) merge(modelAgg[mi], tl_modelAgg[t][mi]);
    }

    auto pct = [](unsigned long long a, unsigned long long b) -> double {
      return (b == 0) ? 0.0 : 100.0 * (double)a / (double)b;
    };
    // recall denominator: each query contributes min(K, candidates_per_query)
    // true neighbors inside the evaluated candidate universe.
    const unsigned long long Kc = std::min<unsigned>(K, candidates_per_query);
    const unsigned long long recall_denom = (unsigned long long)Q * Kc;
    const unsigned long long eval_pairs = (unsigned long long)Q * (unsigned long long)candidates_per_query;

    auto print_full = [&](const char* name, const KmAgg& A) {
      cout << "\n===== KNN METRICS [" << name << "] ("
           << (candidate_limited ? "CANDIDATE-LIMITED" : "FULL-CANDIDATE")
           << ": " << Q << " queries x " << candidates_per_query
           << " candidates, K=" << K << ") =====\n";
      cout << "TopK definite_in=" << A.definite_in
           << " definite_out=" << A.definite_out
           << " unresolved=" << A.unresolved << "\n";
      // Recall@K = fraction of all true top-K memberships we certified IN.
      // Numerator is definite_in_correct (certified AND truly in top-K), so
      // this is bounded by 100% even if bounds were invalid. On a valid
      // metric, precision=1.0 so definite_in_correct == definite_in.
      cout << "Recall@K (definite_in_correct / (Q*K)) = "
           << pct(A.definite_in_correct, recall_denom)
           << "%  (correct_in=" << A.definite_in_correct
           << " / Q*K=" << recall_denom << ")\n";
      // Certified-out coverage: fraction of (query,candidate) pairs proven
      // NOT in top-K, out of all evaluated non-self pairs.
      cout << "CertifiedOut% (definite_out / (Q*candidates_per_query)) = "
           << pct(A.definite_out, eval_pairs)
           << "%  (definite_out=" << A.definite_out
           << " / " << eval_pairs << ")\n";
      cout << "Definite_in precision=";
      if (A.definite_in_total == 0) cout << 0;
      else cout << std::setprecision(6)
                << (double)A.definite_in_correct / (double)A.definite_in_total
                << " (correct=" << A.definite_in_correct
                << "/" << A.definite_in_total << ")";
      cout << "\n";
    };
    print_full("Tri", triAgg);
    for (size_t mi = 0; mi < M; ++mi) print_full(models[mi].name.c_str(), modelAgg[mi]);
    print_full("Combined(max)", combAgg);
    if (do_splub) print_full("SPLUB", splubAgg);
    if (use_knn_cache) {
      cout << "\n[KNN-CACHE SUMMARY] UB hits=" << knn_cache_stats.ub_hits.load()
           << " misses=" << knn_cache_stats.ub_misses.load()
           << " saved=" << knn_cache_stats.ub_saved.load()
           << " | LB hits=" << knn_cache_stats.lb_hits.load()
           << " misses=" << knn_cache_stats.lb_misses.load()
           << " saved=" << knn_cache_stats.lb_saved.load() << "\n";
    }
    return;
  }

  // -------------------------------------------------------------------------
  // SAMPLED MODE (default; London/Porto): unchanged from before.
  // -------------------------------------------------------------------------

  std::mt19937 rng(seed);
  std::uniform_int_distribution<unsigned> unif_node(0, nodes - 1);

  vector<KmAgg> modelAgg(models.size());
  KmAgg triAgg, combAgg, splubAgg;

  vector<double> exact(nodes, 0.0), UB(nodes, 0.0), LB_tri(nodes, 0.0), LB_comb(nodes, 0.0);
  vector<vector<double>> modelLB(models.size(), vector<double>(nodes, 0.0));
  vector<double> LB_splub(do_splub ? nodes : 0, 0.0);

  auto pct = [](unsigned long long a, unsigned long long b) -> double {
    if (b == 0) return 0.0;
    return 100.0 * (double)a / (double)b;
  };

  for (unsigned qi = 0; qi < num_queries; ++qi) {
    unsigned src = unif_node(rng);

    knn_cache_load_or_compute_ub_row(adj_lst, nodes, src, exact, UB,
                                   knn_bound_cache_dir, knn_bound_cache_key,
                                   use_knn_cache ? &knn_cache_stats : nullptr);

    knn_cache_load_or_compute_lb_row(LB_tri,
        [&](unsigned u, unsigned v){ return tri.lookup(u, v); },
        nodes, src, "Tri", knn_bound_cache_dir, knn_bound_cache_key,
        use_knn_cache ? &knn_cache_stats : nullptr);
    for (size_t mi = 0; mi < models.size(); ++mi) {
      knn_cache_load_or_compute_lb_row(modelLB[mi], models[mi].lookup,
          nodes, src, models[mi].name, knn_bound_cache_dir, knn_bound_cache_key,
          use_knn_cache ? &knn_cache_stats : nullptr);
    }

    for (unsigned v = 0; v < nodes; ++v) {
      double x = LB_tri[v];
      for (size_t mi = 0; mi < models.size(); ++mi) x = std::max(x, modelLB[mi][v]);
      LB_comb[v] = x;
    }

    if (do_splub) splub_lb_row(adj_lst, nodes, SP_all, src, LB_splub);

    auto sample_pair = [&]() -> pair<unsigned, unsigned> {
      unsigned a = unif_node(rng);
      unsigned b = unif_node(rng);
      while (b == a) b = unif_node(rng);
      return {a, b};
    };

    unsigned long long tri_resolved = 0, tri_tried = 0;
    for (unsigned s = 0; s < pair_samples_per_query; ++s) {
      auto [a, b] = sample_pair();
      if (a == src || b == src) { --s; continue; }
      tri_tried++;
      if (km_pair_resolved(LB_tri, UB, a, b)) tri_resolved++;
    }
    triAgg.pair_tried += tri_tried;
    triAgg.pair_resolved += tri_resolved;

    for (size_t mi = 0; mi < models.size(); ++mi) {
      unsigned long long resolved = 0, tried = 0;
      for (unsigned s = 0; s < pair_samples_per_query; ++s) {
        auto [a, b] = sample_pair();
        if (a == src || b == src) { --s; continue; }
        tried++;
        if (km_pair_resolved(modelLB[mi], UB, a, b)) resolved++;
      }
      modelAgg[mi].pair_tried += tried;
      modelAgg[mi].pair_resolved += resolved;
    }

    unsigned long long comb_resolved = 0, comb_tried = 0;
    for (unsigned s = 0; s < pair_samples_per_query; ++s) {
      auto [a, b] = sample_pair();
      if (a == src || b == src) { --s; continue; }
      comb_tried++;
      if (km_pair_resolved(LB_comb, UB, a, b)) comb_resolved++;
    }
    combAgg.pair_tried += comb_tried;
    combAgg.pair_resolved += comb_resolved;

    if (do_splub) {
      unsigned long long splub_resolved = 0, splub_tried = 0;
      for (unsigned s = 0; s < pair_samples_per_query; ++s) {
        auto [a, b] = sample_pair();
        if (a == src || b == src) { --s; continue; }
        splub_tried++;
        if (km_pair_resolved(LB_splub, UB, a, b)) splub_resolved++;
      }
      splubAgg.pair_tried += splub_tried;
      splubAgg.pair_resolved += splub_resolved;
    }

    if (candidate_limit > 0 && candidate_limit < nodes - 1) {
      vector<unsigned> cand = km_candidate_pool(nodes, src, candidate_limit, seed + 99173ULL);
      km_eval_topk_fast_candidates(LB_tri, UB, exact, src, K, cand, triAgg);
      for (size_t mi = 0; mi < models.size(); ++mi) {
        km_eval_topk_fast_candidates(modelLB[mi], UB, exact, src, K, cand, modelAgg[mi]);
      }
      km_eval_topk_fast_candidates(LB_comb, UB, exact, src, K, cand, combAgg);
      if (do_splub)
        km_eval_topk_fast_candidates(LB_splub, UB, exact, src, K, cand, splubAgg);
    } else {
      km_eval_topk(LB_tri, UB, exact, src, K, triAgg);
      for (size_t mi = 0; mi < models.size(); ++mi) {
        km_eval_topk(modelLB[mi], UB, exact, src, K, modelAgg[mi]);
      }
      km_eval_topk(LB_comb, UB, exact, src, K, combAgg);
      if (do_splub)
        km_eval_topk(LB_splub, UB, exact, src, K, splubAgg);
    }

    cout << "\n[KNN-METRICS] Query " << (qi + 1) << "/" << num_queries << " src=" << src << "\n";
    cout << "  PairResolved (Tri) = " << pct(tri_resolved, tri_tried) << "%\n";
    for (size_t mi = 0; mi < models.size(); ++mi) {
      unsigned long long r = modelAgg[mi].pair_resolved;
      unsigned long long t = modelAgg[mi].pair_tried;
      (void)r; (void)t;
    }
  }

  auto print_summary = [&](const char* name, const KmAgg& A) {
    cout << "\n===== KNN METRICS [" << name << "] =====\n";
    cout << "Queries: " << num_queries << "  K=" << K
         << "  PairSamplesPerQuery=" << pair_samples_per_query << "\n";
    cout << "PairResolved%: " << pct(A.pair_resolved, A.pair_tried)
         << " (resolved=" << A.pair_resolved << " tried=" << A.pair_tried << ")\n";
    cout << "TopK definite_in=" << A.definite_in
         << " definite_out=" << A.definite_out
         << " unresolved=" << A.unresolved << "\n";
    cout << "Definite_in precision=";
    if (A.definite_in_total == 0) cout << 0;
    else cout << std::setprecision(6)
              << (double)A.definite_in_correct / (double)A.definite_in_total
              << " (correct=" << A.definite_in_correct << "/" << A.definite_in_total << ")";
    cout << "\n";
  };

  print_summary("Tri", triAgg);
  for (size_t mi = 0; mi < models.size(); ++mi) {
    print_summary(models[mi].name.c_str(), modelAgg[mi]);
  }
  print_summary("Combined(max)", combAgg);
  if (do_splub) print_summary("SPLUB", splubAgg);
  if (use_knn_cache) {
    cout << "\n[KNN-CACHE SUMMARY] UB hits=" << knn_cache_stats.ub_hits.load()
         << " misses=" << knn_cache_stats.ub_misses.load()
         << " saved=" << knn_cache_stats.ub_saved.load()
         << " | LB hits=" << knn_cache_stats.lb_hits.load()
         << " misses=" << knn_cache_stats.lb_misses.load()
         << " saved=" << knn_cache_stats.lb_saved.load() << "\n";
  }
}

// ============================================================================
// Single-run ablation component report (R3.D7). Self-contained: measures
// per-component BUILD and QUERY times (SPT, RMQ, HCA) and, for the SAME fixed
// query set, the lower-bound QUALITY (ARLB) with each component optimized vs
// unoptimized -- all in ONE process. Invoked from the top of
// multiple_sample_exec when --ablation-report is passed.
//
// Component semantics (precise, so the table is not misleading):
//   * SPT (Dijkstra trees): foundational. EVERY bound term -- even the
//     edge-folding-only no-SPT bound -- reads path lengths from the
//     landmark-root SPTs, so SPT construction is NOT removed by the no-SPT
//     quality ablation. We report its build cost as a fixed/foundational cost.
//     The no-SPT *quality* ablation drops the tree-path tightening terms at
//     query time (lower ARLB), which is what the ARLB_opt vs ARLB_unopt column
//     captures.
//   * RMQ (Euler-tour / sparse-table LCA preprocessing): needed ONLY by the
//     SPT-tightening terms' find_LCA. Optimized pays a one-time build and O(1)
//     per query; unoptimized (--ablation-no-rmq) pays ZERO build but a naive
//     O(depth) parent-walk per query. RMQ is loss-less: the bound VALUE is
//     identical with RMQ on or off (find_LCA == find_LCA_naive in value), so
//     RMQ's ARLB_opt == ARLB_unopt by construction; only time differs.
//   * HCA (§4.2.1 endpoint scoring): a SELECTION-time optimization. With HCA
//     on, landmark selection scores candidates with the HCA bound; with HCA
//     off it scores with the plain edge+SPT bound, producing a DIFFERENT
//     landmark set. We therefore time the two selections and evaluate both
//     selected models over the same fixed query set for ARLB_opt vs ARLB_unopt.
static void run_ablation_report(
    unsigned int init,
    vector<list<pair<unsigned int, double>> *> *adj_lst,
    map<pair<unsigned int, unsigned int>, double> *known_edges,
    unsigned int nodes,
    unsigned int k,
    unsigned int sampling_size,
    SplubCache *splub_cache,
    unsigned int n_queries) {

  auto OPT = [&](unsigned int a, unsigned int b) -> double {
    return splub_cache ? splub_cached(*splub_cache, a, b) : splub(adj_lst, a, b);
  };

  cout << "\n\n================ ABLATION COMPONENT REPORT ================\n";
  cout << "Config: n=" << nodes << " k=" << k
       << " sampling_size=" << sampling_size
       << " eval_queries(requested)=" << n_queries
       << " seed=" << init << "\n";
  cout << "(Build times = process-CPU us; selection scoring sums CPU across threads.\n"
       << " Query times are single-threaded, repeated for noise stability.)\n";

  const bool small_graph = (nodes < 2000);  // Exact selection is only feasible here.
  if (!small_graph)
    cout << "[ablation-report] n>=2000: skipping the two Exact methods "
            "(infeasible at this scale); reporting Sampling variants only.\n";

  // ---- Fixed evaluation query set (missing edges, OPT>0), shared by every
  //      method so ARLB and timings are apples-to-apples. ------------------
  cout << "[ablation-report] sampling fixed evaluation query set...\n" << std::flush;
  std::mt19937 rng(init ^ 0x5DEECE66Du);
  std::uniform_int_distribution<unsigned int> uni(0, nodes - 1);
  struct EvalQ { unsigned int u, v; double opt; };
  std::vector<EvalQ> queries;
  queries.reserve(n_queries);
  const uint64_t max_attempts = (uint64_t)n_queries * 200ULL + 10000ULL;
  uint64_t attempts = 0;
  while (queries.size() < n_queries && attempts < max_attempts) {
    ++attempts;
    unsigned int u = uni(rng), v = uni(rng);
    if (u == v) continue;
    auto key = make_pair(min(u, v), max(u, v));
    if (known_edges->find(key) != known_edges->end()) continue;
    double opt = OPT(key.first, key.second);
    if (opt <= 0.0) continue;
    queries.push_back({key.first, key.second, opt});
  }
  const size_t Q = queries.size();
  cout << "[ablation-report] eval queries collected=" << Q
       << " (attempts=" << attempts << ")\n" << std::flush;
  if (Q == 0) { cout << "[ablation-report] no valid OPT>0 queries; aborting.\n"; return; }

  cout << "\nColumns: Opt/Unopt ARLB | Select(s) SPT(s) LIFT(s) RMQ(s) Other(s) TotalBuild(s)"
       << " | EdgeFold(ms) LCA/RMQ(ms) SPTmath(ms) Other(ms) TotalQuery(ms)\n";

  // ---- Per-method measurement. select() runs the method's landmark
  //      selection on a fresh model; everything else is identical. ---------
  auto measure_method = [&](const std::string &name,
                            std::function<void(ELM_SPTree*)> select) {
    cout << "\n[ablation-report] >>> method: " << name << " ...\n" << std::flush;
    ELM_SPTree *m = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
    m->reset_build_timers();
    dijkstra_timing_reset();
    uint64_t s0 = cpu_time_us();
    select(m);                                   // selection (builds its SPTs via ensure_spt)
    uint64_t s1 = cpu_time_us();
    const uint64_t sel_cpu = s1 - s0;
    const uint64_t tree_sel = dijkstra_tree_us();
    const uint64_t lift_sel = dijkstra_lift_us();
    const uint64_t rmq_sel  = dijkstra_rmq_us();
    uint64_t p0 = cpu_time_us();
    m->preprocess_landmark_roots();
    uint64_t p1 = cpu_time_us();
    const uint64_t prep_cpu = p1 - p0;
    const uint64_t spt_tree_us = dijkstra_tree_us();
    const uint64_t lift_us     = dijkstra_lift_us();
    const uint64_t rmq_us      = dijkstra_rmq_us() + m->rmq_build_us();
    const long long hca_ll = (long long)sel_cpu - (long long)(tree_sel + lift_sel + rmq_sel);
    const uint64_t select_us = hca_ll > 0 ? (uint64_t)hca_ll : 0;
    const uint64_t total_build = sel_cpu + prep_cpu;
    const long long other_b_ll =
        (long long)total_build - (long long)(select_us + spt_tree_us + lift_us + rmq_us);
    const uint64_t other_build_us = other_b_ll > 0 ? (uint64_t)other_b_ll : 0;

    // Query passes with reps auto-scaled so the full pass runs ~150 ms.
    auto run_reps = [&](std::function<double(const EvalQ&)> fn, unsigned reps, double *arlb)->uint64_t {
      double last = 0.0; uint64_t t0 = cpu_time_us();
      for (unsigned r = 0; r < reps; ++r) { double s = 0.0; for (auto &q : queries) s += fn(q)/q.opt; last = s; }
      uint64_t dt = cpu_time_us() - t0; if (arlb) *arlb = last/(double)Q; return dt;
    };
    auto full_fn = [&](const EvalQ&q){ return m->lookup_ablate(q.u,q.v,false,false); };
    auto fold_fn = [&](const EvalQ&q){ return m->lookup_ablate(q.u,q.v,true,false); };
    auto nolca_fn= [&](const EvalQ&q){ return m->lookup_ablate(q.u,q.v,false,false,true); };
    auto naive_fn= [&](const EvalQ&q){ return m->lookup_ablate(q.u,q.v,false,true); };
    auto skel_fn = [&](const EvalQ&q){ return m->lookup_skeleton(q.u,q.v); };
    uint64_t cal = run_reps(full_fn, 1, nullptr);
    unsigned reps = 1;
    if (cal < 150000ULL) reps = (unsigned)std::min<uint64_t>(500ULL,
                                std::max<uint64_t>(1ULL, 150000ULL/std::max<uint64_t>(1ULL,cal)));
    double arlb_full=0, arlb_fold=0, arlb_naive=0;
    uint64_t T_skel=run_reps(skel_fn,reps,nullptr), T_fold=run_reps(fold_fn,reps,&arlb_fold),
             T_nolca=run_reps(nolca_fn,reps,nullptr), T_full=run_reps(full_fn,reps,&arlb_full),
             T_naive=run_reps(naive_fn,reps,&arlb_naive);
    const double t_full=(double)T_full/reps, t_fold=(double)T_fold/reps,
                 t_nolca=(double)T_nolca/reps, t_naive=(double)T_naive/reps, t_skel=(double)T_skel/reps;
    auto posd=[](double x){return x>0?x:0.0;};
    const double q_edgefold=posd(t_fold-t_skel), q_sptmath=posd(t_nolca-t_fold),
                 q_lca=posd(t_full-t_nolca), q_total=t_full;
    const double q_other=posd(q_total-(q_edgefold+q_sptmath+q_lca));
    auto S=[](double us){return us/1e6;}; auto MS=[](double us){return us/1e3;};
    auto fs=[&](double x,int p){std::ostringstream o;o<<std::fixed<<std::setprecision(p)<<x;return o.str();};

    // Human-readable row.
    cout << "  " << std::left << std::setw(20) << name
         << "  Opt=" << fs(arlb_full,4) << " Unopt=" << fs(arlb_fold,4)
         << " | Sel=" << fs(S(select_us),3) << " SPT=" << fs(S(spt_tree_us),3)
         << " LIFT=" << fs(S(lift_us),3) << " RMQ=" << fs(S(rmq_us),3)
         << " Oth=" << fs(S(other_build_us),3) << " Tot=" << fs(S(total_build),3)
         << " | EF=" << fs(MS(q_edgefold),3) << " LCA=" << fs(MS(q_lca),3)
         << " SPTm=" << fs(MS(q_sptmath),3) << " Oth=" << fs(MS(q_other),3)
         << " Tot=" << fs(MS(q_total),3)
         << "   [reps=" << reps << "; build " << (select_us+spt_tree_us+lift_us+rmq_us+other_build_us)
         << "=" << total_build << " us; query "
         << (uint64_t)llround(q_edgefold+q_sptmath+q_lca+q_other) << "=" << (uint64_t)llround(q_total) << " us]\n";

    // Machine-parseable row.
    cout << "ABLATION_ROW method=" << name
         << " n=" << nodes << " k=" << k << " sampling=" << sampling_size << " n_queries=" << Q
         << " opt_arlb=" << fs(arlb_full,6) << " unopt_arlb=" << fs(arlb_fold,6)
         << " select_hca_s=" << fs(S(select_us),6) << " spt_s=" << fs(S(spt_tree_us),6)
         << " lift_s=" << fs(S(lift_us),6) << " rmq_s=" << fs(S(rmq_us),6)
         << " other_build_s=" << fs(S(other_build_us),6) << " total_build_s=" << fs(S(total_build),6)
         << " edgefold_ms=" << fs(MS(q_edgefold),6) << " lca_rmq_ms=" << fs(MS(q_lca),6)
         << " spt_math_ms=" << fs(MS(q_sptmath),6) << " other_query_ms=" << fs(MS(q_other),6)
         << " total_query_ms=" << fs(MS(q_total),6) << "\n";
    delete m;
  };

  // ---- The four methods (Exact only on small graphs). --------------------
  const bool rel_default = (SGEL_OBJECTIVE_MODE != "raw");
  (void)rel_default;
  if (small_graph) {
    measure_method("ELM_Exact",        [&](ELM_SPTree*m){ m->get_exact_landmarks(); });
    measure_method("ELM_Exact_RelObj", [&](ELM_SPTree*m){ m->get_exact_landmarks_relative_parallel(); });
  }
  measure_method("ELM_Sampling",        [&](ELM_SPTree*m){ m->get_hca_landmarks_parallel(sampling_size, init+73, false); });
  measure_method("ELM_Sampling_RelObj", [&](ELM_SPTree*m){ m->get_hca_landmarks_parallel(sampling_size, init+73, true); });

  cout << "================ END ABLATION COMPONENT REPORT ================\n\n" << std::flush;
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
  unsigned int knn_pair_samples,
  bool knn_full_mode,
  unsigned int knn_candidate_limit,
  SplubCache *splub_cache = nullptr,
  ELM_SPTree *cached_elm_exact = nullptr,
  ELM_SPTree *cached_elm_exact_rel = nullptr,
  SPT_StochasticGreedy *cached_spt_exact = nullptr,
  EdgeLandMark *cached_h1 = nullptr,
  EdgeLandMark *cached_h2 = nullptr,
  SimpleGraph *cached_sg = nullptr,
  unsigned int stoch_cand_cap = 1000,
  // v15.1 / v16.2: one-time prebuild costs for the four shared methods.
  // When the cached pointer is reused inside this function the per-iter
  // build time is genuinely zero; reporting that 0 in the BUILD SUMMARY
  // makes shared methods look "free" to a reader, which is misleading.
  // We instead surface the original prebuild cost in build_us so the
  // table is self-contained. The prebuild_us=<n> lines emitted above the
  // sweep loop remain the authoritative one-time record; this just makes
  // BUILD SUMMARY agree with them.
  uint64_t cached_elm_exact_prebuild_us = 0,
  uint64_t cached_elm_exact_rel_prebuild_us = 0,
  uint64_t cached_spt_exact_prebuild_us = 0,
  uint64_t cached_h1_prebuild_us = 0,
  uint64_t cached_h2_prebuild_us = 0) {

    // Single-run ablation component report short-circuits the normal pipeline.
    // It produces the per-component (SPT/RMQ/HCA) build+query timing breakdown
    // and the optimized-vs-unoptimized ARLB for one fixed query set, then
    // returns so this (D,K,sampling) configuration emits only the report.
    if (ABLATION_REPORT) {
      run_ablation_report(init, adj_lst, known_edges, nodes, k, sampling_size,
                          splub_cache, ABLATION_REPORT_QUERIES);
      return;
    }

    vector<CurveModel> models;
    vector<function<void()>> cleanup;

    uint64_t t0, t1;
    size_t rss;

    // ---- Persistent landmark cache (saves the selection step across runs) ----
    // Filename key: dataset, seed, n, density, K, sampling_size (per the user's
    // spec). The metadata header additionally pins the objective mode + no-hca
    // flag and no-spt flag, so a file is rejected (and rebuilt) if those differ --
    // they change which landmarks get selected. Query-time-only ablations such
    // as no-rmq do not affect selection, but we still include no-rmq in the KNN
    // bound-cache key to keep runtime comparisons isolated.
    auto dens_str = [](double d) {
      std::ostringstream o; o.setf(std::ios::fixed); o.precision(4); o << d; return o.str();
    };
    std::string cache_key =
        CACHE_DATASET + "_s" + std::to_string(init) + "_n" + std::to_string(nodes) +
        "_d" + dens_str(CACHE_DENSITY) + "_k" + std::to_string(k) +
        "_ss" + std::to_string(sampling_size);
    std::string cache_meta =
        "obj=" + SGEL_OBJECTIVE_MODE + ";nohca=" + (ABLATION_NO_HCA ? "1" : "0") +
        ";nospt=" + (ABLATION_NO_SPT ? "1" : "0") +
        ";ds=" + CACHE_DATASET + ";n=" + std::to_string(nodes) +
        ";k=" + std::to_string(k) + ";ss=" + std::to_string(sampling_size);
    auto cache_tag_sanitize = [](std::string s) {
      for (char &c : s) if (c == ' ' || c == '/' || c == '\\') c = '_';
      return s;
    };
    auto cache_path_for = [&](const std::string &method_tag) {
      return CACHE_DIR + "/" + cache_key + "__" + cache_tag_sanitize(method_tag) + ".lm";
    };
    // Try to load `method_tag`'s landmarks for `elm`; on success rebuild the
    // SP-trees/RMQ and return true (selection skipped). On miss return false.
    auto cache_try_load = [&](ELM_SPTree *elm, const std::string &method_tag) -> bool {
      if (CACHE_DIR.empty()) return false;
      std::string p = cache_path_for(method_tag);
      if (elm->load_landmarks(p, cache_meta)) {
        elm->preprocess_landmark_roots();
        cout << "[cache] HIT  " << method_tag << " <- " << p
             << " (landmark selection skipped)" << endl;
        return true;
      }
      return false;
    };
    auto cache_store = [&](ELM_SPTree *elm, const std::string &method_tag) {
      if (CACHE_DIR.empty()) return;
      try { std::filesystem::create_directories(CACHE_DIR); } catch (...) {}
      std::string p = cache_path_for(method_tag);
      try { std::filesystem::create_directories(std::filesystem::path(p).parent_path()); } catch (...) {}
      if (elm->save_landmarks(p, cache_meta))
        cout << "[cache] SAVE " << method_tag << " -> " << p << endl;
      else
        cout << "[cache] WARN could not write " << p
             << " (cache directory creation/write failed)" << endl;
    };

    cout << "[Building ELM exact]" << endl;
    ELM_SPTree *elm_exact = nullptr;
    bool owned_elm_exact = false;
    if (EXACT) {
      // GEL output is independent of `sampling_size`, so reuse the cached
      // version if the caller built it once across the sweep.
      if (cached_elm_exact) {
        elm_exact = cached_elm_exact;
        owned_elm_exact = false;
        // Report the one-time prebuild cost rather than the per-iter 0.
        t0 = 0;
        t1 = cached_elm_exact_prebuild_us;
        rss = get_peak_rss_kb();
        cout << "[Building ELM exact] reused cached GEL output" << endl;
      } else {
        t0 = cpu_time_us();
        elm_exact = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
        elm_exact->set_ablation_no_hca(ABLATION_NO_HCA);
        elm_exact->set_ablation_no_rmq(ABLATION_NO_RMQ);
        elm_exact->set_ablation_no_spt(ABLATION_NO_SPT);
        if (!cache_try_load(elm_exact, "ELM Exact")) {
          elm_exact->get_exact_landmarks();
          elm_exact->preprocess_landmark_roots();
          cache_store(elm_exact, "ELM Exact");
        }
        t1 = cpu_time_us();
        rss = get_peak_rss_kb();
        owned_elm_exact = true;
      }
      models.push_back({
        "ELM Exact",
        [elm_exact](unsigned int u, unsigned int v) { return elm_exact->lookup(u, v); },
        [elm_exact](unsigned int u, unsigned int v) { return elm_exact->lookup_multiple(u, v); },
        [elm_exact]() { return elm_exact->_sizeof(); },
        t1 - t0,
        rss
      });
      if (owned_elm_exact) {
        cleanup.push_back([elm_exact]() { delete elm_exact; });
      }
    }

    cout << "[Building ELM Exact RelObj]" << endl;
    ELM_SPTree *elm_exact_rel = nullptr;
    bool owned_elm_exact_rel = false;
    if (EXACT) {
      if (cached_elm_exact_rel) {
        elm_exact_rel = cached_elm_exact_rel;
        owned_elm_exact_rel = false;
        t0 = 0;
        t1 = cached_elm_exact_rel_prebuild_us;
        rss = get_peak_rss_kb();
        cout << "[Building ELM Exact RelObj] reused cached exact-relative output" << endl;
      } else {
        t0 = cpu_time_us();
        elm_exact_rel = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
        elm_exact_rel->set_ablation_no_hca(ABLATION_NO_HCA);
        elm_exact_rel->set_ablation_no_rmq(ABLATION_NO_RMQ);
        elm_exact_rel->set_ablation_no_spt(ABLATION_NO_SPT);
        if (!cache_try_load(elm_exact_rel, "ELM Exact RelObj")) {
          elm_exact_rel->get_exact_landmarks_relative_parallel();
          elm_exact_rel->preprocess_landmark_roots();
          cache_store(elm_exact_rel, "ELM Exact RelObj");
        }
        t1 = cpu_time_us();
        rss = get_peak_rss_kb();
        owned_elm_exact_rel = true;
      }
      models.push_back({
        "ELM Exact RelObj",
        [elm_exact_rel](unsigned int u, unsigned int v) { return elm_exact_rel->lookup(u, v); },
        [elm_exact_rel](unsigned int u, unsigned int v) { return elm_exact_rel->lookup_multiple(u, v); },
        [elm_exact_rel]() { return elm_exact_rel->_sizeof(); },
        t1 - t0,
        rss
      });
      if (owned_elm_exact_rel) {
        cleanup.push_back([elm_exact_rel]() { delete elm_exact_rel; });
      }
    }

    // cout << "[Building ELM sampling]" << endl;
    // ELM_SPTree *elm_sampling = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
    // t0 = cpu_time_us();
    // elm_sampling->get_sampling_landmarks();
    // elm_sampling->preprocess_landmark_roots();
    // t1 = cpu_time_us();
    
    auto build_sgel_hca_model = [&](const std::string &model_name,
                                    bool relative_objective,
                                    unsigned int seed_offset) {
      cout << "[Building " << model_name << " - HCA/SGEL objective="
           << (relative_objective ? "relative_ARLB" : "raw_L1") << "]" << endl;
      ELM_SPTree *elm_sampling = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
      // Forward ablation flags (R3.D7).
      elm_sampling->set_ablation_no_hca(ABLATION_NO_HCA);
      elm_sampling->set_ablation_no_rmq(ABLATION_NO_RMQ);
      elm_sampling->set_ablation_no_spt(ABLATION_NO_SPT);
      uint64_t local_t0 = cpu_time_us();

      // Cache: skip the (expensive) HCA selection if a matching landmark file
      // exists; otherwise select, then persist for next time.
      if (!cache_try_load(elm_sampling, model_name)) {
        // Use paper-style SGEL/HCA path with OpenMP parallelization over
        // candidate edges. relative_objective=true aligns selection with the
        // reported ARLB metric by weighting each training query by 1/OPT.
        elm_sampling->get_hca_landmarks_parallel(sampling_size, init + seed_offset, relative_objective);
        elm_sampling->preprocess_landmark_roots();
        cache_store(elm_sampling, model_name);
      }
      uint64_t local_t1 = cpu_time_us();

      size_t local_rss = get_peak_rss_kb();
      models.push_back({
        model_name,
        [elm_sampling](unsigned int u, unsigned int v) { return elm_sampling->lookup(u, v); },
        [elm_sampling](unsigned int u, unsigned int v) { return elm_sampling->lookup_multiple(u, v); },
        [elm_sampling]() { return elm_sampling->_sizeof(); },
        local_t1 - local_t0,
        local_rss
      });
      cleanup.push_back([elm_sampling]() { delete elm_sampling; });
    };

    if (SGEL_OBJECTIVE_MODE == "raw" || SGEL_OBJECTIVE_MODE == "both") {
      build_sgel_hca_model("ELM Sampling", false, 73);
    }
    if (SGEL_OBJECTIVE_MODE == "relative" || SGEL_OBJECTIVE_MODE == "both") {
      // Same query-sampling seed as raw objective so the ablation differs only
      // in the objective, not in the training-query set.
      build_sgel_hca_model("ELM Sampling RelObj", true, 73);
    }

    // cout << "[Building ELM stochastic]" << endl;
    // ELM_SPTree *elm_stochastic = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
    // t0 = cpu_time_us();
    // elm_stochastic->get_stochastic_landmarks(STOCH_EPS, sampling_size, init + 7919);
    // elm_stochastic->preprocess_landmark_roots();
    // t1 = cpu_time_us();
    // rss = get_peak_rss_kb();
    // models.push_back({
    //   "ELM Stochastic",
    //   [elm_stochastic](unsigned int u, unsigned int v) { return elm_stochastic->lookup(u, v); },
    //   [elm_stochastic](unsigned int u, unsigned int v) { return elm_stochastic->lookup_multiple(u, v); },
    //   [elm_stochastic]() { return elm_stochastic->_sizeof(); },
    //   t1 - t0,
    //   rss
    // });
    // cleanup.push_back([elm_stochastic]() { delete elm_stochastic; });
    
    // ELM Stochastic: Mirzasoleiman-style candidate sampling. The original
    // code skipped this for LARGE_GRAPH because the unbounded n/k * log(1/eps)
    // candidate sample triggers an SPT+RMQ build per sampled edge (each ~0.5-2s
    // for n=50000). The `stoch_cand_cap` parameter (default 1000) bounds the
    // per-iteration sample size to keep wall clock predictable. Below this
    // cap, the Mirzasoleiman (1 - 1/e - eps) guarantee is preserved.
    {
      cout << "[Building ELM stochastic] (cand_cap=" << stoch_cand_cap << ")" << endl;
      ELM_SPTree *elm_stochastic = new ELM_SPTree(adj_lst, nodes, k, sampling_size);
      elm_stochastic->set_ablation_no_hca(ABLATION_NO_HCA);
      elm_stochastic->set_ablation_no_rmq(ABLATION_NO_RMQ);
      elm_stochastic->set_ablation_no_spt(ABLATION_NO_SPT);
      t0 = cpu_time_us();
      if (!cache_try_load(elm_stochastic, "ELM Stochastic")) {
        elm_stochastic->get_stochastic_landmarks(STOCH_EPS, sampling_size, init + 7919, stoch_cand_cap);
        elm_stochastic->preprocess_landmark_roots();
        cache_store(elm_stochastic, "ELM Stochastic");
      }
      t1 = cpu_time_us();
      rss = get_peak_rss_kb();
      models.push_back({
        "ELM Stochastic",
        [elm_stochastic](unsigned int u, unsigned int v) { return elm_stochastic->lookup(u, v); },
        [elm_stochastic](unsigned int u, unsigned int v) { return elm_stochastic->lookup_multiple(u, v); },
        [elm_stochastic]() { return elm_stochastic->_sizeof(); },
        t1 - t0,
        rss
      });
      cleanup.push_back([elm_stochastic]() { delete elm_stochastic; });
    }

    // Build (or reuse) SimpleGraph for the SPT-based variants.
    SimpleGraph local_sg(0);
    SimpleGraph *sg_ptr = nullptr;
    if (cached_sg) {
      sg_ptr = cached_sg;
    } else {
      local_sg = build_simple_graph(adj_lst, nodes);
      sg_ptr = &local_sg;
    }
    SimpleGraph &sg = *sg_ptr;

    cout << "[Building SPT exact]" << endl;
    SPT_StochasticGreedy *spt_exact = nullptr;
    bool owned_spt_exact = false;
    if (EXACT) {
      // SPT exact (root selection) is also independent of sampling_size; reuse
      // the cached instance across sweep iterations if available.
      if (cached_spt_exact) {
        spt_exact = cached_spt_exact;
        owned_spt_exact = false;
        t0 = 0;
        t1 = cached_spt_exact_prebuild_us;
        rss = get_peak_rss_kb();
        cout << "[Building SPT exact] reused cached SPT-exact output" << endl;
      } else {
        t0 = cpu_time_us();
        spt_exact = new SPT_StochasticGreedy(sg, k, sampling_size, init + 13);
        spt_exact->get_exact_roots();
        t1 = cpu_time_us();
        rss = get_peak_rss_kb();
        owned_spt_exact = true;
      }
      models.push_back({
        "SPT Exact",
        [spt_exact](unsigned int u, unsigned int v) { return spt_exact->lookup(u, v); },
        [spt_exact](unsigned int u, unsigned int v) { return spt_exact->lookup_multiple(u, v); },
        [spt_exact]() { return spt_exact->_sizeof(); },
        t1 - t0,
        rss
      });
      if (owned_spt_exact) {
        cleanup.push_back([spt_exact]() { delete spt_exact; });
      }
    }

    // cout << "[Building SPT stochastic]" << endl;
    // SPT_StochasticGreedy *spt_stochastic = new SPT_StochasticGreedy(sg, k, sampling_size, init + 29);
    // t0 = cpu_time_us();
    // spt_stochastic->get_stochastic_roots(STOCH_EPS);
    // t1 = cpu_time_us();
    // rss = get_peak_rss_kb();
    // models.push_back({
    //   "SPT Stochastic",
    //   [spt_stochastic](unsigned int u, unsigned int v) { return spt_stochastic->lookup(u, v); },
    //   [spt_stochastic](unsigned int u, unsigned int v) { return spt_stochastic->lookup_multiple(u, v); },
    //   [spt_stochastic]() { return spt_stochastic->_sizeof(); },
    //   t1 - t0,
    //   rss
    // });
    // cleanup.push_back([spt_stochastic]() { delete spt_stochastic; });
    
    // SPT Stochastic: same Mirzasoleiman candidate cap rationale as ELM Stochastic.
    {
      cout << "[Building SPT stochastic] (cand_cap=" << stoch_cand_cap << ")" << endl;
      SPT_StochasticGreedy *spt_stochastic = new SPT_StochasticGreedy(sg, k, sampling_size, init + 29);
      t0 = cpu_time_us();
      spt_stochastic->get_stochastic_roots(STOCH_EPS, stoch_cand_cap);
      t1 = cpu_time_us();
      rss = get_peak_rss_kb();
      models.push_back({
        "SPT Stochastic",
        [spt_stochastic](unsigned int u, unsigned int v) { return spt_stochastic->lookup(u, v); },
        [spt_stochastic](unsigned int u, unsigned int v) { return spt_stochastic->lookup_multiple(u, v); },
        [spt_stochastic]() { return spt_stochastic->_sizeof(); },
        t1 - t0,
        rss
      });
      cleanup.push_back([spt_stochastic]() { delete spt_stochastic; });
    }

    cout << "[Building H1]" << endl;
    EdgeLandMark *elm_heuristic1 = nullptr;
    bool owned_h1 = false;
    if (cached_h1) {
      elm_heuristic1 = cached_h1;
      t0 = 0;
      t1 = cached_h1_prebuild_us;
      rss = get_peak_rss_kb();
      cout << "[Building H1] reused cached H1 (sample-size-invariant)" << endl;
    } else {
      // Lazy SP path: large_edge_heuristic doesn't read sp_vector during
      // selection, so we skip the eager all-pairs Dijkstra (40 GB at n=50000)
      // and populate only the landmark-endpoint SP rows afterwards.
      elm_heuristic1 = new EdgeLandMark(adj_lst, nodes, k, sampling_size,
                                        /*precompute_all_pairs=*/false);
      t0 = cpu_time_us();
      elm_heuristic1->large_edge_heuristic();
      elm_heuristic1->ensure_sp_for_landmarks();
      t1 = cpu_time_us();
      rss = get_peak_rss_kb();
      owned_h1 = true;
    }
    models.push_back({
      "ELM H1",
      [elm_heuristic1](unsigned int u, unsigned int v) { return elm_heuristic1->lookup(u, v); },
      [elm_heuristic1](unsigned int u, unsigned int v) { return elm_heuristic1->lookup_multiple(u, v); },
      [elm_heuristic1]() { return elm_heuristic1->_sizeof(); },
      t1 - t0,
      rss
    });
    if (owned_h1) {
      cleanup.push_back([elm_heuristic1]() { delete elm_heuristic1; });
    }

    cout << "[Building H2]" << endl;
    EdgeLandMark *elm_heuristic2 = nullptr;
    bool owned_h2 = false;
    if (cached_h2) {
      elm_heuristic2 = cached_h2;
      t0 = 0;
      t1 = cached_h2_prebuild_us;
      rss = get_peak_rss_kb();
      cout << "[Building H2] reused cached H2 (sample-size-invariant)" << endl;
    } else {
      // H2 reads sp_vector[last_landmark_*] during selection, but the
      // per-landmark-endpoint lazy path inside far_away_heuristic populates
      // those rows on the fly, so we can also skip the eager all-pairs
      // precompute here. Saves another ~40 GB.
      elm_heuristic2 = new EdgeLandMark(adj_lst, nodes, k, sampling_size,
                                        /*precompute_all_pairs=*/false);
      t0 = cpu_time_us();
      elm_heuristic2->far_away_heuristic();
      elm_heuristic2->ensure_sp_for_landmarks();
      t1 = cpu_time_us();
      rss = get_peak_rss_kb();
      owned_h2 = true;
    }
    models.push_back({
      "ELM H2",
      [elm_heuristic2](unsigned int u, unsigned int v) { return elm_heuristic2->lookup(u, v); },
      [elm_heuristic2](unsigned int u, unsigned int v) { return elm_heuristic2->lookup_multiple(u, v); },
      [elm_heuristic2]() { return elm_heuristic2->_sizeof(); },
      t1 - t0,
      rss
    });
    if (owned_h2) {
      cleanup.push_back([elm_heuristic2]() { delete elm_heuristic2; });
    }

    cout << "\n==================== BUILD SUMMARY ====================\n";
    for (const auto &m : models) {
      cout << m.name << " build_us=" << m.build_time_us
           << " peak_rss_kb=" << m.peak_rss_kb
           << " size_bytes=" << m.size_of() << endl;
    }
    cout << "=======================================================\n";

    if (knn_queries > 0 || knn_full_mode) {
      std::string knn_cache_key =
          cache_key + "_obj" + cache_sanitize_token(SGEL_OBJECTIVE_MODE) +
          "_nohca" + std::string(ABLATION_NO_HCA ? "1" : "0") +
          "_normq" + std::string(ABLATION_NO_RMQ ? "1" : "0") +
          "_nospt" + std::string(ABLATION_NO_SPT ? "1" : "0");
      run_knn_metrics_combined(adj_lst, tri, models, nodes, knn_k,
                               knn_queries,
                               knn_pair_samples, init + 4049, knn_full_mode,
                               knn_candidate_limit, KNN_CERT_SELECT,
                               KNN_BOUND_CACHE_DIR, knn_cache_key);
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
      for (unsigned int a = 0; a < nodes; ++a) {
        for (unsigned int b = a + 1; b < nodes; ++b) {
          if (adj_mat->at(a)->at(b) >= -0.1) {
            graph_weight_orig += lb->at(a)->at(b);
          }
        }
      }
    }

    cout << " =====>>> LB begins:\n";

    double relative = 0.;
    unsigned int relative_count = 0;
    double rmse_tri = 0.0;
    unsigned long long sw_saved = 0, tri_saved = 0;

    vector<vector<double>> rmse(models.size(), vector<double>(k, 0.0));
    vector<vector<unsigned long long>> saved(models.size(), vector<unsigned long long>(k, 0));
    vector<uint64_t> duration_model_rt(models.size(), 0);

    // Sum Relative Gain (SRG), per Section 3.3.1 of the paper.
    //   SRG = sum over edges e with exact_LB(e) > 0 of (approx_LB(e) / exact_LB(e)).
    // Reviewer R3.D8 also asks for the per-edge "Sum Relative Error":
    //   sum_rel_err(e) = 1 - approx_LB(e) / exact_LB(e),
    // whose total over the sample (or full edge set in small-graph mode) is
    // exactly N - SRG, where N is the count of edges with exact_LB > 0.
    // We accumulate SRG and N per (model, landmark-prefix); the derived
    // sum-relative-error and per-edge mean are computed at report time.
    vector<vector<double>> srg(models.size(), vector<double>(k, 0.0));
    vector<vector<unsigned long long>> srg_n(models.size(), vector<unsigned long long>(k, 0));
    double srg_tri = 0.0;
    unsigned long long srg_n_tri = 0;

    // v22 terminology/denominator fix for the Zoom-discussed metric.
    // Q is the set of evaluated unknown-query edges.
    //   small/exact datasets: Q = all unknown pairs in the complete graph.
    //   large datasets:       Q = the RMSE sampled missing-edge attempts.
    // We report the denominator explicitly as quality_query_attempts and
    // exclude only queries whose exact/SPLUB lower bound OPT is zero.
    // The average relative lower bound is:
    //   ARLB = sum_{e in Q, OPT(e)>0} approx_LB(e)/OPT(e)
    //          / |{e in Q : OPT(e)>0}|.
    unsigned long long quality_query_attempts = 0;
    unsigned long long zero_opt_queries = 0;

    // SGEL-vs-SPLUB comparison plumbing (Reviewer R1.W2 / R3.D6 / R3.D8).
    // SPLUB is the exact LB reference; the metrics above measure how close
    // each approximate model gets to SPLUB. We additionally track:
    //   - splub_saved: ERC count if we'd used the exact SPLUB LB. This is
    //     the *maximum* ERC any lower-bound-based method can achieve on the
    //     same sample; SGEL/Tri/etc. are bounded above by this number.
    //   - duration_splub_cached_us: total time spent in cached splub_cached()
    //     calls during the RMSE loop. With v6 pre-warm this is O(m) per call.
    //   - splub_raw_us / splub_raw_calls: time and count of *uncached* raw
    //     splub(adj_lst, a, b) calls in the benchmark loop below. This is
    //     the honest O(m + n log n) per-query cost the paper compares
    //     against (SPLUB has no preprocessing per Augustine et al. [3]).
    unsigned long long splub_saved = 0;
    unsigned long long splub_saved_n = 0;
    uint64_t duration_splub_cached_us = 0;
    uint64_t splub_raw_us = 0;
    unsigned int splub_raw_calls = 0;

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

    if (!LARGE_GRAPH) {
      for (unsigned int i = 0; i < nodes; ++i) {
        for (unsigned int j = i + 1; j < nodes; ++j) {
          if (adj_mat->at(i)->at(j) < -0.1) {
            ++quality_query_attempts;
            // v24: Use SPLUB as the exact/OPT reference for the reviewer-facing
            // ARLB/RMSE/AvgRelErr metrics.  Sasha-Wang is sensitive to the
            // initial unknown-UB value; London contains distances around 83,
            // so an old default UB=1.0 made the matrix OPT invalid.  SPLUB
            // computes the exact lower bound directly from the known graph for
            // this queried missing edge.
            uint64_t sp_q0 = cpu_time_us();
            const double exact_lb_ij = splub_cache
                ? splub_cached(*splub_cache, i, j)
                : splub(adj_lst, i, j);
            uint64_t sp_q1 = cpu_time_us();
            duration_splub_cached_us += (sp_q1 - sp_q0);
            if (exact_lb_ij <= 0.0) {
              ++zero_opt_queries;
              continue;
            }
            for (size_t mi = 0; mi < models.size(); ++mi) {
              vector<double> *curve = models[mi].lookup_multiple(i, j);
              unsigned upto = std::min<unsigned>(k, (unsigned)curve->size());
              for (unsigned int idx = 0; idx < upto; ++idx) {
                const double approx = curve->at(idx);
                rmse[mi][idx] += (approx - exact_lb_ij) * (approx - exact_lb_ij);
                saved[mi][idx] += std::distance(sorted_ub.begin(),
                                  lower_bound(sorted_ub.begin(), sorted_ub.end(), approx));
                // SRG term: approx / exact. N - SRG is reported at the end.
                srg[mi][idx] += approx / exact_lb_ij;
                srg_n[mi][idx] += 1;
              }
              delete curve;
            }

            uint64_t q0 = cpu_time_us();
            double tri_val_ij = tri.lookup(i, j);
            uint64_t q1 = cpu_time_us();
            duration_tri_rt += (q1 - q0);

            for (size_t mi = 0; mi < models.size(); ++mi) {
              q0 = cpu_time_us();
              volatile double dummy = models[mi].lookup(i, j);
              q1 = cpu_time_us();
              duration_model_rt[mi] += (q1 - q0);
              (void)dummy;
            }

            rmse_tri += ((tri_val_ij - exact_lb_ij) * (tri_val_ij - exact_lb_ij));
            tri_saved += std::distance(sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(), tri_val_ij));
            sw_saved += std::distance(sorted_ub.begin(), lower_bound(sorted_ub.begin(), sorted_ub.end(), exact_lb_ij));
            srg_tri += tri_val_ij / exact_lb_ij;
            srg_n_tri += 1;
            if (lb_elm != nullptr) {
              relative += 1 - (lb_elm->at(i)->at(j)) / exact_lb_ij;
              relative_count++;
            }
            total_lb_sw += exact_lb_ij;
          }
        }
      }
    } else {
      // NOTE: in v6/v7/v8 we added a parallel pre-warm of the SPLUB cache
      // (`splub_cache->prewarm_all_parallel(nodes)`) on the theory that
      // building all n Dijkstras up-front in parallel would be cheaper than
      // letting them happen on-demand during the sequential RMSE/ERC loops.
      //
      // On the NJIT cluster at n=50000 this turned out to be wrong: the
      // observed pre-warm took 5193 s = 86 min for 50000 Dijkstras at 16
      // threads. Each Dijkstra had a ~25 MB working set (adj list + heap
      // + 400 KB dist_vec + 400 KB max_edge_vec); 16 simultaneous threads
      // produced ~400 MB of contending working sets, far exceeding the L3
      // cache, so parallel scaling collapsed and wall time was close to
      // the sequential 50000 * 0.1 s = 5000 s.
      //
      // Worse, pre-warming the *entire* n sources is wasteful: at smoke
      // (RMSE 2000 + ERC 2000 = 8000 picks) only ~7400 unique sources will
      // ever be queried. Sequential demand-driven caching during the RMSE
      // and ERC loops touches just those ~7400 sources at ~0.1 s each
      // (no parallel contention from a single thread), or ~12 min total.
      //
      // We therefore leave the SplubCache in place (still shared across
      // sweep iterations, still amortizes Dijkstras across both RMSE and
      // ERC via the v6 ERC cache fix) but do NOT pre-warm. Each first
      // call to splub_cache->ensure(src) pays the Dijkstra; subsequent
      // calls for the same src are O(1).
      //
      // If you ever want pre-warming back (e.g. a hardware with high
      // per-core memory bandwidth, or you actually need all n sources),
      // restore `splub_cache->prewarm_all_parallel(nodes);` here. The
      // method itself is unchanged in GraphUtils.{h,cpp}.

      // Raw SPLUB per-query benchmark (Reviewer R1.W2 / R3.D6 / R3.D8).
      // The paper claims SPLUB has "no preprocessing cost" and runs each
      // query in O(m + n log n). To make a fair per-query comparison
      // against SGEL (which runs in O(k) per query but pays heavy
      // preprocessing up front), we time a fixed batch of *uncached* raw
      // splub() calls -- which run two fresh Dijkstras every call. The
      // result is the honest SPLUB per-query cost.
      //
      // NB: do NOT filter out v <= 0 returns. We want to measure SPLUB's
      // time regardless of whether it returns a tight bound. The previous
      // version's filter (`if (v <= 0) continue;`) rejected every sample
      // pair on the sparse SIFT graph because most random pairs have no
      // tight triangle-inequality bound at this density.
      {
        const unsigned int BENCH_N = 50;
        for (unsigned int s = 0; s < BENCH_N; ++s) {
          auto [a, b] = sample_missing_edge(nodes, known_edges);
          uint64_t bq0 = cpu_time_us();
          double v = splub(adj_lst, a, b);
          uint64_t bq1 = cpu_time_us();
          (void)v;
          splub_raw_us += (bq1 - bq0);
          ++splub_raw_calls;
        }
        cout << "[SPLUB raw benchmark] " << splub_raw_calls
             << " calls, total " << splub_raw_us << " us, avg "
             << (splub_raw_calls ? (double)splub_raw_us / splub_raw_calls : 0.0)
             << " us/call" << endl;
      }

      double tri_lb = 0.0;
      cout << "[RMSE phase start] samples=" << RMSE_SAMPLES << std::endl;
      cout.flush();
      auto rmse_start = std::chrono::steady_clock::now();
      const unsigned int RMSE_PROGRESS_EVERY = std::max(100u, RMSE_SAMPLES / 20);
      unsigned long long rmse_accepted = 0, rmse_rejected = 0;

      for (unsigned int s = 0; s < RMSE_SAMPLES; ++s) {
        if (s > 0 && s % RMSE_PROGRESS_EVERY == 0) {
          auto now = std::chrono::steady_clock::now();
          double elapsed = std::chrono::duration<double>(now - rmse_start).count();
          cout << "[RMSE progress] " << s << "/" << RMSE_SAMPLES
               << "  elapsed=" << elapsed << "s"
               << "  accepted=" << rmse_accepted
               << "  rejected=" << rmse_rejected << std::endl;
          cout.flush();
        }
        auto [a, b] = sample_missing_edge(nodes, known_edges);
        ++quality_query_attempts;
        // Use the shared cache if the caller passed one. Otherwise fall back
        // to the original splub() (which spends Dijkstra time per call).
        uint64_t sp_q0 = cpu_time_us();
        double lb_val = splub_cache ? splub_cached(*splub_cache, a, b)
                                    : splub(adj_lst, a, b);
        uint64_t sp_q1 = cpu_time_us();
        duration_splub_cached_us += (sp_q1 - sp_q0);
        // CRITICAL FIX (v11): the original code did `if (lb_val <= 0) { --s;
        // continue; }` -- a reject-and-retry pattern intended to ensure we
        // accumulate exactly RMSE_SAMPLES samples with non-zero SPLUB LB.
        //
        // On sparse-by-density / small-world graphs (e.g. SIFT-HNSW layer 0
        // where most random pairs have SPLUB_LB = 0 because uniform edge
        // weights + short paths make the triangle inequality rarely bind),
        // the accept rate can be ~1%. Combined with the v6 SplubCache, each
        // REJECTED sample still triggers a Dijkstra for any uncached
        // endpoint. Across ~50000 attempts to gather 500 accepted samples,
        // that's ~1.7 hours of serial Dijkstra wall time.
        //
        // The fix is to count rejections but NOT retry. RMSE_SAMPLES becomes
        // a bound on attempts, not on accepted samples. srg_n[mi][idx] (and
        // rmse_accepted) reports the actual valid sample count so the user
        // can see the accept rate. RMSE / SRG / SumRelErr / AvgRelErr are
        // then computed only over the accepted subset (denominator = N).
        // If the accept rate is too low for meaningful statistics, the user
        // can bump --rmse-samples to get more valid samples without paying
        // an unbounded retry cost.
        if (lb_val <= 0) { ++rmse_rejected; ++zero_opt_queries; continue; }
        ++rmse_accepted;

        for (size_t mi = 0; mi < models.size(); ++mi) {
          vector<double> *curve = models[mi].lookup_multiple(a, b);
          unsigned upto = std::min<unsigned>(k, (unsigned)curve->size());
          for (unsigned int idx = 0; idx < upto; ++idx) {
            const double approx = curve->at(idx);
            rmse[mi][idx] += (approx - lb_val) * (approx - lb_val);
            // SRG term per (model, idx): approx / SPLUB_exact_LB. The reported
            // error N - SRG = sum over sampled edges of (1 - approx/exact),
            // which is the Sum Relative Error in Reviewer R3.D8.
            srg[mi][idx] += approx / lb_val;
            srg_n[mi][idx] += 1;
          }
          delete curve;
        }

        uint64_t q0 = cpu_time_us();
        tri_lb = tri.lookup(a, b);
        uint64_t q1 = cpu_time_us();
        duration_tri_rt += (q1 - q0);
        rmse_tri += (tri_lb - lb_val) * (tri_lb - lb_val);
        srg_tri += tri_lb / lb_val;
        srg_n_tri += 1;

        // In sparse-graph mode, lb_elm is null (no n*n matrix); skip the
        // edge-relative metric.
        if (lb_elm != nullptr) {
          relative += 1.0 - (lb_elm->at(a)->at(b)) / lb_val;
          relative_count++;
        }
        total_lb_sw += lb_val;

        for (size_t mi = 0; mi < models.size(); ++mi) {
          q0 = cpu_time_us();
          volatile double dummy = models[mi].lookup(a, b);
          q1 = cpu_time_us();
          duration_model_rt[mi] += (q1 - q0);
          (void)dummy;
        }
      }
      {
        auto rmse_done = std::chrono::steady_clock::now();
        double rmse_elapsed = std::chrono::duration<double>(rmse_done - rmse_start).count();
        cout << "[RMSE phase done] " << RMSE_SAMPLES << " attempted in "
             << rmse_elapsed << "s"
             << "  accepted=" << rmse_accepted
             << "  rejected=" << rmse_rejected
             << "  accept_rate="
             << (rmse_accepted + rmse_rejected > 0
                  ? (100.0 * rmse_accepted / (rmse_accepted + rmse_rejected))
                  : 0.0)
             << "%" << std::endl;
        cout.flush();
      }

      cout << "[ERC phase start] samples=" << ERC_SAMPLES << std::endl;
      cout.flush();
      auto erc_start = std::chrono::steady_clock::now();
      const unsigned int ERC_PROGRESS_EVERY = std::max(100u, ERC_SAMPLES / 20);

      for (int res_i = 0; res_i < (int)ERC_SAMPLES; ++res_i) {
        if (res_i > 0 && res_i % (int)ERC_PROGRESS_EVERY == 0) {
          auto now = std::chrono::steady_clock::now();
          double elapsed = std::chrono::duration<double>(now - erc_start).count();
          cout << "[ERC progress] " << res_i << "/" << ERC_SAMPLES
               << "  elapsed=" << elapsed << "s" << std::endl;
          cout.flush();
        }
        unsigned int u1 = rand() % nodes;
        unsigned int v1 = rand() % nodes;
        while (u1 == v1) v1 = rand() % nodes;

        unsigned int u2 = rand() % nodes;
        unsigned int v2 = rand() % nodes;
        while (u2 == v2) v2 = rand() % nodes;

        while (u1 == u2 && v1 == v2) {
          u2 = rand() % nodes;
          v2 = rand() % nodes;
          while (u2 == v2) v2 = rand() % nodes;
        }

        auto key_u1v1 = make_pair(min(u1, v1), max(u1, v1));
        auto key_u2v2 = make_pair(min(u2, v2), max(u2, v2));
        u1 = key_u1v1.first;
        v1 = key_u1v1.second;
        u2 = key_u2v2.first;
        v2 = key_u2v2.second;

        vector<vector<double>*> curves1(models.size(), nullptr), curves2(models.size(), nullptr);
        for (size_t mi = 0; mi < models.size(); ++mi) {
          curves1[mi] = models[mi].lookup_multiple(u1, v1);
          curves2[mi] = models[mi].lookup_multiple(u2, v2);
        }

        // Use the shared SplubCache when available so we don't run a fresh
        // Dijkstra on every ERC sample. ERC_SAMPLES is typically thousands;
        // at n=50000 each Dijkstra is ~200ms, so uncached this loop alone
        // costs 2 * ERC_SAMPLES * 0.2s = many minutes -> hours. With the
        // cache, repeat sources are O(1) and unique sources still pay the
        // Dijkstra once across both RMSE and ERC phases.
        double ub_u1v1, ub_u2v2;
        if (splub_cache) {
          ub_u1v1 = splub_cache->ensure(u1)->at(v1);
          ub_u2v2 = splub_cache->ensure(u2)->at(v2);
        } else {
          auto sp1 = DijkstraELM(adj_lst, u1);
          auto sp2 = DijkstraELM(adj_lst, u2);
          ub_u1v1 = sp1.first->at(v1);
          ub_u2v2 = sp2.first->at(v2);
          delete sp1.first; if (sp1.second) delete sp1.second;
          delete sp2.first; if (sp2.second) delete sp2.second;
        }

        for (unsigned int idx = 0; idx < k; ++idx) {
          for (size_t mi = 0; mi < models.size(); ++mi) {
            if (idx < curves1[mi]->size() && idx < curves2[mi]->size()) {
              if (ub_u1v1 < curves2[mi]->at(idx) || ub_u2v2 < curves1[mi]->at(idx)) {
                saved[mi][idx] += 1;
              }
            }
          }
        }

        // SPLUB's ERC count (Reviewer R1.W2 / R3.D6). SPLUB returns the
        // exact/tightest lower bound, so this is the upper envelope for any
        // lower-bound-based resolver on the same sampled pair-comparisons.
        // v34 correctness fix: evaluate SPLUB for EVERY ERC sample, not only
        // when all endpoints are already cache hits. The old cache-hit-only
        // shortcut produced an incomplete/biased SPLUB ERC denominator and
        // printed 0/0 in exact/full-quality mode. Using splub_cached keeps
        // the value exact while amortizing Dijkstra rows across samples.
        double splub_lb_u1v1 = splub_cache ? splub_cached(*splub_cache, u1, v1)
                                          : splub(adj_lst, u1, v1);
        double splub_lb_u2v2 = splub_cache ? splub_cached(*splub_cache, u2, v2)
                                          : splub(adj_lst, u2, v2);
        if (ub_u1v1 < splub_lb_u2v2 || ub_u2v2 < splub_lb_u1v1) splub_saved += 1;
        ++splub_saved_n;

        uint64_t q0 = cpu_time_us();
        double tri_lb_u1v1 = tri.lookup(u1, v1);
        double tri_lb_u2v2 = tri.lookup(u2, v2);
        uint64_t q1 = cpu_time_us();
        duration_tri_rt += (q1 - q0);
        if (ub_u1v1 < tri_lb_u2v2 || ub_u2v2 < tri_lb_u1v1) tri_saved += 1;

        for (size_t mi = 0; mi < models.size(); ++mi) {
          delete curves1[mi];
          delete curves2[mi];
        }
      }
      {
        auto erc_done = std::chrono::steady_clock::now();
        double erc_elapsed = std::chrono::duration<double>(erc_done - erc_start).count();
        cout << "[ERC phase done] " << ERC_SAMPLES << " samples in "
             << erc_elapsed << "s" << std::endl;
        cout.flush();
      }
    }

    // v34: In exact/full-quality mode (London/Porto legacy matrix mode), the
    // large-graph ERC loop above is not executed. Still compute the same
    // sampled SPLUB ERC measure so the printed SPLUB resolver count is real
    // instead of "not evaluated". We use the exact shortest-path UB matrix
    // already computed for small graphs and exact SPLUB lower bounds.
    if (!LARGE_GRAPH && ERC_SAMPLES > 0) {
      for (unsigned int res_i = 0; res_i < ERC_SAMPLES; ++res_i) {
        unsigned int u1 = rand() % nodes;
        unsigned int v1 = rand() % nodes;
        while (u1 == v1) v1 = rand() % nodes;

        unsigned int u2 = rand() % nodes;
        unsigned int v2 = rand() % nodes;
        while (u2 == v2) v2 = rand() % nodes;

        auto key_u1v1 = make_pair(min(u1, v1), max(u1, v1));
        auto key_u2v2 = make_pair(min(u2, v2), max(u2, v2));
        while (key_u1v1 == key_u2v2) {
          u2 = rand() % nodes;
          v2 = rand() % nodes;
          while (u2 == v2) v2 = rand() % nodes;
          key_u2v2 = make_pair(min(u2, v2), max(u2, v2));
        }
        u1 = key_u1v1.first;
        v1 = key_u1v1.second;
        u2 = key_u2v2.first;
        v2 = key_u2v2.second;

        const double ub_u1v1 = ub->at(u1)->at(v1);
        const double ub_u2v2 = ub->at(u2)->at(v2);
        const double splub_lb_u1v1 = splub_cache ? splub_cached(*splub_cache, u1, v1)
                                                 : splub(adj_lst, u1, v1);
        const double splub_lb_u2v2 = splub_cache ? splub_cached(*splub_cache, u2, v2)
                                                 : splub(adj_lst, u2, v2);
        if (ub_u1v1 < splub_lb_u2v2 || ub_u2v2 < splub_lb_u1v1) splub_saved += 1;
        ++splub_saved_n;
      }
    }

    // ---------------------------------------------------------------------
    // Raw (UNCACHED) SPLUB per-query benchmark.
    //
    // SPLUB is used in two places with DIFFERENT cache policies:
    //   (1) Computing the OPT reference for ARLB/RMSE quality. The cached
    //       path (splub_cached) is used there purely for speed; the LB
    //       *value* returned is identical to uncached splub(), so quality
    //       numbers are unaffected by caching.
    //   (2) Timing the SPLUB per-query cost (this block). Here caching MUST
    //       NOT be used. The cache stores one Dijkstra per source node, so
    //       across thousands of sampled queries that reuse only n unique
    //       endpoints the two Dijkstras (SPLUB's dominant O(m + n log n)
    //       term) amortize to ~0, leaving just the O(m) edge scan. That
    //       badly understates the true per-query cost (e.g. on London the
    //       cached number is ~11 us/query vs the honest ~131 us/query).
    //
    // We therefore time a batch of uncached splub() calls -- two fresh
    // Dijkstras + the known-edge scan each -- the honest O(m + n log n)
    // per-query cost the paper compares against SGEL's O(k) query.
    //
    // The LARGE_GRAPH branch above already ran this benchmark inside its
    // sampling path; the EXACT/small-graph branch (London, Porto) did not,
    // so historically those runs reported ONLY the misleading cached time.
    // Gating on splub_raw_calls==0 runs it exactly once: here for small
    // graphs, in the branch above for large ones (never double-counted).
    if (splub_raw_calls == 0 && known_edges != nullptr) {
      const unsigned int BENCH_N = 200;  // small graphs are cheap; take more for a stable mean
      for (unsigned int s = 0; s < BENCH_N; ++s) {
        auto [a, b] = sample_missing_edge(nodes, known_edges);
        uint64_t bq0 = cpu_time_us();
        double v = splub(adj_lst, a, b);  // UNCACHED on purpose
        uint64_t bq1 = cpu_time_us();
        (void)v;
        splub_raw_us += (bq1 - bq0);
        ++splub_raw_calls;
      }
      cout << "[SPLUB raw benchmark - exact mode] " << splub_raw_calls
           << " calls, total " << splub_raw_us << " us, avg "
           << (splub_raw_calls ? (double)splub_raw_us / splub_raw_calls : 0.0)
           << " us/call" << endl;
    }

    cout << "Query time tri : " << duration_tri_rt << endl;
    for (size_t mi = 0; mi < models.size(); ++mi) {
      cout << "Query time " << models[mi].name << " : " << duration_model_rt[mi] << endl;
    }

    cout << "Duration SW LB: "  << duration_lb_sw << endl;
    cout << "Total Original Graph weight : " << graph_weight_orig << endl;
    cout << "Total SW LB : " << total_lb_sw << endl;
    cout << "Deprecated lb_elm AvgRelErr on edge (ignore for paper) : "
         << (relative_count ? relative / relative_count : 0.0) << endl;

    // v22: Per (model, landmark-prefix) reporting for the Zoom-discussed
    // relative lower-bound quality metric.  We keep the legacy variable name
    // srg internally, but report it as SumRelLB / ARLB.
    //
    //   SumRelLB = sum_{e in Q, OPT(e)>0} method_LB(e) / OPT(e)
    //   ARLB     = SumRelLB / valid_queries
    //   SumRelErr = sum_{e in Q, OPT(e)>0} (1 - method_LB(e) / OPT(e))
    //             = valid_queries - SumRelLB
    //   AvgRelErr = 1 - ARLB
    //
    // Q is the evaluated unknown-query set, not the set of known/input edges.
    // zero-OPT queries are excluded from the ratio and counted separately.
    const unsigned long long unknown_edges_total =
        (((unsigned long long)nodes * ((unsigned long long)nodes - 1) / 2)
         - (unsigned long long)known_edges->size());
    const unsigned long long N_total = quality_query_attempts;
    cout << "[QUALITY denominator] unknown_edges_total=" << unknown_edges_total
         << "  quality_query_attempts=" << quality_query_attempts
         << "  zero_opt_queries=" << zero_opt_queries
         << "  valid_quality_queries="
         << ((quality_query_attempts >= zero_opt_queries)
              ? (quality_query_attempts - zero_opt_queries) : 0)
         << endl;

    for (size_t mi = 0; mi < models.size(); ++mi) {
      for (unsigned int idx = 0; idx < k; ++idx) {
        const unsigned long long N_valid = srg_n[mi][idx];
        const double rmse_val = (N_valid > 0) ? sqrt(rmse[mi][idx] / (double)N_valid) : 0.0;
        const double srg_val = srg[mi][idx];
        const double arlb_val = (N_valid > 0) ? srg_val / (double)N_valid : 0.0;
        const double sum_rel_err = (N_valid > 0) ? ((double)N_valid - srg_val) : 0.0;
        const double avg_rel_err = (N_valid > 0) ? (1.0 - arlb_val) : 0.0;
        cout << models[mi].name << " RMSE[" << idx << "] : " << rmse_val << endl;
        cout << models[mi].name << " Saved[" << idx << "] : " << saved[mi][idx] << endl;
        cout << models[mi].name << " SumRelLB[" << idx << "] : " << srg_val
             << "  ARLB[" << idx << "] : " << arlb_val
             << "  SumRelErr[" << idx << "] : " << sum_rel_err
             << "  AvgRelErr[" << idx << "] : " << avg_rel_err
             << "  zero_opt_queries[" << idx << "] : " << zero_opt_queries
             << "  (valid_queries=" << N_valid
             << ", query_attempts=" << N_total
             << ", unknown_edges_total=" << unknown_edges_total << ")" << endl;
      }
    }

    {
      const unsigned long long N_valid = srg_n_tri;
      const double tri_rmse_val = (N_valid > 0) ? sqrt(rmse_tri / (double)N_valid) : 0.0;
      const double tri_arlb = (N_valid > 0) ? srg_tri / (double)N_valid : 0.0;
      const double tri_sum_rel_err = (N_valid > 0) ? ((double)N_valid - srg_tri) : 0.0;
      const double tri_avg_rel_err = (N_valid > 0) ? (1.0 - tri_arlb) : 0.0;
      cout << "TriSearch RMSE : " << tri_rmse_val << endl;
      cout << "Saved TriSearch : " << tri_saved << endl;
      cout << "Saved SW : " << sw_saved << endl;
      cout << "TriSearch SumRelLB : " << srg_tri
           << "  ARLB : " << tri_arlb
           << "  SumRelErr : " << tri_sum_rel_err
           << "  AvgRelErr : " << tri_avg_rel_err
           << "  zero_opt_queries : " << zero_opt_queries
           << "  (valid_queries=" << N_valid
           << ", query_attempts=" << N_total
           << ", unknown_edges_total=" << unknown_edges_total << ")" << endl;
      cout << "Size of TriSearch " << tri._sizeof() << endl;
    }

    // v25: explicit SPLUB quality row. SPLUB is the exact reference used as
    // OPT, so it must appear in result tables with perfect quality and with
    // its measured per-query time. This directly addresses reviewer requests
    // to include SPLUB and report the quality-time tradeoff.
    {
      const unsigned long long N_valid = srg_n_tri;
      const double splub_sum_rel_lb = (double)N_valid;
      const double splub_arlb = (N_valid > 0) ? 1.0 : 0.0;
      const double splub_sum_rel_err = 0.0;
      const double splub_avg_rel_err = 0.0;
      cout << "SPLUB RMSE : 0" << endl;
      if (!LARGE_GRAPH) {
        cout << "Saved SPLUB : " << sw_saved
             << " (full exact-quality resolver count over all evaluated unknown-edge queries)"
             << endl;
      } else {
        cout << "Saved SPLUB : not_computed_in_sampled_large_graph_mode"
             << " (use Saved SPLUB (ERC, exact LB) below for sampled SPLUB resolution)"
             << endl;
      }
      cout << "SPLUB SumRelLB : " << splub_sum_rel_lb
           << "  ARLB : " << splub_arlb
           << "  SumRelErr : " << splub_sum_rel_err
           << "  AvgRelErr : " << splub_avg_rel_err
           << "  zero_opt_queries : " << zero_opt_queries
           << "  (valid_queries=" << N_valid
           << ", query_attempts=" << N_total
           << ", unknown_edges_total=" << unknown_edges_total << ")" << endl;
      // Cached SPLUB timing is accumulated for EVERY query attempt (the
      // duration is added before the zero-OPT `exact_lb_ij <= 0` skip), so the
      // honest divisor and call count is quality_query_attempts, NOT N_valid.
      // For London the gap is negligible (~6 of 12080 are zero-OPT); for SIFT
      // it is large (~99% zero-OPT), where dividing by N_valid would inflate
      // the per-call time by ~N_total/N_valid (often ~100x).
      cout << "SPLUB QUALITY_TIME cached_per_query_us : "
           << ((quality_query_attempts > 0)
                 ? (double)duration_splub_cached_us / (double)quality_query_attempts
                 : 0.0)
           << "  total_us : " << duration_splub_cached_us
           << "  calls : " << quality_query_attempts << endl;
      cout << "SPLUB QUALITY_TIME raw_per_query_us : "
           << ((splub_raw_calls > 0) ? (double)splub_raw_us / (double)splub_raw_calls : 0.0)
           << "  total_us : " << splub_raw_us
           << "  calls : " << splub_raw_calls << endl;
      // Headline: the honest SPLUB per-query cost is the RAW (uncached) one.
      // The cached number is the amortized cost when answering many queries
      // and pre-paying n Dijkstras; it must NOT be quoted as the per-query
      // cost. Emit one explicit line so downstream tooling/readers use raw.
      cout << "SPLUB PER_QUERY_US (headline, uncached) : "
           << ((splub_raw_calls > 0) ? (double)splub_raw_us / (double)splub_raw_calls : 0.0)
           << "  | amortized_cached_us : "
           << ((quality_query_attempts > 0)
                 ? (double)duration_splub_cached_us / (double)quality_query_attempts
                 : 0.0)
           << endl;
    }

    // SGEL-vs-SPLUB comparison block (Reviewer R1.W2 / R3.D6 / R3.D8).
    // SPLUB serves as the exact baseline. For RMSE and SumRelLB/ARLB/
    // AvgRelErr, the per-model lines above measure each method against
    // SPLUB even in small-graph mode (v24). Here we add:
    //   - SPLUB's own ERC count (maximum any LB-based method can achieve)
    //   - SPLUB's per-query time, both raw (one-shot, no preprocessing)
    //     and cached (amortized after pre-warm). The raw number is what
    //     the paper's complexity claim O(m + n log n) refers to.
    {
      cout << "Saved SPLUB (ERC, exact LB) : " << splub_saved
           << " / " << splub_saved_n
           << " (sampled pair-comparisons; exact SPLUB LB evaluated for every sample)"
           << endl;
      cout << "SPLUB ERC% (exact LB) : "
           << ((splub_saved_n > 0) ? 100.0 * (double)splub_saved / (double)splub_saved_n : 0.0)
           << "%  (resolved=" << splub_saved << " / tried=" << splub_saved_n
           << "; exact-LB upper envelope for edge resolution)" << endl;
      cout << "SPLUB cached per-query us : "
           << ((quality_query_attempts > 0)
                 ? (double)duration_splub_cached_us / (double)quality_query_attempts
                 : 0.0)
           << "  (total " << duration_splub_cached_us << " us over "
           << quality_query_attempts << " cached calls)" << endl;
      cout << "SPLUB raw per-query us : "
           << ((splub_raw_calls > 0) ? (double)splub_raw_us / (double)splub_raw_calls : 0.0)
           << "  (total " << splub_raw_us << " us over "
           << splub_raw_calls << " uncached benchmark calls)" << endl;
    }

    // Downstream tasks (Reviewer R3.D2): MST via Prim's + PAM K-Medoid.
    // Runs only when --downstream is set AND we have a true distance matrix
    // (legacy mode, not --sparse-graph). The lookup matrices for each method
    // are materialized inside run_downstream() at O(n^2 * k) per method; at
    // n <= 4900 (London + Porto) this is cheap. For SIFT (n=50000, sparse
    // mode) distance == nullptr and run_downstream() returns immediately
    // with a skip notice -- per the user's scope decision (Q4), downstream
    // is not in scope for the n=50000 high-dim runs.
    if (RUN_DOWNSTREAM) {
      vector<downstream::LookupFn> model_lookups;
      vector<string> model_names;
      model_lookups.reserve(models.size() + 1);
      model_names.reserve(models.size() + 1);

      // v25: SPLUB as an explicit downstream bound-provider baseline.
      // The lookup is exact but expensive; by default it is included only
      // for small matrix-mode datasets (London). Users may raise
      // --downstream-splub-max-n to include it elsewhere.
      if (splub_cache != nullptr && nodes <= DOWNSTREAM_SPLUB_MAX_N) {
        cout << "[Downstream SPLUB] enabled n=" << nodes
             << " <= cap=" << DOWNSTREAM_SPLUB_MAX_N << endl;
        model_lookups.push_back([splub_cache](unsigned u, unsigned v) {
          return splub_cached(*splub_cache, u, v);
        });
        model_names.push_back("SPLUB");
      } else {
        cout << "[Downstream SPLUB] skipped n=" << nodes
             << " cap=" << DOWNSTREAM_SPLUB_MAX_N
             << " (raise --downstream-splub-max-n to include exact SPLUB downstream)"
             << endl;
      }

      for (auto &m : models) {
        model_lookups.push_back(m.lookup);
        model_names.push_back(m.name);
      }
      downstream::run_downstream(distance, adj_mat,
                                 tri, model_lookups, model_names,
                                 nodes, CLUSTER_KS, init,
                                 DOWNSTREAM_LB_AS_DISTANCE);
    }

    for (auto &fn : cleanup) fn();
}

int main(int argc, char **argv) {
  // Make stdout/stderr line-buffered so every "\n" triggers a flush. Without
  // this, a multi-hour build can produce no visible output until the very
  // end (cout is fully-buffered when stdout is a file/pipe like SLURM .out
  // redirection). If the process is killed by walltime before the buffer
  // fills, all progress is lost. v12: force line-buffering so cluster
  // timeouts still leave a useful breadcrumb trail of which phase was
  // running.
  setvbuf(stdout, NULL, _IOLBF, 0);
  setvbuf(stderr, NULL, _IOLBF, 0);

  unsigned int init = time(NULL);
  unsigned int nodes = 1000, k = 10, sampling_size = 2000;
  double prob = 0.05;

  cout << "This program can be run with following options" << endl;
  cout << "./edge_landmarks <random seed> <number of nodes> <edge density> <size of landmarks> <number of samples> <filename> [--knn q pair_samples topk] [--knn-all num_queries topk] [--knn-pool num_queries candidates topk] [--knn-cert-select num_queries topk]" << endl;
  cout << "  --knn q pair_samples topk      : sampled KNN (q query nodes, pair_samples pairs/query, topk=K); add --knn-candidates C to restrict TopK universe" << endl;
  cout << "  --knn-all num_queries topk     : KNN over num_queries query nodes, each vs ALL candidates (no pair sampling), exhaustive TopK=K. Pass num_queries>=N for all nodes." << endl;
  cout << "  --knn-candidates C            : restrict KNN evaluation universe to C seeded candidates/query (works with --knn-all or --knn-cert-select)" << endl;
  cout << "  --knn-pool num_queries C topk : shorthand for --knn-all num_queries topk --knn-candidates C" << endl;
  cout << "  --knn-cert-select num_queries topk : SIFT/high-D diagnostic: score every node by decreasing-LB/increasing-UB separation and evaluate top num_queries selected nodes" << endl;
  cout << "  --cache-dir DIR               : persist selected landmark sets across runs" << endl;
  cout << "  --knn-bound-cache-dir DIR     : persist KNN UB rows and per-method LB rows; later runs reuse existing query rows and append missing rows" << endl;

  if (argc > 1) init = stoi(argv[1]);
  if (argc > 2) nodes = stoi(argv[2]);
  if (argc > 3) prob = stof(argv[3]);
  if (argc > 4) k = stoi(argv[4]);
  if (argc > 5) sampling_size = stoi(argv[5]);

  // Cache-key inputs: dataset basename + density. (seed/n/K/sampling_size are
  // available where the key is assembled inside multiple_sample_exec.)
  CACHE_DENSITY = prob;
  if (argc > 6 && string(argv[6]).rfind("--", 0) != 0) {
    string fp = argv[6];
    size_t slash = fp.find_last_of("/\\");
    CACHE_DATASET = (slash == string::npos) ? fp : fp.substr(slash + 1);
  }

  bool LARGE_GRAPH = (nodes >= 2000);
  EXACT = !LARGE_GRAPH;

  unsigned int knn_queries = 50;
  unsigned int knn_pair_samples = 2000;
  unsigned int knn_k = 50;
  bool knn_full_mode = false;  // SIFT: all nodes as queries, no sampling
  unsigned int knn_candidate_limit = 0; // 0 means all non-src candidates

  for (int ai = 6; ai < argc; ++ai) {
    string arg = string(argv[ai]);
    if (arg == "--knn") {
      if (ai + 2 >= argc) { cerr << "Bad --knn usage\n"; exit(2); }
      knn_queries = (unsigned)stoul(argv[ai + 1]);
      knn_pair_samples = (unsigned)stoul(argv[ai + 2]);
      if (ai + 3 < argc) {
        string maybe = argv[ai + 3];
        if (!maybe.empty() && maybe[0] != '-') {
          knn_k = (unsigned)stoul(maybe);
          ai += 3;
        } else ai += 2;
      } else ai += 2;
    } else if (arg == "--knn-all") {
      // KNN over a USER-PROVIDED number of query nodes, each evaluated against
      // ALL N candidate nodes (no pair sampling). Exhaustive TopK per query.
      //   Usage: --knn-all <num_queries> <topk>
      //   - num_queries query nodes are sampled (seeded, distinct); pass a
      //     value >= N to use every node.
      //   - topk = K.
      // Backward-compatible: if only one numeric arg is given it is treated as
      // <topk> and num_queries keeps its current value (knn_queries).
      knn_full_mode = true;
      bool got_q = false;
      if (ai + 1 < argc) {
        string a1 = argv[ai + 1];
        if (!a1.empty() && a1[0] != '-') {
          // Look ahead for a second numeric -> (num_queries, topk)
          if (ai + 2 < argc) {
            string a2 = argv[ai + 2];
            if (!a2.empty() && a2[0] != '-') {
              knn_queries = (unsigned)stoul(a1);
              knn_k = (unsigned)stoul(a2);
              ai += 2;
              got_q = true;
            }
          }
          if (!got_q) {  // single arg = topk
            knn_k = (unsigned)stoul(a1);
            ai += 1;
          }
        }
      }
    } else if (arg == "--knn-candidates") {
      if (ai + 1 >= argc) { cerr << "Bad --knn-candidates usage\n"; exit(2); }
      knn_candidate_limit = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--knn-pool") {
      // Candidate-pool KNN: num_queries query nodes, C candidates/query, topk=K.
      if (ai + 3 >= argc) { cerr << "Bad --knn-pool usage: --knn-pool num_queries candidates topk\n"; exit(2); }
      knn_full_mode = true;
      knn_queries = (unsigned)stoul(argv[ai + 1]);
      knn_candidate_limit = (unsigned)stoul(argv[ai + 2]);
      knn_k = (unsigned)stoul(argv[ai + 3]);
      ai += 3;
    } else if (arg == "--knn-cert-select") {
      // Certifiable-query KNN: scan all nodes, score each possible query by
      // decreasing-LB / increasing-UB separation, then evaluate the top
      // num_queries selected query nodes with exhaustive TopK certification.
      // Candidate universe is all nodes unless --knn-candidates C is also set.
      if (ai + 2 >= argc) { cerr << "Bad --knn-cert-select usage: --knn-cert-select num_queries topk\n"; exit(2); }
      knn_full_mode = true;
      KNN_CERT_SELECT = true;
      knn_queries = (unsigned)stoul(argv[ai + 1]);
      knn_k = (unsigned)stoul(argv[ai + 2]);
      ai += 2;
    } else if (arg == "--rmse-samples") {
      if (ai + 1 >= argc) { cerr << "Bad --rmse-samples usage\n"; exit(2); }
      RMSE_SAMPLES = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--erc-samples") {
      if (ai + 1 >= argc) { cerr << "Bad --erc-samples usage\n"; exit(2); }
      ERC_SAMPLES = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--stoch-eps") {
      if (ai + 1 >= argc) { cerr << "Bad --stoch-eps usage\n"; exit(2); }
      STOCH_EPS = stod(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--stoch-cand-cap") {
      if (ai + 1 >= argc) { cerr << "Bad --stoch-cand-cap usage\n"; exit(2); }
      STOCH_CAND_CAP = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--ablation-no-hca") {
      ABLATION_NO_HCA = true;
    } else if (arg == "--ablation-report") {
      ABLATION_REPORT = true;
    } else if (arg == "--ablation-report-queries") {
      if (ai + 1 >= argc) { cerr << "Bad --ablation-report-queries usage\n"; exit(2); }
      ABLATION_REPORT_QUERIES = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--ablation-no-rmq") {
      ABLATION_NO_RMQ = true;
    } else if (arg == "--ablation-no-spt" || arg == "--ablation-edge-only") {
      ABLATION_NO_SPT = true;
    } else if (arg == "--ablation-binary-lifting") {
      cerr << "ERROR: --ablation-binary-lifting has been removed. "
           << "Use --ablation-no-spt for the meaningful edge-only/no-SPT quality ablation.\n";
      exit(2);
    } else if (arg == "--sgel-objective") {
      if (ai + 1 >= argc) { cerr << "Bad --sgel-objective usage\n"; exit(2); }
      SGEL_OBJECTIVE_MODE = string(argv[ai + 1]);
      if (SGEL_OBJECTIVE_MODE != "raw" && SGEL_OBJECTIVE_MODE != "relative" &&
          SGEL_OBJECTIVE_MODE != "both") {
        cerr << "--sgel-objective must be one of: raw, relative, both\n";
        exit(2);
      }
      ai += 1;
    } else if (arg == "--sweep") {
      // Comma-separated list of sample sizes for the sweep.
      // Example: --sweep 50,100,200,300
      if (ai + 1 >= argc) { cerr << "Bad --sweep usage\n"; exit(2); }
      string s = argv[ai + 1];
      SWEEP_SAMPLE_SIZES.clear();
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(',', pos);
        string tok = s.substr(pos, comma == string::npos ? string::npos : comma - pos);
        if (!tok.empty()) SWEEP_SAMPLE_SIZES.push_back((unsigned)stoul(tok));
        if (comma == string::npos) break;
        pos = comma + 1;
      }
      if (SWEEP_SAMPLE_SIZES.empty()) {
        cerr << "Empty --sweep list\n"; exit(2);
      }
      ai += 1;
    } else if (arg == "--downstream") {
      // Enable MST + PAM K-Medoid downstream evaluation (Reviewer R3.D2).
      // Requires legacy matrix mode (i.e. NOT --sparse-graph); auto-skips
      // otherwise. Runs once per sweep iteration on the materialized LB
      // matrix for each method (TriSearch + each CurveModel).
      RUN_DOWNSTREAM = true;
    } else if (arg == "--cache-dir") {
      if (ai + 1 >= argc) { cerr << "Bad --cache-dir usage\n"; exit(2); }
      CACHE_DIR = argv[ai + 1];
      ai += 1;
    } else if (arg == "--knn-bound-cache-dir" || arg == "--bound-cache-dir") {
      if (ai + 1 >= argc) { cerr << "Bad --knn-bound-cache-dir usage\n"; exit(2); }
      KNN_BOUND_CACHE_DIR = argv[ai + 1];
      ai += 1;
    } else if (arg == "--downstream-lb-as-distance") {
      // Enable the NON-PAPER quality-degradation metrics (LB used as the
      // distance, yielding a different MST/PAM per method). Off by default.
      DOWNSTREAM_LB_AS_DISTANCE = true;
    } else if (arg == "--downstream-splub-max-n") {
      // v25: include SPLUB as an explicit downstream baseline only when
      // n <= this cap. Default is 500, which includes London but skips
      // Porto by default because all-pairs SPLUB materialization is costly.
      if (ai + 1 >= argc) { cerr << "Bad --downstream-splub-max-n usage\n"; exit(2); }
      DOWNSTREAM_SPLUB_MAX_N = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--knn-splub") {
      // Force the exact-SPLUB baseline row in the KNN / ERC metrics on.
      KNN_SPLUB_ENABLE = true;
    } else if (arg == "--no-knn-splub") {
      // Disable the SPLUB baseline row in the KNN / ERC metrics.
      KNN_SPLUB_ENABLE = false;
    } else if (arg == "--knn-splub-max-n") {
      // Only emit the SPLUB KNN/ERC row when n <= this cap (default 2000).
      // SPLUB needs all-pairs shortest paths, so raising this for Porto/SIFT
      // can be very expensive. Use 0 with --knn-splub to force regardless.
      if (ai + 1 >= argc) { cerr << "Bad --knn-splub-max-n usage\n"; exit(2); }
      KNN_SPLUB_MAX_N = (unsigned)stoul(argv[ai + 1]);
      if (KNN_SPLUB_MAX_N == 0) KNN_SPLUB_MAX_N = std::numeric_limits<unsigned>::max();
      ai += 1;
    } else if (arg == "--cluster-k") {
      // Comma-separated list of PAM cluster counts. Default {3, 5, 10}.
      // Example: --cluster-k 3,5,10
      if (ai + 1 >= argc) { cerr << "Bad --cluster-k usage\n"; exit(2); }
      string s = argv[ai + 1];
      CLUSTER_KS.clear();
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(',', pos);
        string tok = s.substr(pos, comma == string::npos ? string::npos : comma - pos);
        if (!tok.empty()) CLUSTER_KS.push_back((unsigned)stoul(tok));
        if (comma == string::npos) break;
        pos = comma + 1;
      }
      if (CLUSTER_KS.empty()) {
        cerr << "Empty --cluster-k list\n"; exit(2);
      }
      ai += 1;
    } else if (arg == "--subsample-nodes") {
      // v16: node-subsample the loaded sparse graph to this many vertices
      // before any preprocessing. Only meaningful in --sparse-graph mode;
      // ignored (with a warning at startup) otherwise.
      if (ai + 1 >= argc) { cerr << "Bad --subsample-nodes usage\n"; exit(2); }
      SUBSAMPLE_NODES = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--subsample-seed") {
      // v16: RNG seed for the node subsample. Independent of the binary's
      // main `init` seed so the subsample can vary while everything else
      // stays fixed -- this is what drives the multi-seed variance run.
      if (ai + 1 >= argc) { cerr << "Bad --subsample-seed usage\n"; exit(2); }
      SUBSAMPLE_SEED = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--subsample-max-retries") {
      // v16: cap on retry-with-seed+1 attempts when the induced subgraph
      // fails connectivity. Default 50; raise if working with very sparse
      // base graphs. Only meaningful with --subsample-strict-connectivity.
      if (ai + 1 >= argc) { cerr << "Bad --subsample-max-retries usage\n"; exit(2); }
      SUBSAMPLE_MAX_RETRIES = (unsigned)stoul(argv[ai + 1]);
      ai += 1;
    } else if (arg == "--subsample-strict-connectivity") {
      // v16.1: opt back into the v16 retry-then-hard-error policy. With
      // this flag, the binary refuses to accept any induced subgraph that
      // has more than one component, retrying with seed+1 up to
      // --subsample-max-retries times. Without it (default), the LCC
      // fallback kicks in: drop fringe nodes outside the giant component,
      // accept the result. Use this flag only if the experimental setup
      // requires the exact requested n; on the production SIFT n=100k,M=24
      // graph this flag will fail (consistent fringe structure leaves 10-25
      // components in every random 50k sample).
      SUBSAMPLE_STRICT_CONNECTIVITY = true;
    }
  }

  // Sparse-graph mode (for very large n where the full n*n matrix is
  // infeasible): scan separately from the main CLI loop so it can be set
  // independently of positional args.
  bool sparse_graph_mode = false;
  for (int ai = 6; ai < argc; ++ai) {
    if (string(argv[ai]) == "--sparse-graph") { sparse_graph_mode = true; break; }
  }

  cout << "knn_queries = " << knn_queries
       << ", knn_pair_samples = " << knn_pair_samples
       << ", knn_k = " << knn_k
       << ", knn_full_mode = " << (knn_full_mode ? "ON (all nodes, no sampling)" : "off")
       << endl;
  cout << "RMSE_SAMPLES=" << RMSE_SAMPLES
       << " ERC_SAMPLES=" << ERC_SAMPLES
       << " STOCH_EPS=" << STOCH_EPS
       << " STOCH_CAND_CAP=" << STOCH_CAND_CAP << endl;
  cout << "Ablation flags:"
       << " no_hca=" << (ABLATION_NO_HCA ? "ON" : "off")
       << " no_rmq=" << (ABLATION_NO_RMQ ? "ON" : "off")
       << " no_spt=" << (ABLATION_NO_SPT ? "ON (edge-only)" : "off")
       << " report=" << (ABLATION_REPORT ? "ON" : "off")
       << (ABLATION_REPORT ? (" (queries=" + std::to_string(ABLATION_REPORT_QUERIES) + ")") : std::string())
       << endl;
  cout << "SGEL objective mode=" << SGEL_OBJECTIVE_MODE
       << "  (raw=L1, relative=ARLB-aligned, both=ablation)" << endl;
  if (!CACHE_DIR.empty()) cout << "landmark_cache_dir=" << CACHE_DIR << endl;
  cout << "KNN certifiable-query selection=" << (KNN_CERT_SELECT ? "ON" : "off") << endl;
  if (!KNN_BOUND_CACHE_DIR.empty()) cout << "knn_bound_cache_dir=" << KNN_BOUND_CACHE_DIR << endl;
  cout << "DOWNSTREAM_SPLUB_MAX_N=" << DOWNSTREAM_SPLUB_MAX_N
       << "  (SPLUB downstream baseline enabled only when n <= cap)" << endl;
  cout << "Sweep sample sizes:";
  for (auto x : SWEEP_SAMPLE_SIZES) cout << " " << x;
  cout << endl;
  if (SUBSAMPLE_NODES > 0) {
    cout << "Subsample: nodes=" << SUBSAMPLE_NODES
         << " seed=" << SUBSAMPLE_SEED
         << " max_retries=" << SUBSAMPLE_MAX_RETRIES
         << " strict_connectivity="
         << (SUBSAMPLE_STRICT_CONNECTIVITY ? "ON" : "off (LCC fallback)")
         << endl;
    if (!sparse_graph_mode) {
      cerr << "[Subsample] WARNING: --subsample-nodes is only meaningful "
           << "in --sparse-graph mode; ignoring it for this legacy-matrix "
           << "run." << endl;
    }
  }

  srand(init);

  vector<vector<double> *> *distance;
  Graph g;
  vector<list<pair<unsigned int, double>> *> *sparse_adj_lst = nullptr;
  if (sparse_graph_mode) {
    // ---- Sparse-graph mode for very large n ----
    // The file contains ONLY the known edges (no full n*n matrix).
    // Required for n in the tens of thousands where the all-pairs matrix
    // is multi-GB. The Python preprocessing script writes this format
    // when invoked with --sparse-density.
    if (argc <= 6 || string(argv[6]).rfind("--", 0) == 0) {
      cerr << "ERROR: --sparse-graph requires a file path as argv[6]\n";
      exit(2);
    }
    cout << "[Sparse-graph mode] loading adjacency list from "
         << argv[6] << " (no full distance matrix)\n";
    sparse_adj_lst = get_adj_list_file(const_cast<char*>(argv[6]));
    // CRITICAL: in --sparse-graph mode the input file is authoritative for
    // n (its first line). The CLI nodes arg is overridden here so users can
    // pass `0` for nodes and let the file dictate. v12 silently relied on
    // the CLI matching the file; v13 catches the mismatch.
    if (sparse_adj_lst && sparse_adj_lst->size() != nodes) {
      unsigned int file_n = (unsigned int)sparse_adj_lst->size();
      cout << "[Sparse-graph mode] CLI nodes=" << nodes
           << " overridden by file nodes=" << file_n << "\n";
      nodes = file_n;
      // Re-evaluate the LARGE_GRAPH / EXACT switches now that we know the
      // real n. Without this, passing `0` for the CLI nodes leaves
      // LARGE_GRAPH=false and EXACT=true, which routes large inputs through
      // the O(k n^2 m) GEL build that's infeasible above n~2000.
      LARGE_GRAPH = (nodes >= 2000);
      EXACT = !LARGE_GRAPH;
      cout << "[Sparse-graph mode] LARGE_GRAPH=" << (LARGE_GRAPH ? "true" : "false")
           << " EXACT=" << (EXACT ? "true" : "false") << "\n";
    }

    // Node subsampling (v16). When --subsample-nodes N is set and N < nodes,
    // induce a random N-vertex subgraph BEFORE edge-density sub-sampling
    // and before any preprocessing. The induced graph replaces sparse_adj_lst
    // and `nodes` is updated; downstream code sees the smaller graph as if
    // it had been loaded that way from disk. Connectivity is enforced via
    // BFS; failures retry with seed+1 up to SUBSAMPLE_MAX_RETRIES times,
    // then hard-error (no silent fallback to largest component, since that
    // would invalidate cross-seed variance reporting by changing `n`).
    if (SUBSAMPLE_NODES > 0) {
      if (SUBSAMPLE_NODES >= nodes) {
        cout << "[Subsample] WARNING: --subsample-nodes=" << SUBSAMPLE_NODES
             << " >= nodes=" << nodes
             << "; no subsampling performed (using full graph)\n";
      } else {
        auto sres = subsample::induce_node_subgraph(
            sparse_adj_lst, nodes, SUBSAMPLE_NODES,
            SUBSAMPLE_SEED, SUBSAMPLE_MAX_RETRIES,
            SUBSAMPLE_STRICT_CONNECTIVITY);
        if (!sres.success) {
          cerr << "[Subsample] aborting due to persistent connectivity "
               << "failures\n";
          exit(3);
        }
        // Free the original full graph and swap in the induced one.
        clean_up_adj_list(sparse_adj_lst);
        sparse_adj_lst = sres.new_adj_lst;
        nodes = sres.n_new;
        // Re-evaluate mode switches in case the new n crosses the threshold
        // (e.g. subsampling 50k from 100k keeps LARGE_GRAPH=true; subsampling
        // 500 from a 100k graph would flip it).
        LARGE_GRAPH = (nodes >= 2000);
        EXACT = !LARGE_GRAPH;
        cout << "[Subsample] post-induce: nodes=" << nodes
             << " m=" << sres.m_new
             << " seed_used=" << sres.seed_used
             << " retries=" << sres.retries_consumed
             << " LARGE_GRAPH=" << (LARGE_GRAPH ? "true" : "false")
             << " EXACT=" << (EXACT ? "true" : "false") << "\n";
      }
    }

    // Density sub-sampling (v13). When the 3rd CLI arg `prob` is in (0, 1),
    // randomly keep that fraction of the loaded edges. This is how London
    // / Porto density-sweep experiments (d = 0.02, 0.10, 0.15, ...) are
    // realized: the input file holds the FULL graph, and we sub-sample at
    // load time to the requested density. Passing `prob = 0.0` or `>= 1.0`
    // disables sub-sampling and keeps every edge (current SIFT behavior).
    //
    // Strategy: compute a spanning tree (Kruskal/union-find on edges in
    // input order, NOT MST -- weight order is irrelevant for connectivity)
    // and ALWAYS keep its n-1 edges. Then randomly add non-spanning edges
    // until the target count is reached. Guarantees the sub-sampled graph
    // is connected as long as the input graph is connected.
    if (prob > 0.0 && prob < 1.0) {
      size_t total_directed = 0;
      for (unsigned int i = 0; i < nodes; ++i) total_directed += sparse_adj_lst->at(i)->size();
      size_t total_undirected = total_directed / 2;
      size_t target_undirected = (size_t)std::round((double)total_undirected * prob);
      if (target_undirected > total_undirected) target_undirected = total_undirected;

      std::vector<std::tuple<unsigned int, unsigned int, double>> all_edges;
      all_edges.reserve(total_undirected);
      for (unsigned int i = 0; i < nodes; ++i) {
        for (auto &e : *sparse_adj_lst->at(i)) {
          if (i < e.first) all_edges.emplace_back(i, e.first, e.second);
        }
      }

      // Spanning-tree anchor: union-find. Iterate edges in random order so
      // the resulting tree is unbiased w.r.t. edge weight.
      std::mt19937_64 rng((uint64_t)init);
      std::vector<size_t> order(all_edges.size());
      for (size_t i = 0; i < order.size(); ++i) order[i] = i;
      std::shuffle(order.begin(), order.end(), rng);

      std::vector<int> uf_parent(nodes), uf_rank(nodes, 0);
      for (unsigned int i = 0; i < nodes; ++i) uf_parent[i] = (int)i;
      auto uf_find = [&](int x) {
        while (uf_parent[x] != x) { uf_parent[x] = uf_parent[uf_parent[x]]; x = uf_parent[x]; }
        return x;
      };
      auto uf_union = [&](int a, int b) {
        a = uf_find(a); b = uf_find(b);
        if (a == b) return false;
        if (uf_rank[a] < uf_rank[b]) std::swap(a, b);
        uf_parent[b] = a;
        if (uf_rank[a] == uf_rank[b]) uf_rank[a]++;
        return true;
      };

      std::vector<bool> in_kept(all_edges.size(), false);
      std::vector<size_t> spanning_idx;
      std::vector<size_t> non_spanning_idx;
      spanning_idx.reserve(nodes);
      non_spanning_idx.reserve(all_edges.size() - nodes + 1);
      for (size_t oi : order) {
        auto [u, v, w] = all_edges[oi];
        if (uf_union((int)u, (int)v)) {
          spanning_idx.push_back(oi);
          in_kept[oi] = true;
        } else {
          non_spanning_idx.push_back(oi);
        }
      }

      size_t spanning_kept = spanning_idx.size();
      if (target_undirected < spanning_kept) {
        cout << "[Sparse-graph mode] WARNING: target density "
             << prob << " (" << target_undirected << " edges) is below the "
             << "spanning-tree floor (" << spanning_kept << " edges). "
             << "Bumping target to keep the graph connected.\n";
        target_undirected = spanning_kept;
      }
      size_t extras_needed = target_undirected - spanning_kept;
      if (extras_needed > non_spanning_idx.size()) extras_needed = non_spanning_idx.size();
      // Shuffle non-spanning indices to pick a uniform random subset.
      std::shuffle(non_spanning_idx.begin(), non_spanning_idx.end(), rng);
      for (size_t i = 0; i < extras_needed; ++i) {
        in_kept[non_spanning_idx[i]] = true;
      }

      size_t kept_count = 0;
      for (bool b : in_kept) if (b) kept_count++;
      cout << "[Sparse-graph mode] density sub-sampling: " << total_undirected
           << " edges -> " << kept_count << " edges (target prob=" << prob
           << ", spanning=" << spanning_kept << ", extras=" << extras_needed << ")\n";

      // Rebuild adj_lst with only kept edges.
      for (unsigned int i = 0; i < nodes; ++i) sparse_adj_lst->at(i)->clear();
      for (size_t i = 0; i < all_edges.size(); ++i) {
        if (in_kept[i]) {
          auto [u, v, w] = all_edges[i];
          sparse_adj_lst->at(u)->emplace_back(v, w);
          sparse_adj_lst->at(v)->emplace_back(u, w);
        }
      }
    }

    // Build a minimal Boost Graph from a flat (u, v) edge list. Using the
    // adjacency_list constructor that takes an edge-iterator pair avoids
    // having to deal with vertex_descriptor types directly.
    std::vector<std::pair<int,int>> edge_list;
    edge_list.reserve(64);
    for (unsigned int i = 0; i < nodes; ++i) {
      for (auto &e : *sparse_adj_lst->at(i)) {
        if (i < e.first) edge_list.emplace_back((int)i, (int)e.first);
      }
    }
    g = Graph(edge_list.begin(), edge_list.end(), nodes);
    distance = nullptr;  // No n*n matrix in sparse mode.
  } else if (argc > 6 && string(argv[6]).rfind("--", 0) != 0) {
    distance = get_adj_matrix_file(argv[6], 0.);
    // CRITICAL (v13): override CLI nodes with the file's n. The matrix
    // loader reads n from the file's first line and sizes the matrix
    // correctly, but the outer `nodes` variable is unchanged unless we
    // do this. Without the override, downstream code (ER graph
    // generation, all preprocessing) sees the CLI nodes value (often 0
    // in our slurms), producing an empty graph.
    if (distance && distance->size() != nodes) {
      unsigned int file_n = (unsigned int)distance->size();
      cout << "[Legacy matrix mode] CLI nodes=" << nodes
           << " overridden by file nodes=" << file_n << "\n";
      nodes = file_n;
      // Same fix as sparse mode: re-evaluate LARGE_GRAPH / EXACT after the
      // override so n=4900 (Porto) doesn't accidentally trigger the
      // infeasible O(k n^2 m) GEL build path.
      LARGE_GRAPH = (nodes >= 2000);
      EXACT = !LARGE_GRAPH;
      cout << "[Legacy matrix mode] LARGE_GRAPH=" << (LARGE_GRAPH ? "true" : "false")
           << " EXACT=" << (EXACT ? "true" : "false") << "\n";
    }
    if (RENYI_ERDOS) {
      boost::minstd_rand gen(init);
      g = Graph(ERGen(gen, nodes, prob), ERGen(), nodes);
    } else {
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

  boost::property_map<Graph, unsigned int VertexProperties::*>::type id =
      get(&VertexProperties::index, g);
  boost::graph_traits<Graph>::vertex_iterator vi, viend;

  int vert_num = 0;
  for (tie(vi, viend) = vertices(g); vi != viend; ++vi) id[*vi] = vert_num++;

  boost::graph_traits<Graph>::vertex_iterator i, end;
  boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;

  int total = 0;
  for (tie(i, end) = vertices(g); i != end; ++i) {
    int count = 0;
    for (tie(ei, edge_end) = out_edges(*i, g); ei != edge_end; ++ei) ++count;
    total += count;
  }

  cout << "Total nodes : " << nodes << endl;
  cout << "Total edges : " << total / 2. << endl;
  cout << "Components " << check_connected(g, id) << endl;

  vector<list<pair<unsigned int, double>> *> *adj_lst;
  vector<vector<double> *> *lb_elm = nullptr;
  if (sparse_graph_mode) {
    // Sparse mode: reuse the adj_lst loaded directly from the file.
    // No lb_elm matrix (n*n would be infeasible); the sum-relative-error
    // code paths that need it are gated below.
    adj_lst = sparse_adj_lst;
    LARGE_GRAPH = true;  // force the no-Sasha-Wang / sampled-SPLUB path
  } else {
    adj_lst = get_adjacency_list(g, id, distance);
    lb_elm = get_adjacency_matrix(g, id, distance, (double)0.);
  }
  vector<vector<double> *> *adj_mat = nullptr;
  vector<vector<double> *> *lb = nullptr;
  vector<vector<double> *> *ub = nullptr;

  if (!LARGE_GRAPH) {
    adj_mat = get_adjacency_matrix(g, id, distance, (double)-1.);
    lb = get_adjacency_matrix(g, id, distance, (double)0.);
    // v24: Unknown distances in legacy matrix mode are not always in [0,1].
    // London has one outlier node with pairwise distances around 83.  The old
    // default unknown UB=1.0 was invalid and contaminated Sasha-Wang bounds.
    // Use a valid global upper cap from the full metric matrix instead.
    const double default_unknown_ub = matrix_max_offdiag(distance);
    cout << "[matrix-mode] default_unknown_ub=" << default_unknown_ub
         << " (max off-diagonal input distance; replaces old 1.0 cap)" << endl;
    ub = get_adjacency_matrix(g, id, distance, default_unknown_ub);
  }

  map<pair<unsigned int, unsigned int>, double> *known_edges = convert_adjList_to_knownEdges(adj_lst);

  cout << " =====>> Building TriSearch " << endl;
  uint64_t start_tri_time = cpu_time_us();
  TriSearch tri = TriSearch(adj_lst, nodes);
  size_t tri_mem = get_peak_rss_kb();
  uint64_t stop_tri_end = cpu_time_us();
  cout << "TriSearch Time : " << (stop_tri_end - start_tri_time) << endl;
  cout << "PeakRSS(KB) after TriSearch:    " << tri_mem << "\n";

  // Build expensive sample-size-INDEPENDENT artifacts once and reuse them
  // across the sweep below:
  //  - SplubCache amortizes Dijkstra costs across thousands of ground-truth
  //    SPLUB calls.
  //  - GEL (ELM Exact) and SPT exact root selection don't depend on
  //    sampling_size, so rebuilding them per sweep iteration is pure waste.
  //  - H1 (kLE) and H2 (kFE) are pure functions of (graph, k); they don't
  //    depend on sampling_size either. Each rebuild costs O(k) Dijkstras
  //    with the new lazy SP path (~hundreds of MB and minutes instead of
  //    ~40 GB and hours with the old eager all-pairs precompute).
  //  - SimpleGraph is identical across sweep iterations.
  SplubCache shared_splub(adj_lst);
  ELM_SPTree *shared_elm_exact = nullptr;
  ELM_SPTree *shared_elm_exact_rel = nullptr;
  SPT_StochasticGreedy *shared_spt_exact = nullptr;
  SimpleGraph shared_sg = build_simple_graph(adj_lst, nodes);

  // Helpers for persistent landmark cache at the shared-prebuild level.  The
  // cache inside multiple_sample_exec skips SGEL/stochastic selection, while
  // exact GEL/GEL-RelObj are normally prebuilt before the sweep.  These helpers
  // make the first cache level complete for exact GEL as well.
  auto dens_str_main = [](double d) {
    std::ostringstream o; o.setf(std::ios::fixed); o.precision(4); o << d; return o.str();
  };
  auto lm_key_for = [&](unsigned int ss) {
    return CACHE_DATASET + "_s" + std::to_string(init) + "_n" + std::to_string(nodes) +
           "_d" + dens_str_main(CACHE_DENSITY) + "_k" + std::to_string(k) +
           "_ss" + std::to_string(ss);
  };
  auto lm_meta_for = [&](unsigned int ss) {
    return "obj=" + SGEL_OBJECTIVE_MODE + ";nohca=" + (ABLATION_NO_HCA ? "1" : "0") +
           ";nospt=" + (ABLATION_NO_SPT ? "1" : "0") +
           ";ds=" + CACHE_DATASET + ";n=" + std::to_string(nodes) +
           ";k=" + std::to_string(k) + ";ss=" + std::to_string(ss);
  };
  auto lm_path_for = [&](const std::string &method_tag, unsigned int ss) {
    return CACHE_DIR + "/" + lm_key_for(ss) + "__" + cache_sanitize_token(method_tag) + ".lm";
  };
  auto lm_try_load_shared = [&](ELM_SPTree *elm, const std::string &method_tag, unsigned int ss) -> bool {
    if (CACHE_DIR.empty()) return false;
    const std::string p = lm_path_for(method_tag, ss);
    if (elm->load_landmarks(p, lm_meta_for(ss))) {
      elm->preprocess_landmark_roots();
      cout << "[cache] HIT  " << method_tag << " <- " << p
           << " (shared landmark selection skipped)" << endl;
      return true;
    }
    return false;
  };
  auto lm_store_shared = [&](ELM_SPTree *elm, const std::string &method_tag, unsigned int ss) {
    if (CACHE_DIR.empty()) return;
    try { std::filesystem::create_directories(CACHE_DIR); } catch (...) {}
    const std::string p = lm_path_for(method_tag, ss);
    try { std::filesystem::create_directories(std::filesystem::path(p).parent_path()); } catch (...) {}
    if (elm->save_landmarks(p, lm_meta_for(ss)))
      cout << "[cache] SAVE " << method_tag << " -> " << p << endl;
    else
      cout << "[cache] WARN could not write " << p << " (cache directory creation/write failed)" << endl;
  };

  // Capture one-time prebuild costs so the BUILD SUMMARY block emitted
  // inside multiple_sample_exec can report a meaningful build_us for these
  // shared methods (instead of 0, which is technically correct per-iter but
  // misleading when read on its own).
  uint64_t shared_elm_exact_prebuild_us = 0;
  uint64_t shared_elm_exact_rel_prebuild_us = 0;
  uint64_t shared_spt_exact_prebuild_us = 0;
  uint64_t shared_h1_prebuild_us = 0;
  uint64_t shared_h2_prebuild_us = 0;
  if (EXACT) {
    cout << "[Pre-building shared GEL once] " << endl;
    {
      auto _t0 = cpu_time_us();
      shared_elm_exact = new ELM_SPTree(adj_lst, nodes, k, /*sampling_size unused*/ 0);
      shared_elm_exact->set_ablation_no_hca(ABLATION_NO_HCA);
      shared_elm_exact->set_ablation_no_rmq(ABLATION_NO_RMQ);
      shared_elm_exact->set_ablation_no_spt(ABLATION_NO_SPT);
      if (!lm_try_load_shared(shared_elm_exact, "ELM Exact", /*ss*/0)) {
        shared_elm_exact->get_exact_landmarks();
        shared_elm_exact->preprocess_landmark_roots();
        lm_store_shared(shared_elm_exact, "ELM Exact", /*ss*/0);
      }
      shared_elm_exact_prebuild_us = cpu_time_us() - _t0;
      cout << "ELM Exact prebuild_us=" << shared_elm_exact_prebuild_us << endl;
    }
    cout << "[Pre-building shared ELM Exact RelObj once] " << endl;
    {
      auto _t0 = cpu_time_us();
      shared_elm_exact_rel = new ELM_SPTree(adj_lst, nodes, k, /*sampling_size unused*/ 0);
      shared_elm_exact_rel->set_ablation_no_hca(ABLATION_NO_HCA);
      shared_elm_exact_rel->set_ablation_no_rmq(ABLATION_NO_RMQ);
      shared_elm_exact_rel->set_ablation_no_spt(ABLATION_NO_SPT);
      if (!lm_try_load_shared(shared_elm_exact_rel, "ELM Exact RelObj", /*ss*/0)) {
        shared_elm_exact_rel->get_exact_landmarks_relative_parallel();
        shared_elm_exact_rel->preprocess_landmark_roots();
        lm_store_shared(shared_elm_exact_rel, "ELM Exact RelObj", /*ss*/0);
      }
      shared_elm_exact_rel_prebuild_us = cpu_time_us() - _t0;
      cout << "ELM Exact RelObj prebuild_us=" << shared_elm_exact_rel_prebuild_us << endl;
    }
    cout << "[Pre-building shared SPT exact once] " << endl;
    {
      auto _t0 = cpu_time_us();
      shared_spt_exact = new SPT_StochasticGreedy(shared_sg, k, /*samples unused for exact*/ 0, init + 13);
      shared_spt_exact->get_exact_roots();
      shared_spt_exact_prebuild_us = cpu_time_us() - _t0;
      cout << "SPT Exact prebuild_us=" << shared_spt_exact_prebuild_us << endl;
    }
  }
  cout << "[Pre-building shared H1 once (lazy SP path)] " << endl;
  EdgeLandMark *shared_h1 = nullptr;
  {
    auto _t0 = cpu_time_us();
    shared_h1 = new EdgeLandMark(adj_lst, nodes, k, /*samples unused for H1*/ 0,
                                 /*precompute_all_pairs=*/false);
    shared_h1->large_edge_heuristic();
    shared_h1->ensure_sp_for_landmarks();
    shared_h1_prebuild_us = cpu_time_us() - _t0;
    cout << "ELM H1 prebuild_us=" << shared_h1_prebuild_us << endl;
  }
  cout << "[Pre-building shared H2 once (lazy SP path)] " << endl;
  EdgeLandMark *shared_h2 = nullptr;
  {
    auto _t0 = cpu_time_us();
    shared_h2 = new EdgeLandMark(adj_lst, nodes, k, /*samples unused for H2*/ 0,
                                 /*precompute_all_pairs=*/false);
    shared_h2->far_away_heuristic();
    // Defensive: far_away_heuristic() already materializes SP rows for every
    // selected landmark endpoint via its per-iter parallel pair-builder, so
    // this call is a no-op after the v17 fix (returns early on empty
    // to_build). Kept in case a future code path produces an H2 landmark set
    // without filling rows.
    shared_h2->ensure_sp_for_landmarks();
    shared_h2_prebuild_us = cpu_time_us() - _t0;
    cout << "ELM H2 prebuild_us=" << shared_h2_prebuild_us << endl;
  }
  cout << "PeakRSS(KB) after shared builds:    " << get_peak_rss_kb() << "\n";

  for (size_t idx = 0; idx < SWEEP_SAMPLE_SIZES.size(); ++idx) {
    sampling_size = SWEEP_SAMPLE_SIZES[idx];
    srand(init);
    cout << "\n\n================ RUN FOR sampling_size=" << sampling_size << " ================\n";
    multiple_sample_exec(init, LARGE_GRAPH,
      adj_lst, adj_mat, lb, ub, lb_elm, tri, known_edges, distance,
      nodes, k, sampling_size, knn_k, knn_queries, knn_pair_samples, knn_full_mode,
      knn_candidate_limit,
      &shared_splub, shared_elm_exact, shared_elm_exact_rel, shared_spt_exact,
      shared_h1, shared_h2, &shared_sg, STOCH_CAND_CAP,
      shared_elm_exact_prebuild_us, shared_elm_exact_rel_prebuild_us, shared_spt_exact_prebuild_us,
      shared_h1_prebuild_us, shared_h2_prebuild_us);
  }

  delete shared_elm_exact;
  delete shared_elm_exact_rel;
  delete shared_spt_exact;
  delete shared_h1;
  delete shared_h2;

  delete known_edges;
  clean_up_adj_list(adj_lst);
  if (lb_elm != nullptr) clean_up_adj_matrix(lb_elm);

  if (!LARGE_GRAPH) {
    clean_up_adj_matrix(adj_mat);
    clean_up_adj_matrix(lb);
    clean_up_adj_matrix(ub);
  }

  // for (unsigned int ii = 0; ii < nodes; ++ii) delete distance->at(ii);
  // delete distance;
  // return 0;
  
  if (distance != nullptr) {
      for (unsigned int ii = 0; ii < nodes; ++ii) {
        delete distance->at(ii);
      }
      delete distance;
    }
    return 0;
}
