// KnnMetrics.h
#pragma once

#include <list>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

// Project headers you already have
#include "Dijkstra.h"     // DijkstraELM(adj_lst, src) -> pair<vector<double>*, vector<double>*>
#include "ELM_SPTree.h"   // lookup(u,v), lookup_multiple(u,v)
#include "EdgeLandMark.h" // lookup(u,v), lookup_multiple(u,v)
#include "TriSearch.h"    // lookup(u,v)

static inline bool km_isfinite(double x) { return std::isfinite(x); }

struct KmAgg {
  unsigned long long definite_in = 0;
  unsigned long long definite_out = 0;
  unsigned long long unresolved = 0;

  unsigned long long definite_in_total = 0;
  unsigned long long definite_in_correct = 0;

  unsigned long long pair_resolved = 0;
  unsigned long long pair_tried = 0;

  void reset() { *this = KmAgg(); }
};

static inline std::vector<unsigned> km_truth_topk(
    const std::vector<double>& exact,
    unsigned src,
    unsigned K)
{
  const unsigned n = (unsigned)exact.size();
  std::vector<unsigned> idx;
  idx.reserve(n > 0 ? n - 1 : 0);
  for (unsigned v = 0; v < n; ++v) {
    if (v == src) continue;
    idx.push_back(v);
  }

  auto cmp = [&](unsigned a, unsigned b) {
    double da = exact[a], db = exact[b];
    if (!km_isfinite(da)) da = std::numeric_limits<double>::infinity();
    if (!km_isfinite(db)) db = std::numeric_limits<double>::infinity();
    if (da != db) return da < db;
    return a < b; // stable tie-break
  };

  if (K > idx.size()) K = (unsigned)idx.size();
  std::nth_element(idx.begin(), idx.begin() + K, idx.end(), cmp);
  idx.resize(K);
  std::sort(idx.begin(), idx.end(), cmp);
  return idx;
}

// Return {min, median, kth} for a vector over v!=src, with sanitization.
// "kth" here means K-th smallest (1-indexed), i.e., order-stat at K, not (K-1).
static inline std::tuple<double, double, double> km_diag_min_med_kth(
    const std::vector<double>& X,
    unsigned src,
    unsigned K)
{
  std::vector<double> vals;
  vals.reserve(X.size() ? X.size() - 1 : 0);

  for (unsigned v = 0; v < X.size(); ++v) {
    if (v == src) continue;
    double x = X[v];
    if (!km_isfinite(x) || x < 0) x = 0.0;
    vals.push_back(x);
  }
  if (vals.empty()) return {0.0, 0.0, 0.0};

  std::sort(vals.begin(), vals.end());

  double mn = vals.front();
  double md = vals[vals.size() / 2];
  if (K < 1) K = 1;
  if (K > vals.size()) K = (unsigned)vals.size();
  double kth = vals[K - 1];
  return {mn, md, kth};
}

// ===============================
// Top-K certification (final)
// ===============================
//
// OUT certificate for v:
//   at least K nodes are definitely closer than v
//   "definitely closer"  means UB(u) < LB(v)   (strict)
//
// IN certificate for v:
//   at most K-1 nodes can possibly be closer than v
//   "possibly closer" means LB(u) < UB(v)      (strict)
//
// Strict < is important to avoid ties/zeros destroying IN-certification.
//
static inline void km_eval_topk(
    const std::vector<double>& LB,
    const std::vector<double>& UB,
    const std::vector<double>& exact,
    unsigned src,
    unsigned K,
    KmAgg& agg)
{
  if (LB.size() < 2) return;
  const unsigned n = (unsigned)LB.size();
  K = std::min<unsigned>(K, n - 1);

  // Ground truth top-K
  std::vector<unsigned> truth = km_truth_topk(exact, src, K);
  std::unordered_set<unsigned> truth_set(truth.begin(), truth.end());

  for (unsigned v = 0; v < n; ++v) {
    if (v == src) continue;

    double lbv = LB[v], ubv = UB[v];
    if (!km_isfinite(lbv) || lbv < 0) lbv = 0.0;
    if (!km_isfinite(ubv)) ubv = std::numeric_limits<double>::infinity();

    // ---- definite OUT: at least K nodes are definitely closer than v
    // definitely closer means UB(u) < LB(v)
    unsigned definitely_closer = 0;
    for (unsigned u = 0; u < n; ++u) {
      if (u == src || u == v) continue;
      double ubu = UB[u];
      if (!km_isfinite(ubu)) ubu = std::numeric_limits<double>::infinity();
      if (ubu < lbv) {
        definitely_closer++;
        if (definitely_closer >= K) break;
      }
    }
    bool out_cert = (definitely_closer >= K);

    // ---- definite IN: at most K-1 nodes can possibly be closer than v
    // possibly closer means LB(u) < UB(v)    (strict!)
    unsigned possibly_closer = 0;
    for (unsigned u = 0; u < n; ++u) {
      if (u == src || u == v) continue;
      double lbu = LB[u];
      if (!km_isfinite(lbu) || lbu < 0) lbu = 0.0;
      if (lbu < ubv) {
        possibly_closer++;
        if (possibly_closer >= K) break; // we only care if it reaches K
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

// ===============================
// Pair-resolution metric
// ===============================
// A pair (a,b) is "resolved" if you can certify an ordering:
//   UB(a) < LB(b)  OR  UB(b) < LB(a)
static inline bool km_pair_resolved(
    const std::vector<double>& LB,
    const std::vector<double>& UB,
    unsigned a,
    unsigned b)
{
  double lba = LB[a], lbb = LB[b];
  double uba = UB[a], ubb = UB[b];

  if (!km_isfinite(lba) || lba < 0) lba = 0.0;
  if (!km_isfinite(lbb) || lbb < 0) lbb = 0.0;
  if (!km_isfinite(uba)) uba = std::numeric_limits<double>::infinity();
  if (!km_isfinite(ubb)) ubb = std::numeric_limits<double>::infinity();

  return (uba < lbb) || (ubb < lba);
}

// Build LB vector for a query source using a bound structure that supports lookup(src,v).
template <class BoundDS>
static inline void km_fill_lb_from_lookup(
    std::vector<double>& LB,
    BoundDS& ds,
    unsigned src)
{
  const unsigned n = (unsigned)LB.size();
  for (unsigned v = 0; v < n; ++v) {
    if (v == src) { LB[v] = 0.0; continue; }
    double x = ds.lookup(src, v);
    if (!km_isfinite(x) || x < 0) x = 0.0;
    LB[v] = x;
  }
}

// ELM_SPTree has same lookup signature, so above works.
// If you ever want lookup_multiple for speed, you can add a batched version later.

// ===============================
// Main driver
// ===============================
static inline void run_knn_metrics(
    std::vector<std::list<std::pair<unsigned, double>>*>* adj_lst,
    TriSearch& tri,
    ELM_SPTree* elm_exact,            // can be nullptr
    ELM_SPTree* elm_sampling,         // can be nullptr
    EdgeLandMark* elm_h1,             // can be nullptr
    EdgeLandMark* elm_h2,             // can be nullptr
    unsigned nodes,
    unsigned K,
    unsigned num_queries,
    unsigned pair_samples_per_query,
    unsigned seed)
{
  if (num_queries == 0) return;
  if (nodes < 2) return;

  std::mt19937 rng(seed);
  std::uniform_int_distribution<unsigned> unif_node(0, nodes - 1);

  KmAgg triAgg, exactAgg, samplingAgg, h1Agg, h2Agg, combAgg;

  // Buffers
  std::vector<double> exact(nodes, 0.0);
  std::vector<double> UB(nodes, 0.0);

  std::vector<double> LB_tri(nodes, 0.0);
  std::vector<double> LB_exact(nodes, 0.0);
  std::vector<double> LB_sampling(nodes, 0.0);
  std::vector<double> LB_h1(nodes, 0.0);
  std::vector<double> LB_h2(nodes, 0.0);
  std::vector<double> LB_comb(nodes, 0.0);

  for (unsigned qi = 0; qi < num_queries; ++qi) {
    unsigned src = unif_node(rng);

    // Exact shortest paths from src (acts as UB as well, safe & tight)
    auto sp = DijkstraELM(adj_lst, src);
    assert(sp.first && "DijkstraELM returned null dist vector");
    exact.assign(sp.first->begin(), sp.first->end());
    delete sp.first;
    if (sp.second) delete sp.second; // if DijkstraELM returns parents/preds in second

    UB = exact;
    for (unsigned v = 0; v < nodes; ++v) {
      if (!km_isfinite(UB[v])) UB[v] = std::numeric_limits<double>::infinity();
      if (UB[v] < 0) UB[v] = 0.0;
    }

    // Fill lower bounds
    km_fill_lb_from_lookup(LB_tri, tri, src);

    if (elm_exact) km_fill_lb_from_lookup(LB_exact, *elm_exact, src);
    else std::fill(LB_exact.begin(), LB_exact.end(), 0.0);

    if (elm_sampling) km_fill_lb_from_lookup(LB_sampling, *elm_sampling, src);
    else std::fill(LB_sampling.begin(), LB_sampling.end(), 0.0);

    if (elm_h1) km_fill_lb_from_lookup(LB_h1, *elm_h1, src);
    else std::fill(LB_h1.begin(), LB_h1.end(), 0.0);

    if (elm_h2) km_fill_lb_from_lookup(LB_h2, *elm_h2, src);
    else std::fill(LB_h2.begin(), LB_h2.end(), 0.0);

    // Combined LB = max across available LBs
    for (unsigned v = 0; v < nodes; ++v) {
      double x = LB_tri[v];
      if (elm_exact)    x = std::max(x, LB_exact[v]);
      if (elm_sampling) x = std::max(x, LB_sampling[v]);
      if (elm_h1)       x = std::max(x, LB_h1[v]);
      if (elm_h2)       x = std::max(x, LB_h2[v]);
      LB_comb[v] = x;
    }

    // Pair-resolution sampling
    auto sample_pair = [&]() -> std::pair<unsigned, unsigned> {
      unsigned a = unif_node(rng);
      unsigned b = unif_node(rng);
      while (b == a) b = unif_node(rng);
      return {a, b};
    };

    unsigned long long tri_resolved = 0, tri_tried = 0;
    for (unsigned s = 0; s < pair_samples_per_query; ++s) {
      auto [a, b] = sample_pair();
      if (a == src || b == src) { --s; continue; } // keep it comparable across runs
      tri_tried++;
      if (km_pair_resolved(LB_tri, UB, a, b)) tri_resolved++;
    }
    triAgg.pair_tried += tri_tried;
    triAgg.pair_resolved += tri_resolved;

    unsigned long long exact_resolved = 0, exact_tried = 0;
    if (elm_exact) {
      for (unsigned s = 0; s < pair_samples_per_query; ++s) {
        auto [a, b] = sample_pair();
        if (a == src || b == src) { --s; continue; }
        exact_tried++;
        if (km_pair_resolved(LB_exact, UB, a, b)) exact_resolved++;
      }
      exactAgg.pair_tried += exact_tried;
      exactAgg.pair_resolved += exact_resolved;
    }

    unsigned long long sampling_resolved = 0, sampling_tried = 0;
    if (elm_sampling) {
      for (unsigned s = 0; s < pair_samples_per_query; ++s) {
        auto [a, b] = sample_pair();
        if (a == src || b == src) { --s; continue; }
        sampling_tried++;
        if (km_pair_resolved(LB_sampling, UB, a, b)) sampling_resolved++;
      }
      samplingAgg.pair_tried += sampling_tried;
      samplingAgg.pair_resolved += sampling_resolved;
    }

    unsigned long long h1_resolved = 0, h1_tried = 0;
    if (elm_h1) {
      for (unsigned s = 0; s < pair_samples_per_query; ++s) {
        auto [a, b] = sample_pair();
        if (a == src || b == src) { --s; continue; }
        h1_tried++;
        if (km_pair_resolved(LB_h1, UB, a, b)) h1_resolved++;
      }
      h1Agg.pair_tried += h1_tried;
      h1Agg.pair_resolved += h1_resolved;
    }

    unsigned long long h2_resolved = 0, h2_tried = 0;
    if (elm_h2) {
      for (unsigned s = 0; s < pair_samples_per_query; ++s) {
        auto [a, b] = sample_pair();
        if (a == src || b == src) { --s; continue; }
        h2_tried++;
        if (km_pair_resolved(LB_h2, UB, a, b)) h2_resolved++;
      }
      h2Agg.pair_tried += h2_tried;
      h2Agg.pair_resolved += h2_resolved;
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

    // Top-K certification counts
    km_eval_topk(LB_tri,  UB, exact, src, K, triAgg);
    if (elm_exact)     km_eval_topk(LB_exact, UB, exact, src, K, exactAgg);
    if (elm_sampling)  km_eval_topk(LB_sampling, UB, exact, src, K, samplingAgg);
    if (elm_h1)        km_eval_topk(LB_h1,  UB, exact, src, K, h1Agg);
    if (elm_h2)        km_eval_topk(LB_h2,  UB, exact, src, K, h2Agg);
    km_eval_topk(LB_comb, UB, exact, src, K, combAgg);

    // Per-query printing (matches your style)
    auto pct = [](unsigned long long a, unsigned long long b) -> double {
      if (b == 0) return 0.0;
      return 100.0 * (double)a / (double)b;
    };

    std::cout << "\n[KNN-METRICS] Query " << (qi + 1) << "/" << num_queries
              << " src=" << src << "\n";

    std::cout << "  PairResolved (Tri)  = " << pct(tri_resolved, tri_tried) << "%\n";
    std::cout << "  PairResolved (ELM_EXACT)  = " << pct(exact_resolved, exact_tried) << "%\n";
    std::cout << "  PairResolved (ELM_SAMPLING)  = " << pct(sampling_resolved, sampling_tried) << "%\n";
    std::cout << "  PairResolved (H1_ELM)  = " << pct(h1_resolved, h1_tried) << "%\n";
    std::cout << "  PairResolved (H2_ELM)  = " << pct(h2_resolved, h2_tried) << "%\n";
    std::cout << "  TopK(Tri): in=" << (triAgg.definite_in)
              << " out=" << (triAgg.definite_out)
              << " unk=" << (triAgg.unresolved)
              << " in_precision=";
    if (triAgg.definite_in_total == 0) std::cout << 0;
    else std::cout << std::setprecision(6)
                   << (double)triAgg.definite_in_correct / (double)triAgg.definite_in_total;
    std::cout << "\n";

    // Diagnostics (min/med/kth)
    auto print_diag = [&](const char* name, const std::vector<double>& LBv) {
      auto [lmn, lmd, lk] = km_diag_min_med_kth(LBv, src, K);
      auto [umn, umd, uk] = km_diag_min_med_kth(UB,  src, K);
      std::cout << "  DIAG[" << name << "] "
                << "LB(min/med/kth)=" << lmn << "/" << lmd << "/" << lk << "  "
                << "UB(min/med/kth)=" << umn << "/" << umd << "/" << uk << "\n";
    };

    // print_diag("Tri",  LB_tri);
    // if (elm_exact) print_diag("ELM",  LB_exact);
    // if (elm_sampling) print_diag("ELM",  LB_sampling);
    // if (elm_h1)       print_diag("H1",   LB_h1);
    // if (elm_h2)       print_diag("H2",   LB_h2);
    // print_diag("Comb", LB_comb);
  }

  auto pct = [](unsigned long long a, unsigned long long b) -> double {
    if (b == 0) return 0.0;
    return 100.0 * (double)a / (double)b;
  };

  auto print_summary = [&](const char* name, const KmAgg& A) {
    std::cout << "\n===== KNN METRICS [" << name << "] =====\n";
    std::cout << "Queries: " << num_queries << "  K=" << K
              << "  PairSamplesPerQuery=" << pair_samples_per_query << "\n";
    std::cout << "PairResolved%: " << pct(A.pair_resolved, A.pair_tried)
              << " (resolved=" << A.pair_resolved << " tried=" << A.pair_tried << ")\n";
    std::cout << "TopK definite_in=" << A.definite_in
              << " definite_out=" << A.definite_out
              << " unresolved=" << A.unresolved << "\n";
    std::cout << "Definite_in precision=";
    if (A.definite_in_total == 0) std::cout << 0;
    else std::cout << std::setprecision(6)
                   << (double)A.definite_in_correct / (double)A.definite_in_total
                   << " (correct=" << A.definite_in_correct << "/" << A.definite_in_total << ")";
    std::cout << "\n";
  };

  print_summary("Tri", triAgg);
  if (elm_exact)     print_summary("ELM_exact", exactAgg);
  if (elm_sampling)  print_summary("ELM_sampling", samplingAgg);
  if (elm_h1)        print_summary("ELM_h1", h1Agg);
  if (elm_h2)        print_summary("ELM_h2", h2Agg);
  print_summary("Combined(max)", combAgg);
}


