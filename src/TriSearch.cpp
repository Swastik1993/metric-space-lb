#include "TriSearch.h"
#include <limits>

TriSearch::TriSearch(vector<list<pair<unsigned int, double>> *> *adj_list,
                     unsigned int n_nodes) {
  no_nodes = n_nodes;
  ord_adj_list = new vector<map<unsigned int, double> *>();
  for (unsigned int i = 0; i < no_nodes; ++i) {
    ord_adj_list->push_back(new map<unsigned int, double>());
  }
  for (unsigned int i = 0; i < no_nodes; ++i) {
    for (auto itr = adj_list->at(i)->begin(); itr != adj_list->at(i)->end();
         ++itr) {
      unsigned int j = itr->first;
      double dist = itr->second;
      if (ord_adj_list->at(i)->find(j) == ord_adj_list->at(i)->end()) {
        ord_adj_list->at(i)->insert({j, dist});
        ord_adj_list->at(j)->insert({i, dist});
      }
    }
  }
}

// long TriSearch::_sizeof() {
//   long size_val = sizeof(no_nodes);
//   for (unsigned int i = 0; i < ord_adj_list->size(); ++i) {
//     size_val += sizeof(*ord_adj_list->at(i));
//   }
//   size_val += sizeof(*ord_adj_list);
//   return size_val;
// }

long TriSearch::_sizeof() {
  long s = 0;

  s += sizeof(*this);              // includes pointers + no_nodes
  s += sizeof(*ord_adj_list);      // vector object
  s += ord_adj_list->capacity() * sizeof(map<unsigned int,double>*); // pointer array

  // each map allocation + nodes
  for (unsigned int i = 0; i < ord_adj_list->size(); ++i) {
    auto *mp = ord_adj_list->at(i);
    if (!mp) continue;

    s += sizeof(*mp); // map header

    // estimate per RB-tree node:
    // value = pair<const unsigned int, double>
    // overhead = parent/left/right pointers + color (implementation dependent)
    size_t kv = sizeof(std::pair<const unsigned int, double>);
    size_t node_overhead = 3*sizeof(void*) + sizeof(unsigned int);

    s += (long)(mp->size() * (kv + node_overhead));
  }

  return s;
}

TriSearch::~TriSearch() {
  for (unsigned int i = 0; i < no_nodes; ++i) {
    delete ord_adj_list->at(i);
  }
  delete ord_adj_list;
}

double TriSearch::lookup(unsigned int u, unsigned int v) {
  double lb = 0.0;
  auto itr_u = ord_adj_list->at(u)->begin();
  auto itr_v = ord_adj_list->at(v)->begin();
  if (u == v) {
    return lb;
  }
  while (itr_u != ord_adj_list->at(u)->end() &&
         itr_v != ord_adj_list->at(v)->end()) {
    if (itr_u->first == itr_v->first) {
      lb = max(lb, abs(itr_u->second - itr_v->second));
      ++itr_u;
      ++itr_v;
    } else if (itr_u->first < itr_v->first) {
      ++itr_u;
    } else {
      ++itr_v;
    }
  }
  return lb;
}

double TriSearch::lookup_ub(unsigned int u, unsigned int v) {
  // Trivial UB on a metric is +infinity (or shortest-path if known); a
  // hard-coded 1.0 silently underbounds whenever the graph has any edge with
  // weight > 1. Initialize to +inf so that, if no common neighbor is found,
  // the caller gets a clearly non-informative value rather than a wrong one.
  double ub = std::numeric_limits<double>::infinity();
  if (u == v) {
    return 0.;
  }
  auto itr_u = ord_adj_list->at(u)->begin();
  auto itr_v = ord_adj_list->at(v)->begin();
  while (itr_u != ord_adj_list->at(u)->end() &&
         itr_v != ord_adj_list->at(v)->end()) {
    // BUGFIX: previously compared itr_u->first == itr_v->second (a key vs. a
    // weight), which was always false on real data. The intent is to detect
    // a common neighbor w (== itr_u->first == itr_v->first) and combine the
    // two known distances d(u,w) + d(w,v) as a triangle-inequality UB.
    if (itr_u->first == itr_v->first) {
      ub = std::min(ub, itr_u->second + itr_v->second);
      ++itr_u;
      ++itr_v;
    } else if (itr_u->first < itr_v->first) {
      ++itr_u;
    } else {
      ++itr_v;
    }
  }
  return ub;
}

// ============================================================================
// Optimization (#13): contiguous sorted-vector mirror.
// ============================================================================
void TriSearch::ensure_sorted() {
  if (sorted_built_) return;
  sorted_adj_.assign(no_nodes, {});
  for (unsigned int i = 0; i < no_nodes; ++i) {
    auto *mp = ord_adj_list->at(i);
    if (!mp) continue;
    sorted_adj_[i].reserve(mp->size());
    for (auto &kv : *mp) sorted_adj_[i].push_back({kv.first, kv.second});
    // Already sorted because std::map ordering is by key.
  }
  sorted_built_ = true;
}

double TriSearch::lookup_fast(unsigned int u, unsigned int v) {
  if (u == v) return 0.0;
  if (!sorted_built_) ensure_sorted();
  const auto &au = sorted_adj_[u];
  const auto &av = sorted_adj_[v];
  size_t i = 0, j = 0;
  double lb = 0.0;
  while (i < au.size() && j < av.size()) {
    if (au[i].first == av[j].first) {
      double diff = au[i].second - av[j].second;
      if (diff < 0) diff = -diff;
      if (diff > lb) lb = diff;
      ++i; ++j;
    } else if (au[i].first < av[j].first) {
      ++i;
    } else {
      ++j;
    }
  }
  return lb;
}
