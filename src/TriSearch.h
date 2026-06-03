#ifndef TRISEARCH_H
#define TRISEARCH_H

#include <algorithm>
#include <cmath>
#include <list>
#include <map>
#include <utility>
#include <vector>


using namespace std;

class TriSearch {
private:
  unsigned int no_nodes;
  vector<map<unsigned int, double> *> *ord_adj_list;
  // Optimization (#13): contiguous sorted-vector mirror of ord_adj_list.
  // Built lazily on first lookup; eliminates the std::map node-chasing
  // overhead in the inner two-pointer walk. Old ord_adj_list is preserved
  // for backward compatibility (size accounting, anyone holding a pointer).
  vector<vector<pair<unsigned int, double>>> sorted_adj_;
  bool sorted_built_ = false;
  void ensure_sorted();

public:
  TriSearch(vector<list<pair<unsigned int, double>> *> *adj_list,
            unsigned int n_nodes);
  ~TriSearch();
  long _sizeof();
  double lookup(unsigned int u, unsigned int v);
  double lookup_ub(unsigned int u, unsigned int v);
  // Faster two-pointer lookup using the contiguous mirror. Same answer as
  // lookup(); call this from hot paths.
  double lookup_fast(unsigned int u, unsigned int v);
};

#endif // TRISEARCH_H
