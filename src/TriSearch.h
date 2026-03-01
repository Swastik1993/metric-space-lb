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

public:
  TriSearch(vector<list<pair<unsigned int, double>> *> *adj_list,
            unsigned int n_nodes);
  ~TriSearch();
  long _sizeof();
  double lookup(unsigned int u, unsigned int v);
  double lookup_ub(unsigned int u, unsigned int v);
};

#endif // TRISEARCH_H
