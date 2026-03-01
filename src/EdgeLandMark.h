#ifndef EDGELANDMARK_H
#define EDGELANDMARK_H

#include "Dijkstra.h"
#include "GraphDefinitions.h"
#include <algorithm>
#include <boost/heap/pairing_heap.hpp>
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <utility>
#include <vector>


using namespace std;
using namespace boost::heap;

class EdgeLandMark {
private:
  vector<list<pair<unsigned int, double>> *> *adj_list;
  unsigned int no_nodes;
  unsigned int no_landmarks;
  unsigned int no_samples;
  vector<vector<double> *> *sp_vector;
  vector<vector<double> *> *le_vector;
  map<pair<unsigned int, unsigned int>, double> *landmarks;
  map<pair<unsigned int, unsigned int>, double> *edges_map;

public:
  EdgeLandMark(vector<list<pair<unsigned int, double>> *> *adj_list,
               unsigned int n_nodes, unsigned int k,
               unsigned int sampling_size);
  ~EdgeLandMark();
  long _sizeof();
  void store_map();
  void large_edge_heuristic();
  void far_away_heuristic();
  void get_exact_landmarks();
  void clean_unwanted_shortest_paths();
  void get_landmarks();
  double lookup(unsigned int u, unsigned int v);
  double lookup_ub(unsigned int u, unsigned int v);
  vector<double> *lookup_multiple(unsigned int u, unsigned int v);
};

#endif // EDGELANDMARK_H
