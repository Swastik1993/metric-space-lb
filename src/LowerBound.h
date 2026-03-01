#ifndef LOWERBOUND_H
#define LOWERBOUND_H

#include "GraphDefinitions.h"
#include <algorithm>
#include <boost/heap/pairing_heap.hpp>
#include <cmath>
#include <list>
#include <utility>
#include <vector>


using namespace std;
using namespace boost::heap;

void compute_lb(vector<list<pair<unsigned int, double>> *> *adj_lst,
                vector<vector<double> *> *adj_mat,
                vector<vector<double> *> *lb_mat);

void SashaWang(vector<vector<double> *> *lb, vector<vector<double> *> *ub);

#endif // LOWERBOUND_H
