#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include "GraphDefinitions.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/random/linear_congruential.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


using namespace std;

int check_connected(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id);
unsigned long count_better(vector<double> *lower_bounds,
                           vector<double> *upper_bounds);
vector<list<pair<unsigned int, double>> *> *get_adjacency_list(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id,
    vector<vector<double> *> *dist);
vector<vector<double> *> *get_adjacency_matrix(
    Graph g,
    boost::property_map<Graph, unsigned int VertexProperties::*>::type id,
    vector<vector<double> *> *dist, double default_missing = -1.);
vector<vector<double> *> *distance_matrix(unsigned int nodes, unsigned int dims,
                                          unsigned p = 2);
void clean_up_adj_matrix(vector<vector<double> *> *adj_matrix);
void clean_up_adj_list(vector<list<pair<unsigned int, double>> *> *adj_lst);
map<pair<unsigned int, unsigned int>, double> *convert_adjList_to_knownEdges(
    vector<list<pair<unsigned int, double>> *> *adj_lst);
vector<list<pair<unsigned int, double>> *> *get_adj_list_file(char *filename);
vector<vector<double> *> *get_adj_matrix_file(char *filename,
                                              double default_val = -1.);
Graph waxman_graph_distances(unsigned int n, vector<vector<double> *> *distance,
                             double beta = 0.4, double alpha = 0.1);
pair<Graph, vector<vector<double> *> *>
waxman_graph(unsigned int n, double beta = 0.4, double alpha = 0.1,
             unsigned int p = 1, unsigned int dims = 3);

#endif // GRAPHUTILS_H
