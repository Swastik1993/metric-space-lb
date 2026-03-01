#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "ShortestPathTree.h"
#include <list>
#include <utility>
#include <vector>


using namespace std;

vector<shortest_path_tree *> *
Dijkstra(vector<list<pair<unsigned int, double>> *> *adj_lst,
         unsigned int source, bool lifting = true);

pair<vector<double> *, vector<double> *>
DijkstraELM(vector<list<pair<unsigned int, double>> *> *adj_lst,
            unsigned int source);

#endif // DIJKSTRA_H
