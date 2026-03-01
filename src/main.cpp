#include "Dijkstra.h"
#include "ELM_SPTree.h"
#include "EdgeLandMark.h"
#include "GraphDefinitions.h"
#include "GraphUtils.h"
#include "LowerBound.h"
#include "ShortestPathTree.h"
#include "TriSearch.h"
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>


using namespace std;

const int DEBUG = 0;
const bool EXACT = false;
const bool RENYI_ERDOS = false;
const bool WAXMAN = !RENYI_ERDOS;

int main(int argc, char **argv) {
vector<list<pair<unsigned int, double>> *> *adj_list = new vector<list<pair<unsigned int, double>> *>();

// o1 - 0
list<pair<unsigned int, double>> *list1 = new list<pair<unsigned int, double>>();
list1->push_back(make_pair(1, 0.2)); // inserts 99 before second element
adj_list->push_back(list1);

// o2 - 1
list1 = new list<pair<unsigned int, double>>();
list1->push_back(make_pair(0, 0.2));
list1->push_back(make_pair(3, 0.8));
list1->push_back(make_pair(5, 0.1));
adj_list->push_back(list1);

// o3 - 2
list1 = new list<pair<unsigned int, double>>();
list1->push_back(make_pair(3, 0.3));
list1->push_back(make_pair(4, 0.3));
list1->push_back(make_pair(5, 0.5));
adj_list->push_back(list1);

// o4 - 3
list1 = new list<pair<unsigned int, double>>();
list1->push_back(make_pair(1, 0.8));
list1->push_back(make_pair(2, 0.3));
list1->push_back(make_pair(4, 0.2));
adj_list->push_back(list1);

// o5 - 4
list1 = new list<pair<unsigned int, double>>();
list1->push_back(make_pair(2, 0.3));
list1->push_back(make_pair(3, 0.2));
adj_list->push_back(list1);

// o6 - 5
list1 = new list<pair<unsigned int, double>>();
list1->push_back(make_pair(1, 0.1));
list1->push_back(make_pair(2, 0.5));
list1->push_back(make_pair(6, 0.2));
adj_list->push_back(list1);

// o7 - 6
list1 = new list<pair<unsigned int, double>>();
list1->push_back(make_pair(5, 0.2));
adj_list->push_back(list1);

// ELM_SPTree(vector<list<pair<unsigned int, double>> *> *adj_list, unsigned int n_nodes, unsigned int k, unsigned int sampling_size);
ELM_SPTree* spt = new ELM_SPTree(adj_list, 7, 1, 5);
//spt->get_sampling_landmarks();
spt->get_exact_landmarks();

delete adj_list;
}

