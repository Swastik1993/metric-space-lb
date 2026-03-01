#ifndef GRAPH_DEFINITIONS_H
#define GRAPH_DEFINITIONS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <utility>


using namespace std;

struct min_heap_edge {
  unsigned int u, v;
  double dist;
  min_heap_edge(unsigned int u, unsigned int v, double dist) {
    this->u = u;
    this->v = v;
    this->dist = dist;
  }

  bool operator<(min_heap_edge const &that) const {
    return this->dist > that.dist;
  }
};

struct VertexProperties {
  unsigned int index;
};

class path {
private:
  pair<unsigned int, unsigned int> edge;
  double distance;

public:
  path(unsigned int end_point1, unsigned int end_point2, double dist) {
    if (end_point1 < end_point2) {
      edge = make_pair(end_point1, end_point2);
    } else {
      edge = make_pair(end_point2, end_point1);
    }
    distance = dist;
    return;
  }

  virtual ~path() {}

  double get_distance() const { return distance; }

  pair<unsigned int, unsigned int> get_edge() const { return edge; }

  void set_distance(double dist) { distance = dist; }

  bool operator<(const path &pathObj) const {
    return distance < pathObj.distance;
  }
};

typedef boost::adjacency_list<boost::vecS, boost::listS, boost::undirectedS,
                              VertexProperties>
    Graph;
typedef boost::erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
typedef boost::adjacency_list_traits<boost::vecS, boost::listS,
                                     boost::undirectedS>
    Traits;

#endif // GRAPH_DEFINITIONS_H
