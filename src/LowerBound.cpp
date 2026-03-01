#include "LowerBound.h"

void compute_lb(vector<list<pair<unsigned int, double>> *> *adj_lst,
                vector<vector<double> *> *adj_mat,
                vector<vector<double> *> *lb_mat) {
  pairing_heap<path> H;
  unsigned int nodes = adj_lst->size();
  vector<bool> unknown(nodes, true);
  vector<vector<pairing_heap<path>::handle_type>> handles;
  for (unsigned int i = 0; i < nodes; ++i) {
    vector<pairing_heap<path>::handle_type> v;
    for (unsigned int j = 0; j < nodes; ++j) {
      unknown.at(j) = true;
      v.push_back((pairing_heap<path>::handle_type)NULL);
    }
    unknown.at(i) = false;
    for (list<pair<unsigned int, double>>::iterator it =
             adj_lst->at(i)->begin();
         it != adj_lst->at(i)->end(); ++it) {
      unsigned int end = it->first;
      unknown.at(end) = false;
    }
    for (unsigned int j = 0; j < i; ++j) {
      v.at(j) = handles.at(j).at(i);
    }
    for (unsigned int j = i + 1; j < nodes; ++j) {
      if (unknown.at(j)) {
        v.at(j) = H.push(path(i, j, 0.));
      }
    }
    handles.push_back(v);
  }
  // All paths of length 2
  for (unsigned int i = 0; i < nodes; ++i) {
    auto end = adj_lst->at(i)->end();
    for (auto it = adj_lst->at(i)->begin(); it != end; ++it) {
      for (auto it_fwd = std::next(it); it_fwd != end; ++it_fwd) {
        unsigned int end1 = it->first, end2 = it_fwd->first;
        if (adj_mat->at(end1)->at(end2) > -0.1) {
          continue;
        }
        double val1 = adj_mat->at(i)->at(end1), val2 = adj_mat->at(i)->at(end2);
        double lb_val = abs(val1 - val2);
        if (handles.at(end1).at(end2).node_->value.get_distance() < lb_val) {
          H.increase(handles.at(end1).at(end2), path(end1, end2, lb_val));
        }
      }
    }
  }
  // Empty heap
  while (H.size() > 0) {
    auto edge = H.top().get_edge();
    double dist = max(H.top().get_distance(), 0.);
    H.pop();
    if (lb_mat->at(edge.first)->at(edge.second) < dist) {
      lb_mat->at(edge.first)->at(edge.second) = dist;
      lb_mat->at(edge.second)->at(edge.first) = dist;
    }
    // Add all neighbours
    if (dist <= 0.) {
      continue;
    }
    for (auto it = adj_lst->at(edge.first)->begin();
         it != adj_lst->at(edge.first)->end(); ++it) {
      // Update lengths of edges on heap
      unsigned int end_point = it->first;
      unsigned int min_end = min(end_point, edge.second),
                   max_end = max(end_point, edge.second);
      if (adj_mat->at(min_end)->at(max_end) < 0.) {
        double val = dist - it->second;
        if (val > handles.at(min_end).at(max_end).node_->value.get_distance()) {
          H.increase(handles.at(min_end).at(max_end),
                     path(min_end, max_end, val));
        }
      }
    }
    for (auto it = adj_lst->at(edge.second)->begin();
         it != adj_lst->at(edge.second)->end(); ++it) {
      // Update lengths of edges on heap
      unsigned int end_point = it->first;
      unsigned int min_end = min(end_point, edge.first),
                   max_end = max(end_point, edge.first);
      if (adj_mat->at(min_end)->at(max_end) < 0.) {
        double val = dist - it->second;
        if (val > handles.at(min_end).at(max_end).node_->value.get_distance()) {
          H.increase(handles.at(min_end).at(max_end),
                     path(min_end, max_end, val));
        }
      }
    }
  }
}

void SashaWang(vector<vector<double> *> *lb, vector<vector<double> *> *ub) {
  unsigned int nodes = lb->size();
  for (unsigned int k = 0; k < nodes; ++k) {
    for (unsigned int j = 0; j < nodes; ++j) {
      for (unsigned int i = 0; i < nodes; ++i) {
        double max_lb = max(lb->at(i)->at(k) - ub->at(k)->at(j),
                            lb->at(k)->at(j) - ub->at(i)->at(k));
        lb->at(i)->at(j) = max(lb->at(i)->at(j), max_lb);
        ub->at(i)->at(j) =
            min(ub->at(i)->at(j), ub->at(i)->at(k) + ub->at(k)->at(j));
      }
    }
  }
}
