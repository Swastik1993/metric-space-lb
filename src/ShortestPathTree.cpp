#include "ShortestPathTree.h"

// Reconstruct the tree-edge weight to a node's parent.
// In an SPT produced by Dijkstra: dist[node] = dist[parent] + w(parent,node)
static inline double parent_edge_weight(const shortest_path_tree *node) {
  if (!node || !node->parent) return 0.0;
  double w = node->path_length - node->parent->path_length;
  // numerical safety (should never be negative in exact arithmetic)
  return (w < 0.0 ? 0.0 : w);
}

// ==============================
// Utilities: RMQ Sparse Table
// ==============================
static void build_rmq(RMQLCA &rmq) {
  int m = (int)rmq.euler.size();
  rmq.lg.assign(m + 1, 0);
  for (int i = 2; i <= m; ++i) {
    rmq.lg[i] = rmq.lg[i / 2] + 1;
  }
  int K = rmq.lg[m] + 1;
  rmq.st.assign(K, vector<int>(m));
  for (int i = 0; i < m; ++i) {
    rmq.st[0][i] = i; // store position in euler[]
  }
  for (int k = 1; k < K; ++k) {
    int len = 1 << k;
    int half = len >> 1;
    for (int i = 0; i + len <= m; ++i) {
      int p1 = rmq.st[k - 1][i];
      int p2 = rmq.st[k - 1][i + half];
      rmq.st[k][i] = (rmq.depth[p1] <= rmq.depth[p2]) ? p1 : p2;
    }
  }
}

static int rmq_query_pos(const RMQLCA &rmq, int l, int r) {
  if (l > r) std::swap(l, r);
  int len = r - l + 1;
  int k = rmq.lg[len];
  int p1 = rmq.st[k][l];
  int p2 = rmq.st[k][r - (1 << k) + 1];
  return (rmq.depth[p1] <= rmq.depth[p2]) ? p1 : p2;
}

// ==============================
// Cartesian tree DSU utilities
// ==============================
struct DSU {
  vector<int> p;
  vector<int> r;
  vector<CTNode *> comp_root; // root CTNode for each component

  explicit DSU(int n) : p(n), r(n, 0), comp_root(n, nullptr) {
    for (int i = 0; i < n; ++i) p[i] = i;
  }

  int find(int x) {
    if (p[x] == x) return x;
    p[x] = find(p[x]);
    return p[x];
  }

  int unite(int a, int b, CTNode *new_root) {
    a = find(a);
    b = find(b);
    if (a == b) return a;
    if (r[a] < r[b]) std::swap(a, b);
    p[b] = a;
    if (r[a] == r[b]) r[a]++;
    comp_root[a] = new_root;
    return a;
  }
};

CartesianPreproc::~CartesianPreproc() {
  for (CTNode *n : nodes) {
    delete n;
  }
  nodes.clear();
}

shortest_path_tree::~shortest_path_tree() {
  if (jump_pointers) {
    delete jump_pointers;
    jump_pointers = nullptr;
  }
  if (children) {
    delete children;
    children = nullptr;
  }
  // rmq_vertex/cartesian are owned by the root node only.
  if (rmq_vertex) {
    delete rmq_vertex;
    rmq_vertex = nullptr;
  }
  if (cartesian) {
    delete cartesian;
    cartesian = nullptr;
  }
}

// Build Euler tour RMQ for vertex-LCA on the SPT
static RMQLCA *build_vertex_rmq(shortest_path_tree *root,
                               vector<shortest_path_tree *> *sp_tree) {
  const int n = (int)sp_tree->size();
  RMQLCA *rmq = new RMQLCA();
  rmq->first.assign(n, -1);

  // Iterative DFS to avoid recursion depth issues
  struct Frame {
    shortest_path_tree *node;
    int child_idx;
  };

  vector<int> depth_stack;
  std::stack<Frame> st;
  st.push({root, -1});

  while (!st.empty()) {
    Frame &f = st.top();
    shortest_path_tree *u = f.node;

    if (f.child_idx == -1) {
      // entering node
      int u_id = (int)u->id;
      if (rmq->first[u_id] == -1) {
        rmq->first[u_id] = (int)rmq->euler.size();
      }
      rmq->euler.push_back(u_id);
      rmq->depth.push_back((int)u->depth);
      f.child_idx = 0;
    }

    if (u->children == nullptr || f.child_idx >= (int)u->children->size()) {
      st.pop();
      if (!st.empty()) {
        // on return to parent, add parent again
        shortest_path_tree *p = st.top().node;
        rmq->euler.push_back((int)p->id);
        rmq->depth.push_back((int)p->depth);
      }
      continue;
    }

    shortest_path_tree *child = u->children->at(f.child_idx);
    f.child_idx++;
    st.push({child, -1});
  }

  build_rmq(*rmq);
  return rmq;
}

// Build DSU-based Cartesian tree over SPT edges; leaves are vertices.
static CartesianPreproc *build_cartesian_dsu(shortest_path_tree *root,
                                            vector<shortest_path_tree *> *sp_tree) {
  const int n = (int)sp_tree->size();
  auto *cp = new CartesianPreproc();
  cp->leaf.assign(n, -1);

  // Create leaves
  cp->nodes.reserve(2 * n);
  for (int v = 0; v < n; ++v) {
    CTNode *leaf = new CTNode();
    leaf->idx = (int)cp->nodes.size();
    leaf->weight = 0.0;
    cp->nodes.push_back(leaf);
    cp->leaf[v] = leaf->idx;
  }

  // Collect SPT edges (parent-child)
  struct Edge {
    int u, v;
    double w;
  };
  vector<Edge> edges;
  edges.reserve(n - 1);
  for (int v = 0; v < n; ++v) {
    shortest_path_tree *node = sp_tree->at(v);
    if (!node || node->parent == nullptr) continue;
    edges.push_back({(int)node->id, (int)node->parent->id, parent_edge_weight(node)});
  }
  std::sort(edges.begin(), edges.end(), [](const Edge &a, const Edge &b) {
    if (a.w != b.w) return a.w < b.w;
    if (a.u != b.u) return a.u < b.u;
    return a.v < b.v;
  });

  DSU dsu(n);
  // initially, each component's Cartesian root is its leaf
  for (int i = 0; i < n; ++i) {
    dsu.comp_root[i] = cp->nodes[cp->leaf[i]];
  }

  for (const auto &e : edges) {
    int cu = dsu.find(e.u);
    int cv = dsu.find(e.v);
    if (cu == cv) continue;

    CTNode *R1 = dsu.comp_root[cu];
    CTNode *R2 = dsu.comp_root[cv];

    CTNode *edge_node = new CTNode();
    edge_node->idx = (int)cp->nodes.size();
    edge_node->weight = e.w;
    edge_node->left = R1;
    edge_node->right = R2;
    cp->nodes.push_back(edge_node);

    int new_rep = dsu.unite(cu, cv, edge_node);
    // If union by rank swapped reps, ensure comp_root points correctly
    dsu.comp_root[new_rep] = edge_node;
  }

  // Determine root of final component
  int rep0 = dsu.find((int)root->id);
  cp->root = dsu.comp_root[rep0]->idx;

  // Build RMQ-LCA on Cartesian tree nodes
  // Euler tour over CTNode indices
  cp->rmq.first.assign((int)cp->nodes.size(), -1);

  struct CFrame {
    CTNode *node;
    int state; // 0 enter, 1 after left, 2 after right
  };
  std::stack<CFrame> st;
  st.push({cp->nodes[cp->root], 0});

  // depth in CT is implicit via traversal stack size
  vector<int> depth_stack;

  while (!st.empty()) {
    CFrame &f = st.top();
    CTNode *u = f.node;

    if (f.state == 0) {
      if (cp->rmq.first[u->idx] == -1) {
        cp->rmq.first[u->idx] = (int)cp->rmq.euler.size();
      }
      cp->rmq.euler.push_back(u->idx);
      // Use current stack depth as depth
      cp->rmq.depth.push_back((int)st.size() - 1);
      f.state = 1;
      if (u->left) {
        st.push({u->left, 0});
      }
      continue;
    }
    if (f.state == 1) {
      f.state = 2;
      if (u->right) {
        st.push({u->right, 0});
      }
      continue;
    }
    // leaving
    st.pop();
    if (!st.empty()) {
      CTNode *p = st.top().node;
      cp->rmq.euler.push_back(p->idx);
      cp->rmq.depth.push_back((int)st.size() - 1);
    }
  }

  build_rmq(cp->rmq);
  return cp;
}

void preprocess_spt_queries(vector<shortest_path_tree *> *sp_tree) {
  if (!sp_tree || sp_tree->empty()) return;

  // root is the only node with depth==0 for this tree
  shortest_path_tree *root = nullptr;
  for (auto *n : *sp_tree) {
    if (n && n->depth == 0) {
      root = n;
      break;
    }
  }
  if (!root) return;

  // avoid double-preprocessing
  if (root->rmq_vertex == nullptr) {
    root->rmq_vertex = build_vertex_rmq(root, sp_tree);
  }
  if (root->cartesian == nullptr) {
    root->cartesian = build_cartesian_dsu(root, sp_tree);
  }
}

// ==============================
// Existing binary-lifting (unchanged)
// ==============================
void _update_pointers(shortest_path_tree *node) {
  unsigned int level = 0;
  pair<shortest_path_tree *, double> fwd_node = node->jump_pointers->at(0);
  double max_edge = fwd_node.second;
  while (fwd_node.first->jump_pointers != nullptr &&
         fwd_node.first->jump_pointers->size() > level) {
    pair<shortest_path_tree *, double> tmp_node =
        fwd_node.first->jump_pointers->at(level);
    max_edge = max(max_edge, tmp_node.second);
    node->jump_pointers->push_back(make_pair(tmp_node.first, max_edge));
    fwd_node = tmp_node;
    ++level;
  }
}

void binary_lifting(shortest_path_tree *root) {
  vector<shortest_path_tree *> queue_sptree;
  if (root->children == nullptr)
    return;
  for (vector<shortest_path_tree *>::iterator root_it = root->children->begin();
       root_it != root->children->end(); ++root_it) {
    shortest_path_tree *node = *root_it;
    node->jump_pointers = new vector<pair<shortest_path_tree *, double>>();
    node->jump_pointers->push_back(make_pair(root, parent_edge_weight(node)));
    queue_sptree.push_back(node);
  }

  while (queue_sptree.size() > 0) {
    shortest_path_tree *node = queue_sptree[0];
    queue_sptree.erase(queue_sptree.begin());
    if (node->children == nullptr) {
      continue;
    }
    for (vector<shortest_path_tree *>::iterator node_it =
             node->children->begin();
         node_it != node->children->end(); ++node_it) {
      shortest_path_tree *child_node = *node_it;
      child_node->jump_pointers =
          new vector<pair<shortest_path_tree *, double>>();
      child_node->jump_pointers->push_back(
          make_pair(node, parent_edge_weight(child_node)));
      _update_pointers(child_node);
      queue_sptree.push_back(child_node);
    }
  }
}

static unsigned int get_exponent(unsigned int depth) {
  unsigned int exponent = static_cast<unsigned int>(floor(log2(depth)));
  if (pow(2, exponent + 1) <= depth || pow(2, exponent) > depth) {
    cout << "Some error occured!" << endl
         << "Check _adjust_height function" << endl;
  }
  return exponent;
}

static pair<shortest_path_tree *, double> _adjust_height(shortest_path_tree *node,
                                                        unsigned int up_depth) {
  if (up_depth == 1) {
    return make_pair(node->jump_pointers->at(0).first, parent_edge_weight(node));
  }

  double max_edge = 0.;
  shortest_path_tree *up_node = node;

  while (up_depth > 0) {
    unsigned int exponent = get_exponent(up_depth);

    pair<shortest_path_tree *, double> temp_node =
        up_node->jump_pointers->at(exponent);
    max_edge = max(max_edge, temp_node.second);
    up_node = temp_node.first;
    up_depth -= pow(2, exponent);
  }
  return make_pair(up_node, max_edge);
}

// ==============================
// New RMQ + DSU Cartesian max-edge query
// ==============================
pair<shortest_path_tree *, double>
find_LCA(vector<shortest_path_tree *> *sp_tree, unsigned int index_i,
         unsigned int index_j) {
  if (!sp_tree) return make_pair((shortest_path_tree *)nullptr, 0.0);
  if (index_i >= sp_tree->size() || index_j >= sp_tree->size()) {
    return make_pair((shortest_path_tree *)nullptr, 0.0);
  }

  shortest_path_tree *node_i = sp_tree->at(index_i);
  shortest_path_tree *node_j = sp_tree->at(index_j);
  if (!node_i || !node_j) {
    return make_pair((shortest_path_tree *)nullptr, 0.0);
  }

  shortest_path_tree *root = node_i->root;
  if (!root) root = node_j->root;
  if (!root) {
    // should never happen if Dijkstra set root pointers
    return make_pair((shortest_path_tree *)nullptr, 0.0);
  }

  if (root->rmq_vertex == nullptr || root->cartesian == nullptr) {
    cerr << "[WARN] preprocess triggered in find_LCA for root " << root->id << "\n";
    preprocess_spt_queries(sp_tree);
  }

  // Vertex-LCA via RMQ on SPT Euler tour
  RMQLCA &v_rmq = *root->rmq_vertex;
  int fi = v_rmq.first[(int)index_i];
  int fj = v_rmq.first[(int)index_j];
  if (fi == -1 || fj == -1) {
    return make_pair((shortest_path_tree *)nullptr, 0.0);
  }
  int pos = rmq_query_pos(v_rmq, fi, fj);
  int lca_id = v_rmq.euler[pos];
  shortest_path_tree *lca_node = sp_tree->at((unsigned int)lca_id);

  // Max-edge via LCA in DSU-built Cartesian tree
  CartesianPreproc &cp = *root->cartesian;
  int leaf_i = cp.leaf[(int)index_i];
  int leaf_j = cp.leaf[(int)index_j];
  int cfi = cp.rmq.first[leaf_i];
  int cfj = cp.rmq.first[leaf_j];
  int cpos = rmq_query_pos(cp.rmq, cfi, cfj);
  int ct_lca_idx = cp.rmq.euler[cpos];
  double max_edge = cp.nodes[ct_lca_idx]->weight;

  return make_pair(lca_node, max_edge);
}

unsigned int find_HCA(vector<shortest_path_tree*> *sp_tree_i,
                      vector<shortest_path_tree*> *sp_tree_j,
                      unsigned int index,
                      vector<unsigned int>& memoized)
{
    const unsigned int UNK = numeric_limits<unsigned int>::max();

    if (index >= memoized.size()) return index; // defensive
    if (memoized[index] != UNK) return memoized[index];

    unsigned int index_i = index;
    unsigned int index_j = index;

    vector<unsigned int> path_common;
    path_common.reserve(32);

    while (true) {
        // If they are not the same vertex id anymore, we stop.
        if (index_i != index_j) break;

        // record this common vertex
        path_common.push_back(index_i);

        // Try to go one step up in BOTH trees.
        shortest_path_tree* ni = sp_tree_i->at(index_i);
        shortest_path_tree* nj = sp_tree_j->at(index_j);

        // If either is null, or either has no parent jump pointer => can't go up further.
        if (!ni || !nj) break;
        if (!ni->jump_pointers || ni->jump_pointers->empty()) break;
        if (!nj->jump_pointers || nj->jump_pointers->empty()) break;

        unsigned int next_i = ni->jump_pointers->at(0).first->id;
        unsigned int next_j = nj->jump_pointers->at(0).first->id;

        // If either "parent" is itself (shouldn't happen, but prevents infinite loops)
        if (next_i == index_i || next_j == index_j) break;

        index_i = next_i;
        index_j = next_j;
    }

    // If we never recorded anything (shouldn't, because first check index_i==index_j),
    // fall back to returning index itself.
    unsigned int answer = path_common.empty() ? index : path_common.back();

    // memoize all nodes on the common prefix
    for (unsigned int v : path_common) {
        if (v < memoized.size()) memoized[v] = answer;
    }
    if (index < memoized.size()) memoized[index] = answer;

    return answer;
}

// // ==============================
// // Existing HCA (changed)
// // ==============================
// unsigned int find_HCA(vector<shortest_path_tree *> *sp_tree_i,
//                       vector<shortest_path_tree *> *sp_tree_j,
//                       unsigned int index) {
//   if (!sp_tree_i || !sp_tree_j) return index;
//   if (index >= sp_tree_i->size() || index >= sp_tree_j->size()) return index;

//   unsigned int node_index = index;

//   shortest_path_tree *node_i = sp_tree_i->at(node_index);
//   shortest_path_tree *node_j = sp_tree_j->at(node_index);
//   if (!node_i || !node_j) return node_index;

//   // Repeatedly jump upward from `node_index` as long as BOTH SPTs share
//   // the same ancestor at some 2^exp distance (common suffix of paths).
//   while (true) {
//     unsigned int depth = min(node_i->depth, node_j->depth);
//     if (depth == 0) return node_index;

//     int exp = static_cast<int>(get_exponent(depth));
//     bool moved = false;

//     for (; exp >= 0; --exp) {
//       if (!node_i->jump_pointers || !node_j->jump_pointers) continue;
//       if (static_cast<size_t>(exp) >= node_i->jump_pointers->size()) continue;
//       if (static_cast<size_t>(exp) >= node_j->jump_pointers->size()) continue;

//       shortest_path_tree *anc_i = node_i->jump_pointers->at(exp).first;
//       shortest_path_tree *anc_j = node_j->jump_pointers->at(exp).first;
//       if (!anc_i || !anc_j) continue;

//       if (anc_i->id == anc_j->id) {
//         node_index = anc_i->id;
//         node_i = sp_tree_i->at(node_index);
//         node_j = sp_tree_j->at(node_index);
//         moved = true;
//         break; // restart from the new node_index (may climb further)
//       }
//     }

//     if (!moved) return node_index;
//   }
// }
