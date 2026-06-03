// spt.h

#pragma once
#include <vector>
#include <utility>
#include <unordered_map>
#include <random>
#include <cstddef>
#include <limits>

// ---------------- Graph ----------------

struct SimpleGraph {
    int n;
    std::vector<std::vector<std::pair<int,double>>> adj;

    SimpleGraph(int n = 0) : n(n), adj(n) {}

    void addEdge(int u, int v, double w) {
        adj[u].push_back({v,w});
        adj[v].push_back({u,w});  // undirected
    }
};

enum class RefineMethod {
    JumpPointers = 1,
    Ladders      = 2,
    Cartesian    = 3,   // Method 3: Cartesian tree over free tree
    HCA          = 4
};

// ---------------- SPT ----------------

class SPT {
public:
    // ---- Ladder for method 2 ----
    struct Ladder {
        std::vector<int>    nodes;   // nodes in this ladder (top..bottom)
        std::vector<double> edgeW;   // edgeW[i] = w(parent(nodes[i]) -> nodes[i]); edgeW[0]=0
        int logLen = 0;
        std::vector<std::vector<double>> rmq;  // RMQ over edgeW

        void buildRMQ();
        double rangeMax(int l, int r) const;   // max on [l, r] in edgeW
    };

    int n;
    int root;

    // Tree structure
    std::vector<int>    parent;
    std::vector<int>    depth;
    std::vector<double> dist;         // distance from current root in *this* rooted tree
    std::vector<double> parentEdgeW;  // weight(parent[v], v), 0 for root

    // Euler tour + RMQ for LCA (methods 1/2 and for distances in stitched LB)
    std::vector<int> euler;
    std::vector<int> depthEuler;
    std::vector<int> firstOcc;
    int logLen = 0;
    std::vector<std::vector<int>> rmq;   // indices into euler[]

    // Jump pointers (method 1)
    int maxLogUp = 0;
    std::vector<std::vector<int>>    up;
    std::vector<std::vector<double>> upMax;

    // Ladders (method 2)
    std::vector<Ladder> ladders;
    std::vector<int> ladderIdOf;
    std::vector<int> posInLadder;

    // Remaining distance from node to deepest leaf in this rooted SPT (for HCA)
    std::vector<double> downMax;

    // ----- ctor -----
    SPT(int n, int root);

    // Build SPT from graph using Dijkstra from this->root,
    // then build all auxiliary structures.
    void buildFromGraph(const SimpleGraph &g);

    // Variant: build only the auxiliary structures needed by `method`.
    // Saves preprocessing time and memory since each query uses exactly
    // one refinement method. Falls back to the full build for HCA, which
    // requires lifting + Euler-RMQ + downMax. Other methods are leaner.
    void buildFromGraph(const SimpleGraph &g, RefineMethod method);

    // Re-root the existing SPT (seen as an unrooted tree) at newRoot.
    // Uses only (parent, parentEdgeW) of "base" and rebuilds all structures.
    static SPT reRootFrom(const SPT &base, int newRoot);

    // Selective re-root: only builds structures needed by `method`.
    static SPT reRootFrom(const SPT &base, int newRoot, RefineMethod method);

    // Max edge weight on the simple path u–v in this rooted tree
    // (for methods 1 and 2; NOT used by Cartesian)
    double maxEdgeOnPath(int u, int v) const;

    // Combined LB term from this root using the chosen refinement method
    double stitchedLB(int i, int j, RefineMethod method) const;

    // Debug helper
    void printEuler() const;

    // ---- exposed for convenience / testing ----
    int    lca(int u, int v) const;
    double maxEdgeOnPathToAncestor(int v, int anc) const;
    double maxEdgeOnPathToAncestor_ladder(int v, int anc) const;

    // Method 3: “Cartesian” = max edge via Cartesian tree + its own LCA
    double maxEdgeCartesian(int u, int v) const;

private:
    // Helper: recompute all jump / ladder / Euler-RMQ structures and Method 3
    void rebuildStructures();

    // Selective rebuild: only build the substructures needed by `method`.
    void rebuildStructuresFor(RefineMethod method);

    // Iterative Euler tour (drop-in for dfsEuler; safe on long paths). The
    // recursive dfsEuler() is preserved for backward compatibility.
    void buildEulerIterative();

    // Iterative DSU-based Cartesian tree builder (alternative to the
    // recursive O(n^2) buildCartesianRec). Same output structure; faster
    // and stack-safe on degenerate trees.
    void buildCartesianTreeDSU();

    // Euler tour + RMQ for LCA (on the rooted SPT)
    void dfsEuler(int node, int p, int d);
    void buildRMQ();

    // Method 1
    void buildJumpPointers();

    // Method 2
    void buildLadders();

    // Common helper (used by methods 1 & 2)
    double maxEdgeToAncestor(int v, int anc, RefineMethod method) const;

    // Build downward "remaining path" (max distance to a leaf)
    void buildDownMax();

    // ================== Method 3: Cartesian tree data ==================
    struct TreeEdge {
        int u, v;
        double w;
    };

    struct CTNode {
        int edgeIndex;   // index into treeEdges
        double w;        // weight of that edge
        int left;        // left child in Cartesian tree (-1 if none)
        int right;       // right child (-1 if none)
    };

    // Free-tree edges derived from parent[] / parentEdgeW[]
    std::vector<TreeEdge> treeEdges;         // size = n-1 in a tree
    std::vector<std::vector<int>> treeAdjByEdge; // node -> list of incident edge indices

    // Cartesian tree over edges
    std::vector<CTNode> ctNodes;
    int ctRoot = -1;

    // Mapping edge index -> Cartesian tree node index
    std::vector<int> edgeToCTNode;

    // Mapping vertex -> an incident edge of minimal weight (for end-edges)
    std::vector<int> vertexToEdge;

    // LCA on Cartesian tree
    std::vector<int> ctEuler;
    std::vector<int> ctDepthEuler;
    std::vector<int> ctFirstOcc;
    int ctLogLen = 0;
    std::vector<std::vector<int>> ctRMQ;

    // Method 3 builders
    void buildCartesianTree();
    int  buildCartesianRec(const std::vector<int> &edgeSet);
    void buildCartesianLCA();
    void dfsCartesianEuler(int node, int depth);
    int  ctLCA(int a, int b) const;

    // Debug
    void printRMQ() const;  // optional; used only if you call it yourself
};







// ---------------- Stochastic greedy wrapper over SPT roots ----------------

class SPT_StochasticGreedy {
public:
    SPT_StochasticGreedy(const SimpleGraph &g,
                         unsigned int no_roots,
                         unsigned int query_samples,
                         unsigned int seed = 73);

    void get_stochastic_roots(double eps,
                              unsigned int cand_cap = std::numeric_limits<unsigned int>::max());
    void get_exact_roots();

    double lookup(unsigned int u, unsigned int v) const;
    std::vector<double>* lookup_multiple(unsigned int u, unsigned int v) const;

    size_t _sizeof() const;

    const std::vector<int>& selected_roots() const { return root_order_; }

private:
    SimpleGraph g_;
    unsigned int n_;
    unsigned int no_roots_;
    unsigned int query_samples_;
    unsigned int seed_;

    std::vector<std::pair<unsigned int, unsigned int>> train_queries_;
    std::vector<double> current_best_;
    std::vector<int> root_order_;

    mutable std::unordered_map<int, SPT*> forest_;

    void build_training_queries(std::mt19937 &rng);
    std::vector<int> sample_candidate_roots(unsigned int count,
                                            std::mt19937 &rng) const;
    SPT* ensure_root(int r) const;
    double marginal_gain_for_root(int r) const;
    static bool contains_root(const std::vector<int> &roots, int r);
};

