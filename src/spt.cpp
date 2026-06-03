// spt.cpp

#include "spt.h"
#include <queue>
#include <limits>
#include <cmath>
#include <iostream>
#include <functional>
#include <algorithm>
#include <unordered_set>

using namespace std;

// ==================== Ladder methods ====================

void SPT::Ladder::buildRMQ() {
    int m = (int)edgeW.size();
    if (m == 0) {
        logLen = 0;
        return;
    }
    logLen = (int)floor(log2(m)) + 1;
    rmq.assign(logLen, vector<double>(m, 0.0));

    for (int i = 0; i < m; ++i) {
        rmq[0][i] = edgeW[i];
    }
    for (int k = 1; k < logLen; ++k) {
        int len  = 1 << k;
        int half = len >> 1;
        for (int i = 0; i + len <= m; ++i) {
            rmq[k][i] = max(rmq[k-1][i], rmq[k-1][i + half]);
        }
    }
}

double SPT::Ladder::rangeMax(int l, int r) const {
    if (l > r) return 0.0;
    int len = r - l + 1;
    int k = (int)floor(log2(len));
    double a = rmq[k][l];
    double b = rmq[k][r - (1 << k) + 1];
    return max(a, b);
}

// ==================== SPT methods ====================

SPT::SPT(int n_, int root_) : n(n_), root(root_) {
    parent.assign(n, -1);
    depth.assign(n, 0);
    dist.assign(n, numeric_limits<double>::infinity());
    parentEdgeW.assign(n, 0.0);
    firstOcc.assign(n, -1);
    logLen = 0;
    maxLogUp = 0;
}

// --------- Build from graph (Dijkstra) ---------

void SPT::buildFromGraph(const SimpleGraph &g) {
    using P = pair<double,int>;
    priority_queue<P, vector<P>, greater<P>> pq;

    // init
    for (int v = 0; v < n; ++v) {
        dist[v] = numeric_limits<double>::infinity();
        parent[v] = -1;
        parentEdgeW[v] = 0.0;
        depth[v] = 0;
    }
    dist[root] = 0.0;
    pq.push({0.0, root});

    while (!pq.empty()) {
        auto [cdist, u] = pq.top();
        pq.pop();
        if (cdist > dist[u]) continue;

        for (auto &e : g.adj[u]) {
            int v = e.first;
            double w = e.second;
            if (cdist + w < dist[v]) {
                dist[v] = cdist + w;
                parent[v] = u;
                parentEdgeW[v] = w;
                depth[v] = (u == -1 ? 0 : depth[u] + 1);
                pq.push({dist[v], v});
            }
        }
    }

    // rebuild all auxiliary structures for this root
    rebuildStructures();
}

// Selective build: same Dijkstra, lighter post-processing.
void SPT::buildFromGraph(const SimpleGraph &g, RefineMethod method) {
    using P = pair<double,int>;
    priority_queue<P, vector<P>, greater<P>> pq;

    for (int v = 0; v < n; ++v) {
        dist[v] = numeric_limits<double>::infinity();
        parent[v] = -1;
        parentEdgeW[v] = 0.0;
        depth[v] = 0;
    }
    dist[root] = 0.0;
    pq.push({0.0, root});

    while (!pq.empty()) {
        auto [cdist, u] = pq.top();
        pq.pop();
        if (cdist > dist[u]) continue;
        for (auto &e : g.adj[u]) {
            int v = e.first;
            double w = e.second;
            if (cdist + w < dist[v]) {
                dist[v] = cdist + w;
                parent[v] = u;
                parentEdgeW[v] = w;
                depth[v] = (u == -1 ? 0 : depth[u] + 1);
                pq.push({dist[v], v});
            }
        }
    }

    rebuildStructuresFor(method);
}

// --------- Re-rooting ---------

SPT SPT::reRootFrom(const SPT &base, int newRoot) {
    SPT t(base.n, newRoot);

    // Build undirected adjacency of the base SPT
    vector<vector<pair<int,double>>> tadj(base.n);
    for (int v = 0; v < base.n; ++v) {
        int p = base.parent[v];
        if (p != -1) {
            double w = base.parentEdgeW[v];
            tadj[p].push_back({v, w});
            tadj[v].push_back({p, w});
        }
    }

    // BFS from newRoot -- using explicit parent in queue
    vector<int> visited(base.n, 0);
    queue<pair<int,int>> q;  // (node, parent)
    q.push({newRoot, -1});
    visited[newRoot] = 1;

    t.dist[newRoot] = 0.0;
    t.parent[newRoot] = -1;
    t.depth[newRoot] = 0;
    t.parentEdgeW[newRoot] = 0.0;

    while (!q.empty()) {
        auto [u, p] = q.front();
        q.pop();

        for (auto &pr : tadj[u]) {
            int v = pr.first;
            double w = pr.second;

            if (v == p) continue;       // use BFS parent, not t.parent[u]

            if (!visited[v]) {
                visited[v] = 1;
                t.parent[v] = u;
                t.parentEdgeW[v] = w;
                t.depth[v] = t.depth[u] + 1;
                t.dist[v] = t.dist[u] + w;
                q.push({v, u});
            }
        }
    }

    t.rebuildStructures();
    return t;
}

// Selective re-root: same BFS, lighter post-processing.
SPT SPT::reRootFrom(const SPT &base, int newRoot, RefineMethod method) {
    SPT t(base.n, newRoot);
    vector<vector<pair<int,double>>> tadj(base.n);
    for (int v = 0; v < base.n; ++v) {
        int p = base.parent[v];
        if (p != -1) {
            double w = base.parentEdgeW[v];
            tadj[p].push_back({v, w});
            tadj[v].push_back({p, w});
        }
    }
    vector<int> visited(base.n, 0);
    queue<pair<int,int>> q;
    q.push({newRoot, -1});
    visited[newRoot] = 1;
    t.dist[newRoot] = 0.0;
    t.parent[newRoot] = -1;
    t.depth[newRoot] = 0;
    t.parentEdgeW[newRoot] = 0.0;
    while (!q.empty()) {
        auto [u, p] = q.front();
        q.pop();
        for (auto &pr : tadj[u]) {
            int v = pr.first;
            double w = pr.second;
            if (v == p) continue;
            if (!visited[v]) {
                visited[v] = 1;
                t.parent[v] = u;
                t.parentEdgeW[v] = w;
                t.depth[v] = t.depth[u] + 1;
                t.dist[v] = t.dist[u] + w;
                q.push({v, u});
            }
        }
    }
    t.rebuildStructuresFor(method);
    return t;
}

// --------- Common rebuild of auxiliary structures ---------

void SPT::rebuildStructures() {
    // clear SPT LCA-related
    euler.clear();
    depthEuler.clear();
    firstOcc.assign(n, -1);

    // 1) jump pointers (method 1)
    buildJumpPointers();

    // 2) ladders (method 2)
    buildLadders();

    // 3) Euler tour + RMQ for LCA on the rooted SPT
    dfsEuler(root, -1, 0);
    buildRMQ();

    // 4) remaining-path DP ("downMax") for HCA extras
    buildDownMax();

    // 5) Cartesian tree (method 3) over the underlying free tree
    buildCartesianTree();
}

// Selective rebuild: build only what `method` consumes. Avoids the cost of
// the unused refinement-method machinery at preprocessing time.
//   - JumpPointers needs: jump pointers + Euler-RMQ for the lca() helper
//     used in stitchedLB.
//   - Ladders needs: ladders + Euler-RMQ for lca().
//   - Cartesian needs: Cartesian tree (its own LCA is internal). Euler-RMQ
//     is also needed by stitchedLB to compute path length via lca().
//   - HCA needs: jump pointers + Euler-RMQ + downMax.
void SPT::rebuildStructuresFor(RefineMethod method) {
    euler.clear();
    depthEuler.clear();
    firstOcc.assign(n, -1);

    // Euler tour + RMQ are used by lca() across all methods (stitchedLB
    // calls lca to compute path length L).
    buildEulerIterative();
    buildRMQ();

    if (method == RefineMethod::JumpPointers || method == RefineMethod::HCA) {
        buildJumpPointers();
    }
    if (method == RefineMethod::Ladders) {
        buildLadders();
    }
    if (method == RefineMethod::Cartesian) {
        buildCartesianTreeDSU();
    } else {
        // Other methods don't need Cartesian; skip the O(n^2) recursive build.
        treeEdges.clear();
        ctNodes.clear();
        ctRoot = -1;
        edgeToCTNode.clear();
        vertexToEdge.clear();
        ctEuler.clear();
        ctDepthEuler.clear();
        ctFirstOcc.clear();
        ctRMQ.clear();
        ctLogLen = 0;
    }
    if (method == RefineMethod::HCA) {
        buildDownMax();
    } else {
        downMax.assign(n, 0.0);
    }
}

// Iterative Euler-tour builder. Same output as recursive dfsEuler(root, -1, 0)
// but with an explicit stack so it can't overflow on long paths.
void SPT::buildEulerIterative() {
    // Construct children list once.
    std::vector<std::vector<int>> children(n);
    for (int v = 0; v < n; ++v) {
        if (parent[v] != -1) children[parent[v]].push_back(v);
    }

    struct Frame { int node; int childIdx; int depth; };
    std::vector<Frame> stack;
    stack.reserve(64);
    stack.push_back({root, 0, 0});

    while (!stack.empty()) {
        Frame &f = stack.back();
        if (f.childIdx == 0) {
            // entering: record first occurrence
            if (firstOcc[f.node] == -1) firstOcc[f.node] = (int)euler.size();
            euler.push_back(f.node);
            depthEuler.push_back(f.depth);
        }
        if (f.childIdx < (int)children[f.node].size()) {
            int c = children[f.node][f.childIdx++];
            stack.push_back({c, 0, f.depth + 1});
        } else {
            stack.pop_back();
            if (!stack.empty()) {
                Frame &p = stack.back();
                euler.push_back(p.node);
                depthEuler.push_back(p.depth);
            }
        }
    }
}

// DSU-based Cartesian tree (max-edge variant). Equivalent to
// ShortestPathTree.cpp::build_cartesian_dsu but adapted to spt.cpp's data
// shapes. Builds in O(n alpha(n) + n log n) instead of O(n^2) worst-case
// from the recursive splitter.
void SPT::buildCartesianTreeDSU() {
    treeEdges.clear();
    treeAdjByEdge.assign(n, {});
    for (int v = 0; v < n; ++v) {
        int p = parent[v];
        if (p != -1) {
            TreeEdge e; e.u = v; e.v = p; e.w = parentEdgeW[v];
            int idx = (int)treeEdges.size();
            treeEdges.push_back(e);
            treeAdjByEdge[e.u].push_back(idx);
            treeAdjByEdge[e.v].push_back(idx);
        }
    }
    int m = (int)treeEdges.size();
    ctNodes.clear();
    ctRoot = -1;
    edgeToCTNode.assign(m, -1);
    if (m == 0) {
        ctEuler.clear();
        ctDepthEuler.clear();
        ctFirstOcc.clear();
        ctRMQ.clear();
        ctLogLen = 0;
        vertexToEdge.assign(n, -1);
        return;
    }

    // Cartesian leaves correspond to vertices for the LCA query path.
    // Use a leaf-per-vertex layout consistent with the existing code path.
    vertexToEdge.assign(n, -1);
    for (int ei = 0; ei < m; ++ei) {
        double w = treeEdges[ei].w;
        int u = treeEdges[ei].u, v = treeEdges[ei].v;
        if (vertexToEdge[u] == -1 || w < treeEdges[vertexToEdge[u]].w)
            vertexToEdge[u] = ei;
        if (vertexToEdge[v] == -1 || w < treeEdges[vertexToEdge[v]].w)
            vertexToEdge[v] = ei;
    }

    // Sort edges by ascending weight; union-find merges them.
    std::vector<int> order(m);
    for (int i = 0; i < m; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return treeEdges[a].w < treeEdges[b].w;
    });

    // Per-vertex DSU. comp_root[parent[v]] tracks the current Cartesian
    // representative of that component.
    std::vector<int> dsu_p(n), dsu_r(n, 0);
    for (int i = 0; i < n; ++i) dsu_p[i] = i;
    std::function<int(int)> find = [&](int x) -> int {
        while (dsu_p[x] != x) { dsu_p[x] = dsu_p[dsu_p[x]]; x = dsu_p[x]; }
        return x;
    };
    std::vector<int> comp_root(n, -1); // CT node id for each component

    for (int ei : order) {
        int u = treeEdges[ei].u, v = treeEdges[ei].v;
        int cu = find(u), cv = find(v);
        if (cu == cv) continue;
        CTNode node;
        node.edgeIndex = ei;
        node.w = treeEdges[ei].w;
        node.left = comp_root[cu]; // -1 if leaf component
        node.right = comp_root[cv];
        int thisId = (int)ctNodes.size();
        ctNodes.push_back(node);
        if (dsu_r[cu] < dsu_r[cv]) std::swap(cu, cv);
        dsu_p[cv] = cu;
        if (dsu_r[cu] == dsu_r[cv]) ++dsu_r[cu];
        comp_root[cu] = thisId;
        edgeToCTNode[ei] = thisId;
    }

    int rep0 = find(0);
    ctRoot = comp_root[rep0];

    buildCartesianLCA();
}

// --------- Euler Tour for SPT LCA ---------

void SPT::dfsEuler(int node, int p, int d) {
    if (firstOcc[node] == -1)
        firstOcc[node] = (int)euler.size();

    euler.push_back(node);
    depthEuler.push_back(d);

    for (int child = 0; child < n; ++child) {
        if (parent[child] == node && child != p) {
            dfsEuler(child, node, d + 1);
            euler.push_back(node);
            depthEuler.push_back(d);
        }
    }
}

void SPT::buildRMQ() {
    int m = (int)depthEuler.size();
    if (m == 0) {
        logLen = 0;
        return;
    }

    logLen = (int)floor(log2(m)) + 1;
    rmq.assign(logLen, vector<int>(m));

    for (int i = 0; i < m; ++i)
        rmq[0][i] = i;

    for (int k = 1; k < logLen; ++k) {
        int len  = 1 << k;
        int half = len >> 1;
        for (int i = 0; i + len <= m; ++i) {
            int i1 = rmq[k-1][i];
            int i2 = rmq[k-1][i + half];
            rmq[k][i] = (depthEuler[i1] < depthEuler[i2] ? i1 : i2);
        }
    }
}

int SPT::lca(int u, int v) const {
    int L = firstOcc[u];
    int R = firstOcc[v];
    if (L > R) std::swap(L, R);

    int length = R - L + 1;
    int k = (int)floor(log2(length));
    int i1 = rmq[k][L];
    int i2 = rmq[k][R - (1 << k) + 1];
    int best = (depthEuler[i1] < depthEuler[i2] ? i1 : i2);
    return euler[best];
}

// --------- Jump pointers (method 1) ---------

void SPT::buildJumpPointers() {
    maxLogUp = (int)floor(log2(std::max(1, n))) + 1;
    up.assign(n, vector<int>(maxLogUp, -1));
    upMax.assign(n, vector<double>(maxLogUp, 0.0));

    for (int v = 0; v < n; ++v) {
        up[v][0]    = parent[v];
        upMax[v][0] = parentEdgeW[v];
    }

    for (int j = 1; j < maxLogUp; ++j) {
        for (int v = 0; v < n; ++v) {
            int mid = up[v][j-1];
            if (mid == -1) {
                up[v][j]    = -1;
                upMax[v][j] = upMax[v][j-1];
            } else {
                up[v][j]    = up[mid][j-1];
                upMax[v][j] = max(upMax[v][j-1], upMax[mid][j-1]);
            }
        }
    }
}

double SPT::maxEdgeOnPathToAncestor(int v, int ancestor) const {
    if (v == ancestor) return 0.0;

    double answer = 0.0;
    int diff = depth[v] - depth[ancestor];
    if (diff < 0) {
        cerr << "Error: ancestor deeper than v in maxEdgeOnPathToAncestor\n";
        return 0.0;
    }
    for (int j = maxLogUp - 1; j >= 0; --j) {
        if ((diff >> j) & 1) {
            answer = max(answer, upMax[v][j]);
            v = up[v][j];
            if (v == -1) break;
        }
    }
    return answer;
}

// --------- Ladder decomposition (method 2) ---------

void SPT::buildLadders() {
    vector<vector<int>> children(n);
    for (int v = 0; v < n; ++v) {
        if (parent[v] != -1) {
            children[parent[v]].push_back(v);
        }
    }

    vector<int> subtreeSize(n, 0), heavyChild(n, -1);

    function<void(int)> dfs1 = [&](int u) {
        subtreeSize[u] = 1;
        int bestSize = 0;
        heavyChild[u] = -1;
        for (int c : children[u]) {
            dfs1(c);
            subtreeSize[u] += subtreeSize[c];
            if (subtreeSize[c] > bestSize) {
                bestSize = subtreeSize[c];
                heavyChild[u] = c;
            }
        }
    };
    dfs1(root);

    ladderIdOf.assign(n, -1);
    posInLadder.assign(n, -1);
    ladders.clear();

    function<void(int,int)> dfs2 = [&](int u, int currentLadder) {
        if (currentLadder == -1) {
            currentLadder = (int)ladders.size();
            ladders.emplace_back();
        }
        Ladder &lad = ladders[currentLadder];
        int pos = (int)lad.nodes.size();
        lad.nodes.push_back(u);
        ladderIdOf[u] = currentLadder;
        posInLadder[u] = pos;

        int hc = heavyChild[u];
        if (hc != -1) dfs2(hc, currentLadder);
        for (int c : children[u]) {
            if (c == hc) continue;
            dfs2(c, -1);
        }
    };
    dfs2(root, -1);

    for (auto &lad : ladders) {
        int m = (int)lad.nodes.size();
        lad.edgeW.assign(m, 0.0);
        for (int i = 1; i < m; ++i) {
            int v = lad.nodes[i];
            lad.edgeW[i] = parentEdgeW[v];
        }
        lad.buildRMQ();
    }
}

// --------- Downward remaining-path DP (for HCA) ---------

void SPT::buildDownMax() {
    // children adjacency derived from parent[]
    std::vector<std::vector<int>> children(n);
    for (int v = 0; v < n; ++v) {
        int p = parent[v];
        if (p != -1) {
            children[p].push_back(v);
        }
    }

    downMax.assign(n, 0.0);

    // post-order DFS to fill downMax[v] = max over children c of (w(v,c) + downMax[c])
    std::function<void(int)> dfs = [&](int u) {
        double best = 0.0;
        for (int c : children[u]) {
            dfs(c);
            // edge weight from u -> c is parentEdgeW[c]
            best = std::max(best, parentEdgeW[c] + downMax[c]);
        }
        downMax[u] = best;
    };
    dfs(root);
}

// --------- Max edge on path u–v (generic: methods 1/2) ---------

double SPT::maxEdgeOnPath(int u, int v) const {
    int c = lca(u, v);
    double m1 = maxEdgeOnPathToAncestor(u, c);
    double m2 = maxEdgeOnPathToAncestor(v, c);
    return max(m1, m2);
}

// --------- Common helper for stitchedLB ---------

double SPT::maxEdgeToAncestor(int v, int anc, RefineMethod method) const {
    if (method == RefineMethod::JumpPointers)
        return maxEdgeOnPathToAncestor(v, anc);
    else if (method == RefineMethod::Ladders)
        return maxEdgeOnPathToAncestor_ladder(v, anc);
    else
        return 0.0;   // Cartesian handled separately
}

// Ladder-based ancestor max (method 2)

double SPT::maxEdgeOnPathToAncestor_ladder(int v, int ancestor) const {
    if (v == ancestor) return 0.0;

    double ans = 0.0;
    int u = v;

    while (u != -1 && ladderIdOf[u] != ladderIdOf[ancestor]) {
        int lid = ladderIdOf[u];
        const Ladder &lad = ladders[lid];

        int posU   = posInLadder[u];
        int posTop = 0;

        if (posU >= posTop + 1) {
            ans = max(ans, lad.rangeMax(posTop + 1, posU));
        }

        int topNode = lad.nodes[posTop];
        if (parent[topNode] != -1) {
            ans = max(ans, parentEdgeW[topNode]);
        }
        u = parent[topNode];
    }

    if (u == -1) return ans;

    int lid = ladderIdOf[u];
    const Ladder &lad = ladders[lid];
    int posU = posInLadder[u];
    int posA = posInLadder[ancestor];

    if (posU > posA) {
        ans = max(ans, lad.rangeMax(posA + 1, posU));
    }
    return ans;
}

// --------- Stitched LB from this root (methods 1–3) ---------

double SPT::stitchedLB(int i, int j, RefineMethod method) const {
    int c = lca(i, j);
    double L = dist[i] + dist[j] - 2.0 * dist[c];

    double M = 0.0;
    if (method == RefineMethod::Cartesian) {
        // method 3: fully independent Cartesian tree query
        M = maxEdgeCartesian(i, j);
    } else {
        // methods 1 & 2: jump pointers / ladders
        double max_i = maxEdgeToAncestor(i, c, method);
        double max_j = maxEdgeToAncestor(j, c, method);
        M = max(max_i, max_j);
    }
    return 2.0 * M - L;
}

// --------- Debug printing ---------

void SPT::printEuler() const {
    cout << "Euler Tour (node, depth):\n";
    for (size_t i = 0; i < euler.size(); ++i) {
        cout << "  idx " << i
             << ": node " << euler[i]
             << ", depth " << depthEuler[i] << "\n";
    }
    cout << "\nFirst occurrence of each node:\n";
    for (int v = 0; v < n; ++v) {
        if (firstOcc[v] != -1)
            cout << "  node " << v << " -> firstOcc idx " << firstOcc[v] << "\n";
    }
    cout << endl;
}

void SPT::printRMQ() const {
    cout << "RMQ sparse table (stores indices into Euler/Depth arrays):\n";
    int m = (int)depthEuler.size();
    for (int k = 0; k < logLen; ++k) {
        cout << "  k = " << k
             << " (interval length = " << (1 << k) << "):\n    ";
        for (int i = 0; i + (1 << k) <= m; ++i) {
            int idx = rmq[k][i];
            cout << "[" << i << "] -> eulerIdx " << idx
                 << " (node " << euler[idx]
                 << ", depth " << depthEuler[idx] << ")  ";
        }
        cout << "\n";
    }
    cout << endl;
}

// ==================== Method 3: Cartesian Tree over free tree ====================

void SPT::buildCartesianTree() {
    // Build undirected free-tree edges from parent[] / parentEdgeW[]
    treeEdges.clear();
    treeAdjByEdge.assign(n, {});

    for (int v = 0; v < n; ++v) {
        int p = parent[v];
        if (p != -1) {
            TreeEdge e;
            e.u = v;
            e.v = p;
            e.w = parentEdgeW[v];
            int idx = (int)treeEdges.size();
            treeEdges.push_back(e);
            treeAdjByEdge[e.u].push_back(idx);
            treeAdjByEdge[e.v].push_back(idx);
        }
    }

    int m = (int)treeEdges.size();
    ctNodes.clear();
    ctRoot = -1;
    edgeToCTNode.assign(m, -1);

    // If there are no edges (single-vertex tree), nothing to do.
    if (m == 0) {
        ctEuler.clear();
        ctDepthEuler.clear();
        ctFirstOcc.clear();
        ctRMQ.clear();
        ctLogLen = 0;
        vertexToEdge.assign(n, -1);
        return;
    }

    // Build vertex -> minimal incident edge mapping
    vertexToEdge.assign(n, -1);
    for (int ei = 0; ei < m; ++ei) {
        double w = treeEdges[ei].w;
        int u = treeEdges[ei].u;
        int v = treeEdges[ei].v;

        if (vertexToEdge[u] == -1 || w < treeEdges[ vertexToEdge[u] ].w)
            vertexToEdge[u] = ei;
        if (vertexToEdge[v] == -1 || w < treeEdges[ vertexToEdge[v] ].w)
            vertexToEdge[v] = ei;
    }

    // Build Cartesian tree recursively over all edges
    std::vector<int> allEdges(m);
    for (int ei = 0; ei < m; ++ei) allEdges[ei] = ei;

    ctRoot = buildCartesianRec(allEdges);

    // Build edge -> Cartesian node mapping
    for (int nodeId = 0; nodeId < (int)ctNodes.size(); ++nodeId) {
        int eIdx = ctNodes[nodeId].edgeIndex;
        edgeToCTNode[eIdx] = nodeId;
    }

    // Build LCA structures for the Cartesian tree itself
    buildCartesianLCA();
}

// Recursive builder: edgeSet is a list of edge indices inside the current component
int SPT::buildCartesianRec(const std::vector<int> &edgeSet) {
    if (edgeSet.empty()) return -1;

    // 1) Find edge with maximum weight
    int bestEdge = edgeSet[0];
    double bestW = treeEdges[bestEdge].w;
    for (int ei : edgeSet) {
        if (treeEdges[ei].w > bestW) {
            bestW = treeEdges[ei].w;
            bestEdge = ei;
        }
    }

    int thisId = (int)ctNodes.size();
    ctNodes.push_back(CTNode{bestEdge, bestW, -1, -1});

    // If only one edge, there are no subtrees below this
    if (edgeSet.size() == 1) return thisId;

    // 2) Split remaining edges into two components after "removing" bestEdge
    int m = (int)treeEdges.size();
    std::vector<char> inSet(m, 0);
    for (int ei : edgeSet) inSet[ei] = 1;
    inSet[bestEdge] = 0; // remove

    // BFS from one endpoint of bestEdge (say u0)
    int u0 = treeEdges[bestEdge].u;
    std::vector<char> visited(n, 0);
    std::vector<int> stack;
    stack.push_back(u0);
    visited[u0] = 1;

    while (!stack.empty()) {
        int x = stack.back();
        stack.pop_back();
        for (int ei : treeAdjByEdge[x]) {
            if (!inSet[ei]) continue;  // not part of this subproblem
            const TreeEdge &e = treeEdges[ei];
            int y = (e.u == x ? e.v : e.u);
            if (!visited[y]) {
                visited[y] = 1;
                stack.push_back(y);
            }
        }
    }

    std::vector<int> leftEdges;
    std::vector<int> rightEdges;
    leftEdges.reserve(edgeSet.size());
    rightEdges.reserve(edgeSet.size());

    for (int ei : edgeSet) {
        if (ei == bestEdge) continue;
        const TreeEdge &e = treeEdges[ei];
        bool a = visited[e.u];
        bool b = visited[e.v];
        if (a && b) {
            leftEdges.push_back(ei);
        } else {
            rightEdges.push_back(ei);
        }
    }

    int leftId  = buildCartesianRec(leftEdges);
    int rightId = buildCartesianRec(rightEdges);
    ctNodes[thisId].left  = leftId;
    ctNodes[thisId].right = rightId;

    return thisId;
}

// --------- Cartesian tree LCA ---------

void SPT::dfsCartesianEuler(int node, int depth) {
    if (node == -1) return;

    if (ctFirstOcc[node] == -1) {
        ctFirstOcc[node] = (int)ctEuler.size();
    }
    ctEuler.push_back(node);
    ctDepthEuler.push_back(depth);

    int L = ctNodes[node].left;
    if (L != -1) {
        dfsCartesianEuler(L, depth + 1);
        ctEuler.push_back(node);
        ctDepthEuler.push_back(depth);
    }
    int R = ctNodes[node].right;
    if (R != -1) {
        dfsCartesianEuler(R, depth + 1);
        ctEuler.push_back(node);
        ctDepthEuler.push_back(depth);
    }
}

void SPT::buildCartesianLCA() {
    ctEuler.clear();
    ctDepthEuler.clear();
    ctFirstOcc.assign(ctNodes.size(), -1);
    ctRMQ.clear();
    ctLogLen = 0;

    if (ctRoot == -1 || ctNodes.empty())
        return;

    dfsCartesianEuler(ctRoot, 0);

    int m = (int)ctDepthEuler.size();
    if (m == 0) return;

    ctLogLen = (int)std::floor(std::log2(m)) + 1;
    ctRMQ.assign(ctLogLen, std::vector<int>(m));
    for (int i = 0; i < m; ++i) {
        ctRMQ[0][i] = i;
    }
    for (int k = 1; k < ctLogLen; ++k) {
        int len  = 1 << k;
        int half = len >> 1;
        for (int i = 0; i + len <= m; ++i) {
            int i1 = ctRMQ[k-1][i];
            int i2 = ctRMQ[k-1][i + half];
            ctRMQ[k][i] = (ctDepthEuler[i1] < ctDepthEuler[i2]) ? i1 : i2;
        }
    }
}

int SPT::ctLCA(int a, int b) const {
    if (a == b) return a;
    int L = ctFirstOcc[a];
    int R = ctFirstOcc[b];
    if (L > R) std::swap(L, R);
    int length = R - L + 1;
    int k = (int)std::floor(std::log2(length));
    int i1 = ctRMQ[k][L];
    int i2 = ctRMQ[k][R - (1 << k) + 1];
    int best = (ctDepthEuler[i1] < ctDepthEuler[i2]) ? i1 : i2;
    return ctEuler[best];
}

// Override the placeholder maxEdgeCartesian with the real one (same signature)
double SPT::maxEdgeCartesian(int u, int v) const {
    // Trivial case
    if (u == v) return 0.0;

    if (ctRoot == -1 || treeEdges.empty()) return 0.0;

    int m = (int)treeEdges.size();
    if ((int)edgeToCTNode.size() != m) return 0.0;
    if ((int)vertexToEdge.size() != n) return 0.0;

    int e1 = (u >= 0 && u < n) ? vertexToEdge[u] : -1;
    int e2 = (v >= 0 && v < n) ? vertexToEdge[v] : -1;
    if (e1 < 0 || e2 < 0) return 0.0;

    int n1 = edgeToCTNode[e1];
    int n2 = edgeToCTNode[e2];
    if (n1 < 0 || n2 < 0) return 0.0;

    int ctnode = ctLCA(n1, n2);
    if (ctnode < 0) return 0.0;

    return ctNodes[ctnode].w;  // max edge weight on path u–v
}








// ==================== Stochastic greedy over SPT roots ====================

SPT_StochasticGreedy::SPT_StochasticGreedy(const SimpleGraph &g,
                                           unsigned int no_roots,
                                           unsigned int query_samples,
                                           unsigned int seed)
    : g_(g),
      n_(g.n),
      no_roots_(no_roots),
      query_samples_(query_samples),
      seed_(seed) {
    std::mt19937 rng(seed_);
    build_training_queries(rng);
    current_best_.assign(train_queries_.size(), 0.0);
}

void SPT_StochasticGreedy::build_training_queries(std::mt19937 &rng) {
    train_queries_.clear();
    if (n_ < 2 || query_samples_ == 0) return;

    std::uniform_int_distribution<unsigned int> uni(0, n_ - 1);

    auto has_edge = [&](unsigned int u, unsigned int v) -> bool {
        for (const auto &pr : g_.adj[u]) {
            if ((unsigned int)pr.first == v) return true;
        }
        return false;
    };

    size_t attempts = 0;
    size_t max_attempts = std::max<size_t>(1000, (size_t)query_samples_ * 50);
    while (train_queries_.size() < query_samples_ && attempts < max_attempts) {
        ++attempts;
        unsigned int u = uni(rng);
        unsigned int v = uni(rng);
        if (u == v) continue;
        if (u > v) std::swap(u, v);
        if (has_edge(u, v)) continue;
        train_queries_.push_back({u, v});
    }

    if (train_queries_.empty()) {
        for (unsigned int u = 0; u < n_; ++u) {
            for (unsigned int v = u + 1; v < n_; ++v) {
                if (!has_edge(u, v)) {
                    train_queries_.push_back({u, v});
                    if (train_queries_.size() >= query_samples_) return;
                }
            }
        }
    }
}

bool SPT_StochasticGreedy::contains_root(const std::vector<int> &roots, int r) {
    for (int x : roots) if (x == r) return true;
    return false;
}

std::vector<int>
SPT_StochasticGreedy::sample_candidate_roots(unsigned int count,
                                             std::mt19937 &rng) const {
    std::vector<int> remaining;
    remaining.reserve(n_);
    for (unsigned int r = 0; r < n_; ++r) {
        if (!contains_root(root_order_, (int)r)) remaining.push_back((int)r);
    }

    std::shuffle(remaining.begin(), remaining.end(), rng);
    if (count < remaining.size()) remaining.resize(count);
    return remaining;
}

SPT* SPT_StochasticGreedy::ensure_root(int r) const {
    auto it = forest_.find(r);
    if (it != forest_.end()) return it->second;

    SPT *t = new SPT((int)n_, r);
    // Use the method-specific build to take the iterative O(n) Euler-tour
    // path (buildEulerIterative) instead of the legacy recursive O(n^2)
    // dfsEuler in rebuildStructures(). All consumers (stitchedLB, kNN,
    // RMSE/ERC eval) call stitchedLB with RefineMethod::JumpPointers.
    // The single-arg buildFromGraph() routes to rebuildStructures() which:
    //   1. uses recursive dfsEuler() (stack overflow on long paths -> SIGSEGV)
    //   2. has O(n^2) child-scan instead of O(n) Euler walk
    //   3. builds Ladders + Cartesian preproc that JumpPointers never reads
    // The method-specific version fixes all three.
    t->buildFromGraph(g_, RefineMethod::JumpPointers);
    forest_[r] = t;
    return t;
}

double SPT_StochasticGreedy::marginal_gain_for_root(int r) const {
    SPT *t = ensure_root(r);
    double gain = 0.0;

    for (size_t qi = 0; qi < train_queries_.size(); ++qi) {
        unsigned int u = train_queries_[qi].first;
        unsigned int v = train_queries_[qi].second;
        double cand = t->stitchedLB((int)u, (int)v, RefineMethod::JumpPointers);
        if (cand > current_best_[qi]) gain += (cand - current_best_[qi]);
    }
    return gain;
}

void SPT_StochasticGreedy::get_stochastic_roots(double eps,
                                                unsigned int cand_cap) {
    if (eps <= 0.0 || eps >= 1.0) eps = 0.1;
    std::mt19937 rng(seed_);
    root_order_.clear();
    current_best_.assign(train_queries_.size(), 0.0);

    for (unsigned int iter = 0; iter < no_roots_; ++iter) {
        unsigned int remaining_budget = no_roots_ - iter;
        unsigned int remaining_candidates = n_ - (unsigned int)root_order_.size();
        if (remaining_candidates == 0) break;

        unsigned int sample_size = (unsigned int)std::ceil(
            ((double)remaining_candidates / (double)std::max(1u, remaining_budget)) *
            std::log(1.0 / eps)
        );
        sample_size = std::max(1u, std::min(sample_size, remaining_candidates));
        // Cap to bound the number of full SPT builds per iteration. Each
        // sampled root triggers one Dijkstra + RMQ + Cartesian preprocess
        // (~0.5-2s for 50000 nodes), so an uncapped sample of e.g. 50000/8
        // would mean ~6000 SPT builds per iteration. Default UINT_MAX
        // preserves original unbounded behavior.
        if (cand_cap > 0 && sample_size > cand_cap) sample_size = cand_cap;

        std::vector<int> candidates = sample_candidate_roots(sample_size, rng);

        // Parallel-build SPTs for sampled candidates that aren't cached.
        // Each thread writes to a different forest_[r], so writes don't
        // alias. forest_ is an unordered_map, so we serialize the insert
        // step through a critical section after construction. Per-tree
        // build is ~95% of the cost; the critical at the end is negligible.
        {
            std::vector<int> to_build;
            to_build.reserve(candidates.size());
            for (int r : candidates) {
                if (forest_.find(r) == forest_.end()) to_build.push_back(r);
            }
            // Build trees into a parallel-local scratch, then insert.
            // Use the method-specific build to take the iterative O(n)
            // Euler-tour path. The single-arg build calls dfsEuler() which
            // is BOTH recursive (segfaults on deep trees in OpenMP worker
            // threads due to limited stack) AND O(n^2) (scans all n entries
            // per visited node to find children). On n=50000 SIFT, a single
            // sequential build took ~41 minutes; in parallel it segfaulted
            // from stack exhaustion. RefineMethod::JumpPointers builds only
            // what every downstream consumer (stitchedLB, marginal_gain,
            // lookup) actually reads.
            std::vector<SPT*> built(to_build.size(), nullptr);
            #pragma omp parallel for schedule(dynamic, 1) if (to_build.size() > 4)
            for (long long i = 0; i < (long long)to_build.size(); ++i) {
                SPT *t = new SPT((int)n_, to_build[i]);
                t->buildFromGraph(g_, RefineMethod::JumpPointers);
                built[i] = t;
            }
            for (size_t i = 0; i < to_build.size(); ++i) {
                forest_[to_build[i]] = built[i];
            }
        }

        // Marginal-gain evaluation parallel across candidates. Each thread
        // accumulates its own best; reduce in a critical section.
        double best_gain = -1.0;
        int best_root = -1;
        {
            #pragma omp parallel
            {
                double local_best_gain = -1.0;
                int local_best_root = -1;
                #pragma omp for schedule(dynamic, 16) nowait
                for (long long ci = 0; ci < (long long)candidates.size(); ++ci) {
                    int r = candidates[ci];
                    double gain = marginal_gain_for_root(r);
                    if (gain > local_best_gain) {
                        local_best_gain = gain;
                        local_best_root = r;
                    }
                }
                #pragma omp critical
                {
                    if (local_best_gain > best_gain) {
                        best_gain = local_best_gain;
                        best_root = local_best_root;
                    }
                }
            }
        }

        if (best_root < 0) break;
        root_order_.push_back(best_root);

        SPT *best_tree = ensure_root(best_root);
        for (size_t qi = 0; qi < train_queries_.size(); ++qi) {
            unsigned int u = train_queries_[qi].first;
            unsigned int v = train_queries_[qi].second;
            double cand = best_tree->stitchedLB((int)u, (int)v, RefineMethod::JumpPointers);
            if (cand > current_best_[qi]) current_best_[qi] = cand;
        }

        // INTRA-ITER CLEANUP: free SPTs of non-winning candidates. Keep
        // only trees rooted at selected landmarks (current and previous).
        // Without this, forest_ grows monotonically and reaches >> 100 GB
        // for production runs at n=50000 (each SPT is ~25-35 MB).
        {
            std::unordered_set<int> keep(root_order_.begin(), root_order_.end());
            for (auto it = forest_.begin(); it != forest_.end(); ) {
                if (keep.find(it->first) == keep.end()) {
                    delete it->second;
                    it = forest_.erase(it);
                } else {
                    ++it;
                }
            }
        }
    }
}

void SPT_StochasticGreedy::get_exact_roots() {
    root_order_.clear();
    current_best_.assign(train_queries_.size(), 0.0);

    for (unsigned int iter = 0; iter < no_roots_; ++iter) {
        double best_gain = -1.0;
        int best_root = -1;
        for (unsigned int r = 0; r < n_; ++r) {
            if (contains_root(root_order_, (int)r)) continue;
            double gain = marginal_gain_for_root((int)r);
            if (gain > best_gain) {
                best_gain = gain;
                best_root = (int)r;
            }
        }

        if (best_root < 0) break;
        root_order_.push_back(best_root);

        SPT *best_tree = ensure_root(best_root);
        for (size_t qi = 0; qi < train_queries_.size(); ++qi) {
            unsigned int u = train_queries_[qi].first;
            unsigned int v = train_queries_[qi].second;
            double cand = best_tree->stitchedLB((int)u, (int)v, RefineMethod::JumpPointers);
            if (cand > current_best_[qi]) current_best_[qi] = cand;
        }
    }
}

double SPT_StochasticGreedy::lookup(unsigned int u, unsigned int v) const {
    double best = 0.0;
    for (int r : root_order_) {
        SPT *t = ensure_root(r);
        best = std::max(best, t->stitchedLB((int)u, (int)v, RefineMethod::JumpPointers));
    }
    return best;
}

std::vector<double>* SPT_StochasticGreedy::lookup_multiple(unsigned int u, unsigned int v) const {
    std::vector<double> *vals = new std::vector<double>();
    vals->reserve(root_order_.size());

    double running_max = 0.0;
    for (int r : root_order_) {
        SPT *t = ensure_root(r);
        running_max = std::max(running_max, t->stitchedLB((int)u, (int)v, RefineMethod::JumpPointers));
        vals->push_back(running_max);
    }
    return vals;
}

size_t SPT_StochasticGreedy::_sizeof() const {
    size_t total = sizeof(*this);
    total += train_queries_.capacity() * sizeof(std::pair<unsigned int, unsigned int>);
    total += current_best_.capacity() * sizeof(double);
    total += root_order_.capacity() * sizeof(int);
    total += forest_.size() * (sizeof(int) + sizeof(SPT*));
    return total;
}


