// spt.cpp

#include "spt.h"
#include <queue>
#include <limits>
#include <cmath>
#include <iostream>
#include <functional>
#include <algorithm>

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

void SPT::buildFromGraph(const Graph &g) {
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






