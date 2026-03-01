// main.cpp

#include "spt.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// small helper to print LB components for methods 1–3
void run_basic_three(const SPT &T2, const SPT &T4, RefineMethod method) {
    int o2 = 2, o3 = 3, o6 = 6, o7 = 7;

    cout << "\n=== SPT rooted at o2 ===\n";
    T2.printEuler();

    cout << "=== SPT rooted at o4 ===\n";
    T4.printEuler();

    double d2_3 = T2.dist[o3];
    double d2_7 = T2.dist[o7];
    double C2   = T2.stitchedLB(o7, o3, method);

    cout << "From root o2:\n";
    cout << "  dist(o2,o3) = " << d2_3 << "\n";
    cout << "  dist(o2,o7) = " << d2_7 << "\n";
    cout << "  stitched LB term C2 = " << C2 << "\n\n";

    double d6_3 = T4.dist[o3];
    double d6_7 = T4.dist[o7];
    double C4   = T4.stitchedLB(o7, o3, method);

    cout << "From root o4:\n";
    cout << "  dist(o4,o3) = " << d6_3 << "\n";
    cout << "  dist(o4,o7) = " << d6_7 << "\n";
    cout << "  stitched LB term C4 = " << C4 << "\n\n";

    // Basic fold term from edge-landmark (o2,o4)
    double d24 = T2.dist[o6]; // equals T4.dist[o2] on this tree
    double C1       = d24 - d2_3 - d6_7;
    double C2_basic = d24 - d2_7 - d6_3;

    cout << "Basic fold terms from (o2,o4):\n";
    cout << "  C1        = " << C1 << "\n";
    cout << "  C2_basic  = " << C2_basic << "\n\n";

    double LB = max({C1, C2_basic, C2, C4});
    cout << "Final combined lower bound LB(o7,o3) = " << LB << "\n";
}

// A simple HCA-style demonstration (still uses base LCA for the splitting node).
double run_hca_demo(const SPT &base,
                    const SPT &T5,
                    const SPT &T7,
                    int o5, int o7)
{
    int h = base.lca(o5, o7);  // “HCA” node in base tree

    std::cout << "HCA-style split point in base tree:\n";
    std::cout << "  HCA(o5,o7) = " << h << "\n\n";

    // Distances from endpoints to h in their own rooted SPTs
    double d5h = T5.dist[h];
    double d7h = T7.dist[h];
    double L   = d5h + d7h;  // total path length between o5 and o7

    std::cout << "Distances:\n";
    std::cout << "  dist(o5,h) = " << d5h << "  (in SPT rooted at o5)\n";
    std::cout << "  dist(o7,h) = " << d7h << "  (in SPT rooted at o7)\n";
    std::cout << "  path length L(o5,o7) via h = " << L << "\n\n";

    // Longest edges on each half-path via jump pointers
    double M5 = T5.maxEdgeOnPathToAncestor(o5, h);
    double M7 = T7.maxEdgeOnPathToAncestor(o7, h);
    double M  = std::max(M5, M7);

    std::cout << "Longest edges on partial paths:\n";
    std::cout << "  max-edge(o5→h) = " << M5 << "\n";
    std::cout << "  max-edge(o7→h) = " << M7 << "\n";
    std::cout << "  M = max(M5, M7)  = " << M << "\n\n";

    // Remaining-path terms at h in each rooted SPT
    double rem5_h = T5.downMax[h];
    double rem7_h = T7.downMax[h];

    std::cout << "Remaining-path lengths at h:\n";
    std::cout << "  remaining_5(h) = " << rem5_h << "\n";
    std::cout << "  remaining_7(h) = " << rem7_h << "\n\n";

    double LB_edge = 2.0 * M - L;
    double LB_rem  = std::fabs(rem5_h - rem7_h);
    double LB      = std::max(0.0, std::max(LB_edge, LB_rem));

    std::cout << "Lower bound using HCA method (option 4):\n";
    std::cout << "  LB_edge  = 2M - L     = " << LB_edge << "\n";
    std::cout << "  LB_rem   = |rem5-rem7| = " << LB_rem  << "\n";
    std::cout << "  LB(o5,o7) = " << LB << "\n";

    return LB;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Build the toy graph from the paper figures.
    // We use vertices 1..7; 0 is unused.
    int n = 8;
    Graph g(n);
    g.addEdge(2,1,0.2);
    g.addEdge(2,4,0.8);
    g.addEdge(2,6,0.1);
    g.addEdge(6,3,0.5);
    g.addEdge(3,5,0.3);
    g.addEdge(6,7,0.2);

    cout << "Select refinement method:\n"
         << " 1 = Jump Pointers (4.4.1)\n"
         << " 2 = Ladders (4.4.2)\n"
         << " 3 = Cartesian (max-edge via Cartesian tree)\n"
         << " 4 = HCA across two SPTs\n"
         << "Choice: ";
    int choice;
    if (!(cin >> choice)) return 0;

    RefineMethod method;
    switch (choice) {
        case 1: method = RefineMethod::JumpPointers; break;
        case 2: method = RefineMethod::Ladders;      break;
        case 3: method = RefineMethod::Cartesian;    break;
        case 4: method = RefineMethod::HCA;          break;
        default: return 0;
    }

    cout << "\nUsing method: ";
    if      (method == RefineMethod::JumpPointers) cout << "Jump Pointers\n";
    else if (method == RefineMethod::Ladders)      cout << "Ladders\n";
    else if (method == RefineMethod::Cartesian)    cout << "Cartesian (max-edge via Cartesian tree)\n";
    else                                           cout << "HCA (demo)\n";

    // ---- Base SPT from o2 via Dijkstra ----
    int o2 = 2, o6 = 6, o5 = 5, o7 = 7;
    SPT base(n, o2);
    base.buildFromGraph(g);   // this is the unique SP-tree of the example

    // ---- Re-root for o2, o6, o5, o7 using the *same* underlying tree ----
    SPT T2 = base;                 // rooted at o2
    SPT T6 = SPT::reRootFrom(base, o6);
    SPT T5 = SPT::reRootFrom(base, o5);
    SPT T7 = SPT::reRootFrom(base, o7);

    if (method == RefineMethod::HCA) {
        std::cout << "\n=== SPT rooted at o5 ===\n";
        T5.printEuler();
        std::cout << "=== SPT rooted at o7 ===\n";
        T7.printEuler();

        run_hca_demo(base, T5, T7, o5, o7);
    } else {
        run_basic_three(T2, T6, method);
    }

    return 0;
}





