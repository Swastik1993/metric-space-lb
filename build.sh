#!/bin/bash -l
# Build the optimized edge_landmarks_sampling binary on the NJIT cluster.
#
# Usage:
#   cd /project/senjutib/sb2785/metric_space/hd_rmq
#   bash build.sh
#
# Output: ./edge_landmarks_sampling

set -Eeuo pipefail

CXX="${CXX:-g++}"

# C++23 is required (uses std::ranges, structured bindings in lambdas, etc).
# -O3 -DNDEBUG are the production flags.
# -fopenmp is REQUIRED to get the parallel HCA + parallel Stochastic SPT
# builds; without it the code falls back to serial paths.
#
# -march=x86-64-v3 sets a portable feature level requiring AVX2 + BMI2 + FMA.
# Every Intel CPU from Haswell (2013+) and every AMD Zen2+ (2019+) supports
# this. Critically, it avoids the build-on-login / run-on-compute mismatch
# that happens with -march=native when the login node has newer instructions
# than the compute nodes (= SIGILL on first AVX-512 instruction).
# -mtune=native still lets the compiler tune scheduling for the build host.
#
# If even x86-64-v3 produces SIGILL on the compute node (= cluster is very
# old), fall back to -march=x86-64-v2 (SSE4.2 only).

CXXFLAGS="-std=c++23 -O3 -DNDEBUG -march=x86-64-v3 -mtune=native -fopenmp"
LDFLAGS="-fopenmp"

SOURCES=(
  Dijkstra.cpp
  ShortestPathTree.cpp
  EdgeLandMark.cpp
  ELM_SPTree.cpp
  spt.cpp
  LowerBound.cpp
  GraphUtils.cpp
  TriSearch.cpp
  fixed_main.cpp
)

echo "[build] using compiler: $($CXX --version | head -1)"
echo "[build] CXXFLAGS: $CXXFLAGS"
echo "[build] LDFLAGS:  $LDFLAGS"

$CXX $CXXFLAGS "${SOURCES[@]}" $LDFLAGS -o edge_landmarks_sampling

echo "[build] OK"
ls -la edge_landmarks_sampling
echo
echo "Quick sanity check (--help-like dump of positional usage):"
./edge_landmarks_sampling 2>&1 | head -5 || true
echo
echo "Done."



