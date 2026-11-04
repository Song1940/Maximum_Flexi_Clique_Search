#pragma once
#include "graph.h"
#include "algo_result.h"

// EBA: exact branch-and-bound algorithm for the maximum flexi-clique
// (Section V of the paper): CP-Branching, Rules 1-3, and the margin pivot,
// started from the FPA incumbent.
//
// seed >= 0: maximum flexi-clique that contains node `seed` (internal index);
//            seed = -1: maximum flexi-clique of G.
// init (optional): a flexi-clique supplied by the caller, used as the initial
//            incumbent when it is valid, larger than the FPA solution, and (in
//            seed mode) contains the seed; otherwise it is ignored.
AlgoResult runEBA(const Graph& G, double tau, int seed = -1,
                  const std::vector<int>* init = nullptr);
