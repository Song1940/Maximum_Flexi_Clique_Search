#pragma once
#include "graph.h"
#include "algo_result.h"

// NPA: the existing peeling heuristic for the maximum flexi-clique (Baseline of
// FPA, Section VI): a k-core seed, then repeated removal of the minimum-degree
// node that is not an articulation point until every node meets the threshold.
AlgoResult runNPA(const Graph& G, double tau);
