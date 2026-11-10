#pragma once
#include "graph.h"
#include "algo_result.h"

// FPA: Flexi-Prune Algorithm (Section IV of the paper).
// Global peeling, then an improvement loop on the lambda_H-core of the current
// solution H (augmentation and slack-aware peeling).
// Returns an empty node list if the graph has no flexi-clique (e.g. no edges).
AlgoResult runFPA(const Graph& G, double tau);
