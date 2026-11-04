#pragma once
#include <vector>

// Result of FPA or EBA.
struct AlgoResult {
    std::vector<int> nodes;      // solution (internal node indices)
    double time_sec = 0.0;       // running time, graph loading excluded
    long long branches = 0;      // EBA: number of search nodes
    double fpa_sec = -1.0;       // EBA: time of its FPA phase within this run (-1: none)
};
