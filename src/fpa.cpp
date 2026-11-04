#include "fpa.h"
#include "peel_split.h"
#include <chrono>
#include <vector>

// ============================================================
// runFPA — Flexi-Prune Algorithm (Algorithm 1 of the paper).
//   1. Global peeling: in every connected component, remove a minimum-degree
//      node until all remaining degrees reach floor(|R|^tau); the largest
//      component of the remainder is the initial solution H (Observation 1).
//   2. Improvement loop, repeated while H grows. A flexi-clique larger than H
//      lies in the lambda_H-core U, lambda_H = floor((|H|+1)^tau)
//      (Observation 2); H is (a) augmented node by node, and (b) replaced by
//      the augmented result of slack-aware peeling of U when that is larger.
// Returns the largest flexi-clique found.
//
// Variants of the FPA ablation (full version), selected at compile time:
// FPA_NO_GLOBAL_PEEL starts the loop from an empty solution, FPA_NO_AUGMENT
// skips the augmentation, and FPA_NO_SLACK_PEEL the slack-aware peeling.
// ============================================================
AlgoResult runFPA(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // 1. Global peeling.
#ifndef FPA_NO_GLOBAL_PEEL
    std::vector<int> best = peelAndSplit(G, tau);
#else
    std::vector<int> best;
#endif
    int best_size = (int)best.size();

    // 2. Improvement loop.
    std::vector<bool> core(G.n);
    std::vector<int>  deg(G.n), stack;
    for (;;) {
        int theta = computeTheta(best_size, tau);
        stack.clear();
        for (int v = 0; v < G.n; ++v) {
            deg[v]  = (int)G.adj[v].size();
            core[v] = deg[v] >= theta;
            if (!core[v]) stack.push_back(v);
        }
        while (!stack.empty()) {
            int v = stack.back(); stack.pop_back();
            for (int u : G.adj[v])
                if (core[u] && --deg[u] < theta) { core[u] = false; stack.push_back(u); }
        }
        bool improved = false;
#ifndef FPA_NO_AUGMENT
        if (augmentFlexi(G, tau, best) > 0) {
            best_size = (int)best.size();
            improved  = true;
        }
#endif
#ifndef FPA_NO_SLACK_PEEL
        {
            std::vector<int> sp = slackPeel(G, tau, best_size, &core);
            if ((int)sp.size() > best_size) {
                best = std::move(sp);
#ifndef FPA_NO_AUGMENT
                augmentFlexi(G, tau, best);
#endif
                best_size = (int)best.size();
                improved  = true;
            }
        }
#endif
        if (!improved) break;
    }

    AlgoResult res;
    res.nodes    = best;
    res.time_sec = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t0).count();
    return res;
}
