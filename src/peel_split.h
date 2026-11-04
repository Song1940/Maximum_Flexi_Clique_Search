#pragma once
#include "graph.h"
#include <vector>

// Global peeling of FPA (Peel in Algorithm 1).
// Within every connected component (restricted to active_mask if given), repeatedly
// remove a minimum-degree node until the remaining set R satisfies the flexi-clique
// degree threshold delta(G[R]) >= floor(|R|^tau). Every connected component of
// G[R] is then a flexi-clique (its degrees are unchanged and its threshold is no
// larger), and the largest one is a candidate. Returns the largest candidate over
// all components (empty if none).
std::vector<int> peelAndSplit(const Graph& G, double tau,
                              const std::vector<bool>* active_mask = nullptr);

// Greedy augmentation (local search). Given a flexi-clique F, repeatedly adds a
// node v outside F such that F ∪ {v} is still a flexi-clique, choosing the
// candidate with the most neighbours in F, until no node can be added. F stays
// connected (v needs at least one neighbour in F). Returns the number of nodes
// added. Cost per added node: O(sum of degrees of F's neighbourhood).
int augmentFlexi(const Graph& G, double tau, std::vector<int>& F);

// Slack-aware thresholded peeling. For each threshold k (from high to low) whose
// size range [s_lo(k), s_hi(k)] = {s : floor(s^tau) = k} can beat |best|, takes the
// connected components of the k-core (within active_mask, if given) and removes
// nodes until the size is at most s_hi(k) while keeping every degree >= k,
// preferring nodes with the fewest tight neighbours (degree exactly k) so that the
// removal does not cascade. A remainder of size >= s_lo(k) is a flexi-clique.
// Returns the largest one found (or an empty vector if none beats |best|).
std::vector<int> slackPeel(const Graph& G, double tau, int best,
                           const std::vector<bool>* active_mask = nullptr);
