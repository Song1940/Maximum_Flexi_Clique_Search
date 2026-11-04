#pragma once
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <functional>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

struct Graph {
    int n = 0;
    long long m = 0;  // number of edges; 64-bit, as the degree sum of SK-2005 exceeds INT_MAX
    std::vector<std::vector<int>> adj;  // sorted adjacency lists

    Graph() = default;
    explicit Graph(int n) : n(n), adj(n) {}

    void addEdge(int u, int v) {
        assert(u >= 0 && u < n && v >= 0 && v < n && u != v);
        adj[u].push_back(v);
        adj[v].push_back(u);
        ++m;
    }

    // Remove duplicate edges and self-loops; sort adjacency lists.
    void finalize();

    bool hasEdge(int u, int v) const {
        return std::binary_search(adj[u].begin(), adj[u].end(), v);
    }

    int deg(int v) const { return (int)adj[v].size(); }

    // Returns true if all active nodes form a single connected component.
    bool isConnected(const std::vector<bool>& active) const;

    // Load an edge-list file (whitespace-separated); re-indexes nodes to [0, n).
    // id_map (optional) receives original id -> internal index.
    static Graph loadFromFile(const std::string& path,
                              std::unordered_map<int,int>* id_map = nullptr);
};

// Flexi threshold floor(k^tau), the single definition used everywhere (feasibility,
// theta, size bounds, Rule 2). tau arrives as a double (e.g. 0.6 is not exactly
// 3/5), so pow() can land just below an integer that the intended value hits
// exactly, e.g. 32^0.6 = 8 evaluates to 7.9999999999999991. Values within a
// relative 1e-12 below an integer are therefore rounded up to it; genuine
// non-integers in our ranges are far further from integers than that.
inline int floorPow(long long k, double tau) {
    if (k <= 0) return 0;
    const long double x = std::pow(static_cast<long double>(k), static_cast<long double>(tau));
    const long double r = std::floor(x * (1.0L + 1e-12L) + 1e-12L);
    return r >= static_cast<long double>(INT_MAX) ? INT_MAX : static_cast<int>(r);
}

// Minimum degree threshold needed to beat best_size: a flexi-clique larger than
// best_size has size >= best_size+1, so every member needs internal degree
// >= floorPow(best_size+1); nodes below it belong to no improving flexi-clique
// (Rule 1 and the global degree prune). Defined through floorPow so that theta,
// feasibility and the size bounds use one and the same threshold function.
inline int computeTheta(int best_size, double tau) {
    return floorPow(static_cast<long long>(best_size) + 1, tau);
}

// Size bound (Lemma 1 of the paper): the largest s >= 0 with floor(s^tau) <= d,
// i.e. the largest flexi-clique size whose threshold a node of degree d can meet.
// Mathematically s = ceil((d+1)^{1/tau}) - 1, which is one less than
// floor((d+1)^{1/tau}) whenever (d+1)^{1/tau} is an integer (always at tau = 0.5).
// It is computed against floorPow() itself, the threshold the solver uses for
// feasibility, so the bound is consistent with the feasibility test by
// construction: floorPow(s, tau) <= d and floorPow(s+1, tau) > d.
// Saturates at INT_MAX (small tau makes 1/tau large).
inline int maxSizeForDegree(long long d, double tau) {
    if (d < 0) return 0;
    long double est = std::floor(std::pow(static_cast<long double>(d) + 1.0L,
                                           1.0L / static_cast<long double>(tau)));
    if (est >= static_cast<long double>(INT_MAX) - 1.0L) return INT_MAX;
    long long s = static_cast<long long>(est);
    while (s > 0 && floorPow(s, tau) > d) --s;
    while (floorPow(s + 1, tau) <= d) ++s;
    return static_cast<int>(s);
}

// Lemma 3 of the paper (Yao and Chang, PVLDB 2021): minimum number of nodes of a
// connected graph with minimum degree k and diameter L.
inline int nFunc(int k, int L) {
    if (L <= 0) return 1;
    if (k <= 0) return L + 1;
    if (L <= 2 || k == 1) return k + L;
    return k + L + 1 + (L / 3) * (k - 2);
}

// Returns true if nodes form a valid tau-flexi-clique in G.
bool isFlexiClique(const Graph& G, const std::vector<int>& nodes, double tau);
