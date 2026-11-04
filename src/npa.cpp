#include "npa.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <queue>
#include <utility>
#include <vector>

namespace {

// Bucket-based k-core decomposition (O(n + m)).
std::vector<int> coreNumbers(const Graph& G) {
    std::vector<int> core(G.n, 0);
    if (G.n == 0) return core;
    std::vector<int> deg(G.n);
    for (int v = 0; v < G.n; ++v) deg[v] = (int)G.adj[v].size();

    int max_d = *std::max_element(deg.begin(), deg.end());

    std::vector<int> bin(max_d + 1, 0);
    for (int v = 0; v < G.n; ++v) bin[deg[v]]++;

    int start = 0;
    std::vector<int> pos(max_d + 1);
    for (int d = 0; d <= max_d; ++d) { pos[d] = start; start += bin[d]; }

    std::vector<int> vert(G.n), pos_v(G.n);
    for (int v = 0; v < G.n; ++v) {
        pos_v[v]         = pos[deg[v]];
        vert[pos_v[v]]   = v;
        pos[deg[v]]++;
    }
    for (int d = max_d; d > 0; --d) pos[d] = pos[d - 1];
    pos[0] = 0;

    for (int i = 0; i < G.n; ++i) {
        int v = vert[i];
        core[v] = deg[v];
        for (int u : G.adj[v]) {
            if (deg[u] > deg[v]) {     // u not processed yet (processed nodes have deg <= deg[v])
                int du = deg[u], pu = pos_v[u], pw = pos[du], w = vert[pw];
                if (u != w) {
                    pos_v[u] = pw; vert[pu] = w;
                    pos_v[w] = pu; vert[pw] = u;
                }
                pos[du]++;
                deg[u]--;
            }
        }
    }
    return core;
}

// Mask of the largest connected component among the active nodes.
std::vector<bool> largestCC(const Graph& G, const std::vector<bool>& active) {
    std::vector<bool> in_lcc(G.n, false), visited(G.n, false);
    int best_size = 0;
    std::vector<int> best_comp;

    for (int s = 0; s < G.n; ++s) {
        if (!active[s] || visited[s]) continue;
        std::vector<int> comp;
        std::queue<int> q;
        q.push(s); visited[s] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int w : G.adj[u])
                if (active[w] && !visited[w]) { visited[w] = true; q.push(w); }
        }
        if ((int)comp.size() > best_size) { best_size = (int)comp.size(); best_comp = comp; }
    }

    for (int v : best_comp) in_lcc[v] = true;
    return in_lcc;
}

} // namespace

// NPA: peeling-based heuristic baseline for flexi-clique search, as in the original
// implementation (Kim et al., CIKM 2024; _archive/NPA.py):
//   1. Seed region: the smallest k >= 2 whose k-core has a largest connected component
//      H_k with floor(|H_k|^tau) <= k (k = k_star + 1 if there is none); the seed is the
//      largest connected component H of the (k-1)-core.
//   2. Repeatedly remove the lowest-degree non-articulation node of H until every node
//      meets floor(|H|^tau). H stays connected, so the remainder is a flexi-clique.
// The peeling works on a local copy of G[H0] for the seed H0, so each step costs
// O(|H0| + |E(H0)|) instead of O(n). Ties are broken by the smallest node id.
AlgoResult runNPA(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();
    AlgoResult res;
    auto finish = [&]() {
        res.time_sec = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t0).count();
        return res;
    };
    if (G.n == 0) return finish();

    // Step 1: seed region.
    std::vector<int> core = coreNumbers(G);
    int k_star = *std::max_element(core.begin(), core.end());
    int k = 2;
    for (; k <= k_star; ++k) {
        std::vector<bool> mask(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= k) mask[v] = true;
        std::vector<bool> lcc = largestCC(G, mask);
        int sz = (int)std::count(lcc.begin(), lcc.end(), true);
        if (floorPow(sz, tau) <= k) break;
    }
    std::vector<bool> mask(G.n, false);
    for (int v = 0; v < G.n; ++v)
        if (core[v] >= k - 1) mask[v] = true;
    std::vector<bool> inH = largestCC(G, mask);

    // Local copy of G[H]; local ids follow node ids, so ties resolve to the smallest id.
    std::vector<int> verts;
    for (int v = 0; v < G.n; ++v)
        if (inH[v]) verts.push_back(v);
    const int h = (int)verts.size();
    if (h == 0) return finish();
    std::vector<int> lid(G.n, -1);
    for (int i = 0; i < h; ++i) lid[verts[i]] = i;
    std::vector<std::vector<int>> ladj(h);
    std::vector<int> deg(h);
    for (int i = 0; i < h; ++i) {
        for (int u : G.adj[verts[i]])
            if (lid[u] >= 0) ladj[i].push_back(lid[u]);
        deg[i] = (int)ladj[i].size();
    }

    // Step 2: peeling.
    std::vector<char> alive(h, 1), ap(h);
    std::vector<int> disc(h), low(h), parent(h);
    std::vector<std::pair<int,int>> stk;      // (node, next index into ladj[node])
    stk.reserve(h);
    int cnt = h;
    bool feasible = false;
    while (cnt > 0) {
        int min_deg = INT_MAX;
        for (int i = 0; i < h; ++i)
            if (alive[i]) min_deg = std::min(min_deg, deg[i]);
        if (min_deg >= floorPow(cnt, tau)) { feasible = true; break; }

        // Articulation points of the alive part (iterative Tarjan).
        for (int i = 0; i < h; ++i)
            if (alive[i]) { disc[i] = -1; ap[i] = 0; }
        int timer = 0;
        for (int s = 0; s < h; ++s) {
            if (!alive[s] || disc[s] != -1) continue;
            disc[s] = low[s] = timer++;
            parent[s] = -1;
            int root_children = 0;
            stk.emplace_back(s, 0);
            while (!stk.empty()) {
                int u = stk.back().first;
                int& i = stk.back().second;   // safe: reserve(h) => no reallocation
                bool descended = false;
                while (i < (int)ladj[u].size()) {
                    int v = ladj[u][i++];
                    if (!alive[v]) continue;
                    if (disc[v] == -1) {            // tree edge
                        disc[v] = low[v] = timer++;
                        parent[v] = u;
                        if (u == s) ++root_children;
                        stk.emplace_back(v, 0);
                        descended = true;
                        break;
                    } else if (v != parent[u]) {    // back edge
                        low[u] = std::min(low[u], disc[v]);
                    }
                }
                if (!descended) {
                    int u_low = low[u];
                    stk.pop_back();
                    if (!stk.empty()) {
                        int p = stk.back().first;
                        low[p] = std::min(low[p], u_low);
                        if (parent[p] != -1 && u_low >= disc[p]) ap[p] = 1;
                    }
                }
            }
            if (root_children > 1) ap[s] = 1;
        }

        // Remove the lowest-degree non-articulation node.
        int best = -1;
        for (int i = 0; i < h; ++i)
            if (alive[i] && !ap[i] && (best < 0 || deg[i] < deg[best])) best = i;
        if (best < 0) break;                       // unreachable for a connected H
        alive[best] = 0;
        --cnt;
        for (int j : ladj[best])
            if (alive[j]) --deg[j];
    }

    if (feasible)
        for (int i = 0; i < h; ++i)
            if (alive[i]) res.nodes.push_back(verts[i]);
    return finish();
}
