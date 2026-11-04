#include "peel_split.h"
#include <algorithm>
#include <climits>
#include <queue>
#include <utility>

std::vector<int> peelAndSplit(const Graph& G, double tau,
                              const std::vector<bool>* active_mask) {
    std::vector<int> best;
    int global_best = 0;

    // Initial connected components, restricted to active_mask if given.
    std::vector<bool> visited(G.n, false);
    if (active_mask) {
        for (int v = 0; v < G.n; ++v)
            if (!(*active_mask)[v]) visited[v] = true;
    }

    std::vector<std::vector<int>> init_components;
    for (int s = 0; s < G.n; ++s) {
        if (visited[s]) continue;
        std::vector<int> comp;
        std::queue<int> q;
        q.push(s); visited[s] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int w : G.adj[u])
                if (!visited[w]) { visited[w] = true; q.push(w); }
        }
        init_components.push_back(std::move(comp));
    }

    std::sort(init_components.begin(), init_components.end(),
              [](const auto& a, const auto& b){ return a.size() > b.size(); });

    // Scratch arrays sized G.n, reused across components.
    std::vector<int>  cur_deg(G.n, 0);
    std::vector<char> in_comp(G.n, 0);
    std::vector<char> cc_vis(G.n, 0);

    for (std::vector<int>& comp : init_components) {
        int comp_sz = (int)comp.size();
        if (comp_sz <= global_best) continue;

        for (int v : comp) { in_comp[v] = 1; cur_deg[v] = 0; }
        for (int v : comp)
            for (int u : G.adj[v])
                if (in_comp[u]) ++cur_deg[v];

        // Any flexi-clique inside comp has size at most the size bound of its
        // maximum degree, so a component that cannot beat the best is skipped.
        int comp_max_deg = 0;
        for (int v : comp) comp_max_deg = std::max(comp_max_deg, cur_deg[v]);
        if (maxSizeForDegree(comp_max_deg, tau) <= global_best) {
            for (int v : comp) in_comp[v] = 0;
            continue;
        }

        // Lazy min-heap on the current degree; an entry whose key exceeds the
        // node's current degree is stale (degrees only decrease while peeling).
        using Entry = std::pair<int,int>;
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;
        for (int v : comp) pq.emplace(cur_deg[v], v);

        int remaining = comp_sz;
        while (remaining > global_best) {
            int min_node = -1;
            while (!pq.empty()) {
                auto [d, v] = pq.top();
                if (!in_comp[v] || d != cur_deg[v]) { pq.pop(); continue; }
                min_node = v;
                break;
            }
            if (min_node == -1) break;

            if (cur_deg[min_node] >= floorPow(remaining, tau)) {
                // The remaining set R meets the threshold. Every component of
                // G[R] keeps its degrees and has a threshold no larger than R's,
                // so each component is itself a flexi-clique; keep the largest.
                std::vector<int> rest;
                rest.reserve(remaining);
                for (int v : comp) if (in_comp[v]) rest.push_back(v);
                std::vector<int> largest;
                for (int s0 : rest) {
                    if (cc_vis[s0]) continue;
                    std::vector<int> sub;
                    std::queue<int> bq;
                    bq.push(s0); cc_vis[s0] = 1;
                    while (!bq.empty()) {
                        int u = bq.front(); bq.pop();
                        sub.push_back(u);
                        for (int w : G.adj[u])
                            if (in_comp[w] && !cc_vis[w]) { cc_vis[w] = 1; bq.push(w); }
                    }
                    if (sub.size() > largest.size()) largest = std::move(sub);
                }
                for (int v : rest) cc_vis[v] = 0;
                if ((int)largest.size() > global_best) {
                    global_best = (int)largest.size();
                    best = std::move(largest);
                }
                break;
            }

            // Peel min_node.
            pq.pop();
            for (int u : G.adj[min_node]) {
                if (in_comp[u]) { --cur_deg[u]; pq.emplace(cur_deg[u], u); }
            }
            in_comp[min_node] = 0;
            --remaining;
        }

        for (int v : comp) in_comp[v] = 0;
    }

    return best;
}

int augmentFlexi(const Graph& G, double tau, std::vector<int>& F) {
    if (F.size() < 2) return 0;
    // Scratch arrays are reused across calls (one graph per process) and only the
    // touched entries are reset, so a call costs O(volume of F's neighbourhood)
    // rather than O(n).
    static std::vector<char> inF;
    static std::vector<int>  degF;           // |N(v) ∩ F| for every touched node v
    if ((int)inF.size() != G.n) { inF.assign(G.n, 0); degF.assign(G.n, 0); }
    std::vector<int> touched;                // nodes whose inF/degF changed
    for (int v : F) { inF[v] = 1; touched.push_back(v); }
    for (int v : F) for (int u : G.adj[v]) { if (degF[u]++ == 0) touched.push_back(u); }
    // Lazy max-heap of frontier candidates keyed by degF; an entry is stale when
    // the node joined F or its degF changed since it was pushed.
    std::priority_queue<std::pair<int,int>> frontier;
    for (int v : F) for (int u : G.adj[v]) if (!inF[u]) frontier.push({degF[u], u});
    int added = 0;
    for (;;) {
        const int k = (int)F.size();
        const int thr_next = floorPow(k + 1, tau);
        // Nodes of F that would fall below the new threshold unless adjacent to v.
        // floor((k+1)^tau) - floor(k^tau) <= 1, so each needs exactly one more edge.
        std::vector<int> tight;
        bool fixable = true;
        for (int u : F) if (degF[u] < thr_next) {
            if (degF[u] + 1 < thr_next) { fixable = false; break; }
            tight.push_back(u);
        }
        if (!fixable) break;
        int best = -1;
        if (!tight.empty()) {
            // v must be adjacent to every tight node: scan the smallest adjacency.
            int pivot = tight[0];
            for (int u : tight) if (G.adj[u].size() < G.adj[pivot].size()) pivot = u;
            int best_deg = -1;
            for (int v : G.adj[pivot]) {
                if (inF[v] || degF[v] < thr_next || degF[v] <= best_deg) continue;
                bool ok = true;
                for (int u : tight) if (!G.hasEdge(u, v)) { ok = false; break; }
                if (ok) { best = v; best_deg = degF[v]; }
            }
        } else {
            while (!frontier.empty()) {
                auto [d, v] = frontier.top();
                if (inF[v] || d != degF[v]) { frontier.pop(); continue; }
                if (d >= thr_next) best = v;
                break;
            }
        }
        if (best < 0) break;
        F.push_back(best); inF[best] = 1; ++added;
        touched.push_back(best);
        for (int u : G.adj[best]) {
            if (degF[u]++ == 0) touched.push_back(u);
            if (!inF[u]) frontier.push({degF[u], u});
        }
    }
    for (int v : touched) { inF[v] = 0; degF[v] = 0; }
    return added;
}

std::vector<int> slackPeel(const Graph& G, double tau, int best,
                           const std::vector<bool>* active_mask) {
    std::vector<int> result;
    // Core numbers of the active subgraph (bucket peeling).
    std::vector<char> act(G.n, 1);
    if (active_mask) for (int v = 0; v < G.n; ++v) act[v] = (*active_mask)[v];
    std::vector<int> core(G.n, 0), deg(G.n, 0);
    int maxd = 0;
    for (int v = 0; v < G.n; ++v) if (act[v]) {
        for (int u : G.adj[v]) deg[v] += act[u];
        maxd = std::max(maxd, deg[v]);
    }
    {
        std::vector<std::vector<int>> bucket(maxd + 1);
        for (int v = 0; v < G.n; ++v) if (act[v]) bucket[deg[v]].push_back(v);
        std::vector<char> done(G.n, 0);
        std::vector<int> d = deg;
        for (int k = 0; k <= maxd; ++k) {
            for (size_t i = 0; i < bucket[k].size(); ++i) {
                int v = bucket[k][i];
                if (done[v] || d[v] != k) continue;
                done[v] = 1; core[v] = k;
                for (int u : G.adj[v]) if (act[u] && !done[u] && d[u] > k) {
                    --d[u]; bucket[std::max(d[u], k)].push_back(u);
                }
            }
        }
    }
    int max_core = 0, active_count = 0;
    for (int v = 0; v < G.n; ++v) if (act[v]) { max_core = std::max(max_core, core[v]); ++active_count; }
    if (best >= active_count) return result;

    std::vector<char> inR(G.n, 0), seen(G.n, 0);
    std::vector<int>  rdeg(G.n, 0), tightcnt(G.n, 0);
    for (int k = max_core; k >= 1; --k) {
        const int s_hi = maxSizeForDegree(k, tau);
        // A saturated bound (INT_MAX) means no graph size reaches threshold k, so
        // every component size is admissible; avoid the signed overflow of +1.
        const int s_lo_m = maxSizeForDegree(k - 1, tau);
        const int s_lo = s_lo_m == INT_MAX ? 0 : s_lo_m + 1;
        const int need = std::max(s_lo, (int)std::max<size_t>(best, result.size()) + 1);
        if (s_hi < need) break;                 // smaller k only gives smaller sizes
        // Components of the k-core.
        std::vector<int> knodes;
        for (int v = 0; v < G.n; ++v) if (act[v] && core[v] >= k) knodes.push_back(v);
        if ((int)knodes.size() < need) continue;
        for (int s0 : knodes) {
            if (seen[s0]) continue;
            std::vector<int> comp{s0}; seen[s0] = 1;
            for (size_t i = 0; i < comp.size(); ++i)
                for (int u : G.adj[comp[i]])
                    if (act[u] && core[u] >= k && !seen[u]) { seen[u] = 1; comp.push_back(u); }
            const int cur_need = std::max(s_lo, (int)std::max<size_t>(best, result.size()) + 1);
            if ((int)comp.size() < cur_need) continue;
            if ((int)comp.size() <= s_hi) {
                // A connected k-core component already within the size range: every
                // node has degree >= k >= floor(|comp|^tau), so it is a flexi-clique.
                if (floorPow((int)comp.size(), tau) <= k) result = comp;
                continue;
            }
            // Peel comp down to <= s_hi keeping degrees >= k.
            for (int v : comp) inR[v] = 1;
            for (int v : comp) { rdeg[v] = 0; for (int u : G.adj[v]) rdeg[v] += inR[u]; }
            for (int v : comp) { tightcnt[v] = 0; for (int u : G.adj[v]) if (inR[u] && rdeg[u] == k) ++tightcnt[v]; }
            using Key = std::pair<std::pair<int,int>, int>;   // ((tight, deg), v)
            std::priority_queue<Key, std::vector<Key>, std::greater<Key>> pq;
            for (int v : comp) pq.push({{tightcnt[v], rdeg[v]}, v});
            int size = (int)comp.size();
            std::vector<int> casc;
            auto removeNode = [&](int v) {
                const bool v_was_tight = (rdeg[v] == k);
                inR[v] = 0; --size;
                for (int u : G.adj[v]) if (inR[u]) {
                    if (v_was_tight) --tightcnt[u];         // u loses v as a tight neighbour
                    if (rdeg[u] == k)                       // u was tight; neighbours lose a tight neighbour
                        for (int w : G.adj[u]) if (inR[w]) { --tightcnt[w]; pq.push({{tightcnt[w], rdeg[w]}, w}); }
                    --rdeg[u];
                    if (rdeg[u] == k)                       // u became tight
                        for (int w : G.adj[u]) if (inR[w]) { ++tightcnt[w]; pq.push({{tightcnt[w], rdeg[w]}, w}); }
                    if (rdeg[u] < k) casc.push_back(u);
                    pq.push({{tightcnt[u], rdeg[u]}, u});
                }
            };
            while (size > s_hi && size >= cur_need) {
                int v = -1;
                while (!pq.empty()) {
                    auto [key, x] = pq.top(); pq.pop();
                    if (inR[x] && key.first == tightcnt[x] && key.second == rdeg[x]) { v = x; break; }
                }
                if (v < 0) break;
                removeNode(v);
                while (!casc.empty()) {                    // restore degrees >= k
                    int u = casc.back(); casc.pop_back();
                    if (inR[u] && rdeg[u] < k) removeNode(u);
                }
            }
            if (size >= cur_need && size <= s_hi) {
                // Keep the largest connected component of the remainder.
                std::vector<int> rest, best_comp;
                for (int v : comp) if (inR[v]) rest.push_back(v);
                for (int v : rest) seen[v] = 2;
                for (int r0 : rest) {
                    if (seen[r0] != 2) continue;
                    std::vector<int> c{r0}; seen[r0] = 3;
                    for (size_t i = 0; i < c.size(); ++i)
                        for (int u : G.adj[c[i]]) if (inR[u] && seen[u] == 2) { seen[u] = 3; c.push_back(u); }
                    if (c.size() > best_comp.size()) best_comp = std::move(c);
                }
                for (int v : rest) seen[v] = 1;
                // Every node has degree >= k; a component C keeps it, and is valid if floor(|C|^tau) <= k.
                if ((int)best_comp.size() >= cur_need && floorPow((int)best_comp.size(), tau) <= k)
                    result = std::move(best_comp);
            }
            for (int v : comp) inR[v] = 0;
        }
        for (int v : knodes) seen[v] = 0;
    }
    return result;
}
