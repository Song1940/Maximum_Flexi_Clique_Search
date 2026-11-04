#include "graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

void Graph::finalize() {
    m = 0;
    for (int v = 0; v < n; ++v) {
        auto& a = adj[v];
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
        auto it = std::find(a.begin(), a.end(), v);
        if (it != a.end()) a.erase(it);
    }
    for (int v = 0; v < n; ++v) m += (long long)adj[v].size();
    m /= 2;
}

bool Graph::isConnected(const std::vector<bool>& active) const {
    int start = -1, cnt = 0;
    for (int v = 0; v < n; ++v)
        if (active[v]) { ++cnt; if (start == -1) start = v; }
    if (cnt <= 1) return true;

    int visited = 0;
    std::vector<bool> seen(n, false);
    std::queue<int> q;
    q.push(start); seen[start] = true; ++visited;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int w : adj[u])
            if (active[w] && !seen[w]) { seen[w] = true; ++visited; q.push(w); }
    }
    return visited == cnt;
}

Graph Graph::loadFromFile(const std::string& path, std::unordered_map<int,int>* id_map_out) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot open: " + path);

    std::unordered_map<int,int> id_map;
    std::vector<std::pair<int,int>> edges;
    int next_id = 0;

    auto getID = [&](int raw) -> int {
        auto [it, inserted] = id_map.emplace(raw, next_id);
        if (inserted) ++next_id;
        return it->second;
    };

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int u, v;
        if (!(ss >> u >> v) || u == v) continue;
        edges.push_back({getID(u), getID(v)});
    }

    Graph G(next_id);
    for (auto [u, v] : edges) G.addEdge(u, v);
    G.finalize();
    if (id_map_out) *id_map_out = std::move(id_map);
    std::cout << "Loaded graph: n=" << G.n << " m=" << G.m << "\n";
    return G;
}

bool isFlexiClique(const Graph& G, const std::vector<int>& nodes, double tau) {
    if (nodes.empty()) return false;
    int k = (int)nodes.size();
    // No singleton special case: a single node has degree 0 < floor(1^tau) = 1,
    // so it is never a flexi-clique (consistent with the solvers' feasibility test).
    int threshold = floorPow(k, tau);

    std::unordered_set<int> nodeSet(nodes.begin(), nodes.end());
    for (int v : nodes) {
        int cnt = 0;
        for (int u : G.adj[v])
            if (nodeSet.count(u)) ++cnt;
        if (cnt < threshold) return false;
    }

    std::vector<bool> active(G.n, false);
    for (int v : nodes) active[v] = true;
    return G.isConnected(active);
}
