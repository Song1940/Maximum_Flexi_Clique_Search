#include "eba.h"
#include "fpa.h"
#include "peel_split.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <queue>
#include <utility>
#include <vector>
#include <cstdio>
#include <iterator>
#ifdef ORDER_RANDOM
#include <cstdlib>
#include <random>
#endif



#if defined(FLAT_SYMSE) && defined(BINARY_BRANCHING)
#error "FLAT_SYMSE and BINARY_BRANCHING are mutually exclusive"
#endif
#if (defined(ORDER_DESC) || defined(ORDER_RANDOM)) && !defined(NO_PIVOT)
#error "ORDER_DESC and ORDER_RANDOM require NO_PIVOT"
#endif
#if defined(ROOT_STATS) && defined(NO_RULE3)
#error "ROOT_STATS reports the Rule 3 upper bound and requires Rule 3"
#endif

#ifdef TRACE_WIDTH

static long long width_hist[256];
static void printWidthHist() {
    std::fprintf(stderr, "WidthHist:");
    for (int w = 0; w < 256; ++w)
        if (width_hist[w]) std::fprintf(stderr, " %d:%lld", w, width_hist[w]);
    std::fprintf(stderr, "\n");
}
#endif

namespace {

#ifndef NO_RULE3
int componentPrune(const Graph& G, std::vector<bool>& active, double tau, int best,
                   int* max_ub = nullptr);
#endif

struct CumState {
    const Graph& G;
    double tau;

    std::vector<int> best;
    int best_size = 0;
    int force_seed = -1;              // >= 0: search only connected sets containing this node

    // Partition / scope state (indexed by vertex).
    std::vector<bool> in_S;       // v ∈ S
    std::vector<bool> in_D;       // v ∈ D (excluded)
    std::vector<bool> in_scope;   // v ∈ R = S ∪ Cr ∪ Cun
    std::vector<int>  local_adj_deg;  // d_a(v) = |adj[v] ∩ R|
    std::vector<int>  deg_in_S;       // |adj[v] ∩ S|

    // Rule 2 (diameter pruning): dist_S[u] = max_{v∈S} sp_{G[H]}(u,v),
    // the longest shortest-path from u to any S vertex within the induced
    // subgraph G[H]. 0 when S=∅; INT_MAX when u is disconnected from S in G[H].
    // Maintained incrementally (one BFS per pushToS, max-combined, undone on pop).
    // NOTE: distances are NOT re-raised when a later scopeRemove shrinks G[H], so a
    // stored finite dist_S[u] may underestimate the true sp in the reduced G[H].
    // This only ever weakens Rule 2 (an underestimate prunes less), never
    // over-prunes, so soundness holds; some disconnection prunes are just missed.
    std::vector<int>  dist_S;
    // BFS scratch reused via touched-list resets (bfs_d kept all-INT_MAX between calls).
    std::vector<int>  bfs_d, bfs_touched, bfs_queue;

    // Scope size |H_B| of the node last refined (cap of the size bounds).
    int node_cap = INT_MAX;

    long long branches = 0;

#ifdef FLAT_SYMSE
    // Scratch for the connectivity test in feasibleS().
    mutable std::vector<char> conn_seen;
    mutable std::vector<int>  conn_stack, conn_touched;
#endif
#ifdef ORDER_RANDOM
    std::mt19937 rng{ []() -> unsigned {
        const char* s = std::getenv("FLEXI_SEED");
        return s ? static_cast<unsigned>(std::strtoul(s, nullptr, 10)) : 0xC0FFEEu;
    }() };
#endif

    CumState(const Graph& G, double tau)
        : G(G), tau(tau),
          in_S(G.n, false), in_D(G.n, false), in_scope(G.n, false),
          local_adj_deg(G.n, 0), deg_in_S(G.n, 0),
          dist_S(G.n, 0), bfs_d(G.n, INT_MAX), size_bound_cache(G.n + 1, -1) {
#ifdef FLAT_SYMSE
        conn_seen.assign(G.n, 0);
#endif
    }

    int theta() const { return computeTheta(best_size, tau); }

    // Size bound (Lemma 1): the largest flexi-clique size a node of
    // adjusted degree d can belong to, i.e. the largest s with floor(s^τ) <= d
    // (= ceil((d+1)^{1/τ}) - 1). Saturates at INT_MAX; every caller either caps it
    // by scope_size or compares it against |S| (both <= n).
    // Memoised per degree (d <= n), since it is queried for every S vertex at
    // every search node.
    mutable std::vector<int> size_bound_cache;   // -1 = not yet computed
    int sizeBound(int d) const {
        if (d < 0) return 0;
        if (d >= (int)size_bound_cache.size()) return maxSizeForDegree(d, tau);
        int& c = size_bound_cache[d];
        if (c < 0) c = maxSizeForDegree(d, tau);
        return c;
    }

    // S is always connected (CP invariant), so only the degree test is needed.
    bool feasibleS(const std::vector<int>& S) const {
        int k = (int)S.size();
        if (k == 0) return false;
        int thr = floorPow(k, tau);
        for (int v : S)
            if (deg_in_S[v] < thr) return false;
#ifdef FLAT_SYMSE
        // Flat Sym-SE does not keep S connected: DFS over G[S] (conn_seen is
        // all-zero between calls).
        conn_stack.assign(1, S[0]);
        conn_touched.assign(1, S[0]);
        conn_seen[S[0]] = 1;
        while (!conn_stack.empty()) {
            int x = conn_stack.back(); conn_stack.pop_back();
            for (int y : G.adj[x])
                if (in_S[y] && !conn_seen[y]) {
                    conn_seen[y] = 1; conn_touched.push_back(y); conn_stack.push_back(y);
                }
        }
        const bool connected = ((int)conn_touched.size() == k);
        for (int x : conn_touched) conn_seen[x] = 0;
        if (!connected) return false;
#endif
        return true;
    }

    // Remove v from scope R; record undo (decremented neighbours + removed v).
    void scopeRemove(int v, std::vector<std::pair<int,int>>& deg_undo,
                     std::vector<int>& removed) {
        in_scope[v] = false;
        removed.push_back(v);
        for (int u : G.adj[v])
            if (in_scope[u]) { deg_undo.push_back({u, local_adj_deg[u]}); --local_adj_deg[u]; }
    }

    // Undo scopeRemove/cascade: restore degrees (reverse) then scope membership.
    void undoScope(std::vector<int>& removed,
                   std::vector<std::pair<int,int>>& deg_undo,
                   bool reset_in_D = false) {
        for (auto it = deg_undo.rbegin(); it != deg_undo.rend(); ++it)
            local_adj_deg[it->first] = it->second;
        for (int v : removed) { in_scope[v] = true; if (reset_in_D) in_D[v] = false; }
    }

    // Move v into S (still in scope); update incremental deg_in_S.
    // deg_in_S[u] = |adj[u] ∩ S| is maintained for EVERY vertex u (not just S
    // members), because buildChildSets uses it to tell Cr (adjacent to S) from
    // Cun. v's own count is already accumulated by earlier pushes of its
    // neighbours, so only its neighbours are bumped here.
    void pushToS(int v, std::vector<int>& S) {
        in_S[v] = true;
        S.push_back(v);
        for (int u : G.adj[v]) ++deg_in_S[u];
    }
    void popFromS(int v, std::vector<int>& S) {
        in_S[v] = false;
        S.pop_back();
        for (int u : G.adj[v]) --deg_in_S[u];
    }

    // Non-neighbours of v inside S: |S| (minus 1 if v∈S) − |adj[v] ∩ S|.
    int nonDegInS(int v, int s_size) const {
        return (s_size - (in_S[v] ? 1 : 0)) - deg_in_S[v];
    }

    int dbarMax(int v, int scope_size) const {
        int da = local_adj_deg[v];
        int smax = std::max(1, std::min(sizeBound(da), scope_size));
        return smax - 1 - floorPow(smax, tau);
    }
    int smaxOf(int v, int scope_size) const {
        int da = local_adj_deg[v];
        return std::max(1, std::min(sizeBound(da), scope_size));
    }

    // Lemma (degree-diameter): a connected graph with min degree k and diameter L
    // has at least n(k,L) nodes. Long-long version of graph.h::nFunc to avoid the
    // (L/3)*(k-2) term overflowing int for large H.
    static long long nFuncLL(int k, int L) {
        if (L <= 0) return 1;
        if (k <= 0) return (long long)L + 1;
        if (L <= 2 || k == 1) return (long long)k + L;
        return (long long)k + L + 1 + (long long)(L / 3) * (k - 2);
    }

    // BFS from w within G[H] (in_scope vertices); update dist_S[u] = max(dist_S[u],
    // sp(w,u)) for every u in `universe` (INT_MAX when w cannot reach u in G[H]).
    // Records each raised value in dist_undo for exact restoration on pop.
    void updateDistOnPush(int w, const std::vector<int>& universe,
                          std::vector<std::pair<int,int>>& dist_undo) {
        bfs_touched.clear();
        bfs_queue.clear();
        bfs_d[w] = 0;
        bfs_touched.push_back(w);
        bfs_queue.push_back(w);
        for (size_t qi = 0; qi < bfs_queue.size(); ++qi) {
            int x = bfs_queue[qi];
            int dx = bfs_d[x];
            for (int y : G.adj[x]) {
                if (in_scope[y] && bfs_d[y] == INT_MAX) {
                    bfs_d[y] = dx + 1;
                    bfs_touched.push_back(y);
                    bfs_queue.push_back(y);
                }
            }
        }
        for (int u : universe) {
            if (!in_scope[u] || in_S[u]) continue;
            int sp = bfs_d[u];                         // INT_MAX if unreachable in G[H]
            if (sp > dist_S[u]) { dist_undo.push_back({u, dist_S[u]}); dist_S[u] = sp; }
        }
        for (int x : bfs_touched) bfs_d[x] = INT_MAX;  // reset scratch
    }

    // Undo dist_S raises recorded at index >= `from` (LIFO).
    void restoreDist(std::vector<std::pair<int,int>>& dist_undo, int from) {
        for (int j = (int)dist_undo.size() - 1; j >= from; --j)
            dist_S[dist_undo[j].first] = dist_undo[j].second;
        dist_undo.resize(from);
    }

    // refine: per-node pruning (Rules 1 and 2) and the saturation test of
    // Section V-C —
    //   Rule 1 (degree): drop candidates with d_a < θ and cascade; prune the node
    //                    if any S vertex has d_a < θ;
    //   Rule 2 (diameter): drop candidates too far from / disconnected from S;
    //   size bound: prune the node if |S| exceeds its size bound;
    //   saturation: prune the node if some S vertex already has more
    //                    non-neighbours in S than any feasible flexi-clique tolerates.
    // Mutates Cr/Cun (drops pruned candidates); records removals for undo.
    // Returns false to prune the whole node.
    bool refine(const std::vector<int>& S,
                std::vector<int>& Cr, std::vector<int>& Cun,
                std::vector<std::pair<int,int>>& deg_undo,
                std::vector<int>& removed) {
#if !defined(NO_RULE1) || !defined(NO_RULE2)
        int th = theta();
#endif

        // Initial removals: Rule 1 (low adjusted degree) + Rule 2 (diameter).
        // Rule 2 (only when S≠∅) discards a candidate u that would force an
        // infeasible subgraph given its distance du = dist_S[u] from S in G[H],
        // via the two conditions of the diameter rule:
        //   Rule 2(i):  n(θ,du) > |H|              (too many nodes needed to reach u)
        //   Rule 2(ii): ⌊n(θ,du)^τ⌋ > min_{v∈S} d_a(v)  (would over-constrain some S vertex)
        // du = INT_MAX (u disconnected from S in G[H]) is always pruned — this is
        // what lets the rule delete unreachable candidates. Both conditions use the
        // node's incoming |H| and min-degree, evaluated before any cascade.
#ifndef NO_RULE2
        const bool Snon = !S.empty();
        const int  scope0 = (int)S.size() + (int)Cr.size() + (int)Cun.size();
        int minAdS = INT_MAX;
        if (Snon) for (int v : S) minAdS = std::min(minAdS, local_adj_deg[v]);
#endif

        std::vector<int> q;
        auto consider = [&](int u) {
            if (!in_scope[u]) return;
#ifdef NO_RULE1
            bool rm = false;
#else
            bool rm = local_adj_deg[u] < th;                       // Rule 1
#endif
#ifndef NO_RULE2
            if (!rm && Snon) {                                     // Rule 2 (diameter)
                int du = dist_S[u];
                if (du == INT_MAX) rm = true;                     // disconnected from S
                else {
                    long long nf = nFuncLL(th, du);
                    if (nf > scope0) rm = true;                                       // Rule 2(i)
                    // Rule 2(ii): floor(nf^τ) > minAdS, with the shared threshold function.
                    else if (floorPow(nf, tau) > minAdS)
                        rm = true;
                }
            }
#endif
            if (rm) q.push_back(u);
        };
        for (int u : Cr)  consider(u);
        for (int u : Cun) consider(u);

        // Cascade (Rule 1 propagation): removing a candidate may drop neighbours
        // below θ, which cascades further.
        for (size_t qi = 0; qi < q.size(); ++qi) {
            int v = q[qi];
            if (!in_scope[v]) continue;          // already removed
            in_scope[v] = false;
            removed.push_back(v);
            for (int u : G.adj[v]) {
                if (!in_scope[u]) continue;
                deg_undo.push_back({u, local_adj_deg[u]});
                --local_adj_deg[u];
#ifndef NO_RULE1
                if (local_adj_deg[u] < th && !in_S[u]) q.push_back(u);   // Rule 1 cascade
#endif
            }
        }
        Cr.erase(std::remove_if(Cr.begin(), Cr.end(),
                                [&](int u){ return !in_scope[u]; }), Cr.end());
        Cun.erase(std::remove_if(Cun.begin(), Cun.end(),
                                 [&](int u){ return !in_scope[u]; }), Cun.end());

        // Rule 1 on S.
#ifndef NO_RULE1
        for (int v : S)
            if (local_adj_deg[v] < th) return false;
#endif

        // Size bound (Lemma 1) from the minimum adjusted degree in S.
        // F_max = sizeBound(min_{v∈S} d_a(v)) is the largest size any flexi
        // containing S can reach. Prune if |S| already exceeds it, OR |S| == F_max
        // but S is not yet a valid flexi (it can no longer grow, so it never will).
        if (!S.empty()) {
            int md = INT_MAX;
            for (int v : S) md = std::min(md, local_adj_deg[v]);
            int fmax = sizeBound(md);
            if ((int)S.size() > fmax) return false;
            if ((int)S.size() == fmax && !feasibleS(S)) return false;
        }

        // Saturation test (Section V-C): a node is infeasible once some S vertex has more
        // non-neighbours inside S than any feasible flexi-clique can tolerate.
        int scope_size = (int)S.size() + (int)Cr.size() + (int)Cun.size();
        node_cap = scope_size;
        int s_size = (int)S.size();
        for (int v : S) {
            if (smaxOf(v, node_cap) < s_size) return false;
            if (nonDegInS(v, s_size) > dbarMax(v, node_cap)) return false;
        }
        return true;
    }

    // Classify the node's candidate universe into the child's (Cr, Cun) by
    // current state: in-S / out-of-scope are skipped; deg_in_S>0 ⇒ reachable.
    void buildChildSets(const std::vector<int>& cand_uni,
                        std::vector<int>& childCr, std::vector<int>& childCun) {
        childCr.clear();
        childCun.clear();
        for (int u : cand_uni) {
            if (in_S[u] || !in_scope[u]) continue;
#ifdef FLAT_SYMSE
            childCr.push_back(u);              // every candidate is branchable
            continue;
#endif
            if (deg_in_S[u] > 0) childCr.push_back(u);
            else                 childCun.push_back(u);
        }
    }

    // Candidate ordering for CP-Branching: the margin pivot (Section V-C).
    //
    // Pick the least-slack node p = argmin_{v∈S∪Cr} margin(v), where
    //   margin(v) = dbar_max(v) − dbar_S(v)   (non-neighbour slack until v breaks).
    // Branch ONLY on p (if p∈Cr) plus the first t = min(margin(p), |N̄(p)|)
    // non-neighbours of p in Cr (ascending margin); every other candidate —
    // notably p's neighbours — is NOT branched individually but left in Cr and
    // deferred to the include-all child and deeper recursion, which Proposition 1
    // permits. Front-loading p's non-neighbours drives dbar_S(p) to dbar_max(p),
    // so the saturation test stops the sibling loop after about margin(p)+1
    // exclusion children.
    std::vector<int> orderCr(const std::vector<int>& S,
                             const std::vector<int>& Cr,
                             const std::vector<int>& Cun, int scope_size) {
        (void)Cun;
#ifdef NO_PIVOT
        (void)S; (void)scope_size;
        std::vector<int> order = Cr;
        sortByOrder(order);
        return order;
#else
        auto margin = [&](int v) {
            return dbarMax(v, scope_size) - nonDegInS(v, (int)S.size());
        };
        int p = -1, pm = INT_MAX;
        for (int v : S)  { int m = margin(v); if (m < pm) { pm = m; p = v; } }
        for (int v : Cr) { int m = margin(v); if (m < pm) { pm = m; p = v; } }

        // N̄(p) = non-neighbours of p inside Cr, ascending margin.
        std::vector<std::pair<int,int>> non_nbr;
        for (int u : Cr)
            if (u != p && !G.hasEdge(p, u)) non_nbr.push_back({margin(u), u});
        std::sort(non_nbr.begin(), non_nbr.end());

        int t = std::min(std::max(pm, 0), (int)non_nbr.size());  // min(margin(p), |N̄(p)|)
        std::vector<int> order;
        order.reserve(t + 1);
        if (!in_S[p]) order.push_back(p);          // p ∈ Cr → branch on p itself first
        for (int i = 0; i < t; ++i) order.push_back(non_nbr[i].second);

        // Empty ordering (t==0 and p∈S): fall back to the full Cr so the node can
        // still make progress.
        if (order.empty()) {
            order = Cr;
            std::sort(order.begin(), order.end(),
                      [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
        }
        return order;
#endif
    }

    // Order of the root seeds (and, with NO_PIVOT, of the branched candidates):
    // ascending adjusted degree by default.
    void sortByOrder(std::vector<int>& vs) {
#if defined(ORDER_RANDOM)
        std::shuffle(vs.begin(), vs.end(), rng);
#elif defined(ORDER_DESC)
        std::sort(vs.begin(), vs.end(),
                  [&](int a, int b){ return local_adj_deg[a] > local_adj_deg[b]; });
#else
        std::sort(vs.begin(), vs.end(),
                  [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
#endif
    }

    void recurse(std::vector<int>& S, std::vector<int> Cr, std::vector<int> Cun) {
        ++branches;

        // Incumbent update (before refine).
        if ((int)S.size() > best_size && feasibleS(S)) {
            best = S;
            best_size = (int)best.size();
        }

        std::vector<std::pair<int,int>> ref_deg_undo;
        std::vector<int> ref_removed;
        if (!refine(S, Cr, Cun, ref_deg_undo, ref_removed)) {
            undoScope(ref_removed, ref_deg_undo);
            return;
        }
        // Size cap of this node (|H_B|); valid for every
        // child, whose scope is contained in H.
        const int cap = node_cap;

        // Sound size bound: the whole reachable universe can't beat the
        // incumbent, so neither can any descendant. (max achievable ≤ |R|.)
        if ((int)S.size() + (int)Cr.size() + (int)Cun.size() <= best_size) {
            undoScope(ref_removed, ref_deg_undo);
            return;
        }

        // Root (S = ∅, Case (ii)): each Cun vertex starts one child, and
        // earlier seeds are excluded so every connected subgraph is enumerated
        // exactly once by its lowest-ranked vertex. (Flat Sym-SE has no root
        // case: it branches on all candidates from S = ∅.)
#ifndef FLAT_SYMSE
        if (S.empty()) {
            std::vector<int> seeds = Cun;       // snapshot of refined Cun
            if (force_seed >= 0) {              // seed mode: every enumerated set contains the seed
                seeds.clear();
                if (in_scope[force_seed]) seeds.push_back(force_seed);
            }
            // Seeds in ascending order of adjusted degree (see sortByOrder). Every
            // connected subgraph is enumerated exactly once (by its first-processed
            // seed), so the order affects efficiency, not correctness.
            sortByOrder(seeds);
            std::vector<std::pair<int,int>> seed_deg_undo;
            std::vector<int> seed_removed;
            std::vector<int> childCr, childCun;
            std::vector<std::pair<int,int>> seed_dist_undo;
#ifndef NO_RULE3
            int pruned_at = best_size;          // incumbent of the last component prune
#endif
            for (int w : seeds) {
#ifndef NO_RULE3
                // The incumbent grew in an earlier seed's subtree: re-apply the
                // component-aware pruning to the remaining scope. Removals share
                // the seed undo log, so they are restored in the right order.
                if (best_size > pruned_at) {
                    std::vector<bool> act(in_scope.begin(), in_scope.end());
                    componentPrune(G, act, tau, best_size);
                    for (int v : seeds)
                        if (in_scope[v] && !act[v]) scopeRemove(v, seed_deg_undo, seed_removed);
                    pruned_at = best_size;
                }
#endif
                if (!in_scope[w]) continue;
                pushToS(w, S);
                updateDistOnPush(w, Cun, seed_dist_undo);
                buildChildSets(Cun, childCr, childCun);
                recurse(S, childCr, childCun);
                restoreDist(seed_dist_undo, 0);
                popFromS(w, S);
                in_D[w] = true;
                scopeRemove(w, seed_deg_undo, seed_removed);
            }
            undoScope(seed_removed, seed_deg_undo, /*reset_in_D=*/true);
            undoScope(ref_removed, ref_deg_undo);
            return;
        }
#endif

        if (Cr.empty()) {                       // leaf
            undoScope(ref_removed, ref_deg_undo);
            return;
        }

        std::vector<int> order = orderCr(S, Cr, Cun, cap);

        std::vector<int> cand_uni;
        cand_uni.reserve(Cr.size() + Cun.size());
        cand_uni.insert(cand_uni.end(), Cr.begin(), Cr.end());
        cand_uni.insert(cand_uni.end(), Cun.begin(), Cun.end());

        // Sibling-loop break (Section V-C): across siblings S_0⊆S_1⊆…⊆S_k, so
        // the partial set only grows. If the current S (= S_i) is already
        // size-infeasible (Rule 2) or saturated, then every later sibling and the
        // include-all child B_k contains it and is infeasible too — stop. Evaluated
        // in the node's (un-excluded) scope, which is stable under prefix
        // accumulation and conservative w.r.t. each child's smaller scope, so the
        // break is sound. This is what bounds cumulative branching width: without
        // it the size bound never bites (each child excludes only one candidate).
        auto siblingInfeasible = [&]() -> bool {
            int s_size = (int)S.size();
            if (s_size == 0) return false;      // only reachable in flat Sym-SE
            int md = INT_MAX;
            for (int v : S) md = std::min(md, local_adj_deg[v]);
            int fmax = sizeBound(md);
            if (s_size > fmax) return true;
            if (s_size == fmax && !feasibleS(S)) return true;
            for (int v : S) {
                if (smaxOf(v, cap) < s_size) return true;
                if (nonDegInS(v, s_size) > dbarMax(v, cap)) return true;
            }
            return false;
        };

        std::vector<int> childCr, childCun;
        std::vector<std::pair<int,int>> dist_undo;   // dist_S raises across the prefix

#ifdef BINARY_BRANCHING
        // Exclude the first candidate of the ordering, then include it; both
        // children are refined by the recursive call.
        {
            int v = order[0];
            std::vector<std::pair<int,int>> ex_deg_undo;
            std::vector<int> ex_removed;
            in_D[v] = true;
            scopeRemove(v, ex_deg_undo, ex_removed);
            buildChildSets(cand_uni, childCr, childCun);
            recurse(S, childCr, childCun);
            undoScope(ex_removed, ex_deg_undo, /*reset_in_D=*/true);

            pushToS(v, S);
            updateDistOnPush(v, cand_uni, dist_undo);
            buildChildSets(cand_uni, childCr, childCun);
            recurse(S, childCr, childCun);
            restoreDist(dist_undo, 0);
            popFromS(v, S);
#ifdef TRACE_WIDTH
            ++width_hist[2];
#endif
        }
        undoScope(ref_removed, ref_deg_undo);
        return;
#endif
        std::vector<int> dist_marks;                 // dist_undo size before each push
        int pushed = 0;
        bool broke = false;
        for (size_t idx = 0; idx < order.size(); ++idx) {
            if (siblingInfeasible()) { broke = true; break; }
            int vi = order[idx];                // S == nodeS ∪ {order[0..idx-1]}
            // Child B_i: exclude v_i only.
            std::vector<std::pair<int,int>> ex_deg_undo;
            std::vector<int> ex_removed;
            in_D[vi] = true;
            scopeRemove(vi, ex_deg_undo, ex_removed);
            buildChildSets(cand_uni, childCr, childCun);
            recurse(S, childCr, childCun);
            undoScope(ex_removed, ex_deg_undo, /*reset_in_D=*/true);
            // Accumulate v_i into S for later siblings; extend dist_S by a BFS from v_i.
            pushToS(vi, S);
            dist_marks.push_back((int)dist_undo.size());
            updateDistOnPush(vi, cand_uni, dist_undo);
            ++pushed;
        }
        // Include-all child B_k: S == nodeS ∪ order, no exclusion (skipped if the
        // sibling break already proved S_k infeasible).
        if (!broke) {
            buildChildSets(cand_uni, childCr, childCun);
            recurse(S, childCr, childCun);
        }
#ifdef TRACE_WIDTH
        ++width_hist[std::min(pushed + (broke ? 0 : 1), 255)];
#endif

        for (int i = pushed - 1; i >= 0; --i) {
            restoreDist(dist_undo, dist_marks[i]);
            popFromS(order[i], S);
        }
        undoScope(ref_removed, ref_deg_undo);
    }
};

// Global degree prune (Rule 1 applied to the whole graph) on the active mask.
#ifndef NO_RULE1
void globalDegreePrune(const Graph& G, std::vector<bool>& active,
                       std::vector<int>& adeg, int th) {
    std::queue<int> q;
    for (int v = 0; v < G.n; ++v)
        if (active[v] && adeg[v] < th) q.push(v);
    while (!q.empty()) {
        int v = q.front(); q.pop();
        if (!active[v]) continue;
        active[v] = false;
        for (int u : G.adj[v])
            if (active[u] && --adeg[u] < th) q.push(u);
    }
}
#endif

#ifndef NO_RULE3
// Component-aware global pruning. A flexi-clique X of size s containing v lies in
// the k-core of G[active], k = floor(s^tau), and, being connected, in the
// component of that core containing v; also s <= sizeBound(k). Hence
//   ub(v) = max_{k <= core(v)} min(sizeBound(k), |component of the k-core containing v|)
// bounds every flexi-clique through v, and v is removed if ub(v) <= best.
// Component sizes of all k-cores are obtained in one pass by adding nodes in
// decreasing core order to a union-find structure. Returns #removed nodes.
int componentPrune(const Graph& G, std::vector<bool>& active, double tau, int best,
                   int* max_ub) {
    const int n = G.n;
    // Core numbers of G[active].
    std::vector<int> deg(n, 0), core(n, -1);
    int maxd = 0;
    for (int v = 0; v < n; ++v) if (active[v]) {
        for (int u : G.adj[v]) deg[v] += active[u];
        maxd = std::max(maxd, deg[v]);
    }
    {
        std::vector<std::vector<int>> bucket(maxd + 1);
        for (int v = 0; v < n; ++v) if (active[v]) bucket[deg[v]].push_back(v);
        std::vector<int> d = deg;
        for (int k = 0; k <= maxd; ++k)
            for (size_t i = 0; i < bucket[k].size(); ++i) {
                int v = bucket[k][i];
                if (core[v] >= 0 || d[v] != k) continue;
                core[v] = k;
                for (int u : G.adj[v]) if (active[u] && core[u] < 0 && d[u] > k) {
                    --d[u]; bucket[std::max(d[u], k)].push_back(u);
                }
            }
    }
    // Union-find over nodes added in decreasing core order; after adding all
    // nodes with core >= k, the sets are the components of the k-core. A merge
    // tree records each event that creates or enlarges a component at level k
    // with value min(sizeBound(k), size); a component unchanged at lower k cannot
    // improve this value (sizeBound is non-decreasing in k). ub(v) is the maximum
    // value on the path from v's leaf to its root, computed top-down at the end.
    std::vector<int> parent(n, -1), csize(n, 0), tnode(n, -1);  // DSU over vertices
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    std::vector<int> tval, tpar;                 // merge-tree nodes: value, parent
    auto newNode = [&](int value) {
        tval.push_back(value); tpar.push_back(-1); return (int)tval.size() - 1;
    };
    std::vector<std::vector<int>> byCore(maxd + 1);
    for (int v = 0; v < n; ++v) if (active[v]) byCore[core[v]].push_back(v);
    std::vector<int> leaf(n, -1);
    for (int k = maxd; k >= 1; --k) {
        const int sb = maxSizeForDegree(k, tau);
        for (int v : byCore[k]) {
            parent[v] = v; csize[v] = 1;
            leaf[v] = tnode[v] = newNode(std::min(sb, 1));
        }
        for (int v : byCore[k])
            for (int u : G.adj[v]) if (active[u] && parent[u] >= 0) {
                int a = find(v), b = find(u);
                if (a == b) continue;
                if (csize[a] < csize[b]) std::swap(a, b);
                parent[b] = a; csize[a] += csize[b];
                int t = newNode(std::min(sb, csize[a]));
                tpar[tnode[a]] = t; tpar[tnode[b]] = t;
                tnode[a] = t;
            }
    }
    // Top-down: a parent is created after its children, so visiting nodes in
    // decreasing index order sees every parent before its children.
    std::vector<int> pathmax(tval.size(), 0);
    for (int t = (int)tval.size() - 1; t >= 0; --t)
        pathmax[t] = std::max(tval[t], tpar[t] >= 0 ? pathmax[tpar[t]] : 0);
    std::vector<int> ub(n, 0);
    for (int v = 0; v < n; ++v) if (leaf[v] >= 0) ub[v] = pathmax[leaf[v]];
    int removed = 0;
    for (int v = 0; v < n; ++v)
        if (active[v] && ub[v] <= best) { active[v] = false; ++removed; }
    if (max_ub) {
        *max_ub = 0;
        for (int v = 0; v < n; ++v) if (active[v]) *max_ub = std::max(*max_ub, ub[v]);
    }
    return removed;
}
#endif

} // namespace

AlgoResult runEBA(const Graph& G, double tau, int seed,
                            const std::vector<int>* init) {
    auto t0 = std::chrono::high_resolution_clock::now();
#ifdef TRACE_WIDTH
    std::fill(std::begin(width_hist), std::end(width_hist), 0LL);
#endif

    if (G.n == 0) {                       // empty graph: empty clique
#ifdef ROOT_STATS
        std::fprintf(stderr, "RootStats: F=0 V1=0 V13=0 UB=0\n");
#endif
#ifdef TRACE_WIDTH
        printWidthHist();
#endif
        AlgoResult res;
        res.time_sec = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t0).count();
        return res;
    }

    CumState st(G, tau);
    st.force_seed = seed;
    std::vector<int> init_sol;
    double fpa_sec = -1.0;
#ifndef NO_INCUMBENT
    // Initial incumbent from FPA (Algorithm 1).
    {
        AlgoResult fr = runFPA(G, tau);
        init_sol = std::move(fr.nodes);
        fpa_sec = fr.time_sec;               // reported so that FPA <= EBA holds per run
    }
#endif
    // Seed mode: the FPA solution is a valid incumbent only if it contains the seed.
    if (seed >= 0 && std::find(init_sol.begin(), init_sol.end(), seed) == init_sol.end())
        init_sol.clear();
    if (!init_sol.empty() && isFlexiClique(G, init_sol, tau)) {
        st.best = init_sol;
        st.best_size = (int)init_sol.size();
    }
    // Caller-supplied incumbent (--init): verified, and kept only if it beats the
    // FPA solution and, in seed mode, contains the seed.
    if (init && !init->empty()) {
        std::vector<int> cand(*init);
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
        const bool has_seed_node = seed < 0 || std::binary_search(cand.begin(), cand.end(), seed);
        const bool in_range = cand.front() >= 0 && cand.back() < G.n;
        if (has_seed_node && in_range && (int)cand.size() > st.best_size
            && isFlexiClique(G, cand, tau)) {
            st.best = cand;
            st.best_size = (int)cand.size();
        }
    }

    // Global pruning (Rule 1 applied to the whole graph) to define V'.
    std::vector<bool> active(G.n, true);
    std::vector<int> adeg(G.n, 0);
    for (int v = 0; v < G.n; ++v) adeg[v] = (int)G.adj[v].size();
#ifndef NO_RULE1
    globalDegreePrune(G, active, adeg, st.theta());
#endif
#ifdef ROOT_STATS
    int rs_V1 = 0;                       // |V'| after the global Rule 1
    for (int v = 0; v < G.n; ++v) rs_V1 += active[v];
    int rs_ub = 0;                       // Rule 3 bound of the last call (on the final V')
#endif
#ifndef NO_RULE3
    // Rule 3 (upper bound ub(v)), alternated with the global Rule 1 until
    // neither removes a node.
#ifdef ROOT_STATS
    while (componentPrune(G, active, tau, st.best_size, &rs_ub) > 0) {
#else
    while (componentPrune(G, active, tau, st.best_size) > 0) {
#endif
        for (int v = 0; v < G.n; ++v) {
            if (!active[v]) continue;
            int d = 0;
            for (int u : G.adj[v]) if (active[u]) ++d;
            adeg[v] = d;
        }
#ifndef NO_RULE1
        globalDegreePrune(G, active, adeg, st.theta());
#endif
    }
#endif

    std::vector<int> V_prime;
    for (int v = 0; v < G.n; ++v) if (active[v]) V_prime.push_back(v);
    // Seed mode: if the global pruning removed the seed, no flexi-clique containing it
    // beats the incumbent (which, if any, contains the seed).
    if (seed >= 0 && (seed >= G.n || !active[seed])) V_prime.clear();

#ifdef ROOT_STATS
    // |F'|, |V'| after Rule 1 (V1) and after Rules 1 and 3 (V13), and the upper
    // bound max(|F'|, max ub(v) over V'); the search is skipped.
    std::fprintf(stderr, "RootStats: F=%d V1=%d V13=%d UB=%d\n", st.best_size, rs_V1,
                 (int)V_prime.size(), std::max(st.best_size, rs_ub));
    V_prime.clear();
#endif

    if (V_prime.empty()) {
        auto t1 = std::chrono::high_resolution_clock::now();
        AlgoResult res;
        res.nodes = st.best;   // only validated solutions enter st.best
        res.time_sec = std::chrono::duration<double>(t1 - t0).count();
        res.fpa_sec = fpa_sec;
#ifdef TRACE_WIDTH
        printWidthHist();
#endif
        return res;
    }

    for (int v : V_prime) {
        st.in_scope[v] = true;
        int d = 0;
        for (int u : G.adj[v]) if (active[u]) ++d;
        st.local_adj_deg[v] = d;
    }

#ifdef FLAT_SYMSE
    std::vector<int> S, Cr = V_prime, Cun;   // every node is branchable from the root
#else
    std::vector<int> S, Cr, Cun = V_prime;
#endif
    st.recurse(S, Cr, Cun);

    auto t1 = std::chrono::high_resolution_clock::now();
    AlgoResult res;
    res.nodes = st.best;   // only validated solutions enter st.best
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    res.fpa_sec = fpa_sec;
    res.branches = st.branches;
#ifdef TRACE_WIDTH
    printWidthHist();
#endif
    return res;
}
