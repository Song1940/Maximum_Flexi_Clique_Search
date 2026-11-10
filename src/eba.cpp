// ============================================================
// EBA: exact maximum flexi-clique (Section V of the paper).
//
// A search node B = (S, Cr, Cun, D) keeps a connected partial set S, the
// reachable candidates Cr (adjacent to S), the unreachable candidates Cun, and
// the excluded nodes D; its scope is H_B = S ∪ Cr ∪ Cun.
//
// CP-Branching (Section V-A). For an ordering <v_1..v_k> of (a subset of) Cr,
// child B_i adds v_1..v_{i-1} to S and excludes v_i; the include-all child
// B_{k+1} adds all of them. Candidates of Cun that become adjacent to S are
// promoted to Cr of the child. The children partition the connected node sets
// represented by B, so every connected node set is enumerated exactly once
// (Proposition 1). It is realised as a single sweep that pushes one node into
// S per step, so the incremental degree bookkeeping applies:
//
//   for v_i in order:                 // S == S_B ∪ {v_1..v_{i-1}}
//       exclude v_i  -> recurse       // child B_i
//       restore v_i; push v_i -> S    // prefix for the next sibling
//   recurse                           // include-all child B_{k+1}
//   pop the accumulated prefix
//
// At the root (S = ∅, Case (ii)) the children start from single nodes, taken
// in the reverse of a degeneracy ordering of the refined scope (rootChildren):
// the scope is peeled by repeatedly removing a node of minimum adjusted
// degree (with the Rule 1 cascade) until it is no larger than the incumbent,
// and the peeled nodes are then expanded in reverse order, each inside the
// scope it was peeled from, so that the densest part of the graph is searched
// first with the fewest candidates. Earlier children are excluded from later
// ones (Case (ii) of the paper).
//
// Pruning (Section V-B), applied in refine() at every search node and, in its
// global form, at the root:
//   Rule 1 (adjusted degree): candidates with d_a < theta_F' are removed with
//     cascading; a node of S below theta_F' prunes the search node.
//   Rule 2 (distance): candidates too far from S (or disconnected from it),
//     using the lower bounds dist_S maintained by one BFS per node added to S.
//   Rule 3 (theta_F'-core components): nodes whose component of the
//     theta_F'-core of the scope has at most |F'| nodes, applied at the root
//     together with the global form of Rule 1 and re-applied whenever the
//     incumbent F' improves, before the next root child is expanded.
//
// Margin pivot (Section V-C). Each node v has a non-neighbour budget
// dbar_max(v) = s_max(v) - 1 - floor(s_max(v)^tau), s_max(v) the size bound of
// Lemma 1 capped by the scope. The node p with the smallest margin
// dbar_max(p) - dbar_S(p) is selected, and only p (if p ∈ Cr) and its first
// margin(p) non-neighbours in Cr are branched on; the other candidates stay
// in the include-all child. A node of S is saturated when its size bound or
// its budget is exceeded; the sibling loop stops at the first saturated
// child. The same test removes candidates (candidate saturation): a candidate
// that is a non-neighbour of a node of S with margin <= 0, or whose own
// non-neighbours in S exceed its budget, lies in no flexi-clique containing S.
//
// Variants of the experiments (Section VI and the technical report), selected
// at compile time; the default build defines none of these macros (see
// CMakeLists.txt for the binaries built from them):
//   NO_RULE1, NO_RULE2, NO_RULE3   disable Rule 1, 2, or 3;
//   NO_INCUMBENT       start from an empty incumbent instead of the FPA solution;
//   NO_PIVOT           branch on all reachable candidates in ascending order of
//                      adjusted degree (no pivot ordering);
//   NO_CANDIDATE_SATURATION
//                      do not remove candidates by the saturation test (the
//                      sibling loop still stops at a saturated partial set);
//                      EBA-Pivot of the paper = NO_PIVOT NO_CANDIDATE_SATURATION;
//   NO_DEGENERACY_ORDER
//                      expand the root children in ascending order of adjusted
//                      degree instead of the reverse degeneracy ordering
//                      (EBA-Degeneracy of the paper);
//   FLAT_SYMSE         flat Sym-SE branching of FastQC: every candidate is
//                      branchable, S may be disconnected, and connectivity is
//                      checked when S is tested as a solution; same rules and
//                      incumbent as EBA, no candidate saturation;
//   BINARY_BRANCHING   include/exclude branching on one reachable candidate,
//                      the first of the pivot ordering; no candidate saturation;
//   ROOT_STATS         print the root pruning statistics and stop before branching;
//   TRACE_WIDTH        print the histogram of the number of children explored
//                      per branching node.
// ============================================================
#include "eba.h"
#include "fpa.h"
#include "peel_split.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iterator>
#include <queue>
#include <utility>
#include <vector>

#if defined(FLAT_SYMSE) && defined(BINARY_BRANCHING)
#error "FLAT_SYMSE and BINARY_BRANCHING are mutually exclusive"
#endif
#if defined(ROOT_STATS) && defined(NO_RULE3)
#error "ROOT_STATS reports the Rule 3 upper bound and requires Rule 3"
#endif

#ifdef TRACE_WIDTH
// width_hist[w]: number of branching nodes that explored w children (the
// exclusion children plus the include-all child, if explored); w >= 255 is
// counted in width_hist[255]. Leaves, pruned nodes and the root are not counted.
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

    // Rule 2 (distance pruning): dist_S[u] = max_{v∈S} sp_{G[H]}(u,v),
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
#ifndef NO_CANDIDATE_SATURATION
    std::vector<int> sat_cnt, sat_tmp;   // candidate saturation scratch (sat_cnt all 0 between uses)
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

    // Lemma 3 (degree-diameter): a connected graph with min degree k and diameter L
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

    // refine: per-node pruning (Rules 1 and 2) and the saturation tests of
    // Section V-C —
    //   Rule 1 (degree): drop candidates with d_a < θ and cascade; prune the node
    //                    if any S vertex has d_a < θ;
    //   Rule 2 (distance): drop candidates too far from / disconnected from S;
    //   candidate saturation: drop candidates that would saturate a node of S or
    //                    are saturated themselves;
    //   size bound: prune the node if |S| exceeds its size bound;
    //   saturation: prune the node if some S vertex already has more
    //                    non-neighbours in S than any feasible flexi-clique tolerates.
    // Mutates Cr/Cun (drops pruned candidates); records removals for undo.
    // Returns false to prune the whole node.
    bool refine(const std::vector<int>& S,
                std::vector<int>& Cr, std::vector<int>& Cun,
                std::vector<std::pair<int,int>>& deg_undo,
                std::vector<int>& removed) {
#if !defined(NO_RULE1) || !defined(NO_RULE2) || !defined(NO_CANDIDATE_SATURATION)
        int th = theta();
#endif

        // Initial removals: Rule 1 (low adjusted degree) + Rule 2 (distance).
        // Rule 2 (only when S≠∅) discards a candidate u that would force an
        // infeasible subgraph given its distance du = dist_S[u] from S in G[H],
        // via the two conditions of the distance rule:
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
            if (!rm && Snon) {                                     // Rule 2 (distance)
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

#ifndef NO_CANDIDATE_SATURATION
        // Candidate saturation (Section V-C): for any flexi-clique X with
        // S ⊆ X ⊆ scope and any v ∈ X, |X| <= s_max(v), so dbar_X(v) <= dbar_max(v)
        // (s - 1 - floor(s^τ) is non-decreasing in s), and dbar_X(v) >= dbar_S(v) for
        // v ∈ S. A candidate u therefore lies in no such X if
        //  (a) u ∉ N(z) for some z ∈ S with margin(z) = dbar_max(z) - dbar_S(z) <= 0
        //      (adding u gives z one non-neighbour more than it can tolerate), or
        //  (b) its own non-neighbours in S already exceed its budget:
        //      |S| - |N(u) ∩ S| > dbar_max(u).
        // Tests (a) and (b) never read the incumbent (they hold for every such X); the
        // Rule 1 cascade they trigger uses θ of the incumbent, as everywhere. Removals
        // shrink the scope and the degrees, which only lowers dbar_max, so the tests
        // are repeated until nothing changes.
        if (!S.empty()) {
            const int s_size = (int)S.size();
            if (sat_cnt.size() != (size_t)G.n) sat_cnt.assign(G.n, 0);
            for (;;) {
                const int scope = s_size + (int)Cr.size() + (int)Cun.size();
                sat_tmp.clear();
                for (int v : S)
                    if (dbarMax(v, scope) - nonDegInS(v, s_size) <= 0) sat_tmp.push_back(v);
                for (int z : sat_tmp)
                    for (int u : G.adj[z]) if (in_scope[u] && !in_S[u]) ++sat_cnt[u];
                std::vector<int> q2;
                auto test = [&](int u) {
                    if (!in_scope[u]) return;
                    if (sat_cnt[u] < (int)sat_tmp.size()                    // (a)
                        || s_size - deg_in_S[u] > dbarMax(u, scope))        // (b)
                        q2.push_back(u);
                };
                for (int u : Cr)  test(u);
                for (int u : Cun) test(u);
                for (int z : sat_tmp)                                       // reset counts
                    for (int u : G.adj[z]) sat_cnt[u] = 0;
                if (q2.empty()) break;
                for (size_t qi = 0; qi < q2.size(); ++qi) {
                    int v = q2[qi];
                    if (!in_scope[v]) continue;
                    in_scope[v] = false;
                    removed.push_back(v);
                    for (int u : G.adj[v]) {
                        if (!in_scope[u]) continue;
                        deg_undo.push_back({u, local_adj_deg[u]});
                        --local_adj_deg[u];
#ifndef NO_RULE1
                        if (local_adj_deg[u] < th && !in_S[u]) q2.push_back(u);   // Rule 1 cascade
#endif
                    }
                }
                Cr.erase(std::remove_if(Cr.begin(), Cr.end(),
                                        [&](int u){ return !in_scope[u]; }), Cr.end());
                Cun.erase(std::remove_if(Cun.begin(), Cun.end(),
                                         [&](int u){ return !in_scope[u]; }), Cun.end());
            }
        }
#endif

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
    // exclusion children. If no candidate is selected (p ∈ S without a usable
    // non-neighbour), all of Cr is branched on in ascending order of adjusted
    // degree.
    std::vector<int> orderCr(const std::vector<int>& S,
                             const std::vector<int>& Cr,
                             const std::vector<int>& Cun, int scope_size) {
        (void)Cun;
#ifdef NO_PIVOT
        (void)S; (void)scope_size;
        std::vector<int> order = Cr;
        sortByDegree(order);
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
            sortByDegree(order);
        }
        return order;
#endif
    }

    // Ascending adjusted degree (the root order without the degeneracy ordering,
    // and the candidate order without the pivot).
    void sortByDegree(std::vector<int>& vs) {
        std::sort(vs.begin(), vs.end(),
                  [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
    }

#if !defined(FLAT_SYMSE) && !defined(NO_DEGENERACY_ORDER)
    // Root of the search (S = ∅, refined scope `uni`, |uni| > |F'|), Case (ii) with
    // the reverse degeneracy ordering.
    //
    // Forward (no search): peel a minimum-adjusted-degree node w_i, remove it and
    // cascade Rule 1 (adjusted degree < θ), until |scope| <= |F'|. Let V_i be the
    // scope just before w_i is removed: V_1 = uni ⊇ V_2 ⊇ ... ⊇ V_m.
    // Backward (i = m..1): restore V_i, search S = {w_i} inside V_i with the
    // unchanged CP recursion.
    //
    // Exactness. Let X be a connected flexi-clique with |X| > |F'| (F' = the
    // incumbent at peel time). Every node of X has >= floor(|X|^τ) >= θ
    // neighbours inside X, so as long as X ⊆ scope no node of X falls below θ
    // and the Rule-1 cascade never removes one. The peel stops only when
    // |scope| <= |F'| < |X|, so some node of X is peeled as a w; the first one,
    // w_i, uniquely owns X: X ⊆ V_i and X contains w_i, so the search of seed w_i
    // (every connected set containing w_i inside V_i) covers X, while for j > i
    // w_i ∉ V_j. A larger incumbent found later may prune X without visiting it,
    // which is safe: only sets that can still beat the incumbent must survive,
    // and Rule 1 / Rule 3 with the current incumbent are sound inside V_i.
    void rootChildren(std::vector<int>& S, const std::vector<int>& uni) {
        const int entry_best = best_size;
        const int th = theta();
        (void)th;                                     // unused under NO_RULE1
        int maxd = 0;
        for (int v : uni) maxd = std::max(maxd, local_adj_deg[v]);
        std::vector<std::vector<int>> bucket(maxd + 1);
        for (int v : uni) bucket[local_adj_deg[v]].push_back(v);

        std::vector<int> chain;                       // w_1..w_m
        std::vector<size_t> rem_mark, deg_mark;       // log sizes before step i
        std::vector<int> premoved;                    // removed nodes (all steps)
        std::vector<std::pair<int,int>> pdeg_undo;    // (node, old degree)
        std::vector<int> q;
        int scope_sz = (int)uni.size();
        int cur = 0;

        // Forward peel. best_size is constant here (no search yet).
        while (scope_sz > best_size) {
            int w = -1;
            while (cur <= maxd) {
                if (bucket[cur].empty()) { ++cur; continue; }
                int v = bucket[cur].back();
                bucket[cur].pop_back();
                if (in_scope[v] && local_adj_deg[v] == cur) { w = v; break; }
            }
            if (w < 0) break;                         // unreachable: scope_sz > 0
            chain.push_back(w);
            rem_mark.push_back(premoved.size());
            deg_mark.push_back(pdeg_undo.size());
            q.assign(1, w);
            for (size_t qi = 0; qi < q.size(); ++qi) {
                int v = q[qi];
                if (!in_scope[v]) continue;
                in_scope[v] = false;
                in_D[v] = true;
                premoved.push_back(v);
                --scope_sz;
                for (int u : G.adj[v]) {
                    if (!in_scope[u]) continue;
                    pdeg_undo.push_back({u, local_adj_deg[u]});
                    int d = --local_adj_deg[u];
#ifndef NO_RULE1
                    if (d < th) { q.push_back(u); continue; }   // Rule 1 cascade
#endif
                    bucket[d].push_back(u);
                    if (d < cur) cur = d;
                }
            }
        }

        // The peel queue is no longer needed; free it before the (long) search.
        std::vector<std::vector<int>>().swap(bucket);
        std::vector<int>().swap(q);

        // Backward search, densest scope first.
        // Rule 3 is re-applied once the incumbent exceeds the one the root scope was
        // refined with. Rule 3 is O(n + m) per call. Backward scopes grow, so its
        // removals cannot be kept across steps; re-apply it only when the incumbent
        // grew since the last application or the scope at least doubled since then.
        // Skipping a call only skips optional pruning, so exactness is unaffected.
        const int peel_best = entry_best;
        int cp_best = peel_best, cp_scope = 0;
        (void)peel_best; (void)cp_best; (void)cp_scope;    // unused under NO_RULE3
        std::vector<int> childCr, childCun;
        std::vector<std::pair<int,int>> seed_dist_undo;
        for (int i = (int)chain.size() - 1; i >= 0; --i) {
            // Restore V_i: degrees (LIFO), then membership of the step's removals.
            for (size_t j = pdeg_undo.size(); j > deg_mark[i]; --j)
                local_adj_deg[pdeg_undo[j - 1].first] = pdeg_undo[j - 1].second;
            pdeg_undo.resize(deg_mark[i]);
            for (size_t j = rem_mark[i]; j < premoved.size(); ++j) {
                in_scope[premoved[j]] = true;
                in_D[premoved[j]] = false;
            }
            scope_sz += (int)(premoved.size() - rem_mark[i]);
            premoved.resize(rem_mark[i]);

            const int w = chain[i];
            if (scope_sz <= best_size) continue;      // V_i cannot beat the incumbent
#ifndef NO_RULE1
            if (local_adj_deg[w] < theta()) continue; // Rule 1 on the seed
#endif
            // The incumbent exceeds the root's: Rule 3 on V_i (undone below).
            std::vector<std::pair<int,int>> c_deg_undo;
            std::vector<int> c_removed;
#ifndef NO_RULE3
            if (best_size > peel_best &&
                (best_size > cp_best || (long long)scope_sz >= 2LL * cp_scope)) {
                cp_best = best_size;
                cp_scope = scope_sz;
                std::vector<bool> act(in_scope.begin(), in_scope.end());
                componentPrune(G, act, tau, best_size);
                for (int v : uni)
                    if (in_scope[v] && !act[v]) scopeRemove(v, c_deg_undo, c_removed);
            }
#endif
            if (in_scope[w] && scope_sz - (int)c_removed.size() > best_size) {
                pushToS(w, S);
                updateDistOnPush(w, uni, seed_dist_undo);
                buildChildSets(uni, childCr, childCun);
                recurse(S, childCr, childCun);
                restoreDist(seed_dist_undo, 0);
                popFromS(w, S);
            }
            undoScope(c_removed, c_deg_undo);
        }
        // All steps restored: scope == uni again, degrees as on entry.
    }
#endif

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

        // Root (S = ∅, Case (ii)): each child starts from one node of Cun, and
        // earlier children are excluded from later ones, so every connected
        // subgraph is enumerated exactly once. (Flat Sym-SE has no root case:
        // it branches on all candidates from S = ∅.)
#ifndef FLAT_SYMSE
        if (S.empty()) {
#ifndef NO_DEGENERACY_ORDER
            if (force_seed < 0) {
                rootChildren(S, Cun);               // restores the scope it changes
                undoScope(ref_removed, ref_deg_undo);
                return;
            }
#endif
            // Seed mode (one child, the seed), or the root children in ascending
            // order of adjusted degree (NO_DEGENERACY_ORDER).
            std::vector<int> seeds = Cun;       // snapshot of refined Cun
            if (force_seed >= 0) {              // seed mode: every enumerated set contains the seed
                seeds.clear();
                if (in_scope[force_seed]) seeds.push_back(force_seed);
            }
            sortByDegree(seeds);
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

        // Sibling-loop break (Section V-C): across siblings S_1⊆S_2⊆…⊆S_{k+1}, so
        // the partial set only grows. If the current S (= S_i) is already
        // size-infeasible or saturated, then every later sibling and the
        // include-all child contains it and is infeasible too — stop. Evaluated
        // in the node's (un-excluded) scope, which is stable under prefix
        // accumulation and conservative w.r.t. each child's smaller scope, so the
        // break is sound.
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
        // Include-all child B_{k+1}: S == nodeS ∪ order, no exclusion (skipped if the
        // sibling break already proved S_{k+1} infeasible).
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
// Rule 3 (theta-core components) in its general form. A flexi-clique X of size
// s containing v lies in the k-core of G[active], k = floor(s^tau), and, being
// connected, in the component of that core containing v; also s <= sizeBound(k).
// Hence
//   ub(v) = max_{k <= core(v)} min(sizeBound(k), |component of the k-core containing v|)
// bounds every flexi-clique through v, and v is removed if ub(v) <= best. For
// k < theta_F' the size bound is at most best, so ub(v) > best exactly when the
// component of the theta_F'-core containing v has more than best nodes, which is
// Rule 3 as stated in the paper; the general form also yields the upper bound
// max_v ub(v) reported by ROOT_STATS. Component sizes of all k-cores are
// obtained in one pass by adding nodes in decreasing core order to a union-find
// structure. Returns #removed nodes.
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
    // Initial incumbent from FPA (Section IV).
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
    // Rule 3 (theta-core components), alternated with the global Rule 1 until
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
