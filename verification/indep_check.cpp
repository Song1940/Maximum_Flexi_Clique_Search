// Independent checker for verification/certify_optimality.py. It shares no code with src/.
//
//   indep_check <edge list> <query file>
//
// Edge list: one edge per line; empty lines and lines starting with '#' are skipped; the first two
// integer tokens are the endpoints, further tokens are ignored; self-loops and duplicate edges are
// dropped; the graph is undirected. Other lines (no two integer tokens separated by spaces or tabs,
// or a token of more than 18 digits) are skipped and counted as malformed, and endpoints outside the
// range of a 32-bit int are counted: flexi reads such lines differently, so the caller requires both
// counts to be zero before it compares the two graphs.
// Query file, one query per line:
//   core <lambda> <f>   connected components of the lambda-core of G: their number, the three largest
//                       sizes, and how many have more than f nodes
//   sol <path>          the node set on the "nodes=" line of a solution file written by
//                       `flexi --out` (original node ids): size, malformed or unknown ids, duplicates,
//                       connectivity, minimum degree inside the set, and the size of the component of
//                       the (minimum degree)-core of G that contains the set
// Output: one line per query, with key=value fields, after two lines describing the graph.
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

typedef int64_t i64;
typedef uint32_t u32;

static bool dense_map = false;
static std::vector<i64> dense_id;                 // original id -> internal id (or -1)
static std::unordered_map<i64, i64> hash_id;      // used when the ids are too sparse for dense_id
static i64 n = 0;
static std::vector<i64> off;                      // CSR offsets, size n + 1
static std::vector<u32> adj;                      // CSR targets: sorted, unique, no self-loops
static std::vector<int> core;

static i64 lookup(i64 raw) {
    if (dense_map) return (raw >= 0 && raw < (i64)dense_id.size()) ? dense_id[raw] : -1;
    auto it = hash_id.find(raw);
    return it == hash_id.end() ? -1 : it->second;
}

static bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r'; }

// Decimal integer token at p (optional sign, at most 18 digits), ending at a space or at e.
static bool parse_int(const char*& p, const char* e, i64& out) {
    while (p < e && is_space(*p)) ++p;
    if (p >= e) return false;
    bool neg = false;
    if (*p == '-' || *p == '+') { neg = (*p == '-'); ++p; }
    const char* digits = p;
    i64 v = 0;
    while (p < e && *p >= '0' && *p <= '9') {
        if (p - digits == 18) return false;
        v = v * 10 + (*p - '0');
        ++p;
    }
    if (p == digits || (p < e && !is_space(*p))) return false;
    out = neg ? -v : v;
    return true;
}

static void load(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    std::vector<i64> raw;                          // endpoints, two per accepted line
    const size_t BUF = 1 << 26;
    std::vector<char> buf(BUF);
    std::string carry;
    i64 lines = 0, comments = 0, malformed = 0, out_of_int = 0, selfloops = 0;
    auto handle = [&](const char* s, const char* e) {
        ++lines;
        const char* p = s;
        while (p < e && is_space(*p)) ++p;
        if (p >= e || *p == '#') { ++comments; return; }
        i64 a, b;
        if (!parse_int(p, e, a) || !parse_int(p, e, b)) { ++malformed; return; }
        if (a < INT_MIN || a > INT_MAX || b < INT_MIN || b > INT_MAX) ++out_of_int;
        if (a == b) { ++selfloops; return; }
        raw.push_back(a); raw.push_back(b);
    };
    size_t got;
    while ((got = fread(buf.data(), 1, BUF, f)) > 0) {
        const char* s = buf.data();
        const char* e = s + got;
        const char* ls = s;
        for (const char* p = s; p < e; ++p) {
            if (*p != '\n') continue;
            if (!carry.empty()) { carry.append(ls, p); handle(carry.data(), carry.data() + carry.size()); carry.clear(); }
            else handle(ls, p);
            ls = p + 1;
        }
        carry.append(ls, e);
    }
    if (!carry.empty()) handle(carry.data(), carry.data() + carry.size());
    fclose(f);

    // internal ids in order of first appearance
    i64 maxid = -1, minid = 0;
    for (i64 x : raw) { maxid = std::max(maxid, x); minid = std::min(minid, x); }
    const i64 E = (i64)raw.size() / 2;
    std::vector<u32> ep(raw.size());
    if (minid >= 0 && maxid <= 4 * E + 1000000) {
        dense_map = true;
        dense_id.assign(maxid + 1, -1);
        for (size_t i = 0; i < raw.size(); ++i) {
            i64& d = dense_id[raw[i]];
            if (d < 0) d = n++;
            ep[i] = (u32)d;
        }
    } else {
        for (size_t i = 0; i < raw.size(); ++i) {
            auto it = hash_id.emplace(raw[i], n);
            if (it.second) ++n;
            ep[i] = (u32)it.first->second;
        }
    }
    if (n >= (i64)UINT32_MAX) { fprintf(stderr, "too many nodes\n"); exit(1); }
    std::vector<i64>().swap(raw);

    off.assign(n + 1, 0);
    for (i64 i = 0; i < E; ++i) { off[ep[2 * i] + 1]++; off[ep[2 * i + 1] + 1]++; }
    for (i64 v = 0; v < n; ++v) off[v + 1] += off[v];
    adj.resize(off[n]);
    {
        std::vector<i64> pos(off.begin(), off.end() - 1);
        for (i64 i = 0; i < E; ++i) {
            u32 a = ep[2 * i], b = ep[2 * i + 1];
            adj[pos[a]++] = b; adj[pos[b]++] = a;
        }
    }
    std::vector<u32>().swap(ep);
    // sort and deduplicate each list, compacting in place (writes stay behind the unread lists)
    i64 w = 0;
    for (i64 v = 0; v < n; ++v) {
        i64 s = off[v], e = off[v + 1];
        std::sort(adj.begin() + s, adj.begin() + e);
        i64 start = w;
        for (i64 i = s; i < e; ++i)
            if (i == s || adj[i] != adj[i - 1]) adj[w++] = adj[i];
        off[v] = start;
    }
    off[n] = w;
    adj.resize(w);
    adj.shrink_to_fit();
    printf("graph lines=%lld comments=%lld malformed=%lld out_of_int=%lld selfloops=%lld n=%lld m=%lld\n",
           (long long)lines, (long long)comments, (long long)malformed, (long long)out_of_int,
           (long long)selfloops, (long long)n, (long long)(w / 2));
    fflush(stdout);
}

// Batagelj-Zaversnik core decomposition, O(n + m).
static void cores() {
    core.assign(n, 0);
    i64 maxd = 0;
    for (i64 v = 0; v < n; ++v) maxd = std::max(maxd, off[v + 1] - off[v]);
    if (maxd > INT_MAX) { fprintf(stderr, "degree too large\n"); exit(1); }
    for (i64 v = 0; v < n; ++v) core[v] = (int)(off[v + 1] - off[v]);
    std::vector<i64> bin(maxd + 1, 0), vert(n), pos(n);
    for (i64 v = 0; v < n; ++v) bin[core[v]]++;
    i64 start = 0;
    for (i64 d = 0; d <= maxd; ++d) { i64 c = bin[d]; bin[d] = start; start += c; }
    for (i64 v = 0; v < n; ++v) { pos[v] = bin[core[v]]; vert[pos[v]] = v; bin[core[v]]++; }
    for (i64 d = maxd; d >= 1; --d) bin[d] = bin[d - 1];
    bin[0] = 0;
    for (i64 i = 0; i < n; ++i) {
        i64 v = vert[i];
        for (i64 j = off[v]; j < off[v + 1]; ++j) {
            i64 u = adj[j];
            if (core[u] > core[v]) {
                int du = core[u];
                i64 pu = pos[u], pw = bin[du], w = vert[pw];
                if (u != w) { pos[u] = pw; vert[pu] = w; pos[w] = pu; vert[pw] = u; }
                bin[du]++;
                core[u]--;
            }
        }
    }
    int kmax = 0;
    for (i64 v = 0; v < n; ++v) kmax = std::max(kmax, core[v]);
    printf("cores maxcore=%d maxdeg=%lld\n", kmax, (long long)maxd);
    fflush(stdout);
}

static void query_core(int lam, i64 f) {
    std::vector<char> seen(n, 0);
    std::vector<u32> q;
    i64 in_core = 0, comps = 0, over = 0;
    std::vector<i64> sizes;
    for (i64 s = 0; s < n; ++s) {
        if (core[s] < lam || seen[s]) continue;
        q.clear(); q.push_back((u32)s); seen[s] = 1;
        for (size_t h = 0; h < q.size(); ++h) {
            u32 v = q[h];
            for (i64 j = off[v]; j < off[v + 1]; ++j) {
                u32 u = adj[j];
                if (core[u] >= lam && !seen[u]) { seen[u] = 1; q.push_back(u); }
            }
        }
        ++comps; in_core += (i64)q.size();
        if ((i64)q.size() > f) ++over;
        sizes.push_back((i64)q.size());
    }
    std::sort(sizes.rbegin(), sizes.rend());
    printf("core lam=%d f=%lld nodes=%lld comps=%lld over_f=%lld top=", lam, (long long)f,
           (long long)in_core, (long long)comps, (long long)over);
    for (size_t i = 0; i < sizes.size() && i < 3; ++i) printf("%s%lld", i ? "," : "", (long long)sizes[i]);
    printf("\n");
    fflush(stdout);
}

static void query_sol(const std::string& path) {
    std::ifstream in(path);
    std::string line, nodes;
    bool found = false;
    while (std::getline(in, line))
        if (line.rfind("nodes=", 0) == 0) { nodes = line.substr(6); found = true; }
    std::vector<i64> ids;
    i64 missing = 0, bad = found ? 0 : 1;
    {
        std::stringstream ss(nodes);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            const char* p = tok.c_str();
            const char* e = p + tok.size();
            i64 r;
            if (!parse_int(p, e, r) || p != e) { ++bad; continue; }
            i64 v = lookup(r);
            if (v < 0) ++missing; else ids.push_back(v);
        }
    }
    std::sort(ids.begin(), ids.end());
    const i64 listed = (i64)ids.size();
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    const i64 dup = listed - (i64)ids.size();
    std::vector<char> in_set(n, 0);
    for (i64 v : ids) in_set[v] = 1;
    i64 mindeg = -1;
    for (i64 v : ids) {
        i64 d = 0;
        for (i64 j = off[v]; j < off[v + 1]; ++j) d += in_set[adj[j]];
        mindeg = (mindeg < 0) ? d : std::min(mindeg, d);
    }
    i64 reached = 0;                               // nodes of the set reachable inside the set
    if (!ids.empty()) {
        std::vector<char> seen(n, 0);
        std::vector<u32> q{(u32)ids[0]};
        seen[ids[0]] = 1;
        for (size_t h = 0; h < q.size(); ++h)
            for (i64 j = off[q[h]]; j < off[q[h] + 1]; ++j) {
                u32 u = adj[j];
                if (in_set[u] && !seen[u]) { seen[u] = 1; q.push_back(u); }
            }
        reached = (i64)q.size();
    }
    i64 kcomp = -1;                                // component of the mindeg-core containing the set
    if (!ids.empty() && mindeg >= 1) {
        std::vector<char> seen(n, 0);
        std::vector<u32> q{(u32)ids[0]};
        seen[ids[0]] = 1;
        for (size_t h = 0; h < q.size(); ++h)
            for (i64 j = off[q[h]]; j < off[q[h] + 1]; ++j) {
                u32 u = adj[j];
                if (core[u] >= mindeg && !seen[u]) { seen[u] = 1; q.push_back(u); }
            }
        kcomp = (i64)q.size();
    }
    printf("sol size=%lld bad=%lld missing=%lld dup=%lld connected=%d mindeg=%lld mindeg_core_comp=%lld\n",
           (long long)ids.size(), (long long)bad, (long long)missing, (long long)dup,
           (int)(!ids.empty() && reached == (i64)ids.size()), (long long)mindeg, (long long)kcomp);
    fflush(stdout);
}

int main(int argc, char** argv) {
    if (argc != 3) { fprintf(stderr, "usage: %s <edge list> <query file>\n", argv[0]); return 1; }
    std::ifstream qf(argv[2]);
    if (!qf) { perror(argv[2]); return 1; }
    load(argv[1]);
    cores();
    std::string line;
    while (std::getline(qf, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        std::string kind;
        ls >> kind;
        int lam; long long f;
        if (kind == "core" && (ls >> lam >> f)) query_core(lam, f);
        else if (kind == "sol" && line.size() > 4) query_sol(line.substr(4));   // the path may contain spaces
        else { fprintf(stderr, "bad query: %s\n", line.c_str()); return 1; }
    }
    if (qf.bad()) { perror(argv[2]); return 1; }
    return 0;
}
