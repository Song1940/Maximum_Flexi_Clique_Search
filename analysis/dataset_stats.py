#!/usr/bin/env python3
"""Statistics of Table II: |V|, |E|, maximum coreness k*, average degree, and global clustering
coefficient CC of edge-list files.

  python3 analysis/dataset_stats.py graphs/FB.txt graphs/AM.txt ... [--csv FILE]

A graph is read as `flexi` reads it: one edge "u v" per line (further columns are ignored), lines
starting with # or % skipped, self-loops dropped, duplicate and reverse edges merged. |V| counts the
nodes with at least one edge and |E| the undirected edges; average degree 2|E| / |V|; k* is the largest
k such that the k-core is nonempty; CC = 3 x #triangles / #connected triples = 6T / sum_v d(v)(d(v) - 1).
Node ids must be integers in [0, 2^32 - 2].

The computation is vectorised with numpy (sorted edge keys, CSR, k-core peeling in batches, triangle
counting over the edges oriented from lower to higher degree). It is slow for the largest graphs:
seconds up to FD, a few minutes for SO and LJ, about 8 minutes for OK, and an hour or more each for
IT and SK (the paper used a parallel C++ implementation of the same computation). Memory is about
40 bytes per edge (40 GB for IT, about 75 GB for SK).
Needs numpy.
"""
import argparse, csv, os, sys, time, warnings

import numpy as np

MASK = np.uint64(0xFFFFFFFF)
SHIFT = np.uint64(32)
MAX_ID = 2 ** 32 - 2
CHUNK = 1 << 24                       # elements per step of the chunked array operations


def parse_block(data):
    """(u, v, number of malformed lines) of the complete lines in the bytes data."""
    lines = data.count(b"\n") + (not data.endswith(b"\n"))
    if b"#" not in data and b"%" not in data:
        # fast path, if every line holds the same number c >= 2 of integer tokens
        b = np.frombuffer(data, dtype=np.uint8)
        ws = (b == 32) | (b == 9) | (b == 10) | (b == 13)
        start = ~ws
        start[1:] &= ws[:-1]                                          # first byte of every token
        per_line = np.bincount(np.cumsum(b == 10)[start], minlength=lines)
        del b, ws, start
        c = int(per_line[0]) if per_line.size else 0
        if c >= 2 and per_line.size == lines and np.all(per_line == c):
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                try:
                    a = np.fromstring(data, dtype=np.int64, sep=" ")
                except (DeprecationWarning, ValueError):
                    a = None
            if a is not None and a.size == c * lines:
                return a[0::c], a[1::c], 0
    us, vs, bad = [], [], 0                                           # general path, line by line
    for line in data.split(b"\n"):
        tok = line.split()
        if not tok or tok[0][:1] in (b"#", b"%"):
            continue
        if len(tok) < 2 or not tok[0].isdigit() or not tok[1].isdigit():
            bad += 1
            continue
        us.append(int(tok[0]))
        vs.append(int(tok[1]))
    return np.array(us, dtype=np.int64), np.array(vs, dtype=np.int64), bad


def read_edges(path, block=1 << 24):
    """Sorted array of the distinct undirected edges, as keys (u << 32) | v with u < v."""
    parts, bad, rest = [], 0, b""
    with open(path, "rb") as f:
        while True:
            buf = f.read(block)
            if buf:
                cut = buf.rfind(b"\n") + 1
                if cut == 0:                                          # no complete line yet
                    rest += buf
                    continue
                data, rest = rest + buf[:cut], buf[cut:]
            else:
                data, rest = rest, b""
            if data.strip():
                u, v, nbad = parse_block(data)
                bad += nbad
                if u.size and (min(u.min(), v.min()) < 0 or max(u.max(), v.max()) > MAX_ID):
                    sys.exit(f"{path}: node ids must be in [0, {MAX_ID}]")
                keep = u != v
                u, v = u[keep], v[keep]
                parts.append((np.minimum(u, v).astype(np.uint64) << SHIFT) | np.maximum(u, v).astype(np.uint64))
            if not buf:
                break
    if bad:
        print(f"  {path}: skipped {bad} lines without two node ids", file=sys.stderr)
    keys = np.concatenate(parts) if parts else np.zeros(0, dtype=np.uint64)
    del parts
    keys.sort()
    if keys.size:
        keys = keys[np.concatenate(([True], keys[1:] != keys[:-1]))]
    return keys


def gather(off, adj, nodes):
    """Concatenated adjacency lists adj[off[x]:off[x+1]] of the given nodes."""
    starts = off[nodes]
    lens = off[nodes + 1] - starts
    total = int(lens.sum())
    if total == 0:
        return adj[:0]
    return adj[np.repeat(starts - (np.cumsum(lens) - lens), lens) + np.arange(total)]


def pieces(weights, budget):
    """Split range(len(weights)) into consecutive ranges [a, b) whose weights sum to about budget or less
    (a single element may exceed it)."""
    cum = np.cumsum(weights)
    if cum.size == 0:
        return []
    cuts = np.searchsorted(cum, np.arange(budget, int(cum[-1]), budget), side="right")
    bounds = np.unique(np.concatenate(([0], cuts, [len(weights)])))
    return list(zip(bounds[:-1].tolist(), bounds[1:].tolist()))


def max_coreness(off, adj, deg, budget=1 << 26):
    """Largest k with a nonempty k-core: for k = 0, 1, ..., remove the nodes of degree <= k, in rounds,
    until every remaining node has degree > k."""
    deg = deg.copy()
    n = len(deg)
    alive = deg > 0
    live = np.flatnonzero(alive)
    k = 0
    while live.size:
        k = max(k, int(deg[live].min()))
        rem = live[deg[live] <= k]
        while rem.size:
            alive[rem] = False
            touched = []
            for a, b in pieces(off[rem + 1] - off[rem], budget):
                nb = gather(off, adj, rem[a:b])
                nb = nb[alive[nb]]
                if nb.size > n // 4:
                    cnt = np.bincount(nb, minlength=n)
                    deg -= cnt
                    touched.append(np.flatnonzero(cnt))
                elif nb.size:
                    x, c = np.unique(nb, return_counts=True)
                    deg[x] -= c
                    touched.append(x)
            t = np.unique(np.concatenate(touched)) if touched else live[:0]
            rem = t[deg[t] <= k]
        live = live[alive[live]]
    return k


def triangles(okeys, ooff, budget=1 << 20, block=1 << 16):
    """Triangles of the graph whose edges u -> w, oriented from lower to higher rank, are the sorted keys
    (u << 32) | w, with CSR offsets ooff: for every edge u -> v and every w with v -> w, test u -> w.
    Every triangle is counted once, from its lowest-ranked node."""
    out = np.empty(len(okeys), dtype=np.uint32)
    for i in range(0, len(okeys), CHUNK):
        out[i:i + CHUNK] = okeys[i:i + CHUNK] & MASK
    dout = np.diff(ooff)
    tri = 0
    for e0 in range(0, len(okeys), block):
        src = (okeys[e0:e0 + block] >> SHIFT).astype(np.int64)
        v = out[e0:e0 + block].astype(np.int64)
        lens = dout[v]
        for a, b in pieces(lens, budget):
            w = gather(ooff, out, v[a:b])
            if w.size == 0:
                continue
            q = (np.repeat(src[a:b], lens[a:b]).astype(np.uint64) << SHIFT) | w.astype(np.uint64)
            region = okeys[ooff[src[a]]:ooff[src[b - 1] + 1]]               # the out-edges of these sources
            pos = np.searchsorted(region, q)
            pos[pos == len(region)] = 0
            tri += int(np.count_nonzero(region[pos] == q))
    return tri


def graph_stats(path):
    keys = read_edges(path)
    m = len(keys)
    if m == 0:
        return dict(V=0, E=0, kmax=0, avg_deg=0.0, triangles=0, cc=0.0)
    hi = (keys & MASK).astype(np.int64)
    n = int(hi.max()) + 1
    deg = np.bincount(hi, minlength=n)
    del hi
    deg += np.bincount((keys >> SHIFT).astype(np.int64), minlength=n)
    nv = int(np.count_nonzero(deg))
    wedges = int((deg * (deg - 1) // 2).sum())

    # undirected CSR with sorted adjacency lists (both directions of every edge), for the k-core
    both = np.empty(2 * m, dtype=np.uint64)
    both[:m] = keys
    np.left_shift(keys, SHIFT, out=both[m:])
    both[m:] |= keys >> SHIFT
    both.sort()
    adj = np.empty(2 * m, dtype=np.uint32)
    for i in range(0, 2 * m, CHUNK):
        adj[i:i + CHUNK] = both[i:i + CHUNK] & MASK
    del both
    off = np.concatenate(([0], np.cumsum(deg)))
    kmax = max_coreness(off, adj, deg)
    del adj

    # orient every edge u - v (u < v) from lower to higher (degree, id) rank, in place
    for i in range(0, m, CHUNK):
        k = keys[i:i + CHUNK]
        u, v = (k >> SHIFT).astype(np.int64), (k & MASK).astype(np.int64)
        rev = deg[u] > deg[v]
        k[rev] = (v[rev].astype(np.uint64) << SHIFT) | u[rev].astype(np.uint64)
    keys.sort()
    ooff = np.searchsorted(keys, np.arange(n + 1, dtype=np.uint64) << SHIFT).astype(np.int64)
    tri = triangles(keys, ooff)
    return dict(V=nv, E=m, kmax=kmax, avg_deg=2 * m / nv, triangles=tri, cc=3 * tri / wedges if wedges else 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="edge-list files")
    ap.add_argument("--csv", help="also write the statistics to this CSV file")
    a = ap.parse_args()
    for path in a.files:
        if not os.path.exists(path):
            sys.exit(f"{path} not found")
    rows = []
    for path in a.files:
        t0 = time.time()
        r = dict(graph=os.path.splitext(os.path.basename(path))[0], **graph_stats(path))
        r["seconds"] = round(time.time() - t0, 1)
        rows.append(r)
        print(f"{r['graph']}: |V|={r['V']} |E|={r['E']} k*={r['kmax']} avg_deg={r['avg_deg']:.4f} "
              f"triangles={r['triangles']} CC={r['cc']:.6f} ({r['seconds']} s)", flush=True)
    print(f"\n{'Graph':<12} {'|V|':>13} {'|E|':>15} {'k*':>7} {'Avg. Deg':>9} {'CC':>6}")
    for r in rows:
        print(f"{r['graph']:<12} {r['V']:>13,} {r['E']:>15,} {r['kmax']:>7,} {r['avg_deg']:>9.2f} {r['cc']:>6.3f}")
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
