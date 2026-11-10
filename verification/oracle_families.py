#!/usr/bin/env python3
"""Compare FPA and EBA with exhaustive enumeration on small graphs of several structures.

Families: Erdos-Renyi, a planted dense block, two dense blocks joined by a few edges, circulant
graphs, core-periphery graphs, and random connected induced subgraphs (BFS balls) of FB, CA, and DB.
Each graph has 8 to 16 nodes (10 to 16 for the balls). For every graph, all node subsets are
enumerated once and, for each size s, the largest minimum degree over connected subsets is kept, so
the optimum is max{s >= 2 : that degree >= floor(s^tau)} for each tau in 0.3..0.7 (thresholds are
computed exactly with integers). The node sets returned by flexi are checked independently.

A run fails if EBA does not return a valid flexi-clique of maximum size, if FPA returns an invalid
set or one larger than the maximum, or if flexi exits with an error. The table reports, per family,
how often FPA is optimal and how often EBA finishes without branching.

  python3 verification/oracle_families.py [--graphs-per-family 500] [--workers 8]

Exit status 1 if a run fails.
"""
import argparse, multiprocessing as mp, os, random, subprocess, sys, tempfile
from collections import Counter, defaultdict
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAUS = ["0.3", "0.4", "0.5", "0.6", "0.7"]
FAMILIES = ["er", "planted", "twoblocks", "circulant", "coreper"]
BALLS = ["FB", "CA", "DB"]


def iroot(x, q):
    lo, hi = 0, 1
    while hi ** q <= x:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if mid ** q <= x else (lo, mid)
    return lo


def floor_pow(s, tau):
    fr = Fraction(tau)
    return iroot(s ** fr.numerator, fr.denominator)


def connected(adj, mask):
    reach = front = mask & -mask
    while front:
        nb, f = 0, front
        while f:
            b = f & -f
            nb |= adj[b.bit_length() - 1]
            f ^= b
        front = nb & mask & ~reach
        reach |= front
    return reach == mask


def optimum(n, edges):
    adj = [0] * n
    for a, b in edges:
        adj[a] |= 1 << b
        adj[b] |= 1 << a
    maxmd = [-1] * (n + 1)
    for S in range(1, 1 << n):
        s = bin(S).count("1")
        if s < 2:
            continue
        md = min(bin(adj[v] & S).count("1") for v in range(n) if S >> v & 1)
        if md >= 1 and md > maxmd[s] and connected(adj, S):
            maxmd[s] = md
    opt = {t: max([s for s in range(2, n + 1) if maxmd[s] >= floor_pow(s, t)], default=0) for t in TAUS}
    return opt, adj


def is_flexi(adj, S, tau):
    if S is None or len(S) < 2 or len(set(S)) != len(S) or any(v >= len(adj) for v in S):
        return False
    mask = sum(1 << v for v in S)
    return connected(adj, mask) and min(bin(adj[v] & mask).count("1") for v in S) >= floor_pow(len(S), tau)


def run(binary, path, algo, tau):
    """Run flexi; return (return code, printed fields, node list or None)."""
    out = path + f".{algo}.{tau}.out"
    r = subprocess.run([binary, "--algo", algo, "--file", path, "--tau", tau, "--out", out],
                       capture_output=True, text=True, timeout=600)
    kv = {}
    for line in r.stdout.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            kv[k.strip()] = v.strip()
    nodes = None
    if os.path.exists(out):
        for line in open(out):
            if line.startswith("nodes="):
                toks = [x for x in line.strip()[6:].split(",") if x]
                nodes = [int(x) for x in toks] if all(x.isdigit() for x in toks) else None
        os.remove(out)
    return r.returncode, kv, nodes


def generate(family, rng):
    n = rng.randint(8, 16)
    E = set()
    def add(a, b):
        if a != b:
            E.add((min(a, b), max(a, b)))
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    if family == "er":
        p = rng.uniform(0.15, 0.9)
        for a, b in pairs:
            if rng.random() < p: add(a, b)
    elif family == "planted":
        k, p_in, p_out = rng.randint(4, n - 2), rng.uniform(0.6, 1.0), rng.uniform(0.05, 0.35)
        for a, b in pairs:
            if rng.random() < (p_in if b < k else p_out): add(a, b)
    elif family == "twoblocks":
        k, p_in = n // 2, rng.uniform(0.5, 1.0)
        for a, b in pairs:
            if (a < k) == (b < k) and rng.random() < p_in: add(a, b)
        for _ in range(rng.randint(1, 3)):
            add(rng.randrange(k), rng.randrange(k, n))
    elif family == "circulant":
        for o in rng.sample(range(1, n // 2 + 1), rng.randint(1, max(1, n // 4))):
            for a in range(n): add(a, (a + o) % n)
        for _ in range(rng.randint(0, 3)):
            add(rng.randrange(n), rng.randrange(n))
    elif family == "coreper":
        k = rng.randint(3, n // 2)
        for a, b in pairs:
            if rng.random() < (0.9 if b < k else 0.4 if a < k else 0.08): add(a, b)
    return n, sorted(E)


def ball(rng, adj, nodes):
    """Random connected induced subgraph (BFS ball) of a real graph, relabelled 0..n-1."""
    target, s = rng.randint(10, 16), rng.choice(nodes)
    seen, queue = [s], [s]
    while queue and len(seen) < target:
        nb = sorted(adj[queue.pop(0)])
        rng.shuffle(nb)
        for w in nb:
            if w not in seen and len(seen) < target:
                seen.append(w); queue.append(w)
    idx = {v: i for i, v in enumerate(seen)}
    return len(seen), sorted({(min(idx[a], idx[b]), max(idx[a], idx[b])) for a in seen for b in adj[a] if b in idx})


def job(args):
    binary, family, n, E = args
    if not E:
        return family, None
    opt, adj = optimum(n, E)
    fd, path = tempfile.mkstemp(suffix=".txt")
    with os.fdopen(fd, "w") as f:
        f.writelines(f"{a} {b}\n" for a, b in E)
    res = []
    try:
        for t in TAUS:
            row = [t, opt[t]]
            for algo in ("fpa", "eba"):
                rc, kv, S = run(binary, path, algo, t)
                ok = (rc == 0 and kv.get("Valid") == "yes" and S is not None and kv.get("Size") == str(len(S))
                      and is_flexi(adj, S, t))
                row += [len(S) if S is not None else -1, ok]
            res.append(row + [kv.get("Branches", "")])        # branches of the EBA run
    finally:
        os.remove(path)
    return family, res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--graphs-per-family", type=int, default=500)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--graphs", default=os.path.join(ROOT, "graphs"), help="edge lists <ABBR>.txt")
    ap.add_argument("--build", default=os.path.join(ROOT, "build"), help="CMake build directory")
    ap.add_argument("--seed", type=int, default=20261104)
    a = ap.parse_args()
    binary = os.path.join(a.build, "flexi")

    rng = random.Random(a.seed)
    jobs = []
    for fam in FAMILIES:
        for _ in range(a.graphs_per_family):
            jobs.append((binary, fam, *generate(fam, random.Random(rng.randrange(10 ** 9)))))
    for name in BALLS:
        adj = defaultdict(set)
        for line in open(os.path.join(a.graphs, f"{name}.txt")):
            t = line.split()
            if line.startswith("#") or len(t) < 2 or t[0] == t[1]:
                continue
            adj[int(t[0])].add(int(t[1])); adj[int(t[1])].add(int(t[0]))
        nodes = sorted(adj)
        for _ in range(a.graphs_per_family):
            jobs.append((binary, "ball_" + name, *ball(random.Random(rng.randrange(10 ** 9)), adj, nodes)))

    stats, fails = defaultdict(Counter), []
    with mp.Pool(a.workers) as pool:
        for fam, res in pool.imap_unordered(job, jobs, chunksize=4):
            st = stats[fam]
            if res is None:
                st["no_edges"] += 1
                continue
            st["graphs"] += 1
            for t, opt, fs, fok, es, eok, br in res:
                st["runs"] += 1
                if not eok or es != opt: fails.append(("EBA", fam, t, opt, es, eok))
                if not fok or fs > opt: fails.append(("FPA", fam, t, opt, fs, fok))
                st["fpa_optimal"] += fok and fs == opt
                st["ratio"] += fs / opt if fok else 0
                st["eba_no_branching"] += eok and br == "0"
    print(f"{'family':<10}{'graphs':>7}{'runs':>7}{'FPA optimal':>14}{'mean FPA/OPT':>14}{'EBA without branching':>23}")
    for fam in FAMILIES + ["ball_" + b for b in BALLS]:
        st = stats[fam]
        if st["runs"]:
            print(f"{fam:<10}{st['graphs']:>7}{st['runs']:>7}{100 * st['fpa_optimal'] / st['runs']:>13.1f}%"
                  f"{100 * st['ratio'] / st['runs']:>13.1f}%{100 * st['eba_no_branching'] / st['runs']:>22.1f}%")
    print(f"failed runs: {len(fails)}")
    for x in fails[:20]:
        print("  ", x)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
