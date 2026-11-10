#!/usr/bin/env python3
"""Validate the flexi binary against exhaustive enumeration on random small graphs.

For each graph, the maximum flexi-clique is found by checking every node subset, and
  * EBA must return a valid flexi-clique of exactly the maximum size,
  * EBA with --seed q must return a valid flexi-clique containing q of the maximum size
    among those containing q,
  * FPA must return a valid flexi-clique (or nothing) no larger than the maximum.

Usage: python3 tests/validate_bruteforce.py [BINARY] [N_GRAPHS] [N_MAX] [WORKERS]
"""
import math, os, random, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations

TAUS = [0.25, 0.5, 0.5, 0.75, 0.2, 0.35, 0.65, 0.8]   # 0.25/0.5/0.75 hit integer (d+1)^(1/tau)
PS = [0.2, 0.35, 0.5, 0.65, 0.8]


def threshold(k, tau):
    # floor(k^tau), with the same 1e-12 tolerance as floorPow in src/graph.h
    return math.floor(k ** tau * (1 + 1e-12) + 1e-12)


def is_flexi(adj, H, tau):
    H = set(H)
    if len(H) < 2:
        return False
    t = threshold(len(H), tau)
    if any(len(adj[v] & H) < t for v in H):
        return False
    start = next(iter(H)); seen = {start}; stack = [start]
    while stack:
        for w in adj[stack.pop()] & H:
            if w not in seen:
                seen.add(w); stack.append(w)
    return seen == H


def brute_max(adj, V, tau, must=None):
    for r in range(len(V), 1, -1):
        for c in combinations(V, r):
            if (must is None or must in c) and is_flexi(adj, c, tau):
                return r
    return 0


def run(binary, args):
    p = subprocess.run([binary] + args, capture_output=True, text=True, timeout=600)
    kv = {l.split(":")[0].strip(): l.split(":", 1)[1].strip() for l in p.stdout.splitlines() if ":" in l}
    if p.returncode != 0 or "Size" not in kv or "Nodes" not in kv:
        raise RuntimeError(f"{args}: exit {p.returncode}, output {p.stdout!r}")
    nodes = [int(x) for x in kv["Nodes"].split(",") if x]
    if len(nodes) != len(set(nodes)) or len(nodes) != int(kv["Size"]):
        raise RuntimeError(f"{args}: duplicate nodes or size mismatch in {nodes}")
    return nodes


def check(args):
    binary, sd, n_max = args
    r = random.Random(sd)
    n, p, tau = r.randint(4, n_max), r.choice(PS), r.choice(TAUS)
    adj = {i: set() for i in range(n)}
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if r.random() < p]
    if not edges:
        return []
    for u, v in edges:
        adj[u].add(v); adj[v].add(u)
    V = sorted({x for e in edges for x in e})
    best = brute_max(adj, V, tau)
    fails = []
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.writelines(f"{u} {v}\n" for u, v in edges)
        path = f.name
    try:
      try:
        base = ["--file", path, "--tau", str(tau)]
        got = run(binary, ["--algo", "eba"] + base)
        if len(got) != best or (best and not is_flexi(adj, got, tau)):
            fails.append(f"sd={sd} n={n} tau={tau} EBA: max={best} got={got}")
        q = r.choice(V)
        best_q = brute_max(adj, V, tau, must=q)
        got = run(binary, ["--algo", "eba", "--seed", str(q)] + base)
        if len(got) != best_q or (best_q and (q not in got or not is_flexi(adj, got, tau))):
            fails.append(f"sd={sd} n={n} tau={tau} EBA --seed {q}: max={best_q} got={got}")
        got = run(binary, ["--algo", "fpa"] + base)
        if len(got) > best or (got and not is_flexi(adj, got, tau)):
            fails.append(f"sd={sd} n={n} tau={tau} FPA: max={best} got={got}")
      except RuntimeError as e:
        fails.append(f"sd={sd} n={n} tau={tau} {e}")
    finally:
        os.remove(path)
    return fails


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    binary = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "build", "flexi")
    n_graphs = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    n_max = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else os.cpu_count()
    fails = []
    with ProcessPoolExecutor(workers) as ex:
        for fl in ex.map(check, [(binary, sd, n_max) for sd in range(n_graphs)], chunksize=8):
            fails.extend(fl)
    for fl in fails[:20]:
        print("FAIL", fl)
    print(f"{n_graphs} random graphs (n <= {n_max}): {len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
