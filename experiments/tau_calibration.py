#!/usr/bin/env python3

import argparse, csv, os, random, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 2027


def load_graph(path):
    """Adjacency sets of an edge list (comment lines # and %, and self-loops, are skipped)."""
    adj = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line[0] in "#%":
                continue
            p = line.split()
            if len(p) < 2:
                continue
            u, v = int(p[0]), int(p[1])
            if u != v:
                adj.setdefault(u, set()).add(v)
                adj.setdefault(v, set()).add(u)
    return adj


def load_communities(path):
    """One community per line, node ids separated by whitespace."""
    with open(path) as f:
        return [[int(x) for x in line.split()] for line in f if line.strip() and line[0] not in "#%"]


def connected(adj, C):
    S, seen, stack = set(C), {C[0]}, [C[0]]
    while stack:
        for y in adj[stack.pop()] & S:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(S)


def min_degree(adj, S):
    """Minimum degree of the subgraph induced by the node set S."""
    return min(len(adj[v] & S) for v in S)


def random_set(adj, q, size, rng):
    """Grow {q} by nodes drawn uniformly from the frontier (the nodes adjacent to the set). The frontier
    is a list with O(1) removal, filled in id order, so the result depends only on the state of rng."""
    S, front, pos = {q}, [], {}

    def push(v):
        for y in sorted(adj[v]):
            if y not in S and y not in pos:
                pos[y] = len(front)
                front.append(y)

    push(q)
    while len(S) < size and front:
        i = rng.randrange(len(front))
        v, last = front[i], front.pop()
        if i < len(front):
            front[i] = last
            pos[last] = i
        del pos[v]
        S.add(v)
        push(v)
    return S


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("networks", nargs="*", default=["AM", "DB"], help="graph abbreviations (default: AM DB)")
    ap.add_argument("--graphs", default=os.path.join(ROOT, "graphs"), help="directory of <ABBR>.txt")
    ap.add_argument("--communities", default=os.path.join(ROOT, "graphs", "communities"),
                    help="directory of <ABBR>.top5000.txt")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "tau_calibration"), help="output directory")
    ap.add_argument("--seeds", type=int, default=5, help="random sets per community (paper: 5)")
    a = ap.parse_args()
    if a.seeds <= 0:
        sys.exit("--seeds must be positive")
    nets = [n.upper() for n in a.networks]
    for net in nets:
        for p in (os.path.join(a.graphs, net + ".txt"), os.path.join(a.communities, net + ".top5000.txt")):
            if not os.path.exists(p):
                sys.exit(f"{p} not found (python3 graphs/prepare_datasets.py {net} --communities {net})")
    os.makedirs(a.out, exist_ok=True)

    for net in nets:
        t0 = time.time()
        adj = load_graph(os.path.join(a.graphs, net + ".txt"))
        comms = load_communities(os.path.join(a.communities, net + ".top5000.txt"))
        restricted = [sorted({v for v in C if v in adj}) for C in comms]
        eligible = [C for C in restricted if len(C) >= 3 and connected(adj, C)]
        if not eligible:
            sys.exit(f"{net}: no community with at least 3 nodes that induces a connected subgraph")
        rng = random.Random(SEED)
        rows = []
        for C in eligible:
            Cs = set(C)
            q = max(C, key=lambda v: (len(adj[v] & Cs), -v))
            rnd = [random_set(adj, q, len(C), rng) for _ in range(a.seeds)]
            assert all(len(s) == len(C) for s in rnd), "a connected community fits in the component of q"
            rows.append([len(C), min_degree(adj, Cs)] + [min_degree(adj, s) for s in rnd])
        path = os.path.join(a.out, f"records_{net}.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["size", "dmin_community"] + [f"dmin_random{i}" for i in range(a.seeds)])
            w.writerows(rows)
        print(f"{net}: {len(comms)} listed communities ({len({tuple(C) for C in restricted})} distinct node sets), "
              f"{len(eligible)} evaluated, {a.seeds} random sets each -> {path} ({time.time() - t0:.0f} s)",
              flush=True)


if __name__ == "__main__":
    main()
