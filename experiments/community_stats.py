#!/usr/bin/env python3

import argparse, csv, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


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


def stats(adj, S):
    """(size, internal edges, density) of the node set S."""
    s = len(S)
    e = sum(len(adj[v] & S) for v in S) // 2
    return s, e, (2 * e / (s * (s - 1)) if s > 1 else "")


def write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["community", "size", "edges", "density"])
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("networks", nargs="*", default=["AM", "DB"], help="graph abbreviations (default: AM DB)")
    ap.add_argument("--graphs", default=os.path.join(ROOT, "datasets", "graphs"), help="directory of <ABBR>.txt")
    ap.add_argument("--communities", default=os.path.join(ROOT, "datasets", "communities"),
                    help="directory of <ABBR>.all.txt and <ABBR>.top5000.txt")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "community_stats"), help="output directory")
    a = ap.parse_args()
    nets = [n.upper() for n in a.networks]
    for net in nets:
        for p in [os.path.join(a.graphs, net + ".txt")] + [os.path.join(a.communities, f"{net}.{k}.txt")
                                                           for k in ("all", "top5000")]:
            if not os.path.exists(p):
                sys.exit(f"{p} not found (python3 datasets/prepare_datasets.py {net} --communities {net})")
    os.makedirs(a.out, exist_ok=True)

    for net in nets:
        adj = load_graph(os.path.join(a.graphs, net + ".txt"))
        comms = load_communities(os.path.join(a.communities, net + ".all.txt"))
        write(os.path.join(a.out, f"{net}.all.csv"),
              [(i,) + stats(adj, {v for v in C if v in adj}) for i, C in enumerate(comms)])
        top = load_communities(os.path.join(a.communities, net + ".top5000.txt"))
        rows, seen = [], set()
        for i, C in enumerate(top):
            S = frozenset(v for v in C if v in adj)
            if len(S) >= 3 and S not in seen:
                seen.add(S)
                rows.append((i,) + stats(adj, S))
        write(os.path.join(a.out, f"{net}.top5000.csv"), rows)
        print(f"{net}: {len(comms)} communities -> {net}.all.csv; {len(rows)} distinct of the {len(top)} "
              f"top-5000 communities -> {net}.top5000.csv (in {a.out})", flush=True)


if __name__ == "__main__":
    main()
