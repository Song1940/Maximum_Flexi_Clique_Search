#!/usr/bin/env python3
"""Certify the optimality of FPA solutions independently of EBA (the instances of EQ1).

For a valid flexi-clique of size f, let lam = floor((f+1)^tau). Every flexi-clique X with |X| > f is
connected, and each of its nodes has at least floor(|X|^tau) >= lam neighbours in G[X]; so X lies in
one connected component of the lam-core of G, which then has more than f nodes. Hence, if every
connected component of the lam-core has at most f nodes, no flexi-clique has more than f nodes.
With incumbent f, this is also the condition under which the global pruning of EBA (Rules 1 and 3)
removes every node, so the certified instances should be exactly those EBA settles without branching.

For each instance the script
  1. runs FPA (build/flexi --algo fpa --out) and keeps its solution, and
  2. checks with verification/indep_check, which shares no code with src/, that
       - the edge list is read the same way by both loaders (every line is empty, a comment, or
         starts with two integer tokens within the range of int) and gives the same n and m,
       - the FPA solution is a flexi-clique of the reported size f: connected, with minimum degree
         at least floor(f^tau) (thresholds are computed exactly with integers),
       - every connected component of the lam-core has at most f nodes.
An instance is certified only if all checks pass. --results adds the EBA outcome of a results file
of experiments/run_experiments.py for comparison: a certified instance on which EBA branched, an
instance EBA solved without branching that is not certified, or an FPA size that differs from the
results file is a failure.

  python3 verification/certify_optimality.py [--datasets FB,GO] [--taus 0.5] [--workers 4]
                                             [--results results/results.csv]

Output: results/certificate/certificate.csv, one row per instance, and a summary on stdout; the FPA
solutions are kept in results/certificate/fpa/. Exit status 1 if a graph or an FPA solution fails a
check.
"""
import argparse, csv, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL = "FB GO FW CA CD RO LU DB AM CS WI LJ SK".split()
TAUS = ["0.5", "0.7", "0.9"]
# dataset names in the result files of the paper's runs
FULL = {"fb-pages-company": "FB", "soc-gplus": "GO", "ia-frwikinews-user-edits": "FW", "CAIDA": "CA",
        "condensed": "CD", "soc-gemsec-RO": "RO", "inf-luxembourg_osm": "LU", "dblp": "DB", "amazon": "AM",
        "co-papers-citeseer": "CS", "web-it-2004": "WI", "soc_lj1": "LJ", "sk-2005": "SK"}
FIELDS = ["dataset", "tau", "fpa_size", "lambda", "lambda_core_nodes", "lambda_core_components",
          "largest_components", "graph_check", "witness_valid", "certified", "structure",
          "eba_status", "eba_size", "eba_branches"]


def iroot(x, q):
    """Largest integer r with r**q <= x, for x >= 0."""
    lo, hi = 0, 1
    while hi ** q <= x:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if mid ** q <= x else (lo, mid)
    return lo


def floor_pow(s, tau):
    """floor(s^tau), exactly; tau is a decimal string."""
    fr = Fraction(tau)
    return iroot(s ** fr.numerator, fr.denominator)


def size_cap(k, tau):
    """Largest size s with floor(s^tau) <= k."""
    lo, hi = 1, 2
    while floor_pow(hi, tau) <= k:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if floor_pow(mid, tau) <= k else (lo, mid)
    return lo


def kv_line(line):
    return dict(x.split("=", 1) for x in line.split()[1:])


def run_fpa(binary, graph, tau, sol):
    """FPA solution file `sol` (reused if present); returns (n, m, size) printed by flexi."""
    log = sol + ".log"
    if not (os.path.exists(sol) and os.path.exists(log)):
        r = subprocess.run([binary, "--algo", "fpa", "--file", graph, "--tau", tau, "--out", sol + ".part"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"flexi failed on {graph} tau={tau}: {r.stderr.strip()[-300:]}")
        with open(log + ".part", "w") as f:
            f.write(r.stdout)
        os.replace(sol + ".part", sol)
        os.replace(log + ".part", log)
    text = open(log).read()
    g = re.search(r"Loaded graph: n=(\d+) m=(\d+)", text)
    s = re.search(r"^Size\s*:\s*(\d+)", text, re.M)
    if not g or not s:
        raise RuntimeError(f"unexpected output in {log}")
    return int(g.group(1)), int(g.group(2)), int(s.group(1))


def check_dataset(checker, graph, items, qfile):
    """items: (tau, f, solution file). Returns (n, m, maxcore, [(core fields, solution fields)])."""
    queries = []
    for tau, f, sol in items:
        queries += [f"core {floor_pow(f + 1, tau)} {f}", f"sol {sol}"]
    with open(qfile, "w") as fq:
        fq.write("\n".join(queries) + "\n")
    r = subprocess.run([checker, graph, qfile], capture_output=True, text=True)
    out = r.stdout.splitlines()
    if r.returncode != 0 or len(out) != 2 + len(queries):
        raise RuntimeError(f"indep_check failed on {graph}: {r.stderr.strip()[-300:]}")
    g = kv_line(out[0])
    maxcore = int(kv_line(out[1])["maxcore"])
    pairs = [(kv_line(out[2 + 2 * i]), kv_line(out[3 + 2 * i])) for i in range(len(items))]
    return g, maxcore, pairs


def read_results(path):
    """EBA and FPA rows per (dataset, tau) from a results file of run_experiments.py or of the paper's runs."""
    rows = {"eba": {}, "fpa": {}}
    for r in csv.DictReader(open(path)):
        if r.get("algo") not in rows or r.get("group", r.get("exp")) not in ("main",):
            continue
        d = FULL.get(r["dataset"], r["dataset"])
        rows[r["algo"]][(d, Fraction(r["tau"]))] = r
    return rows["eba"], rows["fpa"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--datasets", default=",".join(REAL), help="comma-separated (default: the 13 graphs)")
    ap.add_argument("--taus", default=",".join(TAUS), help="comma-separated (default: 0.5, 0.7, 0.9)")
    ap.add_argument("--graphs", default=os.path.join(ROOT, "graphs"), help="edge lists <ABBR>.txt")
    ap.add_argument("--build", default=os.path.join(ROOT, "build"), help="CMake build directory")
    ap.add_argument("--work", default=os.path.join(ROOT, "results", "certificate"), help="output directory")
    ap.add_argument("--results", help="results file whose EBA rows are compared with the certificate")
    ap.add_argument("--workers", type=int, default=1, help="runs in parallel (default: 1)")
    a = ap.parse_args()

    datasets, taus = a.datasets.split(","), a.taus.split(",")
    binary = os.path.join(a.build, "flexi")
    checker = os.path.join(a.build, "indep_check")
    source = os.path.join(ROOT, "verification", "indep_check.cpp")
    if not os.path.exists(checker) or os.path.getmtime(checker) < os.path.getmtime(source):
        subprocess.run(["g++", "-O2", "-std=c++17", "-o", checker, source], check=True)
    os.makedirs(os.path.join(a.work, "fpa"), exist_ok=True)
    graph = {d: os.path.join(a.graphs, f"{d}.txt") for d in datasets}
    sol = {(d, t): os.path.join(a.work, "fpa", f"{d}_{t}.txt") for d in datasets for t in taus}

    with ThreadPoolExecutor(a.workers) as pool:
        fpa = dict(zip(sol, pool.map(lambda k: run_fpa(binary, graph[k[0]], k[1], sol[k]), sol)))
        checks = dict(zip(datasets, pool.map(
            lambda d: check_dataset(checker, graph[d], [(t, fpa[(d, t)][2], sol[(d, t)]) for t in taus],
                                    os.path.join(a.work, f"{d}.queries")), datasets)))
    eba, fpa_res = read_results(a.results) if a.results else ({}, {})

    rows, failures = [], []
    for d in datasets:
        g, maxcore, pairs = checks[d]
        n, m = int(g["n"]), int(g["m"])
        for t, (ck, sl) in zip(taus, pairs):
            fn, fm, f = fpa[(d, t)]
            same_graph = (fn, fm) == (n, m) and g["malformed"] == "0" and g["out_of_int"] == "0"
            witness = (int(sl["size"]) == f and f >= 2 and sl["bad"] == "0" and sl["missing"] == "0"
                       and sl["dup"] == "0" and sl["connected"] == "1" and floor_pow(f, t) <= int(sl["mindeg"]))
            certified = same_graph and witness and ck["over_f"] == "0"
            if not same_graph:
                failures.append(f"{d}: flexi reads n={fn} m={fm}, indep_check n={n} m={m} "
                                f"(malformed lines {g['malformed']}, ids outside int {g['out_of_int']})")
            if not witness:
                failures.append(f"{d} tau={t}: FPA solution is not a flexi-clique of size {f} ({sl})")
            structure = ""
            if certified:
                if sl["mindeg_core_comp"] == sl["size"]:
                    structure = "whole_core_component"      # a whole component of the k-core, k = min degree
                elif f == size_cap(floor_pow(f, t), t):
                    structure = "size_cap"                  # largest size whose threshold is floor(f^tau)
                else:
                    structure = "other"
            e = eba.get((d, Fraction(t)), {})
            rows.append(dict(dataset=d, tau=t, fpa_size=f, **{"lambda": ck["lam"]},
                             lambda_core_nodes=ck["nodes"], lambda_core_components=ck["comps"],
                             largest_components=ck["top"], graph_check="ok" if same_graph else "MISMATCH",
                             witness_valid="yes" if witness else "NO", certified="yes" if certified else "no",
                             structure=structure, eba_status=e.get("status", ""), eba_size=e.get("size", ""),
                             eba_branches=e.get("branches", "")))

    out = os.path.join(a.work, "certificate.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    cert = [(r["dataset"], r["tau"]) for r in rows if r["certified"] == "yes"]
    print(f"{len(rows)} instances; FPA solutions valid: {sum(r['witness_valid'] == 'yes' for r in rows)}; "
          f"certified optimal: {len(cert)}")
    for s in ("whole_core_component", "size_cap", "other"):
        print(f"  structure {s}: {sum(r['structure'] == s for r in rows)}")
    if eba:
        root = [(r["dataset"], r["tau"]) for r in rows if r["eba_status"] == "ok" and r["eba_branches"] == "0"]
        solved = [(r["dataset"], r["tau"]) for r in rows if r["eba_status"] == "ok"]
        print(f"EBA: solved {len(solved)}, without branching {len(root)}; "
              f"same instances as the certificate: {set(root) == set(cert)}")
        # A certified instance on which EBA branched, or an uncertified one on which EBA did not
        # branch, contradicts the correspondence; an instance EBA did not solve is only reported.
        for d, t in sorted(set(cert) - set(root)):
            if (d, t) in solved:
                failures.append(f"{d} tau={t}: certified, but EBA branched")
            else:
                print(f"  certified, EBA did not finish: {d} tau={t}")
        for d, t in sorted(set(root) - set(cert)):
            failures.append(f"{d} tau={t}: EBA did not branch, but the certificate does not hold")
        for r in rows:
            if r["eba_status"] == "ok" and int(r["eba_size"]) < r["fpa_size"]:
                failures.append(f"{r['dataset']} tau={r['tau']}: FPA size {r['fpa_size']} above the EBA optimum {r['eba_size']}")
            k = (r["dataset"], Fraction(r["tau"]))
            if k in fpa_res and str(r["fpa_size"]) != fpa_res[k]["size"]:
                failures.append(f"{r['dataset']} tau={r['tau']}: FPA size {r['fpa_size']} differs from the results file ({fpa_res[k]['size']})")
        print(f"  FPA below the EBA optimum: {[(r['dataset'], r['tau'], r['fpa_size'], r['eba_size']) for r in rows if r['eba_status'] == 'ok' and str(r['fpa_size']) != r['eba_size']]}")
    for x in failures:
        print("FAILED:", x)
    print(f"written to {out}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
