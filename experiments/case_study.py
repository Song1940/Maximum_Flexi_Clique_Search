#!/usr/bin/env python3

import argparse, csv, json, math, os, random, subprocess, sys, threading, time
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLE_SEED = 2027
MIN_SIZE, MAX_SIZE = 12, 50
BASELINES = [("clique", 1), ("kplex", 3), ("quasi", 0.5)]
MODEL_ORDER = ["clique", "kplex", "quasi", "flexi"]          # row order of Table tab:case
FIELDS = ["cid", "query", "csize", "model", "param", "pred_size", "precision", "recall", "f1", "jaccard",
          "timeout", "time", "status", "region_size", "local_size", "local_time", "local_status"]


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
    """One community per line, node ids separated by whitespace or commas."""
    comms = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and line[0] != "#":
                comms.append([int(x) for x in line.replace(",", " ").split()])
    return comms


def query_of(adj, members):
    """Member with the most neighbours in the community; ties -> smallest id."""
    C = set(members)
    return max(members, key=lambda v: (len(adj[v] & C), -v))


def select_communities(adj, comms, n, query=None):
    """[(cid, sorted members)]: the seeded sample of n distinct communities of 12..50 nodes (cid = position
    in the sample), or, with query, every such community whose query node is query (cid = index in the
    sorted list of distinct communities)."""
    pool = sorted({tuple(sorted(x for x in c if x in adj)) for c in comms})
    pool = [c for c in pool if MIN_SIZE <= len(c) <= MAX_SIZE]
    if query is not None:
        return [(i, c) for i, c in enumerate(pool) if query in c and query_of(adj, c) == query]
    if n > len(pool):
        sys.exit(f"--n {n}: only {len(pool)} distinct communities of {MIN_SIZE}-{MAX_SIZE} nodes")
    return list(enumerate(random.Random(SAMPLE_SEED).sample(pool, n)))


def ball(adj, q, hops=2):
    """Sorted node ids within `hops` hops of q."""
    seen, dq = {q}, deque([(q, 0)])
    while dq:
        v, d = dq.popleft()
        if d < hops:
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    dq.append((w, d + 1))
    return sorted(seen)


def connected(adj, P):
    P = set(P)
    start = next(iter(P))
    seen, st = {start}, [start]
    while st:
        for w in adj[st.pop()]:
            if w in P and w not in seen:
                seen.add(w)
                st.append(w)
    return len(seen) == len(P)


def is_model(adj, P, model, param):
    """Exact feasibility of the non-empty set P under the model, on the whole graph."""
    s, P = len(P), set(P)
    d = min(len(adj[v] & P) for v in P)
    if model == "flexi":                       # floor(s^tau) <= d  <=>  s^a < (d+1)^b for tau = a/b
        t = Fraction(str(param))
        return connected(adj, P) and s ** t.numerator < (d + 1) ** t.denominator
    if model == "clique":
        return d == s - 1
    if model == "kplex":
        return d >= s - param
    return d >= math.ceil(Fraction(str(param)) * (s - 1))      # degree-based quasi-clique


def density(adj, P):
    s = len(P)
    if s < 2:
        return None
    P = set(P)
    return sum(len(adj[v] & P) for v in P) / (s * (s - 1))


def prf(pred, truth):
    """Precision, recall, F1 and Jaccard similarity of the answer pred against the community truth."""
    inter = len(pred & truth)
    p = inter / len(pred) if pred else 0.0
    r = inter / len(truth)
    f = 2 * p * r / (p + r) if p + r else 0.0
    j = inter / len(pred | truth) if pred | truth else 0.0
    return p, r, f, j


def write_region(path, adj, R):
    """Edge list of the subgraph induced by R."""
    rs = set(R)
    with open(path, "w") as f:
        for u in R:
            for v in adj[u]:
                if v in rs and u < v:
                    f.write(f"{u} {v}\n")


def run_flexi(binary, graph_file, q, tau, timeout, init=None):
    """(nodes, timed_out, seconds, status) of the maximum flexi-clique containing q."""
    cmd = [binary, "--algo", "eba", "--file", graph_file, "--tau", str(tau), "--seed", str(q)]
    if init:
        cmd += ["--init", ",".join(map(str, sorted(init)))]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return set(), 1, time.time() - t0, "TO"
    except OSError as e:
        return set(), 0, time.time() - t0, f"ERR:{type(e).__name__}"
    kv = {k.strip(): v.strip() for k, v in (l.split(":", 1) for l in r.stdout.splitlines() if ":" in l)}
    if r.returncode != 0 or kv.get("Valid") != "yes" or "Nodes" not in kv:
        return set(), 0, time.time() - t0, f"ERR:{r.returncode}"
    try:
        return {int(x) for x in kv["Nodes"].split(",") if x}, 0, time.time() - t0, "ok"
    except ValueError:
        return set(), 0, time.time() - t0, "ERR:parse"


def run_qcs(binary, graph_file, q, model, param, region, timeout, workdir, tag):
    """(nodes, timed_out, seconds, status) of the maximum clique / k-plex / quasi-clique containing q
    within region (passed to qcs in this order)."""
    qf, of = os.path.join(workdir, f"q_{tag}.txt"), os.path.join(workdir, f"o_{tag}.txt")
    qmodel, qparam = ("kplex", 1) if model == "clique" else (model, param)
    t0 = time.time()
    try:
        with open(qf, "w") as f:
            f.write(f"{q} {qmodel} {qparam} {','.join(map(str, region))}\n")
        subprocess.run([binary, "--file", graph_file, "--queries", qf, "--out", of, "--timeout", str(timeout)],
                       check=True, capture_output=True, timeout=timeout + 600)
        with open(of) as f:
            kv = dict(tok.split("=", 1) for tok in f.readline().split() if "=" in tok)
        if kv.get("valid") != "1" or "nodes" not in kv:
            return set(), 0, time.time() - t0, "ERR:output"
        to = int(kv.get("timeout", "0"))
        return ({int(x) for x in kv["nodes"].split(",") if x}, to, float(kv.get("time", "nan")),
                "TO" if to else "ok")
    except subprocess.TimeoutExpired:
        return set(), 1, time.time() - t0, "TO"
    except subprocess.CalledProcessError as e:
        return set(), 0, time.time() - t0, f"ERR:{e.returncode}"
    except (OSError, ValueError) as e:
        return set(), 0, time.time() - t0, f"ERR:{type(e).__name__}"
    finally:
        for p in (qf, of):
            if os.path.exists(p):
                os.remove(p)


def key(r):
    return int(r["cid"]), r["model"], float(r["param"])


def read_runs(out):
    """Complete runs of an output directory: the rows of raw.csv that also have their line in nodes.jsonl
    (a run interrupted while writing leaves at most a partial last line)."""
    rows, recs = [], {}
    if os.path.exists(os.path.join(out, "nodes.jsonl")):
        with open(os.path.join(out, "nodes.jsonl")) as f:
            for line in f:
                try:
                    d = json.loads(line)
                    recs[key(d)] = d
                except (ValueError, KeyError, TypeError):
                    pass
    if os.path.exists(os.path.join(out, "raw.csv")):
        with open(os.path.join(out, "raw.csv"), newline="") as f:
            for r in csv.DictReader(f):
                try:
                    if key(r) in recs and all(r[k] is not None for k in FIELDS):
                        float(r["jaccard"])
                        rows.append(r)
                except (ValueError, KeyError, TypeError):
                    pass
    return rows, [recs[key(r)] for r in rows]


def summary(out, graph, first=None, adj=None):
    """summary.csv: mean precision, recall, F1, Jaccard, answer size and density per model, and the mean
    size and density of the communities, over the runs of raw.csv (with first, over cids < first)."""
    rows, recs = read_runs(out)
    if first is not None:
        keep = [int(r["cid"]) < first for r in rows]
        rows = [r for r, k in zip(rows, keep) if k]
        recs = [d for d, k in zip(recs, keep) if k]
    if not rows:
        sys.exit(f"no complete runs in {out}")
    adj = adj or load_graph(graph)
    groups = {}
    for r, d in zip(rows, recs):
        groups.setdefault((r["model"], float(r["param"])), []).append((r, d))
    mean = lambda xs: sum(xs) / len(xs) if xs else float("nan")
    table = []
    for (m, p), rs in sorted(groups.items(), key=lambda g: (MODEL_ORDER.index(g[0][0]), g[0][1])):
        pd = [x for x in (density(adj, d["pred"]) for _, d in rs) if x is not None]
        table.append(dict(
            model=m, param=rs[0][0]["param"], n=len(rs),
            **{k: f"{mean([float(r[k]) for r, _ in rs]):.3f}" for k in ("precision", "recall", "f1", "jaccard")},
            pred_size=f"{mean([int(r['pred_size']) for r, _ in rs]):.1f}", pred_density=f"{mean(pd):.3f}",
            csize=f"{mean([len(d['community']) for _, d in rs]):.1f}",
            cdensity=f"{mean([density(adj, d['community']) for _, d in rs]):.3f}",
            timeouts=sum(r["status"] == "TO" for r, _ in rs),
            errors=sum(r["status"].startswith("ERR") for r, _ in rs),
            flagged=sum(r["status"] == "ok<4" for r, _ in rs)))
    path = os.path.join(out, "summary.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    print(f"{'model':<8}{'param':>6}{'n':>5}{'Prec':>8}{'Rec':>7}{'F1':>7}{'Jac':>7}{'|K|':>7}{'dens':>7}"
          f"{'|C|':>7}{'dens':>7}{'TO':>5}{'ERR':>5}")
    for t in table:
        print(f"{t['model']:<8}{t['param']:>6}{t['n']:>5}{t['precision']:>8}{t['recall']:>7}{t['f1']:>7}"
              f"{t['jaccard']:>7}{t['pred_size']:>7}{t['pred_density']:>7}{t['csize']:>7}{t['cdensity']:>7}"
              f"{t['timeouts']:>5}{t['errors']:>5}")
    print("wrote", path)


def run(a):
    for path, what in ((a.flexi, "flexi binary"), (a.qcs, "qcs binary"), (a.graph, "graph"),
                       (a.communities, "communities")):
        if not os.path.exists(path):
            sys.exit(f"{what} not found: {path} (build with cmake, see README.md; "
                     f"python3 datasets/prepare_datasets.py AM --communities AM writes the graph and the communities)")
    taus = list(dict.fromkeys(a.tau))
    if any(not 0 < t < 1 or Fraction(str(t)).denominator > 100 for t in taus):
        sys.exit("--tau: give values in (0,1) with at most two decimals")
    os.makedirs(a.out, exist_ok=True)
    raw_path, nodes_path = os.path.join(a.out, "raw.csv"), os.path.join(a.out, "nodes.jsonl")
    done_rows, done_recs = read_runs(a.out)
    if (os.path.exists(raw_path) or os.path.exists(nodes_path)) and not a.resume:
        sys.exit(f"{a.out} already has results; use --resume to continue them, or another --out")

    adj = load_graph(a.graph)
    comms = select_communities(adj, load_communities(a.communities), a.n, a.query)
    if not comms:
        sys.exit("no community selected" + (f": none of {MIN_SIZE}-{MAX_SIZE} nodes has query node {a.query}"
                                            if a.query is not None else ""))
    print(f"{len(comms)} communities, graph with {len(adj)} nodes", flush=True)
    queries = []
    for cid, members in comms:
        q = query_of(adj, members)
        queries.append((cid, q, set(members), ball(adj, q)))
    # flexi first: its whole-graph runs are the longest
    jobs = [(*x, "flexi", t) for t in taus for x in queries] + [(*x, m, p) for x in queries for m, p in BASELINES]
    done = {key(r) for r in done_rows}
    jobs = [j for j in jobs if (j[0], j[4], float(j[5])) not in done]
    print(f"{len(jobs)} runs ({len(done)} already done), time limit {a.timeout:g} s, {a.workers} in parallel",
          flush=True)

    # rewrite the complete runs (drops a partial last line of an interrupted run), then append
    fout = open(raw_path, "w", newline="")
    w = csv.DictWriter(fout, fieldnames=FIELDS, extrasaction="ignore")
    w.writeheader()
    w.writerows(done_rows)
    fnodes = open(nodes_path, "w")
    for d in done_recs:
        fnodes.write(json.dumps(d) + "\n")
    fout.flush()
    fnodes.flush()
    work_dir = os.path.join(a.out, "work")
    os.makedirs(work_dir, exist_ok=True)
    lock = threading.Lock()

    def work(job):
        cid, q, C, R, model, param = job
        local, lt, lst = None, "", ""
        if model == "flexi":
            rf = os.path.join(work_dir, f"{cid}_{q}_{param}.txt")
            write_region(rf, adj, R)
            local, _, lt, lst = run_flexi(a.flexi, rf, q, param, a.timeout)
            os.remove(rf)
            if lst == "ok" and (q not in local or not is_model(adj, local, model, param)):
                lst = "ERR:invalid"
            if lst != "ok":
                local = None
            P, to, t, st = run_flexi(a.flexi, a.graph, q, param, a.timeout, init=local)
            if st == "ok" and (q not in P or not is_model(adj, P, model, param)):
                P, st = set(), "ERR:invalid"
            if st == "ok" and local and len(P) < len(local):
                P, st = set(), "ERR:shrunk"
            if st != "ok" and local:           # best-known answer: the 2-hop answer
                P = local
        else:
            P, to, t, st = run_qcs(a.qcs, a.graph, q, model, param, R, a.timeout, work_dir, f"{cid}_{model}")
            if P and (q not in P or not is_model(adj, P, model, param)):
                P, st = set(), "ERR:invalid"
            if st == "ok" and model == "kplex" and len(P) < 4:
                st = "ok<4"
        p, r, f1, j = prf(P, C)
        row = dict(cid=cid, query=q, csize=len(C), model=model, param=param, pred_size=len(P),
                   precision=f"{p:.6f}", recall=f"{r:.6f}", f1=f"{f1:.6f}", jaccard=f"{j:.6f}", timeout=to,
                   time=f"{t:.3f}", status=st, region_size=len(R), local_size=len(local) if local else "",
                   local_time=f"{lt:.3f}" if lt != "" else "", local_status=lst)
        with lock:
            w.writerow(row)
            fout.flush()
            fnodes.write(json.dumps(dict(cid=cid, query=q, model=model, param=param, community=sorted(C),
                                         pred=sorted(P), local=sorted(local) if local else None)) + "\n")
            fnodes.flush()
            print(f"  {cid:4d} q={q} {model}({param}) |K|={len(P)} F1={f1:.3f} {st} {t:.1f}s", flush=True)

    with ThreadPoolExecutor(a.workers) as ex:
        list(ex.map(work, jobs))
    fout.close()
    fnodes.close()
    if not os.listdir(work_dir):
        os.rmdir(work_dir)
    summary(a.out, a.graph, adj=adj)


def main():
    ap = argparse.ArgumentParser(description="Case study (EQ7): community recovery on Amazon.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("run", help="sample the communities, run the four models, and summarise",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ps = sub.add_parser("summary", help="rebuild summary.csv from raw.csv and nodes.jsonl",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    for p in (pr, ps):
        p.add_argument("--out", default=os.path.join(ROOT, "results", "case_study"),
                       help="output directory")
        p.add_argument("--graph", default=os.path.join(ROOT, "datasets", "graphs", "AM.txt"), help="edge list")
    pr.add_argument("--communities", default=os.path.join(ROOT, "datasets", "communities", "AM.all.txt"),
                    help="ground-truth communities, one per line")
    pr.add_argument("--flexi", default=os.path.join(ROOT, "build", "flexi"), help="flexi binary")
    pr.add_argument("--qcs", default=os.path.join(ROOT, "build", "qcs"), help="qcs binary")
    pr.add_argument("--n", type=int, default=50, help="number of sampled communities")
    pr.add_argument("--tau", type=float, nargs="+", default=[0.5], help="tau of the flexi-clique")
    pr.add_argument("--timeout", type=float, default=3600, help="time limit of each solver run (seconds)")
    pr.add_argument("--workers", type=int, default=4, help="runs in parallel")
    pr.add_argument("--query", type=int, default=None,
                    help="instead of the sample, the communities whose query node is this node (EQ7 figure: 335413)")
    pr.add_argument("--resume", action="store_true", help="continue the runs already in --out")
    ps.add_argument("--first", type=int, default=None, help="only the communities with cid < FIRST")
    a = ap.parse_args()
    if a.cmd == "run":
        run(a)
    else:
        summary(a.out, a.graph, a.first)


if __name__ == "__main__":
    main()
