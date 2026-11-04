#!/usr/bin/env python3
"""
CP-SAT baseline of the branching-scheme comparison (EQ4): the maximum flexi-clique
with the OR-Tools CP-SAT solver.

Given the same initial incumbent F' as EBA (FPA, from the C++ binary) and the same
global degree pruning (theta = floor((|F'|+1)^tau)), the model searches for a
flexi-clique of size >= |F'|+1 inside the pruned graph G[V'], or proves none
exists (then F' is optimal).

The pruned graph is split into connected components, and each component larger
than the current incumbent is solved separately (largest first).

Model per component (OR-Tools CP-SAT):
  x_v in {0,1}          node v selected
  s = sum_v x_v         in [|F'|+1, |V'|]
  t = floor(s^tau)      via an element constraint over a precomputed table
  x_v = 1  =>  sum_{u in N(v)} x_u >= t                         (degree)
  sum_{u in N(v)} x_u >= t_min x_v                              (redundant)
  single-commodity flow from one root r over selected nodes      (connectivity):
    r = the lowest-index selected node (symmetry breaking via prefix ORs),
    0 <= g_v <= |V'| r_v, f_(u,v) = 0 unless x_u = x_v = 1,
    inflow(v) - outflow(v) = x_v - g_v.
Objective: maximise s.

Usage (requires the Python package ortools; the paper used ortools 9.15 with one worker):
  python3 baselines/cpsat/cpsat_baseline.py GRAPH TAU [--fpa build/flexi] [--time 86400] [--workers 1]
Prints one RESULT line: size, status (OPTIMAL, PROVED_INCUMBENT, PROVED_BY_PRUNING, or
TIMEOUT), validity, total time, algorithm time (graph loading excluded; the time
reported in the paper), CP-SAT solve time, |V'|, number of components, and |F'|.
"""
import argparse, math, os, re, subprocess, sys, tempfile, time
from collections import defaultdict


def load_graph(path):
    """Same re-indexing as Graph::loadFromFile (first appearance, self-loops skipped).
    Returns (n, adjacency sets, original id -> index)."""
    ids, edges = {}, set()
    def gid(x):
        if x not in ids: ids[x] = len(ids)
        return ids[x]
    # Like "ss >> u >> v" into C++ ints: two integers in the int32 range, each read after
    # optional white space.
    num = re.compile(r"\s*([+-]?\d+)")
    with open(path) as f:
        for line in f:
            if not line.strip() or line[0] == '#': continue
            m1 = num.match(line)
            m2 = num.match(line, m1.end()) if m1 else None
            if not m2: continue
            u, v = int(m1.group(1)), int(m2.group(1))
            if not (-2**31 <= u < 2**31 and -2**31 <= v < 2**31): continue
            if u == v: continue
            a, b = gid(u), gid(v)
            edges.add((min(a, b), max(a, b)))
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    return len(ids), adj, ids


def floor_pow(s, tau):          # = graph.h floorPow (shared threshold, tolerant to fp)
    return math.floor(s ** tau * (1 + 1e-12) + 1e-12) if s > 0 else 0


def compute_theta(best, tau):   # = graph.h computeTheta = floorPow(best+1)
    return floor_pow(best + 1, tau)


def is_flexi(adj, nodes, tau):
    S = set(nodes); k = len(S)
    if k < 2: return False
    thr = floor_pow(k, tau)
    if any(len(adj[v] & S) < thr for v in S): return False
    start = next(iter(S)); seen = {start}; st = [start]
    while st:
        x = st.pop()
        for y in adj[x] & S:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == k


def run_fpa(binary, path, tau, ids):
    """Run the C++ FPA; return (nodes as indices of load_graph, FPA time excluding graph loading)."""
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "fpa.out")
        subprocess.run([binary, "--algo", "fpa", "--file", path, "--tau", str(tau), "--out", out],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        nodes, el = [], 0.0
        for line in open(out):
            if line.startswith("nodes="):
                s = line.strip().split("=", 1)[1]
                nodes = [ids[int(x)] for x in s.split(",")] if s else []
            elif line.startswith("time="):
                el = float(line.strip().split("=", 1)[1])
    return nodes, el


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("graph"); ap.add_argument("tau", type=float)
    ap.add_argument("--fpa", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "build", "flexi"),
                    help="the flexi binary, used for the FPA incumbent")
    ap.add_argument("--time", type=float, default=86400.0, help="time limit in seconds")
    ap.add_argument("--workers", type=int, default=1, help="CP-SAT workers (threads)")
    a = ap.parse_args()
    if not (0.0 < a.tau < 1.0):
        print("ERROR: tau must be in (0,1)", file=sys.stderr); sys.exit(2)
    from ortools.sat.python import cp_model

    t_start = time.time()
    n, adj, ids = load_graph(a.graph)
    F, t_fpa = run_fpa(a.fpa, a.graph, a.tau, ids)
    t_after_fpa = time.time()
    if F and not is_flexi(adj, F, a.tau):
        print("ERROR: FPA solution invalid", file=sys.stderr); sys.exit(2)
    best = len(F)

    # Global degree pruning with theta = floor((|F'|+1)^tau) (same as EBA).
    theta = compute_theta(best, a.tau)
    alive = set(range(n)); deg = {v: len(adj[v]) for v in alive}
    stack = [v for v in alive if deg[v] < theta]
    removed = set(stack)
    while stack:
        v = stack.pop()
        for u in adj[v]:
            if u not in removed:
                deg[u] -= 1
                if deg[u] < theta: removed.add(u); stack.append(u)
    Vp = sorted(alive - removed)
    nv = len(Vp)
    t_pre = time.time() - t_start
    if nv <= best:
        algo = t_fpa + (time.time() - t_after_fpa)
        print(f"RESULT size={best} status=PROVED_BY_PRUNING valid={bool(F)} total={t_pre:.4f} algo={algo:.4f} solve=0 Vprime={nv} comps=0 incumbent={best}")
        return

    # Connected components of G[V'], largest first; only those larger than the
    # incumbent can contain a better solution.
    vset = set(Vp); seen = set(); comps = []
    for v in Vp:
        if v in seen: continue
        comp = [v]; seen.add(v); st = [v]
        while st:
            y = st.pop()
            for z in adj[y]:
                if z in vset and z not in seen: seen.add(z); comp.append(z); st.append(z)
        comps.append(sorted(comp))
    comps.sort(key=len, reverse=True)

    best_sol = list(F); status = "PROVED_INCUMBENT"; t_solve = 0.0
    for comp in comps:
        if len(comp) <= len(best_sol): break
        remaining = a.time - (time.time() - t_start)
        if remaining <= 1.0: status = "TIMEOUT"; break
        t1 = time.time()
        st_name, sol = solve_component(cp_model, adj, comp, len(best_sol), a.tau, remaining, a.workers)
        t_solve += time.time() - t1
        if sol is not None:
            if not is_flexi(adj, sol, a.tau):
                print("ERROR: CP-SAT solution invalid", file=sys.stderr); sys.exit(3)
            if len(sol) > len(best_sol): best_sol = sol
        if st_name == "OPTIMAL": status = "OPTIMAL"
        elif st_name != "INFEASIBLE": status = "TIMEOUT"; break
    total = time.time() - t_start
    algo = t_fpa + (time.time() - t_after_fpa)     # excludes graph loading
    valid = is_flexi(adj, best_sol, a.tau) if best_sol else False
    print(f"RESULT size={len(best_sol)} status={status} valid={valid} total={total:.4f} "
          f"algo={algo:.4f} solve={t_solve:.4f} Vprime={nv} comps={len(comps)} incumbent={best}")


def solve_component(cp_model, adj, comp, best, tau, time_limit, workers):
    """Search for a flexi-clique of size > best inside one component.
    Returns (status name, solution or None)."""
    nv = len(comp)
    idx = {v: i for i, v in enumerate(comp)}
    nb = [[idx[u] for u in adj[v] if u in idx] for v in comp]
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(nv)]
    lo = best + 1
    s = m.NewIntVar(lo, nv, "s")
    m.Add(sum(x) == s)
    table = [floor_pow(k, tau) for k in range(lo, nv + 1)]
    t = m.NewIntVar(min(table), max(table), "t")
    sidx = m.NewIntVar(0, nv - lo, "sidx")
    m.Add(sidx == s - lo)
    m.AddElement(sidx, table, t)
    tmin = min(table)
    for i in range(nv):
        m.Add(sum(x[j] for j in nb[i]) >= t).OnlyEnforceIf(x[i])
        m.Add(sum(x[j] for j in nb[i]) >= tmin * x[i])            # redundant
    # Root = lowest-index selected node: p_i = OR_{j<=i} x_j, r_i = p_i - p_{i-1}.
    p = [m.NewBoolVar(f"p{i}") for i in range(nv)]
    r = [m.NewBoolVar(f"r{i}") for i in range(nv)]
    for i in range(nv):
        if i == 0:
            m.Add(p[0] == x[0]); m.Add(r[0] == x[0])
        else:
            m.AddImplication(x[i], p[i])
            m.AddImplication(p[i - 1], p[i])
            m.AddBoolOr([p[i - 1], x[i]]).OnlyEnforceIf(p[i])
            m.Add(r[i] == p[i] - p[i - 1])
        m.AddImplication(r[i], x[i])
    m.AddExactlyOne(r)
    g = [m.NewIntVar(0, nv, f"g{i}") for i in range(nv)]
    inflow = [[] for _ in range(nv)]; outflow = [[] for _ in range(nv)]
    for i in range(nv):
        m.Add(g[i] <= nv * r[i])
        for j in nb[i]:
            f = m.NewIntVar(0, nv - 1, f"f{i}_{j}")
            m.Add(f == 0).OnlyEnforceIf(x[i].Not())
            m.Add(f == 0).OnlyEnforceIf(x[j].Not())
            outflow[i].append(f); inflow[j].append(f)
    for i in range(nv):
        m.Add(sum(inflow[i]) - sum(outflow[i]) == x[i] - g[i])
    m.Maximize(s)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = max(1.0, time_limit)
    solver.parameters.num_workers = workers
    st = solver.Solve(m)
    name = solver.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return name, [comp[i] for i in range(nv) if solver.Value(x[i])]
    return name, None


if __name__ == "__main__":
    main()
