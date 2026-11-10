#!/usr/bin/env python3
"""Tables of the paper and of the technical report, and the numbers quoted in their text.

Input: the results file of experiments/run_experiments.py.

  python3 analysis/tables.py                                  # results/results.csv -> results/tables/
  python3 analysis/tables.py my_results.csv --out my_tables

Writes to the output directory:
  tables.tex    one LaTeX float per table, each preceded by a line "%%% TABLE <label>":
                tab:efficiency (EQ1, Table III), tab:quality (technical report: sizes of NPA and FPA),
                and tab:width (technical report: branching width with and without the pivot);
                they use the macros \\EBA, \\FPA, \\NPA, \\TO and \\MO of the paper
  <label>.csv   the cells of each table as plain text
  stats.txt     the numbers quoted in the text of EQ1-EQ5 (also printed)

Cells: T/O = time limit, M/O = memory limit, -- = not run. Running times are in seconds with three
decimals below 1 s and four significant digits above, as in the paper. A table lists the datasets
that have rows for its algorithms; when a run has several rows (a failed run that was repeated), the
last row counts. Every table and statistic uses the rows of its own group.

Protocol of the paper: one thread per run, a one-day time limit, a 400 GB memory limit, and the
median of three runs for runs under 60 s, measured on an otherwise idle machine; a file written by
run_experiments.py with --workers 1 on an idle machine comes closest.
"""
import argparse, csv, math, os, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "experiments"))
import run_experiments as rx                                  # noqa: E402  (groups, dataset lists)

REAL = rx.REAL
TIME_FLOOR = 0.01                                             # s; ratios of sub-10 ms times are noise
VARIANT = {"eba": r"{\EBA}", "fpa": r"{\FPA}", "npa": r"{\NPA}",
           "eba_no_rule1": r"{\EBA}$-$Rule 1", "eba_no_rule2": r"{\EBA}$-$Rule 2", "eba_no_rule3": r"{\EBA}$-$Rule 3",
           "eba_no_incumbent": r"{\EBA}$-${\FPA}", "eba_no_pivot": r"{\EBA}$-$Pivot",
           "eba_no_degeneracy": r"{\EBA}$-$Degeneracy", "binary": "Binary", "flat_symse": "Flat Sym-SE",
           "cpsat": "CP-SAT"}


def tau_key(t):
    return f"{float(t):.2f}"


def load(path):
    """last[(group, dataset, tau, algo)] = last row (status ok / TO / MO) of the results file."""
    last = {}
    for r in csv.DictReader(open(path)):
        if r["status"] in ("ok", "TO", "MO"):
            last[(r["group"], r["dataset"], tau_key(r["tau"]), r["algo"])] = r
    return last


class Results:
    def __init__(self, path):
        self.last = load(path)

    def row(self, group, d, t, algo):
        return self.last.get((group, d, tau_key(t), algo))

    def ok(self, group, d, t, algo):
        r = self.row(group, d, t, algo)
        return r if r is not None and r["status"] == "ok" else None

    def groups_of(self, group):
        """Datasets of a group that have at least one row."""
        return [d for d in rx.GROUPS[group]["datasets"] or [] if any(k[0] == group and k[1] == d for k in self.last)]


def fmt_time(x, sig=4):
    """Running time as printed in the paper: < 1 s with three decimals, otherwise `sig`
    significant digits; thousands separated by {,} (math mode)."""
    if x < 1:
        s = f"{x:.3f}"
        if float(s) < 1:
            return s
    digits = len(str(int(x)))
    if digits < sig and len(str(int(round(x, sig - digits)))) > digits:
        digits += 1                                           # rounding carries into a new digit
    if digits >= sig:
        return f"{int(round(x)):,}".replace(",", "{,}")
    s = f"{x:.{sig - digits}f}"
    ip, _, fp = s.partition(".")
    return f"{int(ip):,}".replace(",", "{,}") + ("." + fp if fp else "")


def cell_time(r):
    if r is None:
        return "--"
    if r["status"] == "ok":
        return fmt_time(max(float(r["time"]), 0.001))
    return {"TO": r"\TO", "MO": r"\MO"}[r["status"]]


def plain(s):
    return s.replace("{,}", ",").replace(r"\TO", "T/O").replace(r"\MO", "M/O")


def gm(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float("nan")


def time_ratio(r, ref):
    return max(float(r["time"]), TIME_FLOOR) / max(float(ref["time"]), TIME_FLOOR)


# ---------------------------------------------------------------- tables
def tab_efficiency(R):
    """EQ1: running times of NPA, FPA and EBA (rows) on the real-world graphs (columns)."""
    cols = [d for d in REAL if any(R.row("main", d, t, a) for t in rx.TAU_EQ1 for a in ("npa", "fpa", "eba"))]
    if not cols:
        return None, None
    tex = [r"\begin{table*}[t]", r"\centering", r"\caption{EQ1. Efficiency (running time: second).}",
           r"\label{tab:efficiency}", r"\small",
           r"\begin{tabular}{c|c||" + "|".join(["r"] * len(cols)) + "}", r"\hline",
           r"\textbf{Algorithm} & \textbf{$\tau$} & " + " & ".join(rf"\textbf{{{d}}}" for d in cols) + r" \\",
           r"\hline\hline"]
    rows = [["Algorithm", "tau"] + cols]
    for alg in ("npa", "fpa", "eba"):
        for i, t in enumerate(rx.TAU_EQ1):
            vals = [cell_time(R.row("main", d, t, alg)) for d in cols]
            head = rf"\multirow{{3}}{{*}}{{{VARIANT[alg]}}}" if i == 0 else ""
            tex.append(f"{head} & ${t:.1f}$ & " + " & ".join(v if v.startswith(("--", "\\")) else f"${v}$" for v in vals)
                       + (r" \\ \hline" if i == 2 else r" \\ \cline{2-" + str(len(cols) + 2) + "}"))
            rows.append([alg.upper() if i == 0 else "", f"{t:.1f}"] + [plain(v) for v in vals])
    tex += [r"\end{tabular}", r"\end{table*}"]
    return tex, rows


def tab_quality(R):
    """Technical report: solution sizes of NPA and FPA (group main); bold = larger, underlined FPA = tie."""
    cols = [d for d in REAL if all(R.ok("main", d, t, a) for t in rx.TAU_EQ1 for a in ("npa", "fpa"))]
    if not cols:
        return None, None
    size = {(d, t, a): int(R.ok("main", d, t, a)["size"]) for d in cols for t in rx.TAU_EQ1 for a in ("npa", "fpa")}
    for d in cols:
        for t in rx.TAU_EQ1:
            e = R.ok("main", d, t, "eba")
            if e is not None and max(size[(d, t, "npa")], size[(d, t, "fpa")]) > int(e["size"]):
                sys.exit(f"{d} {t}: a heuristic solution is larger than the EBA optimum")
    tex = [r"\begin{table}[h]", r"\centering", r"\small",
           r"\caption{Solution sizes of {\NPA} and {\FPA} on the datasets of the paper at $\tau\in\{0.5,0.7,0.9\}$. "
           r"Bold marks the larger of the two where they differ; an underlined {\FPA} value equals the {\NPA} value.}",
           r"\label{tab:quality}", r"\begin{tabular}{ll" + "r" * len(cols) + "}", r"\toprule",
           r"Alg. & $\tau$ & " + " & ".join(cols) + r" \\", r"\midrule"]
    rows = [["Alg.", "tau"] + cols]
    for ai, alg in enumerate(("npa", "fpa")):
        other = "fpa" if alg == "npa" else "npa"
        for ti, t in enumerate(rx.TAU_EQ1):
            cells = []
            for d in cols:
                mine, oth = size[(d, t, alg)], size[(d, t, other)]
                s = f"{mine:,}"
                if mine > oth:
                    s = r"\textbf{" + s + "}"
                elif mine == oth and alg == "fpa":
                    s = r"\underline{" + s + "}"
                cells.append(s)
            first = rf"\multirow{{3}}{{*}}{{{VARIANT[alg]}}}" if ti == 0 else ""
            tex.append(f"{first} & {t} & " + " & ".join(cells) + r" \\")
            rows.append([alg.upper() if ti == 0 else "", f"{t}"] + [f"{size[(d, t, alg)]:,}" for d in cols])
        if ai == 0:
            tex.append(r"\midrule")
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return tex, rows


def width_fields(r):
    kv = dict(x.split("=", 1) for x in r["extra"].split(";") if "=" in x)
    return int(kv["internal"]), float(kv["mean_w"]), int(kv["max_w"]), int(kv["capped"])


def branching_instances(R):
    """(dataset, tau) of group branching where EBA branches, in dataset and tau order."""
    out = []
    for d in rx.GROUPS["branching"]["datasets"]:
        for t in rx.GROUPS["branching"]["taus"]:
            e = R.ok("branching", d, t, "eba")
            if e is not None and int(e["branches"]) > 0:
                out.append((d, tau_key(t)))
    return out


def tab_width(R):
    """Technical report: nodes, mean and maximum width, and time of EBA and EBA-Pivot where EBA branches."""
    insts = [(d, t) for d, t in branching_instances(R)
             if R.ok("width", d, t, "width") and R.ok("width", d, t, "width_no_pivot") and R.ok("branching", d, t, "eba_no_pivot")]
    if not insts:
        return None, None
    tex = [r"\begin{table}[h]", r"\centering", r"\small",
           r"\caption{Effect of the margin pivot on the instances of the EQ5 datasets on which {\EBA} branches "
           r"(the other instances are certified at the root and are identical under both variants): optimum "
           r"$|F^*|$, number of search-tree nodes, mean and maximum number of children explored per internal "
           r"node (branching width), and running time, for {\EBA} and for {\EBA} without the pivot.}",
           r"\label{tab:width}", r"\begin{tabular}{llr|rrrr|rrrr}", r"\toprule",
           r" & & & \multicolumn{4}{c|}{{\EBA}} & \multicolumn{4}{c}{{\EBA}$-$Pivot} \\",
           r"Dataset & $\tau$ & $|F^*|$ & Nodes & Mean w. & Max w. & Time (s) & Nodes & Mean w. & Max w. & Time (s) \\",
           r"\midrule"]
    rows = [["Dataset", "tau", "F*", "EBA nodes", "mean w", "max w", "time", "-Pivot nodes", "mean w", "max w", "time"]]
    prev, cap = None, False
    for d, t in insts:
        e, p = R.ok("branching", d, t, "eba"), R.ok("branching", d, t, "eba_no_pivot")
        we, wp = R.ok("width", d, t, "width"), R.ok("width", d, t, "width_no_pivot")
        if len({e["size"], p["size"], we["size"], wp["size"]}) != 1:
            sys.exit(f"{d} {t}: sizes differ between the runs")
        if e["branches"] != we["branches"] or p["branches"] != wp["branches"]:
            sys.exit(f"{d} {t}: the width trace and the timed run differ in the number of search nodes")
        def w(r):
            internal, mean, mx, capped = width_fields(r)
            return ("--", "--") if internal == 0 else (f"{mean:.1f}", f"{mx:,}")
        cap |= any(width_fields(r)[3] > 0 for r in (we, wp))
        cells = [d if d != prev else "", t, f"{int(e['size']):,}", f"{int(e['branches']):,}", *w(we),
                 plain(fmt_time(float(e["time"]))), f"{int(p['branches']):,}", *w(wp), plain(fmt_time(float(p["time"])))]
        prev = d
        tex.append(" & ".join(cells) + r" \\")
        rows.append(cells)
    tex += [r"\bottomrule", r"\end{tabular}"]
    if cap:
        tex.append(r"\par\smallskip\footnotesize Widths of $255$ and above are counted as $255$.")
    tex.append(r"\end{table}")
    return tex, rows


# ---------------------------------------------------------------- numbers quoted in the text
def stats(R):
    out = []
    say = out.append

    # EQ1
    say("## EQ1 (group main)")
    for d in REAL:
        e = [R.row("main", d, t, "eba") for t in rx.TAU_EQ1]
        if any(e):
            say(f"  {d}: EBA " + ", ".join(f"{t}: {plain(cell_time(r))}" for t, r in zip(rx.TAU_EQ1, e)))
    faster = sum(float(R.ok("main", d, t, "npa")["time"]) < float(R.ok("main", d, t, "fpa")["time"])
                 for d in REAL for t in rx.TAU_EQ1 if R.ok("main", d, t, "npa") and R.ok("main", d, t, "fpa"))
    both = sum(1 for d in REAL for t in rx.TAU_EQ1 if R.ok("main", d, t, "npa") and R.ok("main", d, t, "fpa"))
    say(f"  NPA faster than FPA on {faster} of {both} instances where both finished")

    # EQ2
    say("## EQ2 (group synthetic): remaining nodes after Rules 1 and 3 vs EBA time")
    pts = []
    for d in rx.LFR:
        for t in rx.TAU_LFR:
            s, e = R.ok("synthetic", d, t, "root_stats"), R.ok("synthetic", d, t, "eba")
            if s and e:
                kv = dict(x.split("=") for x in s["extra"].split(";"))
                n = int(kv["V13"])
                if n > 0:
                    pts.append((n, float(e["time"])))
    if len(pts) >= 3:
        x = [math.log10(n) for n, _ in pts]; y = [math.log10(t) for _, t in pts]
        mx, my = sum(x) / len(x), sum(y) / len(y)
        b = sum((p - mx) * (q - my) for p, q in zip(x, y)) / sum((p - mx) ** 2 for p in x)
        a0 = my - b * mx
        r2 = 1 - sum((q - a0 - b * p) ** 2 for p, q in zip(x, y)) / sum((q - my) ** 2 for q in y)
        within = sum(abs(q - a0 - b * p) <= 0.5 for p, q in zip(x, y))
        say(f"  {len(pts)} instances; fit log10(time) = {a0:.3f} + {b:.3f} log10(nodes), R^2 = {r2:.3f}; "
            f"{within} within half a decade of the fit; remaining nodes {min(n for n, _ in pts):,}-{max(n for n, _ in pts):,}")
    else:
        say("  not enough instances")

    # EQ3
    say("## EQ3 (group effectiveness): FPA size / optimum")
    ratios, misses = [], []
    for d in rx.GROUPS["effectiveness"]["datasets"]:
        for t in rx.TAU_FIG:
            f, e = R.ok("effectiveness", d, t, "fpa"), R.ok("effectiveness", d, t, "eba")
            if f and e:
                q = int(f["size"]) / int(e["size"])
                if q > 1:
                    sys.exit(f"{d} {t}: FPA larger than the optimum")
                ratios.append((q, d, t))
                if q < 1:
                    misses.append(q)
    if ratios:
        say(f"  {len(ratios)} instances: mean ratio {sum(q for q, *_ in ratios) / len(ratios):.3f}, "
            f"minimum {min(ratios)[0]:.3f} ({min(ratios)[1]} {min(ratios)[2]}), FPA misses the optimum on "
            f"{len(misses)} (mean ratio there {sum(misses) / len(misses) if misses else float('nan'):.3f})")

    # EQ4
    say("## EQ4 (group ablation): variant time / EBA time")
    for v in ("eba_no_rule1", "eba_no_rule2", "eba_no_rule3", "eba_no_incumbent"):
        rat, to, per = [], [], defaultdict(list)
        for d in rx.GROUPS["ablation"]["datasets"]:
            for t in rx.TAU_FIG:
                e, r = R.ok("ablation", d, t, "eba"), R.row("ablation", d, t, v)
                if e is None or r is None:
                    continue
                if r["status"] != "ok":
                    to.append(f"{d} {t}")
                    continue
                q = time_ratio(r, e)
                rat.append((q, d, t, float(r["time"]), float(e["time"])))
                per[d].append(q)
        if rat:
            m = max(rat)
            say(f"  {v}: max slowdown {m[0]:,.0f}x ({m[1]} {m[2]}: {m[3]:.4g} s vs {m[4]:.4g} s), "
                f"failed ({len(to)}): {', '.join(to) or 'none'}; per dataset max: "
                + ", ".join(f"{d} {max(q):,.0f}x" for d, q in per.items()))

    # EQ5
    say("## EQ5 (group branching): time ratios versus EBA")
    insts = [(d, tau_key(t)) for d in rx.GROUPS["branching"]["datasets"] for t in rx.TAU_FIG
             if R.ok("branching", d, t, "eba")]
    br = set(branching_instances(R))
    say(f"  {len(insts)} instances with an EBA run, {len(br)} where EBA branches")
    for v in ("eba_no_pivot", "eba_no_degeneracy", "binary", "flat_symse", "cpsat"):
        root, bran, fail, nodes, faster, faster25 = [], [], [], [], 0, 0
        for d, t in insts:
            e, r = R.ok("branching", d, t, "eba"), R.row("branching", d, t, v)
            if r is None:
                continue
            if r["status"] != "ok":
                fail.append(f"{d} {t}")
                continue
            q = time_ratio(r, e)
            if (d, t) in br:
                bran.append(q)
                faster += q < 1
                faster25 += q < 0.8                       # beyond the noise of sub-second runs
                if v != "cpsat":
                    nodes.append(max(int(r["branches"]), 1) / max(int(e["branches"]), 1))
            else:
                root.append(q)
        say(f"  {v}: root-solved instances max {max(root) if root else float('nan'):.2f}x; branching instances "
            f"geomean {gm(bran):.2f}x, max {max(bran) if bran else float('nan'):,.0f}x, faster than EBA on {faster} of {len(bran)} "
            f"(by more than 25%: {faster25})"
            + (f", search nodes geomean {gm(nodes):.2f}x" if nodes else "")
            + f"; failed ({len(fail)}): {', '.join(fail) or 'none'}")
    narrower = equal = wider = 0
    for d, t in branching_instances(R):
        we, wp = R.ok("width", d, t, "width"), R.ok("width", d, t, "width_no_pivot")
        if we and wp and width_fields(we)[0] and width_fields(wp)[0]:
            a, b = width_fields(we)[2], width_fields(wp)[2]
            narrower += a < b; equal += a == b; wider += a > b
    say(f"  max width (group width): EBA narrower {narrower} / equal {equal} / wider {wider} than EBA-Pivot")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("results", nargs="?", default=os.path.join(ROOT, "results", "results.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "tables"))
    a = ap.parse_args()
    if not os.path.exists(a.results):
        sys.exit(f"{a.results} not found")
    R = Results(a.results)
    os.makedirs(a.out, exist_ok=True)
    tex = []
    for label, fn in (("tab:efficiency", tab_efficiency), ("tab:quality", tab_quality), ("tab:width", tab_width)):
        t, rows = fn(R)
        if t is None:
            print(f"{label}: no data")
            continue
        tex += [f"%%% TABLE {label}"] + t + [""]
        with open(os.path.join(a.out, label.split(":")[1] + ".csv"), "w", newline="") as f:
            csv.writer(f).writerows(rows)
        print(f"{label}: {len(rows) - 1} rows")
    with open(os.path.join(a.out, "tables.tex"), "w") as f:
        f.write("\n".join(tex))
    s = "\n".join(stats(R))
    with open(os.path.join(a.out, "stats.txt"), "w") as f:
        f.write(s + "\n")
    print(s)


if __name__ == "__main__":
    main()
