#!/usr/bin/env python3
"""Figures of Section VI (EQ2-EQ5) from the results file of experiments/run_experiments.py.

  python3 analysis/figures.py                                   # all four -> results/figures/*.{pdf,png}
  python3 analysis/figures.py scalability
  python3 analysis/figures.py branching --datasets FB,GO,FW,CA,LU,SK
  python3 analysis/figures.py ablation --fine
  python3 analysis/figures.py all my_results.csv --out my_figures

  scalability    EQ2 (fig_scalability): EBA running time against the number of nodes left after Rules 1
                 and 3 on the five LFR graphs (group synthetic), one marker per instance; the marker
                 identifies the graph, its shade the value of tau, the dashed line is a least-squares
                 power-law fit. Instances with no remaining node are left out (printed).
  effectiveness  EQ3 (fig_effectiveness): FPA solution size relative to the optimum (group effectiveness).
  ablation       EQ4 (fig_ablation): running time against tau of EBA and EBA without Rule 1, 2, 3, or the
                 FPA incumbent (group ablation), one panel per dataset.
  branching      EQ5 (fig_branching): EBA (CP-Branching), EBA-Pivot, EBA-Degeneracy, Binary, flat Sym-SE
                 and CP-SAT (group branching).

The panels of EQ4 and EQ5 show tau = 0.5, 0.6, ..., 0.9 as in the paper (--fine: every 0.05); a run
that reached the time limit is drawn on the dotted line labelled T/O at the limit recorded in the
results file (the paper: one day), one stopped by the memory limit there too and labelled M/O in the
colour of its curve. --datasets sets the panels (default: those of
the paper). Needs matplotlib; reads the results file with analysis/tables.py.
"""
import argparse, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tables as T                                            # noqa: E402
rx = T.rx

INK, GRID = "#222222", "#dedede"
LIMIT = 86400                                                 # default time limit (one day)
MIN_DECADES = 3
TAU_COARSE = [0.5, 0.6, 0.7, 0.8, 0.9]
# (algo, label, colour, marker, line style)
SERIES = {
    "ablation": [("eba", "EBA", "#e34948", "o", "-"), ("eba_no_rule1", "EBA$-$Rule 1", "#2a78d6", "s", "--"),
                 ("eba_no_rule2", "EBA$-$Rule 2", "#1baf7a", "^", "--"),
                 ("eba_no_rule3", "EBA$-$Rule 3", "#6250d6", "D", "--"),
                 ("eba_no_incumbent", "EBA$-$FPA", "#8c564b", ">", "--")],
    "branching": [("eba", "EBA (CP-Branching)", "#e34948", "o", "-"), ("eba_no_pivot", "EBA$-$Pivot", "#00a3c4", "<", "--"),
                  ("eba_no_degeneracy", "EBA$-$Degeneracy", "#7a7a7a", "*", "--"),
                  ("binary", "Binary", "#262626", "P", "-."), ("flat_symse", "Flat Sym-SE", "#eda100", "v", "-."),
                  ("cpsat", "CP-SAT", "#008300", "h", "-.")],
}
PANELS = {"ablation": "FB FW CA CD CS WI".split(), "branching": "FB GO FW CA LU CS".split()}
EFFECTIVENESS = "FB FW CA CD RO LU".split()
LFR_LABEL = {"lfr_n20000": "L2", "lfr_n40000": "L4", "lfr_n60000": "L6", "lfr_n80000": "L8", "lfr_n100000": "L10"}


def style(size=9):
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Nimbus Roman", "Times New Roman", "STIXGeneral"],
                         "mathtext.fontset": "stix", "font.size": size, "axes.edgecolor": INK,
                         "axes.linewidth": 0.7, "xtick.color": INK, "ytick.color": INK})


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300)
    plt.close(fig)
    print(f"written {os.path.join(out, name)}.{{pdf,png}}")


# ---------------------------------------------------------------- EQ2
def scalability(R, out, **_):
    pts = {}
    for d in rx.LFR:
        for t in rx.TAU_LFR:
            s, e = R.ok("synthetic", d, t, "root_stats"), R.ok("synthetic", d, t, "eba")
            if not (s and e):
                continue
            n = int(dict(x.split("=") for x in s["extra"].split(";"))["V13"])
            if n == 0:
                print(f"  left out (certified at the root): {d} {t}")
                continue
            pts[(d, T.tau_key(t))] = (n, float(e["time"]))
    if len(pts) < 3:
        print("scalability: fewer than three instances, figure not drawn")
        return
    x = [math.log10(n) for n, _ in pts.values()]; y = [math.log10(t) for _, t in pts.values()]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    b = sum((p - mx) * (q - my) for p, q in zip(x, y)) / sum((p - mx) ** 2 for p in x)
    a0 = my - b * mx
    r2 = 1 - sum((q - a0 - b * p) ** 2 for p, q in zip(x, y)) / sum((q - my) ** 2 for q in y)
    print(f"  fit: log10(time) = {a0:.3f} + {b:.3f} log10(nodes), R^2 = {r2:.3f}, {len(pts)} instances")
    style(8)
    taus = [T.tau_key(t) for t in rx.TAU_LFR]
    blues = matplotlib.colormaps["Blues"]
    colour = {t: matplotlib.colors.to_hex(blues(0.25 + 0.7 * i / (len(taus) - 1))) for i, t in enumerate(taus)}
    markers = dict(zip(rx.LFR, ["o", "s", "^", "D", "v"]))
    ms = {"o": 4.5, "s": 4.5, "D": 5.4 / math.sqrt(2), "^": 5.6, "v": 5.6}
    W, H, L, Rm, Tm, Bm = 204, 92, 38, 10, 15, 24
    w, h = L + W + Rm, Bm + H + Tm
    fig = plt.figure(figsize=(w / 72, h / 72))
    ax = fig.add_axes([L / w, Bm / h, W / w, H / h])
    ax.set_xscale("log"); ax.set_yscale("log")
    xmin, xmax = 10 ** math.floor(min(x)), 10 ** math.ceil(max(x))
    ymin, ymax = 10 ** math.floor(min(y)), 10 ** math.ceil(max(y))
    xt = [10.0 ** k for k in range(int(math.log10(xmin)), int(math.log10(xmax)) + 1)]
    yt = [10.0 ** k for k in range(int(math.log10(ymin)), int(math.log10(ymax)) + 1)]
    for v in yt:
        ax.axhline(v, color=GRID, lw=0.5, zorder=0)
    for v in xt:
        ax.axvline(v, color=GRID, lw=0.5, zorder=0)
    xs = [10 ** min(x), 10 ** max(x)]
    ax.plot(xs, [10 ** (a0 + b * math.log10(v)) for v in xs], color="#555555", lw=0.9, ls="--", zorder=2)
    for (d, t), (n, tm) in sorted(pts.items(), key=lambda kv: kv[0][1]):       # darkest (largest tau) on top
        mk = markers[d]
        ax.plot([n], [tm], ls="none", marker=mk, ms=ms[mk], mfc=colour[t], mec="#1a1a1a", mew=0.35, zorder=3)
    for d in rx.LFR:
        mk = markers[d]
        ax.plot([], [], ls="none", marker=mk, ms=ms[mk], mfc="#bdbdbd", mec="#1a1a1a", mew=0.35, label=LFR_LABEL[d])
    kx = [0.64 + 0.032 * i for i in range(len(taus))]                         # tau key, lower right
    for px, t in zip(kx, taus):
        ax.plot([px], [0.13], ls="none", marker="s", ms=ms["s"], mfc=colour[t], mec="#1a1a1a", mew=0.3,
                transform=ax.transAxes, zorder=4)
    ax.text(kx[0] - 0.03, 0.13, r"$\tau$", transform=ax.transAxes, ha="right", va="center", fontsize=7, color=INK)
    ax.text(kx[0], 0.05, taus[0], transform=ax.transAxes, ha="center", va="center", fontsize=6.5, color=INK)
    ax.text(kx[-1], 0.05, taus[-1], transform=ax.transAxes, ha="center", va="center", fontsize=6.5, color=INK)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    pow10 = lambda v: "1" if round(math.log10(v)) == 0 else f"$10^{{{round(math.log10(v))}}}$"
    ax.xaxis.set_major_locator(FixedLocator(xt)); ax.set_xticklabels([pow10(v) for v in xt])
    ax.yaxis.set_major_locator(FixedLocator(yt)); ax.set_yticklabels([pow10(v) for v in yt])
    ax.xaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_locator(NullLocator())
    ax.tick_params(length=2, labelsize=7, pad=1.5)
    ax.set_xlabel("Remaining nodes after pruning", fontsize=8, labelpad=1)
    ax.set_ylabel("Running time (sec.)", fontsize=8, labelpad=2)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=len(rx.LFR), frameon=False, fontsize=7,
              title="LFR", title_fontsize=7, handletextpad=0.05, columnspacing=0.7, borderaxespad=0.2, handlelength=1.0)
    save(fig, out, "fig_scalability")


# ---------------------------------------------------------------- EQ3
def effectiveness(R, out, datasets=None, **_):
    datasets = datasets or EFFECTIVENESS
    taus = rx.TAU_FIG
    style(8)
    ncol = math.ceil(len(datasets) / 2)
    fig, axes = plt.subplots(2, ncol, figsize=(252 / 72, 150 / 72), sharey=True, squeeze=False)
    xs = list(range(len(taus)))
    for ax, d in zip(axes.flat, datasets):
        vals = []
        for t in taus:
            f, e = R.ok("effectiveness", d, t, "fpa"), R.ok("effectiveness", d, t, "eba")
            if f is None or e is None:
                sys.exit(f"effectiveness: no solved FPA and EBA run for {d} {t}")
            if int(f["size"]) > int(e["size"]):
                sys.exit(f"effectiveness: FPA {f['size']} larger than the optimum {e['size']} on {d} {t}")
            vals.append(int(f["size"]) / int(e["size"]))
        ax.bar(xs, vals, width=0.72, color="#aab8e2", edgecolor=INK, linewidth=0.5, zorder=2)
        ax.set_ylim(0, 1.08)
        ax.set_yticks([0, 0.5, 1.0]); ax.set_yticklabels(["0.0", "0.5", "1.0"])
        ax.set_xticks(xs)
        ax.set_xticklabels([f"{t:.1f}" if round(t * 100) % 10 == 0 else "" for t in taus])
        ax.set_xlim(-0.7, len(taus) - 0.3)
        ax.grid(True, axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(length=2, labelsize=8, pad=1.5)
        ax.set_title(d, fontsize=9, pad=2.5)
    for ax in list(axes.flat)[len(datasets):]:
        ax.axis("off")
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.17, top=0.9, wspace=0.12, hspace=0.62)
    fig.text(0.005, 0.53, "Approximation ratio", rotation=90, ha="left", va="center", fontsize=9)
    fig.text(0.56, 0.0, r"$\tau$", ha="center", va="bottom", fontsize=10)
    save(fig, out, "fig_effectiveness")


# ---------------------------------------------------------------- EQ4 / EQ5
CANDIDATES = [(0, 3.5, "center", "bottom"), (3.5, 3.5, "left", "bottom"), (-3.5, 3.5, "right", "bottom"),
              (5, 0, "left", "center"), (-5, 0, "right", "center"),
              (3.5, -3.5, "left", "top"), (-3.5, -3.5, "right", "top"),
              (-12, 3.5, "right", "bottom"), (12, 3.5, "left", "bottom"),
              (-12, -3.5, "right", "top"), (12, -3.5, "left", "top"),
              (-3.5, -14, "right", "top"), (3.5, -14, "left", "top")]


def place_labels(fig, ax, name_txt, labels, points, lines, lim):
    """Failure labels (M/O) near their markers, avoiding the dataset name, markers and lines."""
    r = fig.canvas.get_renderer()
    pt = fig.dpi / 72
    axbb = ax.get_window_extent(r)
    taken = [name_txt.get_window_extent(r)]
    pix = [ax.transData.transform((x, y)) for x, y in points]
    samples = []
    for pl in lines:
        q = [ax.transData.transform(p) for p in pl]
        for (x0, y0), (x1, y1) in zip(q, q[1:]):
            samples += [(x0 + (x1 - x0) * k / 24, y0 + (y1 - y0) * k / 24) for k in range(25)]
    rad = 2.6 * pt
    for xs, text, col in labels:
        best = None
        for xi, x in enumerate(xs):
            px, py = ax.transData.transform((x, lim))
            for ci, (dx, dy, ha, va) in enumerate(CANDIDATES):
                xd, yd = ax.transData.inverted().transform((px + dx * pt, py + dy * pt))
                t = ax.text(xd, yd, text, fontsize=7.5, color=col, ha=ha, va=va, zorder=7)
                bb = t.get_window_extent(r)
                cost = 10 * sum(bb.overlaps(o) for o in taken) + 0.01 * ci + 0.001 * xi
                cost += sum(bb.x0 - rad < qx < bb.x1 + rad and bb.y0 - rad < qy < bb.y1 + rad for qx, qy in pix)
                cost += 0.05 * sum(bb.x0 < qx < bb.x1 and bb.y0 < qy < bb.y1 for qx, qy in samples)
                if not (axbb.x0 <= bb.x0 and bb.x1 <= axbb.x1 and axbb.y0 <= bb.y0 and bb.y1 <= axbb.y1):
                    cost += 100
                if best is None or cost < best[0]:
                    if best is not None:
                        best[1].remove()
                    best = (cost, t)
                else:
                    t.remove()
        taken.append(best[1].get_window_extent(r))


def stack(fig, ax, xs, raw, order, step_pt=1.5, same_dec=0.07):
    """Marker positions with coinciding values stacked upwards so that no marker hides another; the
    lines stay at the true values. Values within same_dec decades form a group; inside a group the
    markers are placed in legend order, each step_pt points above the previous one."""
    pt = fig.dpi / 72
    adj = {a: list(raw[a]) for a in order}
    for i, x in enumerate(xs):
        rest = sorted(order, key=lambda a: raw[a][i])
        while rest:
            base = math.log10(raw[rest[0]][i])
            group = [a for a in rest if math.log10(raw[a][i]) - base <= same_dec]
            rest = [a for a in rest if a not in group]
            for rank, a in enumerate(sorted(group, key=order.index)):
                if rank:
                    px, py = ax.transData.transform((x, raw[a][i]))
                    adj[a][i] = ax.transData.inverted().transform((px, py + rank * step_pt * pt))[1]
    return adj


def value(R, group, d, t, a):
    """seconds, 'TO' or 'MO'; None if the run is missing."""
    r = R.row(group, d, t, a)
    if r is None:
        return None
    if r["status"] == "ok":
        return max(float(r["time"]), 1e-4)
    return r["status"]


def limit_of(R, group, d, taus, series):
    """Time limit at which the failed runs of a panel are drawn: the largest recorded limit of its
    runs (the paper: one day); a panel whose runs used different limits is reported."""
    lims = {int(R.row(group, d, t, a)["timeout"] or LIMIT) for t in taus for a, *_ in series if R.row(group, d, t, a)}
    if len(lims) > 1:
        print(f"  {group}: {d} mixes time limits {sorted(lims)}; failed runs drawn at {max(lims)} s")
    return max(lims) if lims else LIMIT


def panels(R, out, group, datasets=None, fine=False, **_):
    series = SERIES[group]
    taus = rx.TAU_FIG if fine else TAU_COARSE
    datasets = datasets or PANELS[group]
    complete = []
    for d in datasets:
        missing = [(t, a) for t in taus for a, *_ in series if value(R, group, d, t, a) is None]
        if missing:
            print(f"  {group}: {d} left out, no run for " + ", ".join(f"{a} at {t}" for t, a in missing[:4])
                  + (" ..." if len(missing) > 4 else ""))
        else:
            complete.append(d)
    if not complete:
        print(f"{group}: no complete panel, figure not drawn")
        return
    style(9)
    fig, axes = plt.subplots(1, len(complete), figsize=(515.52 / 72, 96 / 72), squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(left=0.072, right=0.995, bottom=0.2, top=0.8, wspace=0.34)
    todo = []
    n = len(series)
    for ax, d in zip(axes, complete):
        vals = {a: [value(R, group, d, t, a) for t in taus] for a, *_ in series}
        lim = limit_of(R, group, d, taus, series)
        ys = [v if isinstance(v, float) else lim for vs in vals.values() for v in vs]
        lo, hi = math.log10(min(ys)), math.log10(max(ys))
        if hi - lo < MIN_DECADES:
            mid = (lo + hi) / 2
            lo, hi = mid - MIN_DECADES / 2, mid + MIN_DECADES / 2
        span = hi - lo
        bottom, top = lo - 0.12 * span, hi + 0.5 * span
        xs = list(taus)
        ax.set_yscale("log")
        ax.set_ylim(10 ** bottom, 10 ** top)
        ax.set_xlim(0.47, max(xs) + 0.03)
        raw = {a: [v if isinstance(v, float) else lim for v in vals[a]] for a, *_ in series}
        adj = stack(fig, ax, xs, raw, [a for a, *_ in series])
        labels, points, lines = [], [], []
        for idx, (a, lab, col, mk, ls) in enumerate(series):
            ax.plot(xs, raw[a], color=col, ls=ls, lw=1.3 if a == "eba" else 0.9, zorder=2 + (n - idx) / 10)
            ax.plot(xs, adj[a], ls="none", marker=mk, ms=5.2 if mk in "*P" else 4.3, mfc=col,
                    mec="#1a1a1a", mew=0.6, zorder=4 + (n - idx) / 10)
            points += list(zip(xs, adj[a]))
            lines.append(list(zip(xs, raw[a])))
            fx = [x for x, v in zip(xs, vals[a]) if v == "MO"]
            if fx:
                labels.append((fx, "M/O", col))
        decades = list(range(math.ceil(bottom), math.floor(top) + 1))
        step = max(1, math.ceil(len(decades) / 4))
        decades = [k for k in decades if k % step == 0]
        has_to = any(v == "TO" for vs in vals.values() for v in vs)
        if has_to:
            decades = [k for k in decades if k < math.log10(lim) - 1]
            ax.axhline(lim, color="#555555", ls=":", lw=0.8, zorder=1)
        ticks = [10.0 ** k for k in decades] + ([lim] if has_to else [])
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.set_yticklabels([f"$10^{{{k}}}$" for k in decades] + (["T/O"] if has_to else []))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_xticks(TAU_COARSE)
        ax.grid(True, color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(length=2, labelsize=8, pad=1.5)
        name_txt = ax.text(0.95, 0.95, d, transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
                           fontweight="bold", zorder=9)
        todo.append((ax, name_txt, labels, points, lines, lim))
    fig.text(0.003, 0.5, "Running time (sec.)", rotation=90, ha="left", va="center", fontsize=9.5)
    fig.text(0.535, 0.0, r"$\tau$", ha="center", va="bottom", fontsize=11)
    handles = [Line2D([], [], color=col, ls=ls, lw=1.3 if a == "eba" else 0.9, marker=mk,
                      ms=5.2 if mk in "*P" else 4.3, mfc=col, mec="#1a1a1a", mew=0.6) for a, lab, col, mk, ls in series]
    fig.legend(handles, [lab for _, lab, *_ in series], loc="upper center", ncol=n, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.535, 1.03), handlelength=2.0, handletextpad=0.5, columnspacing=0.9, borderaxespad=0)
    fig.canvas.draw()
    for ax, name_txt, labels, points, lines, lim in todo:
        place_labels(fig, ax, name_txt, labels, points, lines, lim)
    save(fig, out, f"fig_{group}")


FIGURES = {"scalability": scalability, "effectiveness": effectiveness,
           "ablation": lambda R, out, **kw: panels(R, out, "ablation", **kw),
           "branching": lambda R, out, **kw: panels(R, out, "branching", **kw)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("which", nargs="?", default="all", choices=["all"] + list(FIGURES))
    ap.add_argument("results", nargs="?", default=os.path.join(ROOT, "results", "results.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "figures"))
    ap.add_argument("--datasets", help="panels (comma-separated abbreviations)")
    ap.add_argument("--fine", action="store_true", help="EQ4/EQ5: every 0.05 of tau instead of every 0.1")
    a = ap.parse_args()
    if not os.path.exists(a.results):
        sys.exit(f"{a.results} not found")
    R = T.Results(a.results)
    datasets = a.datasets.split(",") if a.datasets else None
    for name in (list(FIGURES) if a.which == "all" else [a.which]):
        print(f"{name}:")
        FIGURES[name](R, a.out, datasets=datasets, fine=a.fine)


if __name__ == "__main__":
    main()
