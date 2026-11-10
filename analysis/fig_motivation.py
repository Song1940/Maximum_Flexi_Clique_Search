#!/usr/bin/env python3
"""Figure 1 of the paper: density of ground-truth communities against their size, with the minimum
density implied by each cohesive subgraph model, on DBLP (left) and Amazon (right).

  python3 analysis/fig_motivation.py [--results DIR] [--out DIR]
      -> results/figures/fig_motivation_density.pdf (Fig. 1a), fig_motivation_mismatch.pdf (Fig. 1b)

Input: results/community_stats/<ABBR>.{all,top5000}.csv of experiments/community_stats.py.
Minimum density implied at size n (k = 5, gamma = 0.5, tau = 0.5): clique 1; k-plex max(0, n - k) / (n - 1);
gamma-quasi-clique gamma (its density variant); flexi-clique floor(n^tau) / (n - 1) (floor as floorPow in
src/graph.h). The empirical trend is a power law a n^b fitted by least squares on log scale to the
per-size median densities of all SNAP communities (<ABBR>.all.csv; sizes 2..100 with at least 3
communities).
  (a) The points are the distinct top-5000 node sets (<ABBR>.top5000.csv) of sizes 5..100, one point per
      distinct (size, density) pair, shaded like a scatter of opacity 0.4 (k coincident points: opacity
      1 - 0.6^k); the threshold curves start at size 3 (k-plex: 5).
  (b) |threshold density - trend| from size max(3, smallest community size) to 100, so the trend is not
      extrapolated below the data (DBLP communities have at least 6 nodes).
The axes carry no tick labels (the figure shows tendencies only). Needs numpy and matplotlib (3.4 or later).
"""
import argparse, csv, math, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NETS = [("DB", "DBLP"), ("AM", "Amazon")]
START, FLEXI_START, FIT_MIN, SMAX, MIN_PTS = 5, 3, 2, 100, 3
K, GAMMA, TAU = 5, 0.5, 0.5


def floor_pow(n, tau):
    """floor(n^tau) for an array of sizes, with the same 1e-12 tolerance as floorPow in src/graph.h."""
    return np.array([math.floor(x ** tau * (1 + 1e-12) + 1e-12) for x in n.astype(int).tolist()], dtype=float)


MODELS = [  # name, minimum density at size n, first size drawn, colour, line style, width
    ("Clique", lambda n: np.ones_like(n), FLEXI_START, "#7d3c98", ":", 1.2),
    ("$k$-plex", lambda n: np.maximum(n - K, 0) / (n - 1), START, "#1a5276", "-.", 1.1),
    (r"$\gamma$-quasi-clique", lambda n: np.full_like(n, GAMMA), FLEXI_START, "#c0392b", "--", 1.3),
    ("Flexi-clique", lambda n: floor_pow(n, TAU) / (n - 1), FLEXI_START, "#1a7a4a", "-", 1.7),
]
INK, GRIDC, POINT, ALPHA = "#222222", "#dedede", "#777777", 0.4


def load(path):
    """Sizes and densities (NaN for a community of one node) of a community_stats CSV."""
    if not os.path.exists(path):
        sys.exit(f"{path} not found (python3 experiments/community_stats.py)")
    with open(path) as f:
        rows = list(csv.DictReader(f))
    n = np.array([int(r["size"]) for r in rows])
    d = np.array([float(r["density"]) if r["density"] else np.nan for r in rows])
    return n, d


def fit_trend(n, d):
    """Power law a n^b fitted (least squares, log scale) to the per-size median densities of sizes
    FIT_MIN..SMAX with at least MIN_PTS communities."""
    ok = (n >= FIT_MIN) & (n <= SMAX) & ~np.isnan(d)
    n, d = n[ok], d[ok]
    xs, ys = [], []
    for s in np.unique(n):
        v = d[n == s]
        if len(v) >= MIN_PTS and s > 1 and np.median(v) > 0:
            xs.append(float(s))
            ys.append(np.median(v))
    b, loga = np.polyfit(np.log(xs), np.log(ys), 1)
    return np.exp(loga), b


def points(n, d):
    """Distinct (size, density) pairs of sizes START..SMAX, their multiplicities k, and their colours:
    k coincident points of opacity ALPHA composite to opacity 1 - (1 - ALPHA)^k on white."""
    ok = (n >= START) & (n <= SMAX) & ~np.isnan(d) & (d <= 1.0)
    pairs, k = np.unique(np.column_stack([n[ok].astype(float), d[ok]]), axis=0, return_counts=True)
    order = np.argsort(k, kind="stable")                          # darker points drawn last
    pairs, k = pairs[order], k[order]
    rgb = np.array(mcolors.to_rgb(POINT))
    cols = [tuple(1 - (1 - (1 - ALPHA) ** c) * (1 - rgb)) for c in k]   # white * (1 - o) + POINT * o
    return pairs[:, 0], pairs[:, 1], k, cols


def style(ax, name, ytop):
    ax.set_xlim(FLEXI_START - 4, SMAX + 3)                           # curves start a little off the y-axis
    ax.set_ylim(0, ytop)
    ax.set_xticks(np.arange(20, SMAX + 1, 20))
    ax.set_yticks(np.linspace(0, ytop, 5)[1:-1])
    ax.grid(True, color=GRIDC, lw=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelbottom=False, labelleft=False)
    ax.text(0.96, 0.66, name, transform=ax.transAxes, ha="right", va="center", fontsize=9, color=INK)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=os.path.join(ROOT, "results", "community_stats"),
                    help="directory of <ABBR>.all.csv and <ABBR>.top5000.csv")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "figures"), help="output directory")
    a = ap.parse_args()
    data = {net: load(os.path.join(a.results, f"{net}.all.csv")) for net, _ in NETS}
    top = {net: load(os.path.join(a.results, f"{net}.top5000.csv")) for net, _ in NETS}
    fits = {net: fit_trend(*data[net]) for net, _ in NETS}
    for net, name in NETS:
        print(f"{name}: trend density {fits[net][0]:.4f} n^{fits[net][1]:.4f}")
    os.makedirs(a.out, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Nimbus Roman", "Times New Roman", "STIXGeneral"],
                         "mathtext.fontset": "stix", "font.size": 8, "axes.edgecolor": INK,
                         "axes.linewidth": 0.7})

    # (a) density of the communities and minimum density of each model
    fig, axes = plt.subplots(1, 2, figsize=(3.4, 1.32), sharey=True)
    for ax, (net, name) in zip(axes, NETS):
        xs, ys, _, cols = points(*top[net])
        for x, y, c in zip(xs, ys, cols):                            # one marker per point keeps its own shade
            ax.plot([x], [y], ls="none", marker="o", ms=1.5, mfc=c, mec=c, mew=0, zorder=2)
        for lab, f, first, col, ls, lw in MODELS:
            g = np.arange(first, SMAX + 1).astype(float)
            ax.plot(g, f(g), color=col, ls=ls, lw=lw, zorder=3, label=lab)
        style(ax, name, 1.05)
    c1 = points(np.array([10]), np.array([0.5]))[3][0]               # the shade of a single point
    axes[0].plot([], [], ls="none", marker="o", ms=1.5, mfc=c1, mec=c1, mew=0, label="Ground-truth community")
    axes[0].set_ylabel("Density threshold", fontsize=8.5, labelpad=2)
    fig.supxlabel("Community size", fontsize=8.5, y=0.02)
    h, l = axes[0].get_legend_handles_labels()
    order = [len(h) - 1] + list(range(len(h) - 1))                   # the communities first
    fig.legend([h[i] for i in order], [l[i] for i in order], loc="upper center", ncol=5, frameon=False,
               fontsize=5.6, bbox_to_anchor=(0.53, 1.03), handlelength=1.7, columnspacing=0.6, handletextpad=0.3,
               markerscale=2.5)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.12, top=0.86, wspace=0.06)
    fig.savefig(os.path.join(a.out, "fig_motivation_density.pdf"))

    # (b) distance between the minimum density of each model and the trend
    fig, axes = plt.subplots(1, 2, figsize=(3.4, 1.25), sharey=True)
    for ax, (net, name) in zip(axes, NETS):
        ta, tb = fits[net]
        g = np.arange(max(FLEXI_START, int(data[net][0].min())), SMAX + 1).astype(float)
        tr = ta * g ** tb
        for lab, f, first, col, ls, lw in MODELS:
            ax.plot(g, np.abs(f(g) - tr), color=col, ls=ls, lw=lw, label=lab)
        style(ax, name, 1.05)
    axes[0].set_ylabel("Density mismatch", fontsize=8.5, labelpad=2)
    fig.supxlabel("Community size", fontsize=8.5, y=0.02)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.13, top=0.97, wspace=0.06)
    fig.savefig(os.path.join(a.out, "fig_motivation_mismatch.pdf"))
    print(f"wrote fig_motivation_density.pdf and fig_motivation_mismatch.pdf in {a.out}")


if __name__ == "__main__":
    main()
