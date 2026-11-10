#!/usr/bin/env python3
"""Figure of EQ6 (empirical calibration of tau): fraction of ground-truth communities and of randomly
grown connected node sets of the same size that are flexi-cliques, for tau = 0.10, 0.15, ..., 0.90.

  python3 analysis/fig_tau_calibration.py [--results DIR] [--out DIR]
      -> results/figures/fig_tau_calibration.pdf, fig_tau_calibration.txt

Input: results/tau_calibration/records_AM.csv and records_DB.csv of experiments/tau_calibration.py
(per community: its size s, the minimum degree of the community and of each of its random sets). A set
of size s with minimum degree d is a flexi-clique at tau iff d >= floor(s^tau), computed as floorPow in
src/graph.h. The gap at tau is the fraction of communities that are flexi-cliques minus the fraction of
random sets that are.

Sampling (as in the paper): R = 100 samples of M = 1,000 communities drawn uniformly with replacement
from the evaluated communities (numpy default_rng(2027), first Amazon, then DBLP). The plotted fractions
are means over the samples, and the vertical bar marks the tau with the largest mean gap. The text file
lists the curves, the maximising tau and the range of the per-sample maximiser (ties: smallest tau).
Needs numpy (1.22 or later) and matplotlib.
"""
import argparse, csv, math, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NETS = [("AM", "Amazon"), ("DB", "DBLP")]                  # sampling order and panels
TAUS = np.array([i / 20 for i in range(2, 19)])             # 0.10, 0.15, ..., 0.90
R, M, SEED = 100, 1000, 2027
COMM = ("Communities", "#a0522d", "-", 1.3)
RAND = ("Randomly grown sets", "#56606b", "--", 1.1)
INK, GRIDC = "#222222", "#dedede"


def threshold(s, tau):
    """floor(s^tau), with the same 1e-12 tolerance as floorPow in src/graph.h."""
    return math.floor(s ** tau * (1 + 1e-12) + 1e-12)


def load(path):
    """Per community and tau: acc_c = community is a flexi-clique (0/1), acc_r = number of its random sets
    that are (0..k); k = random sets per community."""
    with open(path) as f:
        rows = [list(map(int, r)) for r in list(csv.reader(f))[1:]]
    if not rows:
        sys.exit(f"{path}: no records")
    a = np.array(rows)
    size, dc, dr = a[:, 0], a[:, 1], a[:, 2:]
    th = {s: [threshold(int(s), t) for t in TAUS] for s in np.unique(size)}
    TH = np.array([th[s] for s in size])                                  # communities x taus
    acc_c = (dc[:, None] >= TH).astype(np.int64)
    acc_r = (dr[:, :, None] >= TH[:, None, :]).sum(axis=1)
    return acc_c, acc_r, dr.shape[1]


def calibrate(acc_c, acc_r, k, rng):
    """Mean fractions over R samples of M communities, the tau of the largest mean gap, and the
    maximiser of the gap in every sample (computed on exact integer sums)."""
    idx = rng.integers(0, len(acc_c), size=(R, M))
    comm = acc_c[idx].sum(axis=(0, 1)) / (R * M)
    rand = acc_r[idx].sum(axis=(0, 1)) / (k * R * M)
    per_sample = (k * acc_c - acc_r)[idx].sum(axis=1)                      # R x taus
    j = int(np.argmax(per_sample.sum(axis=0)))
    return comm, rand, j, TAUS[np.argmax(per_sample, axis=1)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=os.path.join(ROOT, "results", "tau_calibration"),
                    help="directory of records_<ABBR>.csv")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "figures"), help="output directory")
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)
    res, lines = {}, [f"EQ6: fraction of flexi-cliques for tau = 0.10, 0.15, ..., 0.90; means over {R} samples "
                      f"of {M} communities drawn with replacement"]
    for net, name in NETS:
        path = os.path.join(a.results, f"records_{net}.csv")
        if not os.path.exists(path):
            sys.exit(f"{path} not found (python3 experiments/tau_calibration.py)")
        acc_c, acc_r, k = load(path)
        comm, rand, j, amax = calibrate(acc_c, acc_r, k, rng)
        res[net] = (comm, rand, j)
        lo, hi = np.quantile(amax, [0.025, 0.975], method="inverted_cdf")
        lines += [f"{name} ({net}): {len(acc_c)} communities, {k} random sets each",
                  "  tau   communities  random sets  gap"]
        lines += [f"  {t:.2f}  {c:.3f}        {r:.3f}        {c - r:+.3f}" for t, c, r in zip(TAUS, comm, rand)]
        lines += [f"  largest mean gap at tau = {TAUS[j]:.2f}: {comm[j] - rand[j]:+.3f} "
                  f"(communities {comm[j]:.3f}, random sets {rand[j]:.3f})",
                  f"  tau of the largest gap per sample: 2.5-97.5% [{lo:.2f}, {hi:.2f}], "
                  f"min-max [{amax.min():.2f}, {amax.max():.2f}]"]
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "fig_tau_calibration.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    plt.rcParams.update({"font.family": "serif", "font.serif": ["Nimbus Roman", "Times New Roman", "STIXGeneral"],
                         "mathtext.fontset": "stix", "font.size": 9, "axes.edgecolor": INK,
                         "axes.linewidth": 0.7, "xtick.color": INK, "ytick.color": INK})
    fig, axes = plt.subplots(1, len(NETS), figsize=(3.5, 1.6), sharey=True)
    for ax, (net, name) in zip(axes, NETS):
        comm, rand, j = res[net]
        tb = TAUS[j]
        for ys, (lab, col, ls, lw) in ((comm, COMM), (rand, RAND)):
            ax.plot(TAUS, ys, color=col, ls=ls, lw=lw, zorder=3, label=lab)
        ax.plot([tb, tb], [rand[j], comm[j]], color=INK, lw=0.9, zorder=4)           # the largest gap
        for y in (rand[j], comm[j]):
            ax.plot([tb - 0.025, tb + 0.025], [y, y], color=INK, lw=0.9, zorder=4)
        ax.text(tb + 0.04, rand[j] + 0.25 * (comm[j] - rand[j]), f"max. gap\n$\\tau={tb:.2f}$",
                fontsize=6.5, color=INK, ha="left", va="center")
        ax.set_xlim(0.1, 0.9)
        ax.set_ylim(0, 1.05)
        ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.grid(True, color=GRIDC, lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(length=2, labelsize=7, pad=1.5)
        ax.text(0.96, 0.95, name, transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK)
        ax.set_xlabel(r"$\tau$", fontsize=8.5, labelpad=0.5)
    axes[0].set_ylabel("Fraction of flexi-cliques", fontsize=8)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False, fontsize=7.5, bbox_to_anchor=(0.55, 1.03),
               handlelength=2.2)
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.21, top=0.85, wspace=0.12)
    fig.savefig(os.path.join(a.out, "fig_tau_calibration.pdf"))
    print(f"wrote {os.path.join(a.out, 'fig_tau_calibration.pdf')}")


if __name__ == "__main__":
    main()
