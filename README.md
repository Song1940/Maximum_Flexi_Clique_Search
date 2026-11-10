# Exact Maximum Flexi-Clique Computation

Source code, datasets, and experiment scripts of the paper

> **Exact Maximum Flexi-Clique Computation under Non-Hereditary Property and Connectivity**
> Song Kim, Kaiqiang Yu, Taejoon Han, Hyewon Kim, Junghoon Kim, Susik Yoon, Jungeun Kim.
> Submitted to IEEE ICDE 2027.

The **technical report** with the complete proofs (hardness for a fixed τ, exactness of EBA) and
additional experiments (solution sizes of NPA and FPA, branching width with and without the pivot)
is [`technical_report.pdf`](technical_report.pdf).

Given a graph *G* and a parameter τ ∈ (0, 1), a node set *H* is a **flexi-clique** if *G*[*H*] is
connected and every node of *H* has at least ⌊|*H*|^τ⌋ neighbours in *H*. The binary `flexi`
computes a maximum flexi-clique with:

| `--algo` | Algorithm | Paper |
|---|---|---|
| `eba` | **EBA** (Exact Branch-and-bound Algorithm): exact search with Connectivity-Preserving Branching (CP-Branching), Rules 1–3, and the margin pivot, starting from the FPA solution | Section V, Algorithm 1 |
| `fpa` | **FPA** (Flexi-Prune Algorithm): scalable heuristic — global peeling, augmentation, and slack-aware peeling | Section IV |
| `npa` | **NPA**: the existing peeling heuristic for flexi-cliques (Kim et al., CIKM 2024), the baseline of FPA | Section VI |

## Contents

| Path | Content |
|---|---|
| `src/` | C++ implementation of EBA, FPA, and NPA, and the variants of the experiments |
| `baselines/` | the other baselines: CP-SAT (`cpsat/`) and the query-centred solver of the case study (`qcs/`) |
| `graphs/` | the nine smaller real-world graphs and the scripts that download the other four and the communities |
| `experiments/` | scripts that run the experiments and write their results to `results/` |
| `analysis/` | scripts that compute the tables and draw the figures from `results/` |
| `tests/` | validation against exhaustive enumeration |
| `verification/` | independent check of the optimality results, and a comparison with exhaustive enumeration on small graphs ([`verification/README.md`](verification/README.md)) |
| `example/` | a small example graph |
| `technical_report.pdf` | technical report: complete proofs and additional experiments |

## Build

Requirements: a C++17 compiler and CMake ≥ 3.15.

```bash
cmake -S . -B build
cmake --build build -j
```

This produces `build/flexi`, `build/qcs`, and the variants in `build/variants/` (see
[Variants](#variants)); `-DFLEXI_VARIANTS=OFF` builds `flexi` and `qcs` only. The code is compiled
with `-O3 -march=native`, as for the paper; `-DFLEXI_NATIVE=OFF` drops `-march=native`.

The scripts need Python 3.8 or later. The experiment and dataset scripts use only the standard
library, the analysis scripts need numpy (1.22 or later) and matplotlib (3.4 or later), and CP-SAT
needs the package `ortools` (the paper used 9.15). The experiment runner needs Linux (`prlimit` and
`/usr/bin/time`). The conversion of the LAW graph needs Java 11 or later.

## Usage

```bash
./build/flexi --algo <eba|fpa|npa> --file <edge list> --tau <t> [--seed <id>] [--init <id,id,...>] [--out <file>]
```

| Option | Meaning |
|---|---|
| `--algo` | `eba` (exact), `fpa` or `npa` (heuristics) |
| `--file` | undirected edge list (see *Input format*) |
| `--tau` | the parameter τ, 0 < τ < 1 |
| `--seed` | (EBA only) compute the maximum flexi-clique that **contains** the given node |
| `--init` | (EBA only) a known flexi-clique, given as comma-separated node ids, used as the initial incumbent when it is valid and larger than the FPA solution |
| `--out` | also write the result to a file (`key=value` lines) |

Example (`example/example.txt`: two 4-cliques joined by two edges, plus a pendant node):

```text
$ ./build/flexi --algo eba --file example/example.txt --tau 0.5
Loaded graph: n=9 m=15
Algorithm : eba
Tau       : 0.5
Size      : 8
Time(s)   : 0.000110099
Valid     : yes
FPA_Time(s): 4.812e-05
Branches  : 0
Nodes     : 1,2,3,4,6,5,7,8
```

* `Size`, `Nodes`: the solution, with the node ids of the input file.
* `Time(s)`: running time, graph loading excluded. For EBA it includes its FPA phase (`FPA_Time(s)`).
* `Valid`: the solution is checked to be a flexi-clique (connected, every degree ≥ ⌊|*H*|^τ⌋).
* `Branches`: number of EBA search nodes (0 when the pruning proves the FPA solution optimal at the root).

Exit status: 0 on success (also when the graph has no flexi-clique, e.g. no edges, in which case the result is empty and `Valid : no`), 1 on invalid arguments or I/O errors, and 2 if a returned solution fails the validity check (this should never happen).

## Input format

A plain-text edge list with one edge `u v` per line, where `u` and `v` are non-negative integer node ids
(32-bit range). Lines that are empty or start with `#` are ignored, and so are self-loops; duplicate
edges are merged. The graph is treated as undirected and simple, and node ids need not be contiguous.

## Datasets

The paper uses 13 real-world graphs from Network Repository, KONECT, SNAP, and LAW, the SNAP
ground-truth communities of Amazon and DBLP, and five LFR graphs (`python3 graphs/prepare_datasets.py --list`
prints the sources and sizes of the real-world graphs).

- The nine smaller graphs (FB, GO, FW, CA, CD, RO, LU, DB, AM) are included in `graphs/`.
- The other four are downloaded from their sources:
  [CS](https://networkrepository.com/co-papers-citeseer.php),
  [WI](https://networkrepository.com/web-it-2004.php) (Network Repository),
  [LJ](https://snap.stanford.edu/data/soc-LiveJournal1.html) (SNAP), and
  [SK](https://law.di.unimi.it/webdata/sk-2005/) (LAW).

```bash
python3 graphs/prepare_datasets.py CS WI LJ                 # -> graphs/<Abbr>.txt
bash graphs/convert_law.sh                                  # SK (Java, about 100 GB of disk)
python3 graphs/prepare_datasets.py --communities AM DB      # communities -> graphs/communities/
```

The files of Network Repository are written in canonical form (every edge once as `u < v`, sorted),
the KONECT and SNAP graphs in the order of the source file; these are exactly the files used in the
experiments (the order of the edges decides ties in the search).

The five LFR graphs of EQ2 (20,000, 40,000, 60,000, 80,000, and 100,000 nodes) were generated with
the LFR benchmark for undirected, unweighted networks (Lancichinetti, Fortunato, and Radicchi, 2008)
with average degree 40, maximum degree 200, mixing parameter 0.4, community sizes 20 to 1,000, and
no overlapping nodes.

```bash
echo <seed> > time_seed.dat
./benchmark -N <nodes> -k 40 -maxk 200 -mu 0.4 -minc 20 -maxc 1000 -on 0 -om 0
mv network.dat graphs/lfr/lfr_n<nodes>.txt                   # read by group synthetic
```

## Reproducing the experiments

All experiments run single-threaded with a time limit of one day and a memory limit of 400 GB, on
two Intel Xeon Gold 6342 CPUs with 1 TB of RAM. `experiments/run_experiments.py` runs the C++
experiments in groups and appends one row per run to `results/results.csv` (it resumes after an
interruption with the same command); `analysis/tables.py` and `analysis/figures.py` turn that file
into the tables and figures. The other experiments have their own scripts.

| Paper | Run | Tables and figures |
|---|---|---|
| Table II (dataset statistics) | — | `python3 analysis/dataset_stats.py graphs/*.txt` |
| EQ1, Table III | `python3 experiments/run_experiments.py --groups main` | `python3 analysis/tables.py` |
| EQ2, Fig. 5 | `--groups synthetic` | `python3 analysis/figures.py scalability` |
| EQ3, Fig. 6 | `--groups effectiveness` | `python3 analysis/figures.py effectiveness` |
| EQ4, Fig. 7 | `--groups ablation` | `python3 analysis/figures.py ablation` |
| EQ5, Fig. 8 | `--groups branching` | `python3 analysis/figures.py branching` |
| EQ6, Fig. 9 | `python3 experiments/tau_calibration.py` | `python3 analysis/fig_tau_calibration.py` |
| EQ7, Table IV and Fig. 10 | `python3 experiments/case_study.py run` | printed, and `results/case_study/summary.csv` |
| Fig. 1 | `python3 experiments/community_stats.py` | `python3 analysis/fig_motivation.py` |
| Technical report: sizes of NPA and FPA, branching width | `--groups main,branching`, then `--groups width` | `python3 analysis/tables.py` |

`run_experiments.py --help` lists the options: `--datasets`, `--taus`, and `--algos` select a
subset, `--timeout` and `--mem-gb` set the limits, and `--workers` runs several runs in parallel.
Running everything takes several CPU-weeks, most of it in runs that reach the one-day limit; for
example, `--groups main --datasets FB,GO,FW --timeout 3600` finishes in minutes. The paper ran
independent runs in parallel and measured the reported times on an otherwise idle machine, taking
the median of three runs below 60 s (CP-SAT and the width traces run once); `--workers 1` corresponds
to the latter. `analysis/tables.py` also prints the numbers quoted in the text of EQ1–EQ5
(`results/tables/stats.txt`). Group `width` takes the instances on which EBA branches from the rows of
group `branching` in the results file, so it runs in a separate invocation after that group.
Group `branching` includes SK (the technical report); the figure of the paper shows the other six graphs.

The case study (EQ7) runs EBA and the baselines for 50 Amazon communities with a time limit of
one hour per run; a flexi-clique run that reaches it keeps the maximum flexi-clique within the
2-hop neighbourhood of the query node, which it started from (`experiments/case_study.py`
documents the protocol). `--query 335413` runs the community shown in Fig. 10.

## Variants

The ablations and comparisons use variants of EBA that are compiled from the same source with the
macros documented in `src/eba.cpp`. They take the options of `flexi` and print their name on the
line `Variant`.

| Binary (`build/variants/`) | Paper |
|---|---|
| `flexi_no_rule1`, `flexi_no_rule2`, `flexi_no_rule3` | EBA−Rule 1, 2, 3 (EQ4) |
| `flexi_no_incumbent` | EBA−FPA: EBA without the FPA incumbent (EQ4) |
| `flexi_no_pivot` | EBA−Pivot: no pivot ordering and no candidate saturation; all reachable candidates are branched on in ascending order of adjusted degree (EQ5) |
| `flexi_no_degeneracy` | EBA−Degeneracy: the root children in ascending order of degree instead of the reverse degeneracy ordering (EQ5) |
| `flexi_binary` | Binary: include/exclude branching over reachable candidates (EQ5) |
| `flexi_flat_symse` | flat Sym-SE: the Sym-SE branching of FastQC on all candidates, connectivity checked at the leaves, with the pruning rules and incumbent of EBA (EQ5) |
| `flexi_root_stats` | nodes left after Rules 1 and 3 with the FPA incumbent, without branching (EQ2) |
| `flexi_width`, `flexi_width_no_pivot` | histogram of the branching width, with and without the pivot (technical report) |

The other baselines are `baselines/cpsat/cpsat_baseline.py` (CP-SAT on a constraint-programming
model, with the FPA incumbent and the global pruning of EBA; EQ5) and `build/qcs` (maximum clique,
*k*-plex, and γ-quasi-clique containing a query node; EQ7).

## Validation

`tests/validate_bruteforce.py` compares the binary with exhaustive enumeration on random small graphs:
EBA must return a valid flexi-clique of maximum size (also with `--seed`), and FPA a valid flexi-clique
no larger than the maximum.

```bash
python3 tests/validate_bruteforce.py build/flexi 1000 14    # 1,000 random graphs with up to 14 nodes
```

`verification/` checks the optimality results of EQ1 without the code of EBA. If every connected
component of the $\lambda$-core, with $\lambda=\lfloor(f+1)^\tau\rfloor$, has at most $f$ nodes, no
flexi-clique is larger than the FPA solution of size $f$; this is also exactly the condition under
which EBA finishes without branching. The directory also compares FPA and EBA with exhaustive
enumeration on eight families of small graphs ([`verification/README.md`](verification/README.md)).

## Source layout

| File | Content |
|---|---|
| `src/eba.cpp` | EBA: CP-Branching with the reverse degeneracy ordering at the root, Rules 1–3, margin pivot and saturation tests (Algorithm 1), and its variants |
| `src/fpa.cpp` | FPA main loop |
| `src/peel_split.cpp` | global peeling (`Peel`), augmentation (`Augment`), slack-aware peeling (`SlackPeel`) |
| `src/npa.cpp` | NPA |
| `src/graph.h`, `src/graph.cpp` | graph container, threshold ⌊k^τ⌋, size bound (Lemma 1), loader, validity check |
| `src/main.cpp` | command-line driver |


