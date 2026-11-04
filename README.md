# Exact Maximum Flexi-Clique Computation

Source code, datasets, and experiment scripts of the paper

> **Exact Maximum Flexi-Clique Computation under Non-Hereditary Feasibility and Connectivity Constraints**
> Song Kim, Kaiqiang Yu, Taejoon Han, Hyewon Kim, Junghoon Kim, Susik Yoon, Jungeun Kim.
> Submitted to IEEE ICDE 2027.

Given a graph *G* and a parameter τ ∈ (0, 1), a node set *H* is a **flexi-clique** if *G*[*H*] is
connected and every node of *H* has at least ⌊|*H*|^τ⌋ neighbours in *H*. The binary `flexi`
computes a maximum flexi-clique with:

| Algorithm | Description | 
|---|---|---|
| `eba` | **EBA** (Exact Branch-and-bound Algorithm): exact search with connectivity-preserving branching (CP-Branching), Rules 1–3, and the margin pivot, starting from the FPA solution | 
| `fpa` | **FPA** (Flexi-Prune Algorithm): scalable heuristic — global peeling, augmentation, and slack-aware peeling |  
| `npa` | **NPA**: the existing peeling heuristic for flexi-cliques (Kim et al., CIKM 2024), the baseline of FPA | 

## Contents

| Path | Content |
|---|---|
| `src/` | C++ implementation of EBA, FPA, and NPA, and the variants of the experiments |
| `baselines/` | the other baselines: CP-SAT (`cpsat/`), the query-centred solver of the case study (`qcs/`), and instructions for IterQC (`iterqc/`) |
| `graphs/` | the six smaller real-world graphs
| `experiments/` | scripts that run the experiments and write their results to `results/` |


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
`/usr/bin/time`). The conversion of the two LAW graphs needs Java 11 or later.

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

## Input format

A plain-text edge list with one edge `u v` per line, where `u` and `v` are non-negative integer node ids
(32-bit range). Lines that are empty or start with `#` are ignored, and so are self-loops; duplicate
edges are merged. The graph is treated as undirected and simple, and node ids need not be contiguous.

## Datasets

The paper uses 12 real-world graphs from KONECT, SNAP, LAW, and Chakraborty et al., the SNAP
ground-truth communities of Amazon and DBLP

- The six smaller graphs (PG, AP, CD, CO, DB, AM) are included in `datasets/graphs/`.
- The five larger graphs are downloaded from their sources:
  [FD](https://snap.stanford.edu/data/sx-stackoverflow.html),
  [SO](https://snap.stanford.edu/data/sx-stackoverflow.html),
  [OK](https://snap.stanford.edu/data/com-Orkut.html),
  [LJ](https://snap.stanford.edu/data/soc-LiveJournal1.html) (SNAP), and
  [IT](https://law.di.unimi.it/webdata/it-2004/),
  [SK](https://law.di.unimi.it/webdata/sk-2005/) (LAW).


