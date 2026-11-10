# Verification of the optimality results

On many instances of EQ1, EBA finishes without branching: the global pruning (Rules 1 and 3 with the
FPA incumbent) removes every node, which proves the FPA solution optimal. The scripts in this
directory check these results without relying on the code of EBA, and compare FPA and EBA with
exhaustive enumeration on small graphs.

## Certificate

Let $F'$ be a flexi-clique of size $f$ and $\lambda=\lfloor (f+1)^\tau\rfloor$. A flexi-clique $X$
with $|X|>f$ is connected, and each of its nodes has at least $\lfloor |X|^\tau\rfloor\ge\lambda$
neighbours in $G[X]$. Hence $X$ lies in one connected component of the $\lambda$-core of $G$, and
that component has more than $f$ nodes. So **if every connected component of the $\lambda$-core
has at most $f$ nodes, no flexi-clique has more than $f$ nodes**, and $F'$ is optimal.

With incumbent $f$, this is exactly the condition under which Rules 1 and 3 remove every node
before branching ($\lambda=\theta_{F'}$ is the threshold of Rule 1, and Rule 3 removes the
components of the $\lambda$-core with at most $f$ nodes). The certified instances must therefore
be exactly those that EBA settles without branching, which `--results` checks.

## Files

| File | Content |
|---|---|
| `indep_check.cpp` | Edge-list loader, core decomposition, components of the $\lambda$-core, and checks of a node set (connectivity, minimum degree). It shares no code with `src/`. |
| `certify_optimality.py` | Runs FPA on every instance and checks four things with `indep_check`: both loaders read the same graph; the FPA solution is a flexi-clique of the reported size; every component of the $\lambda$-core has at most $f$ nodes; and, with `--results`, the outcome agrees with EBA (certified if and only if EBA reports 0 search nodes). |
| `oracle_families.py` | Compares FPA and EBA with exhaustive enumeration on small graphs of eight structures, and checks the returned node sets independently. |
| `data/oracle_families.txt` | Output of `oracle_families.py` with its default settings. |

All thresholds are computed exactly, with $\tau=p/q$ as a fraction and integer roots.

## Running

```bash
python3 verification/certify_optimality.py --workers 4 --results results/results.csv
python3 verification/oracle_families.py --workers 8
```

- `certify_optimality.py` compiles `indep_check` into `build/` and keeps the FPA solutions in
  `results/certificate/fpa/`; it writes `results/certificate/certificate.csv`, one row per instance,
  with the structure of each certified solution: a whole connected component of the $k$-core,
  where $k$ is its minimum degree, or a set whose size is the largest one allowed by its threshold.
- `--datasets` and `--taus` select a subset (default: the 13 graphs of Table II at
  $\tau\in\{0.5,0.7,0.9\}$).
- `--results` takes the results file of `experiments/run_experiments.py --groups main`.
- **SK** needs about 60 GB of memory per process.
- `oracle_families.py` reads **FB**, **CA**, and **DB** from `graphs/` for its subgraph families.
- Both scripts exit with status 1 if a check fails.

## Results

- **Certificate.** With the paper's runs, every FPA solution of EQ1 is a valid flexi-clique of the
  reported size, and the certificate holds exactly on the instances that EBA solves without
  branching, and on no other instance.
- **Exhaustive enumeration (`data/oracle_families.txt`).** 4,000 graphs and 20,000 runs (5 values
  of $\tau$ per graph), with no failed run. EBA always returns a maximum flexi-clique, and FPA
  always returns a valid one.
