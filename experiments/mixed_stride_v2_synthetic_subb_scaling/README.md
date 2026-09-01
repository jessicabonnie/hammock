# Synthetic subB scaling: time, memory, and accuracy

This experiment extends the Figure 3A synthetic benchmark across subB and
both hammock similarity outputs.

## Design

- 10,000 intervals per BED file;
- disjoint query/reference collections with the full N x N cross product;
- N = 2, 4, 8, 16, 32, 64, 128, 256, 512;
- 20 replicates through N=32 and three replicates from N=64 onward;
- corpus seed 20260810, matching the rebuilt Figure 3A experiment;
- p=18 and 16 threads;
- subB = 1, 0.1, 0.01, 0.001 using the public/default mixed-stride v2;
- register-equality and inclusion-exclusion (+IE) output shapes;
- exact BEDTools Jaccard on every generated corpus.

Each replicate generates and pre-sorts one corpus, then rotates BEDTools and
all eight hammock arms through the execution order. This ensures all timing,
peak-RSS, and MAE measurements within a replicate refer to exactly the same
files. BEDTools runs once per replicate, not once per subB or similarity type.

Completed replicate files are atomic and independently resumable under
`results/raw/`, so a timeout at a large N does not discard earlier cells.

## Outputs

- `results/raw/N*_rep*.csv`: one row per tool/shape/subB arm;
- `results/prepared/all_runs.csv`: validated combined replicate rows;
- `results/prepared/summary.csv`: per-cell aggregates;
- `figures/synthetic_subb_scaling.png`: time, peak RSS, and MAE by N, faceted
  by similarity output.

The canonical paper figures are not modified by this experiment.

## Completed run

The 2026-08-28 run used commit `ec2d14a`:

- benchmark job `30328362` completed on `c573` in 1h19m43s;
- finalizer job `30328363` completed in 6s;
- all 112 replicate checkpoints completed, yielding 1,008 validated tool
  rows and 81 aggregate cells;
- the benchmark job's orchestration process reached 599,164 KiB MaxRSS.

At N=512 (262,144 cross-product pairs), BEDTools' median wall time was
637.46s and median peak RSS was 19.95 MiB. Hammock peak RSS was approximately
266.4 MiB for every subB and both outputs because subsampling changes the
number of positions processed, not the fixed p=18 HLL allocation per file.

| subB | RE median (s) | +IE median (s) | RE MAE | +IE MAE |
|---:|---:|---:|---:|---:|
| 1 | 62.63 | 68.60 | 1.640e-1 | 1.197e-3 |
| 0.1 | 25.41 | 31.41 | 1.642e-1 | 1.204e-3 |
| 0.01 | 15.47 | 21.33 | 1.435e-1 | 1.197e-3 |
| 0.001 | 19.71 | 22.06 | 2.477e-2 | 8.563e-4 |

The fastest N=512 cell is therefore register-equality at subB=0.01 (41.2x
faster than BEDTools), not subB=0.001. At the sparsest rate, reduced sketch
work is offset by increased pairwise time. +IE remains close to BEDTools
accuracy across the full grid, whereas register-equality retains its expected
chance-agreement bias except at the sparsest occupancy.
