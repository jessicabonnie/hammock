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
