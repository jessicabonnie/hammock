# Validation record

The exact-feature and independent-rehash validation completed on 2026-08-21
without changing production code, compatibility defaults, or manuscript text.
The isolated production source commit was
`008e9f24653cf4f40d10305c3b56ef29979ea6ae`.

## Execution and gates

- Exact gate: Slurm job `30131672` (array task 11, `k=10,w=30`).
- Remaining exact grid: Slurm array `30131710`.
- Primary p=18 sweep: Slurm array `30131720`, 296 jobs.
- Frozen-leader precision follow-up: Slurm array `30133142`, 24 jobs.
- Exact matrices: 37/37 cells, each with 20 samples and 400 ordered pairs.
- Rehash matrices: 320/320 outputs passed the named-column, finite-value,
  boundedness, symmetry, and unit-diagonal checks; reported clamps were zero.
- Scheduler terminal-state audits found no failed, cancelled, timed-out,
  out-of-memory, or node-failure tasks.

The unique exact-feature leader and the unique eight-seed p=18 robustness
leader were both `k=10,w=30`.

## Historical control

The archived historical grid
`paper/sequence_tissue_clustering/muscle_merged_p18_sensitivity.csv` has SHA-256
`21bc87d2efc92d718df4756d2422598a01edc06189d43fb0a11e12162c0f8164`.
Its historical `k=10,w=30` ARI (`0.9101796407185628`) and NMI
(`0.9609882289630555`) agree with the independently computed exact-feature
scores (up to the final displayed floating-point digit for NMI).

## Precision result

| p | seeds | ARI min | ARI median | ARI max | NMI min | NMI median | NMI max | exact partition |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 8 | 0.641703 | 0.766154 | 0.910180 | 0.886126 | 0.911405 | 0.960988 | 3/8 |
| 18 | 8 | 0.766154 | 0.910180 | 0.910180 | 0.911405 | 0.960988 | 0.960988 | 6/8 |
| 22 | 8 | 0.766154 | 0.910180 | 0.910180 | 0.911405 | 0.960988 | 0.960988 | 5/8 |
| 24 | 8 | 0.766154 | 0.766154 | 0.910180 | 0.911405 | 0.911405 | 0.960988 | 2/8 |

The final classification is **exact but estimator-sensitive**. At p=24,
seeds 1 and 3 recover the exact ARI and partition; the other six prespecified
seeds yield ARI `0.7661538461538462`. This sensitivity persists even though
median exact-Jaccard MAE falls with precision, consistent with the very small
exact ten-cluster cut gap (`4.106835699915767e-06`).

## Independent clustering cross-check

R 4.3.0 independently reconstructed average-linkage, ten-cluster partitions
for exact `k=10,w=30` and seed 42 at p=12,18,22,24. Every case agreed with the
Python/SciPy result for all 400 pairwise co-clustering decisions. The generated
tables and machine-readable comparison are in `r_crosscheck/`.

The 9.3 GiB exact feature binaries, raw matrices, completion metadata, job
manifests, and per-run clustering tables remain in the shared execution
workspace. The curated Git artifact set retains the final aggregate tables,
decision and validation reports, R agreement summary, and figures.
