# Validation record

The exact-feature and independent-rehash validation completed on 2026-08-21,
with the exploratory 101-seed extension completed on 2026-08-24,
without changing production code, compatibility defaults, or manuscript text.
The isolated production source commit was
`008e9f24653cf4f40d10305c3b56ef29979ea6ae`.

## Execution and gates

- Exact gate: Slurm job `30131672` (array task 11, `k=10,w=30`).
- Remaining exact grid: Slurm array `30131710`.
- Primary p=18 sweep: Slurm array `30131720`, 296 jobs.
- Frozen-leader precision follow-up: Slurm array `30133142`, 24 jobs.
- Seed-0 extension gate: Slurm array `30134192`, 4 jobs.
- Remaining seed extension: Slurm array `30134268`, 368 jobs.
- Missing-precision gate: Slurm array `30220609`, 5 jobs.
- Remaining p=13–17,19–21,23 interpolation: Slurm array `30220873`,
  904 fingerprint-locked jobs.
- Exact matrices: 37/37 cells, each with 20 samples and 400 ordered pairs.
- Rehash matrices: 1,601/1,601 outputs passed the named-column, finite-value,
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

| p | exact partition | median ARI | partitions | ranked hierarchies | unranked topologies | median exact MAE |
|---:|---:|---:|---:|---:|---:|---:|
| 12 | 44/101 | 0.766154 | 9 | 101 | 90 | 0.006067553 |
| 13 | 52/101 | 0.910180 | 5 | 101 | 60 | 0.004157942 |
| 14 | 56/101 | 0.910180 | 3 | 98 | 32 | 0.002885017 |
| 15 | 50/101 | 0.766154 | 2 | 82 | 17 | 0.002044305 |
| 16 | 54/101 | 0.910180 | 2 | 36 | 3 | 0.001326408 |
| 17 | 61/101 | 0.910180 | 2 | 18 | 4 | 0.000808426 |
| 18 | 52/101 | 0.910180 | 2 | 7 | 2 | 0.000531063 |
| 19 | 52/101 | 0.910180 | 2 | 4 | 2 | 0.000353679 |
| 20 | 61/101 | 0.910180 | 2 | 4 | 2 | 0.000244677 |
| 21 | 58/101 | 0.910180 | 2 | 4 | 2 | 0.000174753 |
| 22 | 54/101 | 0.910180 | 2 | 4 | 2 | 0.000119898 |
| 23 | 41/101 | 0.766154 | 2 | 3 | 2 | 0.000085429 |
| 24 | 44/101 | 0.766154 | 2 | 2 | 1 | 0.000061938 |

The final classification is **exact but estimator-sensitive**. From p=15
onward, seeds produce only two ten-class partitions and corresponding ARI
values (`0.9101796407185628` or `0.7661538461538462`), but the frequency of the
exact state is nonmonotonic. Paired Cochran's Q detects overall heterogeneity
(Q=21.788, 12 df, p=0.03997), while no adjacent exact McNemar comparison passes
Holm correction (minimum adjusted p=0.2299). Median flips therefore reflect
whether a bistable recovery frequency crosses 50%, not a smooth ARI response.

The coarse two-state result masks strong hierarchical convergence. From p=15
to p=24, distinct ranked hierarchies fall from 82 to 2 and unranked topologies
from 17 to 1, while median cophenetic correlation against exact rises from
0.99757 to 0.999999 and median exact-Jaccard MAE falls from 0.002044 to
0.00006194. In same-seed adjacent comparisons, 70/101 seeds retain their ARI
while changing ranked hierarchy from p=12 to p=13; 68/101 also change unranked
topology. Both counts fall to 1/101 from p=23 to p=24. Thus ARI/NMI conceal
precision-dependent reorganization at low and intermediate precision even
though the ten-class state remains seed-sensitive near the very small exact
cut gap (`4.106835699915767e-06`).

## Eight-class muscle-merged sensitivity

The secondary endpoint collapses `fMuscle_arm`, `fMuscle_back`, and
`fMuscle_leg` into one `fMuscle` label and cuts each same linkage tree into
eight groups. Exact `k=10,w=30` scores are ARI `0.861952` and NMI `0.945740`.
Exact eight-cluster partition recovery rises from 48/101 at p=12 to 60/101,
64/101, 75/101, 86/101, and 95/101 at p=13–17, then reaches 101/101 at every
precision p=18–24. Correspondingly, all 101 seeds have identical eight-class
ARI and NMI from p=18 onward, despite continuing changes in ranked merge order.

Paired Cochran's Q detects the expected overall recovery difference
(Q=393.500, 12 df, p=9.00e-77), while no adjacent exact McNemar comparison
passes Holm correction (minimum adjusted p=0.0886). The exact eight-cluster
cut gap is `0.001442`, roughly 351 times the ten-cluster gap. This sensitivity
analysis therefore localizes the persistent high-precision ten-class ARI split
to merge ranking around the narrow ten-cluster cut: the coarser muscle-merged
partition is fully stable from p=18, and by p=24 the unranked topology is exact
for every seed.

## Independent clustering cross-check

R 4.3.0 independently reconstructed average-linkage, ten-cluster partitions
for exact `k=10,w=30` and seed 42 at p=12,18,22,24. Every case agreed with the
Python/SciPy result for all 400 pairwise co-clustering decisions. The generated
tables and machine-readable comparison are in `r_crosscheck/`.

The 9.3 GiB exact feature binaries, raw matrices, completion metadata, job
manifests, and per-run clustering tables remain in the shared execution
workspace. The curated Git artifact set retains the final aggregate tables,
decision and validation reports, R agreement summary, and figures.
