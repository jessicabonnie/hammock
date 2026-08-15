# The precision frontier on Maurano (Figure 3, Panel B)

Job 29651773, node c718, one exclusive allocation, 2026-08-09. Data:
`docs/data/sweep_precision_maurano_p18_t16.csv` and `..._t8.csv`
(`sweep.py --axis precision --corpus maurano`, 20 real fetal DHS files passed as
both operands, 3 replicates, means).

**400 ordered pairs are timed; 190 unique off-diagonal pairs carry the accuracy.** The 20 self-pairs are
excluded from MAE — both tools return ~1.0 there at zero error, which is free
correctness that dilutes MAE ~5% without measuring anything. They are still part
of the timed workload, which is why the two counts differ.

## Results (t=16; t=8 in the same file, accuracy identical by construction)

| p | 12 | 14 | 16 | 18 | 20 | 22 | 24 |
|---|---|---|---|---|---|---|---|
| IE MAE vs exact bedtools | 8.80e-3 | 4.01e-3 | 2.68e-3 | **1.152e-3** | 5.89e-4 | 2.78e-4 | 1.54e-4 |
| register-equality MAE | 0.1374 | 0.1383 | 0.1380 | 0.1379 | 0.1377 | 0.1378 | 0.1376 |
| wall (s) | 5.42 | 5.46 | 5.70 | 6.13 | 8.17 | 14.16 | 35.61 |
| speedup vs bedtools, t=16 | 1.15× | 1.14× | 1.10× | **1.02×** | 0.76× | 0.44× | 0.18× |
| speedup vs bedtools, t=8 | 1.03× | 1.02× | 0.98× | **0.92×** | 0.73× | 0.47× | 0.21× |
| 9-column block cost | 1.00× | 1.00× | 1.00× | 1.00× | 1.00× | 1.02× | 1.02× |

## Four things worth stating

**1. p=18 reproduces the archived value exactly.** 1.151658×10⁻³ against
1.1517×10⁻³ computed independently from `docs/data/hammock_hll_p18_jaccB.csv`
(renamed to `..._full.csv` in the metrics-column restructure,
`docs/seed-metrics-column-restructure.md` — same file content, tag added).
But note what that does and does not show: `jaccard_ie` is *exactly* symmetric
(0 of 190 pairs asymmetric, exact comparison), so the MAE over 190 unique pairs
is **algebraically guaranteed** to equal the former MAE over 380 ordered rows. It confirms
the corpus, the self-pair filter and the column parse — a count-and-identity
check, not an independent statistical validation. Do not cite it as the latter.

**2. IE MAE tracks HLL theory.** Ratios per two steps of p: 2.19, 1.50, 2.32,
1.95, 2.12, 1.80 against a theoretical 2.00.

**3. Register-equality MAE is flat at ~0.138 for every p from 12 to 24.** The
definitional gap does not shrink with precision, which is the clearest available
demonstration that `jaccard_similarity` is not on bedtools' scale. This is a
better argument for the "read `jaccard_similarity_ie`" rule than anything in the
prose, because it is visible in one row.

**4. hammock is at parity or losing on this corpus, and that is expected.**
Maurano is N=20, far below the N≈64 crossover; Panel A measures 14.75× at N=256
on the same code. Panel B's p=18 point and Panel C's first bar are the same
measurement. Two caveats that both cut against us and must be stated:

- The bedtools leg here also ran at ~1.4× parallel efficiency (see
  `docs/bedtools-parallelism-caveat.md`). A properly parallel baseline would
  make these speedups **worse**, not better. The numbers above are conservative.
- The 9-column block looks free (1.00–1.02×) only because N=20 is
  sketching-dominated. At N=64 it costs 1.59–1.79× on `pair_time`. Do not
  generalise the 1.00× across N.
