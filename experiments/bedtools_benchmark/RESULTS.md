# bedtools_benchmark — May 12 results (subB comparison)

Four runs against `hammock-cpp`, on Rockfish `shared` partition (sr-class
nodes), 16 cpus, 32 GB:

> **Every timing on this page is a 3-column run** — no `jaccard_similarity_ie`,
> no containment/cosketch block. These predate v0.7.0, when the binary emitted
> 3 columns unless `--metrics` was passed; since 0.7.0 the block is the default
> and the harnesses pass `--no-metrics` explicitly to keep new timings on this
> same footing. The block costs a union plus two cardinality estimates per pair,
> which lands entirely in the pairwise phase — 0.61% of wall at N=512, p=14, so
> the effect on these numbers would be about +1.5%. Timings from a run *with*
> the block are not comparable to the tables below.

- **2026-05-10 (morning)** — pre-optimization hammock-cpp, sequential pre-sort.
- **2026-05-10 (evening)** — same sweeps with **parallel** pre-sort and the
  first round of inner-loop optimizations (~3.1× speedup).
- **2026-05-11 (evening)** — second round of optimizations (inlined xxhash,
  removed `thread_local std::string`, incremental ASCII counter). Another
  ~2.1× wall-time win. subB=1.0 throughout.
- **2026-05-12 (afternoon) — current canonical** — same hammock binary as
  May-11 evening plus the new `--subB-method=mixed-stride` default and
  `hash32_short` (commit 9778ef8). **Adds subB ∈ {0.25, 0.1} as new
  hammock variants alongside the subB=1.0 baseline.**

Headline: **mixed-stride subsampling is essentially free in accuracy and
roughly halves hammock wall time per halving of subB.** Subsampling does
not change bedtools' results (it's a hammock-only knob). With subB=0.1,
hammock beats bedtools at **every N≥2** in the t=16 files sweep, with
**52.77× speedup at N=512**. Per-pair Jaccard error is statistically
indistinguishable across subB values — the definitional gap vs bedtools
stays at ~0.164. **Read that as a null result about the gap, not as evidence
that subsampling is unbiased** — see the precision-sweep caveat below.

## Source files

CSVs/text reports live in `results/` (symlinked to
`/vast/blangme2/jbonnie/hammock_claude_experiments/bedtools_benchmark/results/`).
PNGs in `figures/`.

**Current run (canonical, May 12):**

| Sweep         | Job ID      | CSV / report stem                            |
| ------------- | ----------- | -------------------------------------------- |
| precision     | 24057047    | `sweep_precision_20260512_150459`            |
| threads       | 24057048    | `sweep_threads_20260512_150458`              |
| intervals     | 24060777¹   | `sweep_intervals_20260512_160412`            |
| files (t=8)   | 24057050    | `cpp_vs_bedtools_t8_20260512_150458`         |
| files (t=16)  | 24060778¹   | `cpp_vs_bedtools_t16_20260512_160412`        |

¹ Resubmit — original jobs 24057049 (intervals) and 24057051 (files_t16)
were cancelled by a scheduler event at 15:29:19 about 25 minutes in
(simultaneous, different nodes; not a wall-time hit).

**Previous run (May 11, subB=1.0 only):** see git history for the table
of job IDs (`sweep_*_20260511_181918` / `cpp_vs_bedtools_*_20260511_181919`).

**Aug 4 files_t16 rerun (job 29552415, `cpp_vs_bedtools_t16_20260804_172242`,
archived to `docs/data/`).** Same protocol as the May 12 run plus a fourth
hammock arm, `hammock_ie_B` — subB=1.0 run *without* `--no-metrics`, i.e.
emitting `jaccard_similarity_ie` and the containment/cosketch block. It is the
only run on disk that measures the configuration `CLAUDE.md` recommends; the
block costs **1.45% of wall at N=512** (54.15 s vs 53.38 s) and 34.8 MB peak RSS
vs 22.7 MB. Everything else reproduces May 12 to within the drift documented
below: N=512 gives 12.68× / 26.05× / **52.35×** for subB 1.0 / 0.25 / 0.1
against 13.14× / 26.69× / 52.77×. **The tables in this file have not been
restated onto it** — the May 12 run is still canonical here. Whether it stays
canonical is pending the Figure 3 decision in `docs/figure3-candidate-v2.md`.

The precision sweep also dumps `*_pairs.csv` — per-pair bedtools and
hammock jaccards across all 4096 pairs × 5 precisions × 3 subB × 3 runs
= 184,320 rows — for the faceted scatter plot in `figures/`.

## files sweep at t=16 (the headline comparison)

p=14, 10k intervals/file, N up to 512, 3 runs/config:

| N   | bedtools | hammock @ subB=1.0 | @ subB=0.25 | @ subB=0.1 | speedup (bt/hm@0.1) |
| --- | -------- | ------------------ | ----------- | ---------- | ------------------- |
|   2 |   0.80 s |             0.29 s |      0.18 s |     0.11 s |  **7.15×** |
|   4 |   0.49 s |             0.47 s |      0.26 s |     0.17 s |  **2.91×** |
|   8 |   0.62 s |             0.89 s |      0.49 s |     0.30 s |  **2.08×** |
|  16 |   1.10 s |             1.72 s |      0.89 s |     0.51 s |  **2.16×** |
|  32 |   3.09 s |             3.39 s |      1.69 s |     0.92 s |  **3.36×** |
|  64 |  11.04 s |             6.76 s |      3.30 s |     1.73 s |  **6.40×** |
| 128 |  44.36 s |            13.63 s |      6.59 s |     3.28 s | **13.54×** |
| 256 | 176.46 s |            26.80 s |     13.23 s |     6.87 s | **25.70×** |
| 512 | 706.09 s |            53.72 s |     26.45 s |    13.38 s | **52.77×** |

Persistent-crossover summary (smallest N where hammock wins and keeps winning):

| subB | persistent crossover | speedup there | speedup at N=512 |
| ---- | -------------------- | ------------- | ----------------- |
| 1.0  |                 N=64 |         1.63× |             13.14× |
| 0.25 |                 N≤2 |         4.55× |             26.69× |
| 0.1  |                 N≤2 |         7.15× |             52.77× |

For subB=1.0 hammock loses at N∈{4,8,16,32} before winning persistently
at N=64. For subB ∈ {0.25, 0.1} hammock wins at every tested N — at
subB=0.25 the dip is at N=16 (1.24×) and at subB=0.1 the dip is at N=8
(2.08×). Halving subB roughly halves hammock walltime at every N, so
each subB step doubles the speedup over bedtools.

## subB perf scaling (intervals sweep)

16 files, p=14, t=8 — the regime where Mode B is hash-bound:

| N intervals/file | bedtools | subB=1.0 | subB=0.25 | subB=0.1 | hm@0.1 / hm@1.0 |
| ---------------- | -------- | -------- | --------- | -------- | ---------------- |
|        1 000 |   1.39 s |   0.34 s |    0.17 s |   0.09 s | **0.26×** (3.8× faster) |
|       10 000 |   1.05 s |   3.12 s |    1.46 s |   0.67 s | **0.21×** |
|      100 000 |   4.13 s |  30.48 s |   14.08 s |   6.30 s | **0.21×** |
|    1 000 000 |  34.22 s | 302.83 s |  139.63 s |  62.01 s | **0.20×** |

At 1M intervals/file, subB=0.1 brings hammock from 302 s (8.85× slower
than bedtools) down to 62 s (1.81× slower) — and that's at N=16 where
bedtools' O(N²) is trivial (256 pairs). At larger N the hammock advantage
asymptotes regardless.

The factor-2.15 speedup per subB-halving step is consistent across all
sweeps: see the threads / precision / files tables. mixed-stride's per-
position cost is roughly linear in `1/stride`, so we'd expect a clean
≈subB-proportional speedup, which is what we see.

> **Note on subB=0.5.** The 2026-05-11 subB_mixed_stride experiment
> recorded that at **subB=0.5** mixed-stride was actually slightly
> *slower* than subB=1.0 (the per-position gate cost doesn't yet
> amortize). We skipped subB=0.5 in this benchmark; the next informative
> step below 1.0 is subB=0.25, which is firmly in the "real speedup"
> regime (≈2.15× faster than 1.0 across all sweeps).

## threads sweep — subB doesn't change scaling shape

64 files × 10k intervals, p=14:

| t  | bedtools | subB=1.0 | subB=0.25 | subB=0.1 |
| -- | -------- | -------- | --------- | -------- |
|  1 |  58.11 s |  72.52 s |   40.87 s |  16.74 s |
|  2 |  38.20 s |  42.41 s |   21.07 s |   8.81 s |
|  4 |  19.63 s |  24.01 s |   10.88 s |   4.70 s |
|  8 |  10.26 s |  12.45 s |    5.80 s |   2.64 s |
| 16 |   9.90 s |   6.63 s |    3.22 s |   1.60 s |

OpenMP scaling efficiency is preserved across subB values (each line is
near-linear). At t=16 the speedup over bedtools goes 1.49× → 3.07× →
**6.19×** as subB drops from 1.0 → 0.25 → 0.1 (this is the N=64 threads
sweep). Compare with the files sweep at the same (N=64, t=16) point,
which lands at 1.63× / 3.34× / **6.40×** — same shape, slightly larger
numbers because the threads sweep regenerates data per run and absorbs
slightly more wall-time variance.

## precision sweep — accuracy is subB-independent

64 files × 10k intervals, t=8. Wall and MAE per subB:

| p  | wall@1.0 | wall@0.25 | wall@0.1 | MAE vs bedtools (1.0/0.25/0.1) | MAE vs hammock@p=18 (1.0/0.25/0.1) |
| -- | -------- | --------- | -------- | ------------------------------- | ----------------------------------- |
| 10 |  12.33 s |    5.69 s |   2.57 s | 0.1654 / 0.1651 / 0.1642        | 0.0110 / 0.0110 / 0.0110           |
| 12 |  12.36 s |    5.72 s |   2.58 s | 0.1663 / 0.1644 / 0.1638        | 0.0057 / 0.0055 / 0.0056           |
| 14 |  12.47 s |    5.81 s |   2.65 s | 0.1641 / 0.1646 / 0.1638        | 0.0026 / 0.0028 / 0.0028           |
| 16 |  13.53 s |    6.29 s |   2.96 s | 0.1639 / 0.1643 / 0.1644        | 0.0012 / 0.0015 / 0.0015           |
| 18 |  15.95 s |    7.61 s |   3.93 s | 0.1640 / 0.1641 / 0.1642        | 0.0000 / 0.0009 / 0.0010           |

- **MAE vs bedtools** (the bp-set-Jaccard definitional gap) is identical
  to within run noise across all 15 (p, subB) cells. The table is right; the
  original reading of it — "confirming that mixed-stride subsampling does not
  introduce systematic bias" — is not. This column cannot detect that. It is
  dominated by the register-equality chance floor `c`, which is a step
  function of the load factor λ = n/m: flat at 0.1699 for λ ≳ 5 and
  insensitive to everything else. Every file here has the same 10k intervals,
  so λ is uniform across pairs and far above the knee at every p in the
  table; the pairs are also low-J, and `MAE ≈ c·(1 − J̄)` lands just under `c`,
  which is the ~0.164 observed. A subsampling bias of a few 10⁻³ would be
  invisible underneath it. The column that *can* answer the bias question is
  `jaccard_ie_mae_vs_bt`, added to `sweep.py` after these runs; this table
  predates it. See `docs/jaccard-definitional-gap.md`.
- **MAE vs hammock@p=18, subB=1.0** halves per +2 in p (textbook HLL),
  identically for subB=1.0 and within 5×10⁻⁴ for subB ∈ {0.25, 0.1}.
  The tiny residual ~10⁻³ at p=18 for subB<1.0 is the
  different-accepted-position-set noise relative to the subB=1.0
  reference, not HLL noise — both still converge to the same set-Jaccard
  estimate, just from sampling slightly different positions.

The pairs-scatter plot is faceted across the three subB values; all
three panels are visually identical to the eye (same point cluster,
same median gap of 0.164). Subsampling at subB ≥ 0.1 is statistically
indistinguishable from the full-bp sketch **as measured against the p=18
hammock reference** (the `MAE vs hammock@p=18` column, which is a like-for-like
comparison and does carry that information). The bedtools column does not
support the same conclusion, for the reason above.

## Plotting

```bash
ml r/4.3.0
Rscript experiments/bedtools_benchmark/make_graphs.R \
  --files-csv results/cpp_vs_bedtools_t16_20260512_160412.csv \
  --pairs-csv results/sweep_precision_20260512_150459_pairs.csv
```

Produces (multi-subB-aware):

- `*_sketch_compare_split.png` — bedtools vs hammock(sketch+compare),
  one hammock-sketching line per subB. Compare line is plotted once
  (it's subB-agnostic — fixed-cost per pair regardless of sketch
  contents). Persistent-crossover annotation marks the subB=1.0
  crossover; subB<1.0 crossovers are visually obvious from the plot.
- `*_cost_per_pair.png` — wall/N² vs N, one hammock-total and one
  hammock-sketching line per subB.
- `*_jaccard_scatter.png` — facetted by subB, colored by precision.
  Same definitional-gap cluster in every panel.
- `*_jaccard_delta.png` — facetted by subB; per-precision running mean
  of (hammock − bedtools).

## Caveats addressed in earlier runs

- Sort time is captured separately and parallelized (post-May-10-evening).
- `find_hammock_cpp` picks by mtime, not lexicographically (so the
  cp310 wheel is preferred over a stale cp38 one).

## Related docs

- `docs/jaccard-definitional-gap.md` — why hammock jaccard ≠ bedtools jaccard
- `docs/bedtools-parallelism-caveat.md` — the GNU-parallel framing and sort-time fairness
- `docs/mode-b-subsampling-perf.md` — design notes for the inner-loop
  optimizations and the `--subB-method` flag
