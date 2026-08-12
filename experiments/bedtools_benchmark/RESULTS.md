> # ⛔ EVERY BEDTOOLS-RELATIVE NUMBER BELOW IS SUPERSEDED (2026-08-09)
>
> **Do not quote any speedup-over-bedtools figure from this file.** Four defects
> were found in the benchmark harness in one day, **all inflating the bedtools
> baseline, i.e. all in hammock's favour**:
>
> 1. `bedtools.sh` ran three processes per pair instead of one (~2.1×).
> 2. bedtools was pinned to 2.27.1, whose jaccard union is order-dependent.
> 3. `ml` module loads ran inside the timed region.
> 4. **The bedtools module exports `LD_LIBRARY_PATH` with 17 directories, of
>    which bedtools uses 4. The dynamic linker searches the other 13 — on GPFS —
>    on every `execve`, and a pairwise workflow is N² execs.** This alone
>    inflated bedtools 2.4–2.8×.
>
> Measured effect at N=512: bedtools 1978 s → 716 s. **The headline speedup falls
> from 27.61× to roughly 10.2×**, and that replacement is itself provisional —
> it rests on a single job and has not been through the checks below.
>
> Also retracted: the "process creation caps near **123 exec/s** and does not
> scale with cores" mechanism. Measured 16-way: **364 exec/s** clean, 196 slow.
> The `md5sum` control that appeared to confirm it ran in the same polluted
> environment, so it confirmed nothing.
>
> **What survives:** anything hammock-internal — mixed-stride vs hash-threshold,
> the precision sweeps, the fused-pass A/B, Mode D threading. Those divide two
> hammock numbers measured in one run and are unaffected.
>
> Status and worklist: `docs/bedtools-baseline-retraction.md`.

# bedtools_benchmark — subB comparison

Script inventory, flags, and the SLURM wrappers are in
[README.md](README.md); this file is results only.

**Mixed provenance.** The precision, threads, intervals and t=8 sweeps are
from May 12 2026; the t=16 files sweep was re-run on Aug 4 2026 and its
section is the only one restated onto that run. Each section says which.

> ## Every bedtools number in this file is affected by a defect found 2026-08-09
>
> A pairwise bedtools workflow launches one process per pair -- N^2 of them --
> and on these nodes process creation caps near **123 exec/s and does not scale
> with cores**. Every sweep below therefore compares hammock at its stated
> thread count against a bedtools that achieved somewhere between ~1x and ~5.7x,
> unrecorded and varying by node and day. A speedup quoted from this file can be
> inflated by up to ~6x.
>
> Two things were fixed after these runs, so nothing below reproduces exactly:
>
> * `bedtools.sh` now runs **one** process per pair instead of three, which is
>   worth ~2.1x to the baseline. Every bedtools wall time here is from the slow
>   version and is correspondingly too high.
> * `bedtools` is pinned to 2.30.0. The 2.27.1 build used here computes an
>   **order-dependent union**: 93 of 190 unordered Maurano pairs had
>   J(A,B) != J(B,A) in what is supposed to be the exact reference. Magnitude is
>   small (max 1.3e-5) so no conclusion below changes, but do not use these
>   files to argue about differences at the 1e-5 level.
>
> New bedtools rows carry `mean_bedtools_parallel_eff`; these do not. Prefer the
> p=18 re-runs and `docs/bedtools-parallelism-caveat.md` over this file for any
> speed claim.

Four runs against `hammock-cpp`, on Rockfish `shared` partition (sr-class
nodes), 16 cpus, 32 GB:

> **Every timing on this page is a 3-column run** — no `jaccard_similarity_ie`,
> no containment/cosketch block. These predate v0.7.0, when the binary emitted
> 3 columns unless `--metrics` was passed; from 0.7.0 through the pre-restructure
> era the block was the default and the harnesses passed `--no-metrics`
> explicitly to keep new timings on this same footing. Since the metrics-column
> restructure (`docs/seed-metrics-column-restructure.md`) the harnesses pass
> `--register-equality` instead — same reduced-column arm, renamed flag; see
> `benchmark_cpp_vs_bedtools.py`'s `_metrics_off_flag`, which resolves the
> correct flag per binary. The block costs a union plus two cardinality
> estimates per pair,
> which lands entirely in the pairwise phase — 0.61% of wall at N=512, p=14.
> That prediction has since been measured directly: the `hammock_ie_B` arm of
> the Aug 4 files t=16 run costs **+1.45%** of wall at N=512. Timings from a run
> *with* the block are not comparable to the tables below.

- **2026-05-10 (morning)** — pre-optimization hammock-cpp, sequential pre-sort.
- **2026-05-10 (evening)** — same sweeps with **parallel** pre-sort and the
  first round of inner-loop optimizations (~3.1× speedup).
- **2026-05-11 (evening)** — second round of optimizations (inlined xxhash,
  removed `thread_local std::string`, incremental ASCII counter). Another
  ~2.1× wall-time win. subB=1.0 throughout.
- **2026-05-12 (afternoon) — canonical for the precision, threads, intervals
  and t=8 sweeps** — same hammock binary as May-11 evening plus the new
  `--subB-method=mixed-stride` default and `hash32_short` (commit 9778ef8).
  **Adds subB ∈ {0.25, 0.1} as new hammock variants alongside the subB=1.0
  baseline.**
- **2026-08-04 — canonical for the files t=16 sweep only** (job 29552415,
  hammock-cpp 0.7.0). Same protocol, plus a fourth arm `hammock_ie_B` — subB=1.0
  emitting the full similarity block — and microsecond phase timers in place of
  the old integer-millisecond ones. This is the run behind Figure 3 Panel A; see
  `docs/figure3-panel-a-rebuild.md`. **The other four sweeps were not re-run and
  their tables below are still May 12.** Do not read a table here as coming from
  one job unless its section says so.

Headline: **mixed-stride subsampling is essentially free in accuracy and
substantially reduces hammock wall time.** (This previously read "roughly
halves hammock wall time per halving of subB", which the data does not support:
a 10x subsample buys **4.21x** on synthetic at p=14 and **1.83x** on Maurano at
p=18, not 10x. The per-step factor is neither 2 nor stable across corpora.) Subsampling does
not change bedtools' results (it's a hammock-only knob). With subB=0.1,
hammock beats bedtools at **every N≥2** in the t=16 files sweep, with
**52.35× speedup at N=512** (52.77× in the superseded May 12 run). Per-pair Jaccard error is statistically
indistinguishable across subB values — the definitional gap vs bedtools
stays at ~0.164. **Read that as a null result about the gap, not as evidence
that subsampling is unbiased** — see the precision-sweep caveat below.

## Source files

CSVs/text reports live in `results/` (symlinked to
`/vast/blangme2/jbonnie/hammock_claude_experiments/bedtools_benchmark/results/`).
PNGs in `figures/`.

**Current runs:**

| Sweep         | Job ID      | CSV / report stem                            |
| ------------- | ----------- | -------------------------------------------- |
| precision     | 24057047    | `sweep_precision_20260512_150459`            |
| threads       | 24057048    | `sweep_threads_20260512_150458`              |
| intervals     | 24060777¹   | `sweep_intervals_20260512_160412`            |
| files (t=8)   | 24057050    | `cpp_vs_bedtools_t8_20260512_150458`         |
| files (t=16)  | 29552415    | `cpp_vs_bedtools_t16_20260804_172242` ²      |
| files (t=16), superseded | 24060778¹ | `cpp_vs_bedtools_t16_20260512_160412` |

² Aug 4 re-run; the only sweep here not from May 12. Archived to
`docs/data/`. Everything else in this table is May 12.

¹ Resubmit — original jobs 24057049 (intervals) and 24057051 (files_t16)
were cancelled by a scheduler event at 15:29:19 about 25 minutes in
(simultaneous, different nodes; not a wall-time hit).

**Previous run (May 11, subB=1.0 only):** see git history for the table
of job IDs (`sweep_*_20260511_181918` / `cpp_vs_bedtools_*_20260511_181919`).

**The Aug 4 files_t16 re-run** adds a fourth hammock arm, `hammock_ie_B` —
subB=1.0 run *without* `--no-metrics`, i.e. emitting `jaccard_similarity_ie` and
the containment/cosketch block. It is the only run on disk that measures the
configuration `CLAUDE.md` recommends: the block costs **1.45% of wall at N=512**
(54.15 s vs 53.38 s, 12.50× over bedtools against 12.68×) and 34.8 MB peak RSS
against 22.7 MB. Everything else reproduces May 12 inside the drift documented
above — N=512 gives 12.68× / 26.05× / **52.35×** for subB 1.0 / 0.25 / 0.1
against 13.14× / 26.69× / 52.77× — so the qualitative claims are unchanged.
The files t=16 section below has been restated onto it; nothing else has.

The precision sweep also dumps `*_pairs.csv` — per-pair bedtools and
hammock jaccards across all 4096 pairs × 5 precisions × 3 subB × 3 runs
= 184,320 rows — for the faceted scatter plot in `figures/`.

## files sweep at t=16 (the headline comparison)

p=14, 10k intervals/file, N up to 512, 3 runs/config:

| N   | bedtools | hammock @ subB=1.0 | @ subB=0.25 | @ subB=0.1 | speedup (bt/hm@0.1) |
| --- | -------- | ------------------ | ----------- | ---------- | ------------------- |
|   2 |   0.82 s |             0.22 s |      0.11 s |     0.06 s | **12.66×** |
|   4 |   0.52 s |             0.43 s |      0.21 s |     0.11 s |  **4.93×** |
|   8 |   0.64 s |             0.84 s |      0.41 s |     0.21 s |  **3.07×** |
|  16 |   1.12 s |             1.67 s |      0.81 s |     0.40 s |  **2.77×** |
|  32 |   3.04 s |             3.34 s |      1.62 s |     0.80 s |  **3.80×** |
|  64 |  10.69 s |             6.67 s |      3.23 s |     1.60 s |  **6.69×** |
| 128 |  41.10 s |            13.29 s |      6.45 s |     3.19 s | **12.89×** |
| 256 | 166.46 s |            26.62 s |     12.92 s |     6.41 s | **25.98×** |
| 512 | 676.62 s |            53.38 s |     25.98 s |    12.92 s | **52.35×** |

The `hammock_ie_B` arm (subB=1.0, full similarity block) is in the same CSV and
tracks subB=1.0 within 1.5%: 54.15 s at N=512 against 53.38 s.

Persistent-crossover summary (smallest N where hammock wins and keeps winning):

| subB | persistent crossover | speedup there | speedup at N=512 |
| ---- | -------------------- | ------------- | ----------------- |
| 1.0  |                 N=64 |         1.60× |             12.68× |
| 0.25 |                 N≤2 |         7.65× |             26.05× |
| 0.1  |                 N≤2 |        12.66× |             52.35× |

For subB=1.0 hammock loses at N∈{8,16,32} before winning persistently
at N=64. For subB ∈ {0.25, 0.1} hammock wins at every tested N. All three
curves dip at N=16 — 0.67× / 1.37× / 2.77× — which is the worst point for
hammock at every subB: bedtools' O(N²) is still cheap there while hammock is
already paying full sketch construction on 32 files. Halving subB roughly halves hammock walltime at every N, so
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
which lands at 1.60× / 3.31× / **6.69×** — same shape, slightly larger
numbers because the threads sweep regenerates data per run and absorbs
slightly more wall-time variance. Note the two sweeps are now from different
dates (threads May 12, files Aug 4), so a few percent of that gap is
between-day drift rather than the effect being described.

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
  --files-csv results/cpp_vs_bedtools_t16_20260804_172242.csv \
  --pairs-csv results/sweep_precision_20260512_150459_pairs.csv
# (optional) --out-dir <dir>; default is figures/
```

The files CSV above is the Aug 4 run, i.e. the one the files t=16 section is
restated onto. Swap in `cpp_vs_bedtools_t16_20260512_160412.csv` to regenerate
the superseded May 12 panels. `hammock_ie_B` is filtered out of the files
plots by `make_graphs.R`'s `^hammock_cpp_B` tool match, so adding that arm did
not change the figures' content — only the numbers behind them.

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
