# Seed: the performance benchmarks are not cross-run comparable — OPEN

Handoff note for a fresh chat. Written 2026-08-08.

**Status: the defect is established and quantified; the corrective job has RUN
and settled the headline. The eleven open items below are still open.**

Job `29628907` completed 2026-08-09 (sr09, Xeon Gold 6448Y, 16 reserved cores,
59m28s, 140 runs, seed 20260808). Data: `docs/data/fusion_ab_20260808_232422.csv`.

**The control worked.** `post/pre` on `--no-metrics` — byte-identical code —
landed at **0.96–1.02** on `pair_time` and 0.99–1.00 on wall, so the
experiment's own measurement error is ±2–4%. Contrast the cross-run comparison
this replaces, where the same control read 0.65–0.79.

**Settled result** (within-run, paired, one machine, error bar above):

| p | 12 | 14 | 16 | 18 | 20 | 22 | 24 |
|---|---|---|---|---|---|---|---|
| fusion speedup, metrics `pair_time` | 1.60× | 1.63× | 1.66× | 1.67× | 1.72× | 1.77× | 1.78× |
| block cost **before** (metrics÷no_metrics) | 2.90× | 2.65× | 2.67× | 2.65× | 2.69× | 2.98× | 2.97× |
| block cost **after** | 1.79× | 1.62× | 1.61× | 1.59× | 1.64× | 1.67× | 1.67× |

So: **the fusion made the block's pairwise phase 1.6–1.8× faster, rising with
precision, and cut the block's cost multiplier from ≈2.7–3.0× to ≈1.6–1.8×.**

This *supersedes* the "2.09× at p=24" figure derived from the 08-04/08-08
comparison, which was inflated by that run's 1.34× machine factor; the honest
number is 1.78×. It also supersedes the "≈2.5×" pre-fusion multiplier — measured
on one machine with a control, the pre-fusion multiplier is ≈2.7–3.0×.

On `comparison_time` (pair + serial write) the same effect reads 0.90 → 0.56
from p=12 to p=24: at low precision the unchanged `fprintf` loop dominates and
dilutes it, which is open item 3 in a different guise.

Everything here is evidence and open questions. The only *decisions* taken are
the four listed under "What was already changed".

## The finding

**Every performance number in this repo that was obtained by comparing one
benchmark run to another is confounded by an unmeasured machine factor of up to
1.5×, and none of the archived CSVs record enough to detect it.**

This was found while re-measuring the cost of the similarity block after the
pairwise passes were fused (`826d7b4`). The re-measurement compared a
2026-08-04 run against a 2026-08-08 run and reported that the estimator
multiplier "fell from ≈2.5× to ≈1.65×".

The `--no-metrics` arm is the control that breaks it. That arm executes
**byte-identical code** in both binaries — the fused path sits behind
`if (!args.metrics)` in `cpp/app/hammock_cli.cpp`, and `omp_set_num_threads` was
present in both. It should have been 1.00. Measured `new/old` on `pair_time`:

| p | 12 | 14 | 16 | 18 | 20 | 22 | 24 |
|---|---|---|---|---|---|---|---|
| ratio | 0.731 | 0.788 | 0.653 | 0.753 | 0.785 | 0.712 | 0.747 |

The 08-04 machine was **1.27–1.53× slower on unchanged code**. `sacct` shows no
SLURM job at the time of that run — it was taken on a shared dev node with no
reserved cores. An independent check pins the same effect: the serial `fprintf`
write loop costs 0.3036 µs/row in the 08-04 *allocated* files sweep and 0.298
µs/row in the 08-08 allocated run (stable to 1.8%), but **0.427 µs/row** in the
08-04 unallocated run — 43% inflation of a loop that does nothing but format
doubles.

**Consequence:** any absolute `pair_time`, `us_per_pair`, `comparison_time`, or
wall-percentage quoted across two runs is invalid. Deflating the headline p=24
figure by the 1.34× machine factor leaves ≈1.56× of the claimed 2.09×.

## What survives, and why

The **within-run ratio** (metrics ÷ no_metrics) is a paired contrast and is
immune to a uniform machine factor. It reproduces on two independent run-pairs
measured on three machines:

| | pairwise-cost | files sweep |
|---|---|---|
| 08-04 | 2.53 | 2.69 |
| 08-08 | 1.61 | 1.63 |

So "the block got roughly 1.5× cheaper in the pairwise phase" is supported. The
3-significant-figure endpoints are not: these quantities show ±6% machine
sensitivity and 3–17% within-cell CV at the extremes. Quote **≈2.5–2.7 → ≈1.6–1.8**.

## What was already changed (decisions, not open questions)

1. **Provenance in the CSVs** — `hostname`, `cpu_model`, `cpu_count`
   (affinity-aware), `slurm_job_id`, `git_sha`, `binary_version`, `corpus_seed`.
   Added because the 08-04 files sweep records none of it and therefore cannot
   be checked for comparability at all. `cpu_model` is load-bearing:
   `-march=native` is baked in (`CMakeLists.txt:42,47,68`).
2. **`--corpus-seed`** (`benchmark_cpp_vs_bedtools.py`,
   `pairwise_cost_by_precision.py`). Unseeded remains the default because every
   archived CSV was produced that way. `run_i` is mixed into the derived seed in
   the files sweep so replicates still draw *different* corpora — seeding them
   identically would silently redefine `std_wall_time` from "machine noise plus
   corpus variation" to machine noise alone.
3. **The subB sweep's `--metrics` arm now records metrics.** It previously read
   only `jaccard_similarity` and the pair de-duplication rebuilt each row from
   three keys, so `sweep_maurano_ie_*.csv` held jaccards **bit-identical** to
   the no-metrics arm across all 17,100 cells, under a filename containing "ie".
4. **`fusion_ab.py` + `sbatch_fusion_ab.sh`** — the corrective design, below.
5. **`bedtools.sh` is pinned to `bedtools/2.30.0`** (2026-08-09), the load is
   fatal on failure, and the resolved version + path are echoed to stderr. It
   previously read `ml bedtools2 … || ml bedtools … || true`, resolving to
   **2.27.1** (2017) while `estimator_compare.py` hardcoded 2.30.0 — two parts of
   one experiment scored against different bedtools — and the trailing `|| true`
   made a total module failure silent.

   **This is a ground-truth correctness fix, not version hygiene.** Run the full
   400-pair Maurano workload under both builds and 93 pairs disagree. The
   `intersection` and `n_intersections` columns are identical; only the **union**
   differs, and **2.27.1 is order-dependent** — `jaccard -a A -b B` and
   `-a B -b A` return different unions, with the swapped orientation agreeing
   with 2.30.0. So under 2.27.1, **93 of the 190 unordered Maurano pairs have
   J(A,B) ≠ J(B,A)** in what is supposed to be the *exact* reference. Under
   2.30.0: 0.

   That bites hardest exactly where the fixed-corpus sweep operates, since it
   passes one list as both operands and therefore scores both orientations of
   every pair against a reference that was not symmetric while hammock's
   estimate is.

   Magnitude is small — max |Δ| 1.3×10⁻⁵, mean 5.6×10⁻⁷ — so **no archived
   conclusion changes**: against an IE MAE of 1.15×10⁻³ at p=18 the mean shift is
   0.05%, and the p=18 gate still reproduces (1.15166×10⁻³ under 2.30.0 vs
   1.151647×10⁻³ under 2.27.1). But it is **6.4% of the MAE at p=23 in the worst
   single pair**, so do not wave it away at the high-precision end of a frontier
   plot. hammock's own per-pair output is bit-identical across the two runs,
   which is what makes the attribution clean.

   Consequence for archived data: `docs/data/maurano_bedtools_ref.tsv` and every
   bedtools column in the pre-2026-08-09 sweeps carry the asymmetric-union
   artifact. Fine as a gate to ~4 significant figures; do not use them to argue
   about differences at the 10⁻⁵ level.

## The corrective design (job 29628907)

Four arms — `{pre, post} × {--metrics, --no-metrics}` — one seeded corpus, one
allocation, one node, arm order permuted per replicate. All reported quantities
are within-replicate ratios.

Grounded in standard practice for two-build comparisons:

- **Paired + order-randomized.** Measure the difference on identical input;
  randomize order so cache/CPU state does not align with arm. See
  <https://www.bazhenov.me/posts/paired-benchmarking/>.
- **Setup randomization** as the mitigation for measurement bias — Mytkowicz et
  al., ASPLOS 2009, <https://dl.acm.org/doi/10.1145/1508284.1508275>. This repo
  hit precisely their failure mode.
- **Effect size + interval, not point estimates** — *Rigorous Benchmarking in
  Reasonable Time*.
- **Duet** (both binaries concurrent) was considered and rejected:
  <https://arxiv.org/pdf/2001.05811>. It equalizes interference by construction
  and is better where spare cores exist, but both arms size their OpenMP team to
  the whole allocation, so concurrent execution would measure contention.

**The control arm is the point.** `pre`/`post` on `--no-metrics` is identical
code and must land at 1.00; its departure is the error bar on everything else in
the job. Reporting an effect without that calibration is what allowed a machine
factor to be published as a speedup.

Binaries: built from `826d7b4^` (pre) and `HEAD` (post) with the **same**
compiler and cmake flags, in `build/ab/` (`/tmp` is node-local and invisible to
compute nodes). The post-fusion cmake binary's md5 equals the pip-built one, so
the flags agree. Verified: the two binaries emit **byte-identical CSVs at
p=12/16/20/24 on both output shapes** — which tests bit-exactness old-vs-new,
something the 56-CSV gate could not do since its baseline came from the new code.

## What job 29628907 settles, and what it does not

Settles: the fusion effect on `pair_time` and `comparison_time`, free of the
machine factor, with a calibration arm.

**Does not settle — still open:**

1. **The N×p interaction is unmeasured, and a claim was made about it anyway.**
   The two benchmarks form a plus-sign, not a grid: p varies only at N=64, N
   varies only at p=14. `pairwise_cost_by_precision.py`'s docstring asserts
   µs/pair "flattens above N≈16"; the 08-08 files sweep at p=14 contradicts it —
   µs/pair falls 42% from N=16 to N=512 and is still falling, and the
   metrics multiplier **rises monotonically** 1.505 → 1.759 over that range,
   above the top of the "1.61–1.79" band, which is an N=64 band.
   **Needed:** an N ∈ {64, 512} × p ∈ {14, 20, 24} block. Fix the docstring
   either way.
2. **The peak-RSS argument is null by construction — drop it, do not repeat it.**
   Observed peak RSS fits `2N·2^p + T·2^p + baseline` exactly (2309.8 MiB
   predicted vs 2309.734 measured at p=24, N=64, T=16). The `T·2^p` term is the
   per-thread `HLLSketch` allocated inside the *sketching* OpenMP region
   (`cpp/src/processing_modes.cpp:142,249`). It is exactly co-sized with the
   pre-fusion per-pair union scratch and lives in a **different, non-overlapping
   phase**, so `max(RSS)` is identical whether or not the union allocation
   exists. A pre-fusion binary would report the same number. Measuring it needs
   a phase-resolved counter (`getrusage` around the pairwise phase, or a malloc
   hook), not `/usr/bin/time -v`.
3. **The residual `+IE` cost is now mostly `fprintf`, and the paper says
   otherwise.** At N=512: Δpair 0.1422 s vs Δwrite 0.3609 s — **71.7%** of the
   delta is the serial `%.17g` write of six extra columns, up from 46.9%
   pre-fusion (the fusion cut the pair component 2.8× and left write untouched).
   `paper/outline.md` says "a factor of 3.4–3.8 within the comparison phase";
   the new numbers are **2.67–2.89×** for comparison phase (pair+write) and
   **1.63–1.76×** for the pair phase alone. Different denominators — do not swap
   them. If the write loop matters, that is a separate, precision-independent,
   Θ(N²) optimization.
4. **The Maurano leg is underpowered and its arms are blocked, not interleaved.**
   `sbatch_maurano.sh` runs `--corpus maurano` then `--corpus maurano --metrics`
   as two invocations ~12 min apart, so arm is fully confounded with position.
   Predicted effect there is **+0.05% to +0.28%** of wall (the pairwise phase is
   0.073–0.396% of wall on that corpus); observed spread across all 18 cells is
   **−1.62% to +0.65%**, i.e. 6–16× the effect and centred on the wrong sign,
   with several cells "significant" in the wrong direction at 3–9σ. "No
   measurable cost on Maurano" is a foregone conclusion, not a finding.
   **Either** interleave the arms in one loop, **or** stop making the wall-time
   claim and quote the instrumented `Δ(pair+write)/wall` instead.
5. **Job-level n=1.** All three 08-08 benchmarks ran once. The n=5/n=3
   replicates share one allocation and one node microstate, so the
   between-job/between-node variance component — the one every cross-run claim
   depends on — is unestimated. Effective n for "08-04 vs 08-08" is 1 vs 1.
6. **The files sweep reports means, not medians, at n=3.** One cold run poisons
   a cell (bedtools at N=2: 0.772 ± 0.475, CV 61.5%). `pairwise_cost` uses
   medians of 5. The two get quoted side by side with neither labelled.
7. **Arm ordering in the files sweep cannot balance.** `_rotate(arms, run_i)`
   with 4 arms and 3 runs gives summed positions 6 vs 5 for `hammock_ie_B` vs
   `hammock_cpp_B`, i.e. the IE arm sits on average ⅓ of a position later, for
   an effect of ~0.95%. Use `--runs` a multiple of the arm count. Separately,
   **bedtools always runs first**, never rotated, always immediately after the
   parallel pre-sort — its position is confounded with the warmest page cache.
8. **Why bedtools appeared to speed up 3.1–5.1% between the two runs is still
   unexplained.** It is *not* the corpus: simulated over 12 draws, total corpus
   size has CV 0.013% and a pair-sum cost CV 0.035% — 88× too small. It shows in
   **CPU** time (−2.63%), which excludes core-stealing and page cache. Excluded:
   bedtools version (only 2.27.1 available). Remaining: node identity (sr14 →
   sr09) and memory-bandwidth co-tenancy, indistinguishable because the 08-04
   report predates the `cpu_model` field. **Needed:** a fixed CPU-bound
   calibration kernel timed once per job, so every future run carries its own
   machine yardstick.
9. **Provenance holes that cannot be closed retroactively.** Run 2 records
   `git_sha 67d336e-dirty`. Runs 1 and 3 carry no `git_sha` (the column landed
   after them). All three are unseeded. Worst: **both binaries on disk have
   mtime after all three jobs finished**, because the extension was rebuilt for
   an unrelated crash fix at 20:31 — so neither binary on disk is provably the
   one benchmarked. `fusion_ab.py` records a binary md5 per row for this reason;
   consider doing the same everywhere.
10. **`sweep.py` was not updated** — no `--corpus-seed`, no provenance columns,
    though it calls the same generator and its CSVs are what
    `docs/bedtools-parallelism-caveat.md` and `RESULTS.md` cite for cross-run
    comparisons.
11. **Θ(2^p) is stale.** CLAUDE.md quotes "1021× against 1024× predicted"; on
    the 08-08 medians p=14→24 is **967.7×** (−5.5%).

## Where things are

- Harness: `experiments/bedtools_benchmark/fusion_ab.py`,
  `sbatch_fusion_ab.sh`; binaries in `build/ab/bin_pre`, `build/ab/bin_post`
  (gitignored; rebuild from `826d7b4^` and `HEAD` if missing — the recipe is in
  the sbatch header).
- New data (08-08, all allocated): `pairwise_cost_by_precision_20260808_183458.csv`
  (job 29620391, sr09), `cpp_vs_bedtools_t16_20260808_190441.csv` (29620392,
  sr09), `experiments/subB_mixed_stride/results/sweep_maurano_{,ie_,bedtools_}20260808_*.csv`
  (29620619, sr06).
- Archived/suspect: `docs/data/pairwise_cost_by_precision_20260804_164807.csv`
  (**no allocation — do not compare absolutes against it**),
  `docs/data/cpp_vs_bedtools_t16_20260804_172242.csv` (job 29552415, sr14,
  allocated, but no `cpu_model`).

## Documents still carrying superseded numbers

Not yet restated, pending the corrected measurement:
`experiments/bedtools_benchmark/RESULTS.md` (the `+1.45%` / 34.8 vs 22.7 MB /
"3.4–3.8×" block), `docs/figure3-panel-a-rebuild.md`,
`docs/metrics-by-default.md` (both tables, and its "union plus two cardinality
estimates" mechanism sentence), `paper/outline.md` (Figure 3 Panel A prose),
CLAUDE.md's "flat ≈2.5×" paragraph (already flagged superseded in place), and
`hammock-cpp --help` line 54 ("Skips a union plus …"). CLAUDE.md divergence #4's
subB wall medians will also move once the Maurano leg is restated.

**Open decision belonging to the user, not to a future session:** whether Figure
3 Panel B should be redrawn on the `--metrics` arm (the configuration CLAUDE.md
recommends) rather than `--no-metrics`. The data to decide now exists on both
arms; see open question 4 for why the *wall-time* half of that data is not yet
trustworthy.
