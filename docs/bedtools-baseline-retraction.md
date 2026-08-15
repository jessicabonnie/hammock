# The bedtools baseline retraction (2026-08-09)

Status document. Every bedtools-relative performance number in this repo is
suspect until it appears in the "re-measured and gated" list at the bottom.

## What happened

Four defects were found in the benchmark harness in a single day. **All four
inflated the bedtools baseline; all four therefore flattered hammock.** That
one-sidedness is itself the diagnostic — it is what asymmetric scrutiny looks
like from the outside, and it should be assumed to be the cause until shown
otherwise.

| # | defect | effect on bedtools | fixed in |
|---|---|---|---|
| 1 | `bedtools.sh` ran three processes per pair (`bedtools \| cut \| awk`) | ~2.1× slower | `3df44ec` |
| 2 | pinned to bedtools 2.27.1, order-dependent jaccard union (93/190 Maurano pairs asymmetric) | correctness, not speed | `44a244c` |
| 3 | `ml bedtools`, `ml parallel`, two `--version` execs inside the timed region | fixed ~0.7 s/cell | `73c2312` |
| 4 | **module `LD_LIBRARY_PATH` = 17 dirs, 4 used; 13 unused GPFS dirs searched by `ld.so` on every one of N² execs** | **2.4–2.8× slower** | `d93ea20` |

Measured at N=512 (262,144 pairs): bedtools **1978 s → 716 s**. hammock over the
same change reads 71.65 s → 70.33 s, which is the control saying the node is
comparable and only the bedtools arm moved.

**Headline: 27.61× → roughly 10.2×.** The replacement is provisional. It comes
from one job and has not passed the gates below.

## Retracted mechanisms, not just numbers

- **"Process creation caps near 123 exec/s and does not scale with cores."**
  Measured 16-way: **364 exec/s** with a clean environment, 196 with the module
  path. The documented cap was ~3× low.
- **The `md5sum` control that appeared to prove this was tool-independent ran
  inside the same polluted environment.** A control that shares the confound
  confirms nothing. This is the single most important process lesson here: the
  control agreed, and that agreement was taken as verification.
- **The "page-cache warming / bedtools always runs first" explanation** for
  bedtools costing *less* at N=4 than at N=2. Arm rotation was added to fix it.
  Rotation cannot move a fixed floor, and the impossible pattern survived the
  fix — in `docs/data/cpp_vs_bedtools_t16_p18.csv`, which was then published.

## The recurring failure

**Mechanisms were verified; magnitudes were not.** Defect 3 was diagnosed from a
1-pair vs 4-pair test giving a ~0.7 s floor. Nobody asked whether 0.7 s could
account for a 1262 s gap at N=512. It cannot, by three orders of magnitude, and
that arithmetic was available at the moment the claim was made. Defect 4 was
found only because someone later checked a fix against its own predicted effect.

Corollary: **a fix must be checked against the size of the thing it claims to
explain**, not merely shown to be real.

## Gates that must pass before any bedtools-relative number is quoted again

None of these are "be more careful"; all are mechanical.

1. **Impossible-pattern scan** over the CSV: cost must not fall as work rises;
   parallel efficiency must not exceed the core count; no negative times. This
   check, run for the first time during this incident, flagged exactly three
   archived files — all of them this bug at different stages.
2. **Floor vs signal.** Measure each arm's zero-work cost explicitly
   (`harness_floor_ms`) and refuse to publish a ratio for any cell where the
   floor exceeds ~5% of either arm.
3. **Self-ratio control.** Run the baseline twice on identical input with arm
   order swapped; it must read ~1.00 before any cross-tool ratio is trusted.
   This is the design used for the fused-pass A/B, which is the one experiment
   in the repo that was built correctly. Make it the default shape.
4. **Environment capture.** Record `LD_LIBRARY_PATH`, `PATH`, `module list` and
   an env digest per row. Defect 4 would have shown up as a provenance diff.
5. **Unused-search-path check.** Diff `ldd $(command -v bedtools)` against the
   exported `LD_LIBRARY_PATH` and fail on unused directories. Two lines.
6. **Baseline symmetry check.** `bedtools jaccard -a A -b B` must equal
   `-a B -b A`. Would have caught defect 2 immediately.

## Still open, and load-bearing

- **The missing control: an exact, single-process, batch baseline.** hammock's
  advantage is Θ(N) sketching vs Θ(N²) file re-parsing. A batch tool that reads
  each BED once and does N² exact sorted-merge intersections would recover most
  of that **with zero approximation error**. Until it is run, the speedup is
  evidence for *workflow architecture*, not for HyperLogLog — any fixed-size
  sketch would show the same curve. This is the single most important missing
  experiment.
- **Speed and accuracy are never simultaneously measured at one operating
  point.** The speed headline is p=18 on synthetic, where no accuracy against
  bedtools was ever computed (though bedtools runs in that harness and the exact
  answer is a free byproduct). The accuracy headline is p=21 on Maurano, where
  hammock is at parity or slower.
- **The timed configuration is not the recommended one.** Headlines run
  `--no-metrics` (now `--register-equality`/`--re`, renamed by the v0.8.0
  metrics-restructure); `jaccard_similarity_ie`, the column the paper tells
  users to read, costs 8.5% more at N=512 and rising (job 29671317,
  `docs/data/cpp_vs_bedtools_t16_p18.csv`: 78.08 s `hammock_ie_B` vs 70.97 s
  `hammock_cpp_B`, actually a 10.0% delta at that one cell — see
  `paper/outline.md` §4.4 for the "8.5%"/"under 1% up to N=32" framing this
  number feeds; the two figures were never reconciled and that mismatch
  predates Part 9, not introduced by it).

  **Part 9 update, 2026-08-11 (`docs/seed-metrics-column-restructure.md`):**
  that 78.08 s was measured with `--metrics` (the pre-restructure-equivalent
  full 8-column block) — at the time, the *only* way to get
  `jaccard_similarity_ie` at all. Since v0.8.0 the bare/no-flag default gives
  that column directly for less write cost (still pays the same fused
  union-pass compute as `--metrics` — see `run_hammock`'s docstring in
  `benchmark_cpp_vs_bedtools.py` — so this is a write-cost effect only, not a
  compute one). `benchmark_cpp_vs_bedtools.py`'s `--metrics-arm`/`--metrics-all`
  are now retargeted to the bare default (`run_hammock(..., ie_only=True)`
  instead of `metrics=True`) so the *next* rerun of job 29671317's config will
  measure the corrected number automatically.

  **Correction to the paragraph above, same day: this is NOT a harmless
  secondary line.** An earlier draft of this note claimed "Figure 3 itself
  plots the register-equality arm, untouched by this change" — that is
  wrong. Since commit `c024b9e` (2026-08-10, "Simplify Figure 3 to
  jaccard_similarity_ie only"), Panel A's plotted hammock curve, and its
  headline "hammock is 8.4× faster than BEDTools" number (`paper/outline.md`,
  `paper/draft.md`), is computed *exclusively* from the `hammock_ie_B` row —
  `plot_pairwise_scaling.R` hard-fails (`stop()`) if that row is absent. The
  register-equality arm (`hammock_cpp_B`, 70.97 s, "9.2×") moved to
  Supplementary Figure S9 and is no longer what the main text reports. So the
  arm Part 9 retargeted (`hammock_ie_B`, from `--metrics` to the bare
  default) **is the literal source of Figure 3A's headline number**, not a
  side quantity — get this right before deciding whether to skip a
  regeneration.

  Having said that, the *magnitude* is small enough that Part 1's
  "harness-runtime, not published-number" framing holds up on the numbers,
  just not on the mechanism: `hammock_ie_B`'s wall time (78.08 s) is
  `hammock_cpp_B`'s (70.97 s, unaffected by Part 9) plus a ~10.0% overhead —
  see below for the real, measured post-restructure figure.

  **Real dedicated re-run, 2026-08-11/12 (job 29769994, node c700, same CPU
  model — Xeon Gold 6248R — as the archived job 29671317's c529, exclusive
  allocation, same `sbatch_fig3_panelA_v2.sh` config: two passes, N∈{2,4,8,16,32}
  at 20 reps and N∈{64,128,256,512} at 3 reps).** This supersedes the
  estimate that used to sit in this paragraph — user explicitly asked for
  dedicated compute rather than leaving it a projection. Full table
  (`experiments/bedtools_benchmark/results/cpp_vs_bedtools_t16_20260811_231858.csv`
  + `..._232623.csv`, not promoted to `docs/data/` — see below for why),
  monotonic at every N (gate 1 passes):

  | N | bedtools | `hammock_cpp_B` (re) | `hammock_ie_B` (ie, retargeted) |
  |---|---|---|---|
  | 64 | 9.58 s | 7.73 s | 7.80 s |
  | 128 | 37.98 s | 15.74 s | 16.08 s |
  | 256 | 154.48 s | 32.90 s | 34.84 s |
  | 512 | 634.61 s | 71.29 s | 78.13 s |

  **Isolated Part 9 effect at N=512** (same job, same node — the only valid
  same-run comparison): `hammock_ie_B`/`hammock_cpp_B` overhead is
  **9.60%** (78.13/71.29), down from the archived **10.02%**
  (78.08/70.97) — real, measured, not estimated; close to the ~9.0%
  projection this section used to carry. Recomputing the headline holding
  the archived, carefully-corrected bedtools baseline fixed (653.41 s, per
  the user's explicit choice to isolate Part 9's effect rather than adopt a
  wholesale new dataset — see next paragraph): 653.41 / 78.13 ≈ **8.36×**,
  still rounds to the published **8.4×**. No change made to
  `paper/outline.md`/`paper/draft.md` or the actual figure — confirmed
  stable, not merely projected as stable.

  **This run's own bedtools number (634.61 s) is 2.9% below the archived
  653.41 s — do not attribute that to Part 9.** It's a different node
  (c700 vs. c529, same CPU model) on a different day; this repo has
  documented bedtools node-to-node spread up to 28% before (see "Small-N
  cells corrected," above) — a few percent here is unremarkable noise, not
  a regression or an improvement. Folding it into the headline (634.61/78.13
  = 8.12×, "8.1×") would conflate an unrelated bedtools redraw with the
  metrics-restructure's effect — exactly the mistake this document exists to
  warn against. Per explicit user decision (2026-08-12), this run is **not**
  adopted as the new canonical Figure 3A dataset; `docs/data/cpp_vs_bedtools_t16_p18.csv`
  (job 29671317) remains the source of record, and only the isolated,
  bedtools-held-fixed calculation above is used to confirm the headline is
  unaffected. If a genuine reason arises later to re-baseline bedtools
  itself (not a Part 9 concern), that's a separate decision with its own
  writeup, not a byproduct of this one.

  `sweep.py`'s `sweep_precision` retarget was also confirmed this session by
  a real rerun (job 29769995) rather than left as a flagged estimate — see
  the Figure 3 Panel B / Supplementary Fig S8 entry below for that result.
- **`sweep.py`'s `sweep_precision` never rotates bedtools** — it runs first,
  once per replicate, before the per-precision loop. (Corrected: this used to
  say "in all three sweeps." `sweep_threads` and `sweep_intervals` were fixed
  2026-08-09 — bedtools and every hammock arm, including the `--metrics-arm`
  IE arm, now rotate via `_rotate(leg, run_i)`, the same helper and convention
  `benchmark_cpp_vs_bedtools.py` already used. `sweep_precision`'s structure
  is different — bedtools' cost doesn't depend on precision, so it runs once
  outside the swept loop rather than once per point — and rotating it needs
  its own fix, not yet done.)
- **`run_bedtools_maurano.py` bypasses the fixed `run_bedtools` wrapper**, so
  the Maurano baseline re-measured on 2026-08-09 (job 29652709) still carries
  defects 3 and 4.
- **The sort is charged to neither tool**, though bedtools requires it and
  hammock does not — 60% of bedtools' N=2 cell, 0.3% at N=512.
- **Peak RSS for bedtools is a structural undercount** (largest single child,
  not concurrent sum).
- **GIGGLE** — built for many-vs-many interval search — is absent from both the
  benchmark and the related work.

## Re-measured and gated

**Figure 3 Panel A (synthetic N-scaling, p=18, t=16).** Job **29671317**, node
**c529**, `docs/data/cpp_vs_bedtools_t16_p18.csv` (one exclusive allocation,
two passes: N∈{2,4,8,16,32} at 20 replicates, N∈{64,128,256,512} at 3).
N=512: bedtools 653.41 s, hammock 70.97 s, **9.21×**. hammock's own N=512 wall
has moved under 1% across the whole chain of reruns (71.65 s → 70.32 s →
70.97 s) — the control confirming only the bedtools arm's correction is doing
the work. Gate 1 (monotonicity): now passes at **every** N, including N=2→N=4
— see "small-N cells corrected" below for why that's new. Gates 2–6 not yet
run as automated checks (done by hand for this rerun); still worth building as
code.

**Small-N cells corrected (superseding the entry below this one).** This job
supersedes an intermediate rerun, job 29656140/node c531, which fixed the
LD_LIBRARY_PATH bug (moving N=512 from 27.61× to 10.17×) but whose own N=2/N=4
cells were noise-corrupted at n=3 — one bad bedtools replicate made the table
read an impossible cost *drop* from N=2 to N=4 (0.50 s → 0.18 s). A 20-rep,
single-node recheck (job 29670970, node c431) found bedtools winning 20/20 at
N=2, the opposite sign, with a properly monotonic cost curve — but that job's
node differed from 29656140's (c431 vs c531), and their one overlapping cell
(N=32) disagreed by 28%, too much to safely splice. Job 29671317 was designed
as a single-node, two-pass rerun specifically to fix this without introducing
a new node confound. Result: bedtools wins reliably (monotonic, tight spread
— e.g. N=32 spans 2.52–2.63 s across 20 reps, CV ~1%) at every N ≤ 32; the
crossover falls between N=32 (0.67×) and N=64 (1.28×). N=512 also moved
slightly in this rerun (714.77 s → 653.41 s bedtools, 10.17× → 9.21×) — a
normal amount of n=3 run-to-run noise on the unchanged high-N pass, not a
further correction.

**The `hammock total (+IE)` curve is no longer TEMPORARY.** It was added to
Figure 3A pending validation against this job; validation is done. At N=512
it reads **8.4×** faster than bedtools against **9.2×** for the default
column (comparison-phase cost 1.64–1.94× at N≥64, rising to 14.1% of total
wall by N=512).

**Figure 3 Panel B (Maurano, subB sweep, t=8).** `run_bedtools_maurano.py`
ported to call the fixed `run_bedtools()` wrapper instead of `run_with_time`
directly — it had none of the (c)/(d) fixes until this pass.
`docs/data/maurano_bedtools.csv`, `docs/data/maurano_subB_summary.csv`.
bedtools 7.07 s (was 7.34 s under the partial fix, 11.08 s originally).
hammock at subB=1.0 is **1.16× slower** (not 1.12×); subB=0.1 is 1.58× faster
(not 1.64×); subB=0.01 is 2.49× faster (not 2.58×). `docs/figures/maurano_speed_bars.png`
and `maurano_subB_pareto_scatter.png` regenerated from the corrected CSVs;
the bar-chart title and per-bar labels were hardcoded to "faster" and have
been changed to derive from the data's sign (same defect, different file,
as the one already fixed in `plot_pairwise_scaling.R`).

**Supplementary Fig S6 (threading, `paper/figures/threading_supplement.png`).**
Job 29670792, node c531, `docs/data/sweep_threads_p18.csv` (p=18, 64
files/side, 10k intervals, means of 3; replaces the pre-fix `sweep_threads_p18.csv`
generated 17:53, before the LD_LIBRARY_PATH commit at 22:26:45). BEDTools
plateaus by t=8 (4.22× speedup vs its own t=1, 6.4 of 16 cores) and stays flat
through t=48 (4.20×, no further gain, no strong regression); hammock keeps
improving to t=32 (12.11×, 24.0 of 32 cores busy) then eases back slightly at
t=48 (11.59×, 28.1 of 48), consistent with the oversubscription concern in
`docs/seed-mode-d-threading.md`'s sibling note for interval-mode pools.
BEDTools's achieved parallelism plateaus far below its thread count (6.5 of
16-48 cores), matching this doc's "process-creation-bound" finding, though
the mechanism is BEDTools' N² per-pair process model saturating, not a
literal 123 exec/s ceiling (that framing is separately retracted above).

**Supplementary Fig S8 (precision frontier, `paper/figures/precision_frontier.png`).**
Job 29670793, node c432, `docs/data/sweep_precision_maurano_p18_t16.csv` +
`..._t8.csv` (Maurano, 190 unique off-diagonal pairs, subB=1.0; the raw output contains both reciprocal directions; replaces the pre-fix pair
generated 18:00, job 29651773). At the CLI default p=18, t=16: hammock 6.41 s
vs bedtools 4.40 s reference (**0.69×**, i.e. slower — N=20 is below the N≈64
crossover, consistent with Panel B); at t=8, 0.81×. hammock does not beat
BEDTools at any precision in this sweep (0.08×–0.77× across p=12–24) — the
near-parity an earlier, pre-fix draft of this figure showed was itself an
artifact of the same LD_LIBRARY_PATH bug. `jaccard_similarity_ie` MAE falls
monotonically 0.0088 (p=12) → 0.00015 (p=24), a 57× improvement, while
register-equality `jaccard_similarity` MAE stays flat at ≈0.138 across the
whole range (chance-agreement floor, not sampling noise — precision cannot
fix it, hence why the frontier's x-axis is the IE column). Cross-checked
bit-identical against the t=8 CSV's own accuracy columns (max |ΔMAE| = 0, as

**Part 9 confirmation, 2026-08-11/12, job 29769995 (`sbatch_fig3_panelB.sh` —
misleadingly named: it generates this S8 precision-frontier data, not the
main-text Panel B bar chart; that one is a separate arm, see below).**
`sweep.py`'s `sweep_precision` untimed IE pass was retargeted from `--metrics`
to the bare default (Part 9). Re-ran this exact config (Maurano, p=12–24,
t=8 and t=16, 3 reps each, exclusive node c699) as a real dedicated rerun, not
an estimate. Both S8's inputs are unaffected in substance: (1) the
correctness gate — `jaccard_ie_mae_vs_bt` at p=18 must reproduce
`1.1517e-3` — passed exactly (`0.00115165752` at both t=8 and t=16, all 3
reps each), confirming `jaccard_similarity_ie` really is bit-identical
between the retargeted arm and the old `--metrics` arm, as the fused-pass
argument predicted rather than merely asserted. (2) The retargeted arm's own
overhead over the register-equality pass is 0.2–1.25% here (not the ~9-10%
seen at N=512 on synthetic data) — at N=20/400 timed ordered pairs (190 unique
off-diagonal accuracy observations) the sketch-construction
phase (~9 s) completely dwarfs the pairwise union pass (~0.015 s), so there
is essentially nothing for a write-cost saving to be a percentage *of*.
S8's plotted speedup axis reads off the register-equality arm's `wall_time`
(untouched by Part 9 either way), so this job changes nothing about S8's
published numbers — it is a clean, empirical no-op confirmation, not a
number update. `docs/data/sweep_precision_maurano_p18_t16.csv`/`_t8.csv` are
left as-is (new run's numbers agree with them, so overwriting would be
churn, not a correction). New run's own output:
`experiments/bedtools_benchmark/results/sweep_precision_20260811_231858{,_pairs}.csv`
(t=16) and `..._232719{,_pairs}.csv` (t=8), not promoted to `docs/data/`.

**The real main-text Panel B bar chart is a different arm and was already
handled — no rerun was needed for it.** Panel B (`plot_pairwise_scaling.R`'s
`panel_b`) reads `docs/data/maurano_subB_ie_summary.csv`
(`wall_median`/`mae_ie_vs_bedtools`), which comes from
`experiments/subB_mixed_stride/run_sweep.py`'s own `--metrics` arm (via
`sbatch_maurano.sh`), a completely separate script from `sweep.py`/
`sbatch_fig3_panelB.sh` above despite the similar name. That arm's Part 9
retarget was already measured directly (paired, 5 reps, exact p=18/t=8/subB=1
cell): ratio 0.997–1.003, within noise floor — see
`experiments/subB_mixed_stride/RESULTS.md`'s "`--metrics` IE-capture arm
retargeted" entry. So Figure 3's actual "1.87×" Panel B number was already
confirmed unaffected before this dedicated-compute round started; the
dedicated Panel-A/B reruns launched here targeted the headline N-scaling
curve (Panel A, still running as of this writing) and the S8 precision
frontier (this entry, done) — not a third rerun of Panel B, which didn't need
one.

**Original S8 prose below predates the Part 9 confirmation above; numbers
unchanged, cross-checked bit-identical against the t=8 CSV's own accuracy
columns (max |ΔMAE| = 0, as
it must be — accuracy doesn't depend on thread count).

**Supplementary Fig S7 (catalog scale, `paper/figures/largeN_supplement.png`)**
re-rendered against the corrected `cpp_vs_bedtools_t16_p18.csv`, now job
29671317; projected bedtools at N=2048 fell further, 30,093 s/72.1× (original,
pre-LD_LIBRARY_PATH-fix) → 10,463 s/25.1× (job 29656140, before its small-N
fix) → **10,294 s/24.7×** (current). The N=512 overlap check between this
figure's two source CSVs now reads 70.97 s vs 71.73 s (1.1%, was 0.1% against
the superseded job) — still within this cluster's documented node-to-node
spread. No further action needed here.

All three supplementary figures are now on the fixed harness, and Figure 3
Panel A's small-N cells are corrected (see above). **Still open**: `RESULTS.md`
in both experiment directories (prose, not regenerated), and the
bedtools-own-optimum (t=8) reframing that used to accompany the Panel A
table.
