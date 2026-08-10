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
  `--no-metrics`; `jaccard_similarity_ie`, the column the paper tells users to
  read, costs 8.5% more at N=512 and rising.
- **`sweep.py` never rotates bedtools** — it runs first in all three sweeps.
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

**Figure 3 Panel A (synthetic N-scaling, p=18, t=16).** Job 29656140, node
c531, `docs/data/cpp_vs_bedtools_t16_p18.csv`. N=512: bedtools 714.77 s,
hammock 70.32 s, **10.17×** (was wrongly 27.61×). hammock's own N=512 wall is
unchanged by the fix (71.65 s → 70.32 s) — the control confirming only the
bedtools arm moved. Gate 1 (monotonicity): passes at every N except N=2→N=4,
which is genuine n=3 noise on a ~15 ms/pair cell (CV 102%), not a new defect —
see `docs/seed-benchmark-methodology.md` item 6. Gates 2–6 not yet run as
automated checks (done by hand for this rerun); still worth building as code.

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

**Still not re-measured**: Figs S6/S7/S8 (`sweep.py` path), `RESULTS.md` in
both experiment directories, and the bedtools-own-optimum (t=8) reframing
that used to accompany the Panel A table. Do not quote any of these until
they go through `run_bedtools()` and get an entry here.
