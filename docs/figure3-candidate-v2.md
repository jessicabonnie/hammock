# Figure 3 candidate (v2) — for side-by-side review

Nothing in **Panel A** is adopted. `plot_pairwise_scaling.R`'s Panel A code, the
CSV it reads, and the Panel A text in `paper/outline.md` are untouched; the
candidate is a parallel copy so the published figure and it can be opened
together.

Panel B is a separate matter and *was* changed in place, in both scripts: its
bar labels now read `mean |ΔJ|` and the `paper/outline.md:42` caption defines
the quantity. That is a labelling correction to the published figure, unrelated
to the Panel A question, and deliberately applied to v1 as well so the caption
never describes a label only the candidate carries.

| | script | figure | data |
|---|---|---|---|
| published | `paper/pairwise_scaling/plot_pairwise_scaling.R` | `paper/figures/pairwise_scaling.png` | `docs/data/cpp_vs_bedtools_t16_20260512_160412.csv` |
| candidate | `paper/pairwise_scaling/plot_pairwise_scaling_v2.R` | `paper/figures/pairwise_scaling_v2.png` | `docs/data/cpp_vs_bedtools_t16_20260804_172242.csv` |

Panel B is identical between them — its bar labels were relabelled to `mean |ΔJ|` in
**both** scripts at the same time, so the figure decision stays purely about Panel A.

## What changed, and why

**1. A fourth series: `hammock sketch comparison (+IE)`.**

Every timed hammock number the paper has ever published is a `--no-metrics`
run — the three-column output. `CLAUDE.md` tells readers to use
`jaccard_similarity_ie`, which is not in that output. So the figure showed the
cost of a configuration the documentation advises against, and there was no
published measurement of the recommended one.

The new series is the `hammock_ie_B` arm added to
`benchmark_cpp_vs_bedtools.py` in v0.7.0: subB=1.0, emitting
`jaccard_similarity_ie` and the containment/cosketch block. It was measured in
the same job as the series beside it, with arm order rotated by run index, so
the comparison is not confounded with position.

**It costs 1.45% of wall at N=512** (54.15 s against 53.38 s) and 34.8 MB peak
RSS against 22.7 MB. The pairwise *phase* is 2.5–3.8× more expensive; the phase
is 0.4–0.5% of wall, which is why the total barely moves. The `hammock total`
line does not visibly change, and that is the finding.

**2. The secondary axis counts N², not N(N−1)/2.**

The harness builds two **disjoint** sets of N files and both tools run the full
cross product — `benchmark_cpp_vs_bedtools.py` writes `set1_*`/`set2_*`,
`hammock_cli.cpp` loops n×m, `bedtools.sh` is `parallel ::: file1 ::: file2`.
N=512 is **262,144** pairs, not 130,816. `make_graphs.R:174` already had this
right; this figure and `paper/outline.md` did not.

Panel B is **not** affected: its "190 unique pairs" genuinely is N(N−1)/2,
because that run is a single within-set all-pairs comparison of 20 samples.

The `12.7× faster` annotation is a wall-time ratio with no pair-count term, so
it is unaffected by this fix — it moved from 13.1× because the data is a new
run, not because of the correction.

**3. The `pmax(comparison_time, 1e-4)` floor is gone.**

This is the most visible change and it is a bug fix, not a cosmetic one. The
binary used to report integer **milliseconds**, so `mean_comparison_time` was
literally `0.000000` for N ∈ {2,4,8,16} and the floor pinned those four points
to a flat line at 10⁻⁴. **Four of the nine points on that series were a
constant, not a measurement.** v0.7.0 reports microseconds, so the series is now
a real power law across its whole range.

The floor is replaced by a hard error: a zero there now means the CSV predates
the microsecond timers, which is worth failing on rather than papering over.

The published series is ragged for three separable reasons, and only the first
two are fixed by the timer:

- **Four of nine points were the floor**, as above.
- **The other five were quantized to 1 ms per run.** Multiply the archived means
  by the 3 runs and every one is an exact integer number of milliseconds — 19,
  26, 85, 288, 991 — so each mean sat on a ⅓-ms grid. At N=32 that is a 22% grid
  on a 1.5 ms quantity.
- **The archived run was also far noisier than resolution alone explains.**
  Coefficient of variation on `comparison_time` was 7.4% / **59.8%** / 20.2% /
  6.8% / 14.1% at N = 32…512, against 0.6% / 0.4% / 0.5% / 1.4% / **0.13%** in
  the rerun. At N=512 a 1 ms quantum is 0.3% of the value, so it cannot produce
  a 14% spread. The N=64 point is the worst offender: ±60% run-to-run, landing
  at 8.7 ms between neighbours now measured at 1.5 and 17.4 ms.

  Two things changed at once — a different node, and arm order now rotates where
  the old harness always ran subB=1.0 first, immediately after bedtools had
  saturated 16 threads for minutes. Both are consistent with the tightening and
  the archived CSV keeps only aggregates, so they cannot be separated after the
  fact. Do not claim the rotation caused it.

Note the rerun's series still bends at low N, and that is real, not residual
noise: per-pair cost falls 43.3 → 11.1 → 3.6 → … → 1.02 µs from N=2 to N=512,
asymptoting near 1 µs/pair as fixed per-run overhead (file open/close, thread
startup) stops dominating a sub-millisecond phase. It should not be a straight
line.

**4. The y axis (two pre-existing defects, both fixed).**

Both are present in the published figure and neither was introduced by the
changes above, though removing the floor makes them easier to see.

- **Breaks were two decades apart.** The data spans 6.6 decades, and over that
  range `scale_y_log10()`'s default `log_breaks()` picks breaks **100× apart** —
  `1e-4, 0.01, 1, 100`. So consecutive gridlines on the published figure are two
  decades, and the axis appears to step `0.1 → 10 → 1,000`. On a log axis that
  reads as a decade step unless you stop and check the numbers. The candidate
  pins `breaks = 10^(-4:3)`, one decade per gridline.
- **Labels collapsed below 0.1.** `label_number(accuracy = 0.1)` renders every
  break under 0.1 as `"0.0"` — the published figure labels 10⁻⁴ as `0.0`. The
  candidate uses three significant figures in fixed notation.

Together the axis now reads `0.0001 / 0.001 / 0.01 / 0.1 / 1 / 10 / 100 / 1,000`.

## The rerun reproduces the archived run

Job 29552415, 1:01:21 on `sr14`. Same protocol as the archived run: 16 threads,
p=14, 10k intervals/file, N = 2…512, subB ∈ {1.0, 0.25, 0.1}, 3 runs.

| N=512 wall | archived | rerun |
|---|---|---|
| BEDTools | 706.09 s | 676.62 s |
| hammock subB=1.0 | 53.72 s (13.14×) | 53.38 s (12.68×) |
| hammock subB=0.25 | 26.45 s (26.69×) | 25.98 s (26.05×) |
| hammock subB=0.1 | 13.38 s (**52.77×**) | 12.92 s (**52.35×**) |
| hammock subB=1.0 +IE | — | 54.15 s (12.50×) |

Every cell is inside the 1.5–7% between-day drift `RESULTS.md` already
documents, and the qualitative claims are unchanged: hammock's subB=1.0 line
crosses BEDTools at N=64 (1.60×, was 1.63×) and the headline subB=0.1 speedup at
N=512 is 52.35× (was 52.77×).

Full rerun table:

| N | BEDTools | subB=1.0 | subB=0.25 | subB=0.1 | speedup (0.1) |
|---|---|---|---|---|---|
| 2 | 0.82 s | 0.22 s | 0.11 s | 0.06 s | 12.66× |
| 4 | 0.52 s | 0.43 s | 0.21 s | 0.11 s | 4.93× |
| 8 | 0.64 s | 0.84 s | 0.41 s | 0.21 s | 3.07× |
| 16 | 1.12 s | 1.67 s | 0.81 s | 0.40 s | 2.77× |
| 32 | 3.04 s | 3.34 s | 1.62 s | 0.80 s | 3.80× |
| 64 | 10.69 s | 6.67 s | 3.23 s | 1.60 s | 6.69× |
| 128 | 41.10 s | 13.29 s | 6.45 s | 3.19 s | 12.89× |
| 256 | 166.46 s | 26.62 s | 12.92 s | 6.41 s | 25.98× |
| 512 | 676.62 s | 53.38 s | 25.98 s | 12.92 s | **52.35×** |

## If you adopt it

```bash
git mv paper/pairwise_scaling/plot_pairwise_scaling_v2.R \
       paper/pairwise_scaling/plot_pairwise_scaling.R      # overwrite v1
# then in that file: out_png default -> pairwise_scaling.png, drop the v2 header
ml r/4.3.0 && Rscript paper/pairwise_scaling/plot_pairwise_scaling.R
git rm paper/figures/pairwise_scaling_v2.png docs/figure3-candidate-v2.md
```

and the text that then has to move with it:

- `paper/outline.md:42` — the Panel A sentence, N(N−1)/2 → N². **Leave Panel B's
  "190 unique pairs" alone**, it is correct.
- `paper/outline.md:44` — "the correct unique-pair count, N(N−1)/2" is wrong.
- `paper/outline.md` §4.4 — the hand-measured percentages, replaceable with the
  measured table in `docs/metrics-by-default.md`. Note that section does not
  currently disclose that Figure 3's timings come from the standalone C++
  binary rather than the `hammock` CLI a reader would install, and the two have
  different default thread counts (OpenMP-all-cores vs `min(8, cpu_count())`).
- `RESULTS.md:12, 21, 32-44, 57-81, 126-127` — canonical-run framing, the
  headline 52.77×, job provenance, and the files table.
- `docs/paper_outline.md:85-90` — that table's "hammock" column is the
  **subB=0.1** arm (1.7/3/7/13 s), not the subB=1.0 arm Panel A plots. Restating
  it from the wrong column would turn ~50× into ~13× and read as a regression.
  Worth labelling the arm while editing.
- `docs/scripts/synthetic_nscaling.R:3,34` and `docs/paper_outline.md:313` point
  at the 20260512 CSV. They are self-consistent today and only need repointing
  if the new run becomes canonical — in which case regenerate and commit
  `docs/figures/synthetic_nscaling.png` too, since it is embedded at
  `docs/paper_outline.md:102`.

Disclose in any restatement that the timer change alone alters the
`hammock sketch comparison` series at low N, independently of the rerun's noise.

If you keep v1 instead, delete this file, `plot_pairwise_scaling_v2.R` and
`paper/figures/pairwise_scaling_v2.png`. The new CSV is worth keeping either
way — it is the only run containing the `+IE` arm.
