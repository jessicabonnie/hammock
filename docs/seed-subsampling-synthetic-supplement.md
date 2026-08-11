# Seed: synthetic-corpus subB speedup for the mixed-stride paragraph

## Question

`paper/outline.md`/`draft.md` claimed "a 10-fold subsample buys 4.21× on the
synthetic corpus at p=14 but only 1.83× on Maurano at p=18." Neither number
survived a source check (2026-08-11).

## What was wrong

- **`1.83×` (Maurano, p=18) traced to a real file** —
  `docs/data/maurano_subB_summary.csv`, mixed-stride subB 1.0→0.1:
  8.20s → 4.47s = 1.834×. But that file is the **register-equality-era** run,
  superseded by `docs/data/maurano_subB_ie_summary.csv` (the file
  `paper/pairwise_scaling/plot_pairwise_scaling.R` actually reads for
  Figure 3B, corrected baseline, `jaccard_similarity_ie` on): 9.404s → 5.027s
  = **1.871×**. The outline was quoting a different, older run than the one
  behind the figure sitting next to it.
- **`4.21×` (synthetic, p=14) traced to nothing.** Every subB-bearing
  synthetic dataset in the repo was checked directly against its source CSV:

  | source | precision | subB 1.0→0.1 speedup |
  |---|---|---|
  | `experiments/subB_mixed_stride/results/summary_synthetic_20260512_144455.csv` (10k/100k/1M) | p=18 | 3.09× / 3.71× / 3.93× |
  | `experiments/bedtools_benchmark` intervals sweep (1k/10k/100k/1M) | p=14 | 3.97× / 4.67× / 4.83× / 4.88× |
  | `experiments/bedtools_benchmark` threads sweep (t=16) | p=14 | 4.15× |
  | `docs/data/cpp_vs_bedtools_t16_20260804_172242.csv` (job 29552415) | p=14 | 4.13×–4.18× |

  Nothing lands on 4.21×. A broader repo-wide scan of every p=14 cell with a
  subB=1.0→0.1 pair (170 cells across all `experiments/`/`docs/data/` CSVs,
  including one previously-uncited run,
  `experiments/bedtools_benchmark/results/cpp_vs_bedtools_t16_20260808_190441.csv`)
  found a full range of 2.57×–5.11× and a closest single value of 4.1998×
  (N=64, t=16) — close but not exact, and not an obviously "right" cell to
  have quoted. No plausible transposition (2.41×, 1.24×, ...) matches
  anything nearby either. p=14 is also the wrong precision to reach for at
  all: `plot_pairwise_scaling.R`'s own comment says Panel A's synthetic CSV
  "replaces `cpp_vs_bedtools_t16_20260804_172242.csv`, which was p=14 — a
  precision used nowhere else in the paper." p=18 is the CLI default and the
  precision every other synthetic/Maurano number in the paper uses.
- No synthetic-corpus run has ever been taken with `--metrics`
  (`jaccard_similarity_ie`), so there was no IE-era synthetic number to quote
  even at the right precision — `docs/data/cpp_vs_bedtools_t16_p18.csv`
  (Panel A's synthetic source) only has `sub_b=1` rows.

## Fix applied 2026-08-11

- `paper/outline.md` / `draft.md`: dropped the synthetic-corpus half of the
  sentence; corrected Maurano to **1.87×** (matching
  `maurano_subB_ie_summary.csv`, the file Figure 3B is actually built from).
- Submitted job 29756567 (`experiments/subB_mixed_stride/sbatch_synthetic_ie.sh`):
  synthetic corpus, mixed-stride only, subB ∈ {1.0, 0.1, 0.01}, p=18,
  `--metrics`, 5 reps, all three size classes (10k/100k/1M), fixed at 6 files
  per class (bar-chart design, paralleling Figure 3B). **Deliberately
  cancelled mid-run** (~14:17, 35/45 cells done — not a crash, not OOM,
  `sacct` shows a clean external cancel) once the design changed: the user
  asked for a line plot instead (subB as separate lines, N files-per-side on
  the x-axis, paralleling Figure 3A's scaling design rather than Figure 3B's
  bars), so a fixed-N/fixed-corpus sweep no longer answers the right
  question. Partial output is still on disk if the per-size-class numbers are
  ever wanted for something else:
  `experiments/subB_mixed_stride/results/sweep_synthetic_ie_20260811_140716.csv`
  (10k and 100k complete at all 3 subB values; 1M only got 4/5 reps of
  subB=1.0, subB∈{0.1,0.01} never ran). Medians from what did complete: 10k
  3.25×, 100k 3.85× (subB=1.0→0.1) — in the same neighborhood as the old
  register-equality numbers (3.09×/3.71×) and nowhere near 4.21×, which is
  further evidence 4.21× was never a real measurement at the correct
  precision either.
- **Current plan (superseding the above):** an N-sweep line plot, reusing
  `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py`'s existing
  harness (the one behind Figure 3A / `plot_largeN_supplement.R`) rather than
  `subB_mixed_stride/run_sweep.py`. Three hammock lines (subB = 1.0, 0.1,
  0.01), fully measured out to N=2048 (`--no-bedtools`; hammock is
  near-linear at p=18 in this regime, so this is minutes not hours — see
  `plot_largeN_supplement.R`'s header comment for the same argument applied
  to the existing large-N supplement). bedtools appears as a **projected**
  reference line, fit the same way `plot_largeN_supplement.R` does: log-log
  power law on its measured points (`docs/data/cpp_vs_bedtools_t16_p18.csv`,
  N≥32), extrapolated dashed/labeled "projected", fitted exponent printed and
  checked against the theoretical 2.0. Two more reviewers are checking this
  plan (methodology + `--metrics`-per-subB-arm feasibility in the existing
  script) before it's executed.

## Lesson

Same shape as the bedtools-baseline retraction: a plausible-looking number
sitting next to real data was never checked against the CSV it claimed to
summarize. Trace every quoted number to its exact source cell before it goes
in the outline, not just once at figure-build time.
