# subB / subB-method sweep (Mode B)

Compares the three `--subB-method` modes (`hash-threshold`, `mixed-stride`,
`single-hash`) against the no-subsample baseline (`subB=1.0`) on two corpora,
and against `bedtools jaccard` on the real-data corpus.

1. **synthetic** — 18 random BEDs across 3 size classes (10k / 100k / 1M
   intervals).
2. **maurano** — 20 fetal-tissue DNase-seq hotspot BEDs (Maurano et al. 2012,
   hg19), shared with the `maurano_dhs_validation` experiment.

**Results, figures, and recommendation are in [RESULTS.md](RESULTS.md).**

## Design

- **Sweep grid:** all 3 methods × `subB ∈ {1.0, 0.5, 0.25, 0.1, 0.05, 0.01}` ×
  5 replicates per corpus class. Synthetic adds a 3-way size_class facet, so
  it's 3 × 6 × 5 × 3 = 270 runs; Maurano is one class, so 3 × 6 × 5 = 90 runs.
- **Tool:** standalone `hammock-cpp` (Mode B, precision 18, 8 threads). No
  Python in the timing loop. `run_sweep.py` requires the binary to be ≥ 0.7.0
  and passes `--no-metrics`, so timed runs are 3-column
  (`query`/`reference`/`jaccard_similarity`). The archived sweeps in
  `RESULTS.md` were taken with a pre-0.7.0 binary, where that shape was the
  default and the flag did not exist.
- **Accuracy column:** `jaccard_similarity` (register-equality), not
  `jaccard_similarity_ie`. See the caveat box at the top of `RESULTS.md`.
  The IE question is answered separately by `run_ie_subb.py` +
  `analyze_ie_subb.py` (RESULTS.md §4.6): subB does not meaningfully degrade
  `jaccard_similarity_ie` anywhere in 0.01 ≤ subB ≤ 1. Those two scripts measure
  accuracy against exact bedtools truth rather than drift against subB=1.0, and
  write `results/ie_maurano_*.csv` — deliberately not a `sweep_*` name, which the
  `pick_latest` globs in `analyze.R` would adopt.
- **Ground truth for accuracy:** `subB=1.0`, averaged across methods (every
  method should produce the same Jaccard at subB=1 since all points pass the
  gate). Self-consistent; isolates subsample-induced drift from HLL error.
- **No-subsample baseline for speedup:** median wall time across methods at
  `subB=1.0`, per size_class.
- **Bedtools baseline (Maurano):** pairwise `bedtools jaccard` on
  pre-sorted copies of the same BEDs, GNU parallel at 8 threads. Sort time
  is not charged (matches `bedtools_benchmark/` convention).

## Files

- `generate_data.py` — seeded synthetic BED corpus generator.
- `prepare_maurano.sh` — symlinks the 20 Maurano BEDs into `data/maurano/`.
- `run_sweep.py` — sweep driver; `--corpus {synthetic,maurano}` selects
  inputs. Writes `results/sweep_<corpus>_<stamp>.csv`.
- `run_bedtools_maurano.py` — pre-sorts the Maurano BEDs and times
  pairwise `bedtools jaccard` (8-way GNU parallel) as a reference baseline.
  Writes `results/sweep_maurano_bedtools_<stamp>.csv`.
- `analyze.R` — reads the latest CSV for each corpus, writes summary CSVs
  and per-corpus PNG sets to `figures/`. If a bedtools CSV is present,
  overlays it as a reference line on the wall-time and speedup plots.
- `headline_figures.R` — produces the three TL;DR plots
  (`figures/headline_*.png`) referenced from `RESULTS.md`. Run after the
  sweeps and `analyze.R`.
- `pareto_variants.R` — four renderings of the headline Maurano Pareto plot
  (path in subB order / path re-sorted by MAE / points only / points + IQR)
  into a single `figures/pareto_variants.pdf`. Diagnostic for the zigzag
  discussed in `RESULTS.md`.
- `sbatch_sweep.sh` — single-job workflow wrapper.

`data/`, `results/`, and `logs/` are symlinks into `/vast/blangme2/...`.

## How to run

```bash
# synthetic only
python generate_data.py
python run_sweep.py --corpus synthetic

# real data (hammock + bedtools reference)
./prepare_maurano.sh
python run_sweep.py --corpus maurano
python run_bedtools_maurano.py        # bedtools baseline for the speedup plot

# analysis (picks up the latest CSV per corpus + bedtools if present)
ml r/4.3.0 && Rscript analyze.R
ml r/4.3.0 && Rscript headline_figures.R    # publication-style headline PNGs
ml r/4.3.0 && Rscript pareto_variants.R     # figures/pareto_variants.pdf (4 pages)
```

`run_sweep.py` also takes `--size-classes`, `--methods`, `--subB-values`,
`--reps`, `--precision/-p`, `--threads/-t`, `--binary`, `--output`,
`--verbose` and `--smoke` (a one-cell dry run). The three R scripts take no
arguments — each picks the newest matching CSV out of `results/`.
