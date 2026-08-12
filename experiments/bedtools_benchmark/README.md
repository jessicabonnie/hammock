# bedtools_benchmark

Timing and accuracy of `hammock-cpp` Mode B against `bedtools jaccard` on
synthetic BED corpora, plus two follow-on studies that reuse the same harness
(estimator comparison, cost of the similarity block).

**Findings for the main sweeps are in [RESULTS.md](RESULTS.md).** This file is
the inventory: what each script is, how to invoke it, and where its output is
written up. Nothing here restates a result.

## Conventions

- `results/` and `logs/` are symlinks into
  `/vast/blangme2/jbonnie/hammock_claude_experiments/bedtools_benchmark/`.
  Figures live in `figures/`, in-repo.
- BED inputs are generated per run and pre-sorted before the clock starts, so
  sort time is not charged to bedtools. See
  `docs/bedtools-parallelism-caveat.md`.
- **Three metrics shapes** (`docs/seed-metrics-column-restructure.md`, landed
  in this repo's metrics-restructure work): bare/no-flag emits only
  `jaccard_similarity_ie` (tag `_ie`), `--register-equality`/`--re` emits
  `jaccard_similarity` + `register_equality_similarity` (tag `_re`, cheap —
  skips the union pass), `--metrics` emits the full 8-column block (tag
  `_full`). `--no-metrics` is gone, not aliased.
  The timed arms in `RESULTS.md` all pass `--register-equality` (was
  `--no-metrics` pre-restructure); a run without that flag is not comparable
  to those tables. Two scripts deliberately measure the full block instead
  and run *both* shapes: `benchmark_cpp_vs_bedtools.py --metrics-arm` adds a
  fourth arm with `--metrics`, and `pairwise_cost_by_precision.py` exists to
  compare the two arms — that is its whole purpose.
- `benchmark_cpp_vs_bedtools.py` and `sweep.py` probe `--version` and refuse a
  binary older than 0.7.0. `pairwise_cost_by_precision.py` does not.
- The binary is picked by `find_hammock_cpp` (newest by mtime) unless
  `--binary` or `$HAMMOCK_CPP_BIN` is set.

## Scripts

| Script | What it does | Written up in |
| --- | --- | --- |
| `benchmark_cpp_vs_bedtools.py` | The **files axis**: N BED files, all N² pairs, hammock vs bedtools. `--num-files 2,4,…` `--threads` `--precision` `--num-intervals` `--subB-list` `--runs` `--metrics-arm` `--test`. Writes `results/cpp_vs_bedtools_t<T>_<stamp>.{csv,txt}`. | `RESULTS.md` §files sweep |
| `sweep.py` | One-axis sweeps: `--axis precision\|threads\|intervals`, with `--precisions` / `--thread-list` / `--intervals-list` and the same `--subB-list`. Writes `results/sweep_<axis>_<stamp>.csv`; the precision axis also dumps `*_pairs.csv` (per-pair jaccards). | `RESULTS.md` §§precision / threads / intervals |
| `bedtools.sh` | Helper both harnesses shell out to: pairwise `bedtools jaccard` over two file lists, GNU-parallel at N threads. Not run directly. | — |
| `make_graphs.R` | Figures from those CSVs: `--files-csv` and/or `--pairs-csv`, optional `--out-dir`. Needs `ml r/4.3.0`; renders through Cairo. | `RESULTS.md` §Plotting |
| `estimator_compare.py` | Synthesises BED pairs spanning the whole J range and compares `jaccard_similarity` (register-equality) against the inclusion-exclusion estimator, with `bedtools jaccard` as truth. `--tmp-dir` and `--out` are required; `--precisions` `--reps` `--intervals` `--data-seed` `--sketch-seed`. Output on disk: `results/estimator_compare_full.csv`. | `docs/estimator-analysis-findings.md`, `docs/jaccard-definitional-gap.md` |
| `estimator_rank_by_precision.py` | Rank fidelity (Kendall τ) of the two estimators stratified by true J and precision. Reads the CSV above (`--csv`, `--min-n`); prints tables, runs nothing. | `CLAUDE.md` divergence #2, `docs/estimator-analysis-findings.md` §9 |
| `pairwise_cost_by_precision.py` | What the full metrics block (`--metrics`) costs vs. the reduced-column arm (`--register-equality`) as a function of precision (hammock-vs-hammock, no bedtools). `--precisions` `--num-files` `--num-intervals` `--threads` `--runs` `--binary` `--output-dir`. Writes `pairwise_cost_by_precision_<stamp>.csv`; the Aug 4 run is archived at `docs/data/pairwise_cost_by_precision_20260804_164807.csv`. | `docs/metrics-by-default.md`, `CLAUDE.md` build/test notes |

## SLURM wrappers

Each hard-codes its grid and runs on Rockfish `shared` (16 cpus, 32 GB):

- `sbatch_precision.sh` — `sweep.py --axis precision`, p ∈ {10,12,14,16,18}, t=8, 64 files.
- `sbatch_threads.sh` — `sweep.py --axis threads`, t ∈ {1,2,4,8,16}, 64 files.
- `sbatch_intervals.sh` — `sweep.py --axis intervals`, 1k…1M intervals, 16 files.
- `sbatch_files.sh` — files axis at t=8, N ≤ 256.
- `sbatch_files_t16.sh` — files axis at t=16, N ≤ 512. Forwards `"$@"`, which
  is how the Aug 4 run added its `--metrics-arm` leg.
- `sbatch_files_t16_big.sh` — extension to N ∈ {512, 1024}. **No results from
  this one are reported in `RESULTS.md`**; treat it as un-run unless you find
  a matching CSV in `results/`.

## Reference dates

Results in `results/` span 2026-05-10 through 2026-08-04 and were produced by
several different hammock builds. `RESULTS.md` opens with the provenance table
saying which stamp is canonical for which sweep — read it before pairing a CSV
with a claim.
