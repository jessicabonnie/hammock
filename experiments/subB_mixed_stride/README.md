# subB / mixed-stride sweep (Mode B)

Does `--mixed-stride` with `--subB < 1` make hammock Mode B faster than the
no-subsample baseline, and how does Jaccard accuracy degrade as subB drops?

## Design

- **Corpus:** 18 seeded synthetic BEDs — 6 files × 3 size classes (10k, 100k,
  1M intervals). Reproducible from `generate_data.py` (`BASE_SEED=20260511`).
- **Sweep:** `subB ∈ {1.0, 0.5, 0.25, 0.1, 0.05, 0.01}`, with `--mixed-stride`
  always on. 5 replicates per (size_class, subB) cell.
- **Tool:** standalone `hammock-cpp` (Mode B, precision 18, 8 threads). No
  Python in the timing loop.
- **Ground truth for accuracy:** the same binary at `subB=1.0`. Self-
  consistent — this measures *subsample-induced* drift, not absolute HLL
  error.

## Files

- `generate_data.py` — seeded BED corpus generator (writes into `data/`).
- `run_sweep.py` — sweep driver; emits a long-form CSV (one row per pair ×
  subB × rep) to `results/sweep_<timestamp>.csv`.
- `analyze.R` — reads the latest sweep CSV, writes a summary CSV and four
  PNGs into `figures/` via `CairoPNG` (system R has no X11 / native cairo).
- `sbatch_sweep.sh` — single-job workflow: generate BEDs if needed → sweep →
  R analysis.

`data/`, `results/`, and `logs/` are symlinks into `/vast/blangme2/...` so the
heavy artifacts don't live in the repo. Recreate per machine:

```
mkdir -p /vast/blangme2/jbonnie/hammock_claude_experiments/subB_mixed_stride/{data,results,logs}
ln -sfn /vast/.../subB_mixed_stride/data    data
ln -sfn /vast/.../subB_mixed_stride/results results
ln -sfn /vast/.../subB_mixed_stride/logs    logs
```

## How to run

```bash
# one-shot
sbatch experiments/subB_mixed_stride/sbatch_sweep.sh

# or interactively
python experiments/subB_mixed_stride/generate_data.py        # ~once
python experiments/subB_mixed_stride/run_sweep.py            # ≈30 min on the cluster
ml r/4.3.0 && Rscript experiments/subB_mixed_stride/analyze.R
```

## Findings (2026-05-11 run)

| size  | subB  | wall median | speedup vs subB=1.0 | MAE (J) | max |J error| |
|-------|-------|-------------|---------------------|---------|----------------|
| 10k   | 1.00  |   1.14 s    |  1.00×              |  0      |  0             |
| 10k   | 0.50  |   1.33 s    |  0.86×              |  4.4e-4 |  1.8e-3        |
| 10k   | 0.10  |   0.37 s    |  3.09×              |  6.0e-4 |  2.1e-3        |
| 10k   | 0.01  |   0.14 s    |  8.08×              |  2.0e-2 |  2.3e-2        |
| 100k  | 1.00  |   9.92 s    |  1.00×              |  0      |  0             |
| 100k  | 0.10  |   3.01 s    |  3.30×              |  6.1e-4 |  1.6e-3        |
| 100k  | 0.01  |   0.88 s    | 11.22×              |  8.4e-4 |  1.6e-3        |
| 1M    | 1.00  |  95.5 s     |  1.00×              |  0      |  0             |
| 1M    | 0.10  |  28.7 s     |  3.33×              |  1.7e-5 |  6.1e-5        |
| 1M    | 0.01  |   7.09 s    | **13.48×**          |  1.6e-5 |  3.8e-5        |

- **Speedup is non-monotonic.** At `subB=0.5`, `--mixed-stride` is actually
  *slower* than `subB=1.0` (0.77–0.86× across size classes) — the
  mixed-stride bookkeeping overhead outweighs the points saved at that rate.
  Real wins start at `subB ≤ 0.25` and peak at `subB=0.01`.
- **Accuracy is excellent for ≥100k intervals.** MAE stays ≤ 1e-3 all the
  way down to `subB=0.01`. For 1M files the error is essentially noise
  (max ≈ 6e-5).
- **The cliff is small-file + aggressive-subB.** At 10k intervals,
  `subB=0.01` jumps to MAE ≈ 0.02 (max ≈ 0.023). Above `subB=0.05` the 10k
  class is also fine.

**Recommendation:** for large BED corpora, default to `--mixed-stride` with
`subB ∈ [0.05, 0.1]`. Avoid `subB ≥ 0.5` (net slowdown for no accuracy
benefit). Below `subB=0.05`, watch the accuracy on small inputs.

See `figures/{wall_time,speedup,mae,jaccard_error}_vs_subB.png` for the
plots; `results/summary_<timestamp>.csv` for the full numbers.
