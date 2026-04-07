# mode-d-optimality

Evaluates minimizer Mode D accuracy across (k, w) parameter combinations using
synthetic data. Ground truth is **exact k-mer set Jaccard** — not character
similarity or bedtools overlap — which is what the minimizer sketch is actually
estimating.

## Key files

- `config/config.yaml` — k/w sweep ranges, mutation rates, seq length, trials
- `scripts/sweep_kw.py` — runs the sweep, outputs a CSV to `results/`
- `scripts/plot_sweep.py` — bias/MAE heatmaps + sketch-vs-true scatter plots

## Directories

- `data/`    → symlink to `/vast/blangme2/jbonnie/hammock/mode-d-optimality/data`
- `results/` → symlink to `/vast/blangme2/jbonnie/hammock/mode-d-optimality/results`

## Environment

```bash
ml anaconda
conda activate claude-ref-comparison   # same env; shares hammock + pandas/matplotlib
```

## Running

```bash
# From this directory
python scripts/sweep_kw.py                        # uses config/config.yaml
python scripts/sweep_kw.py --trials 5            # quick test run

python scripts/plot_sweep.py results/sweep_*.csv
```

## Design notes

- (k, w) pairs where w < k are automatically skipped — they produce zero
  minimizers per window and silently fall back to full-sequence hashing
- `jaccard_similarity_with_ends` is the primary metric (adds flanking k-mers);
  `jaccard_similarity` is recorded for comparison
- Mutation rate ≠ true k-mer Jaccard (depends on k); true Jaccard is computed
  exactly from k-mer sets for each trial
- Prior dnase1-hypersensitivity sweep (vs bedtools) found k=10 w=20-50 optimal,
  but that used a different ground truth — results here may differ
