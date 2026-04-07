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
ml anaconda && conda activate hammock   # always run this at the start of a terminal session
```

## Running

```bash
# From this directory
python scripts/sweep_kw.py                        # uses config/config.yaml
python scripts/sweep_kw.py --trials 5            # quick test run

python scripts/plot_sweep.py results/sweep_*.csv
```

## Hammock output columns

Always use:
- `hash_with_ends_similarity` — primary metric (minimizers + flanking end k-mers)
- `hash_similarity` — secondary metric (minimizers only)

**Never use `jaccard_similarity`** — legacy column, ignore even if present.

A key goal of this sweep is to determine which of `hash_with_ends_similarity` vs
`hash_similarity` performs better (lower bias/MAE against true k-mer Jaccard) across
comparable (k, w, mutation rate) conditions, to inform which to use as the primary
metric in claude-ref-comparison.

## Design notes

- (k, w) pairs where w < k are automatically skipped — they produce zero
  minimizers per window and silently fall back to full-sequence hashing
- `hash_with_ends_similarity` adds flanking end k-mers to the minimizer sketch;
  `hash_similarity` uses minimizers only
- Mutation rate ≠ true k-mer Jaccard (depends on k); true Jaccard is computed
  exactly from k-mer sets for each trial
- Prior dnase1-hypersensitivity sweep (vs bedtools) found k=10 w=20-50 optimal,
  but that used a different ground truth — results here may differ
