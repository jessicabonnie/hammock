# Experiment Design: Mode D (k, w) Parameter Optimality

**Project:** Determine optimal minimizer parameters (k-mer size k, window size w) for hammock Mode D sequence sketching  
**Status:** Sweep submitted; results pending  
**Last Updated:** 2026-04

---

## Goal

Select (k, w) values for the `claude-ref-comparison` experiments by evaluating how well minimizer sketches approximate exact k-mer Jaccard similarity on synthetic sequence pairs with known ground truth.

---

## Method

For each (k, w, mutation_rate) combination:
1. Generate a random sequence and a mutated copy at the given mutation rate
2. Compute exact k-mer Jaccard (ground truth) from the k-mer sets directly
3. Compute hammock Mode D sketch similarity for the same pair
4. Record error = sketched − true

Summarise bias and RMSE across trials and mutation rates; select (k, w) that minimises error while keeping w ≥ k (w < k produces zero minimizers per window).

Configuration: `config/config.yaml` — k ∈ {3, 5, 7, 10}, w ∈ {5, 10, 15, 20, 40}, precision 24, 20 trials × 6 mutation rates.

---

## Similarity columns

Mode D produces two similarity columns per comparison:

| Column | Description |
|---|---|
| `jaccard_similarity` | Minimizers only |
| `jaccard_similarity_with_ends` | Minimizers + first/last k boundary k-mers |

**`jaccard_similarity_with_ends` is the primary column** for downstream use in `claude-ref-comparison`. However, both columns are recorded in the sweep CSVs, and a goal of this experiment is to determine which column tracks the true k-mer Jaccard more accurately across (k, w) settings. The plot script should report error metrics for both columns so they can be compared directly. The better-performing column should be used as `primary_sim_col` in subsequent experiments.

---

## Output

Results written to `/vast/blangme2/jbonnie/hammock/mode-d-optimality/results/sweep_*.csv`.

Visualize with:
```bash
python scripts/plot_sweep.py results/sweep_*.csv --out results/figs/
```

---

## Next step

Feed the winning (k, w) into `claude-ref-comparison/config/config.yaml` under `kmer_sizes` and `window_sizes` (narrowed to the optimal pair), then run the full validation pipeline.
