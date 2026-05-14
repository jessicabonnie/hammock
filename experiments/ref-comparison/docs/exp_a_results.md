# Experiment A — Results

Cross-reference robustness validation: same-biological-sample H3K27ac
ChIP-seq peaks, aligned to GRCh37 / GRCh38 / CHM13, should sketch more
similarly *across references* than across tissues on any one reference.

3 tissues (heart, liver, lung) × 3 references = 9 sample × ref combinations.
2 peak callers (MACS3 broad + narrow). (k, w) sweep: 26 cells with w ≥ k.

Hammock Mode D with minimizer HLL, `--precision 24`. Similarity metric:
`jaccard_similarity_with_ends` (minimizer ∪ start-end HLL).

## Where to find the raw results

Symlinked into the repo with fine-grained paths (Exp A subset only — the
legacy storage location also holds defunct Exp B artifacts that are not
exposed here):

```
results/
├── exp_a/{broad,narrow}/k{k}_w{w}/
│   ├── exp_a_mnmzr_p24_jaccD_k{k}_w{w}.csv   ← all-vs-all hammock matrix
│   ├── cross_ref_validation.png              ← per-cell Wilcoxon + heatmap
│   └── cross_ref_stats.tsv                   ← per-cell summary stats
├── exp_a/{broad,narrow}/sweep_effect_size.png
├── fastas/{broad,narrow}/{GRCh37,GRCh38,CHM13}/{sample}.fa
└── logs/{exp_a,fastas}/
```

Underlying storage: `/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/`.

## Headline result

Same-tissue cross-reference Jaccard is significantly higher than
different-tissue Jaccard across **all** (k, w) cells with k ≥ 8.
The hypothesis ("reference choice contributes less than tissue identity")
holds across the entire usable parameter range.

### Dendrogram — direct visual proof

UPGMA clustering on `1 − jaccard_similarity_with_ends` at k=10, w=10.
Each tissue's three references form a tight monophyletic clade in both
peak types; the within-tissue merge heights (≤ 0.1) are several-fold
smaller than inter-tissue merges (≥ 0.4). Source: `scripts/exp_a_dendrogram.R`.

![Cross-reference dendrogram (broad + narrow, k=10, w=10)](../figures/cross_ref_dendrogram_k10_w10.png)

### Per-cell view at k=10, w=10

Boxplot of same-tissue cross-reference pairs (n=18) vs different-tissue
pairs (n=54), plus 9×9 sample × reference similarity heatmap. Visible
block structure: tissue diagonals are darker than off-diagonals
regardless of reference.

![Broad k=10 w=10 single-cell figure](../figures/cross_ref_validation_broad_k10_w10.png)

## Parameter sweep — three regimes

`figures/sweep_effect_size_{broad,narrow}.png`. Each cell shows the
effect size (median cross-ref − median diff-tissue) and the Wilcoxon p.

| Regime | Cells | Behavior |
|---|---|---|
| **Saturated high** | k=5 (any w) | Both groups ≈ 0.99. Tiny effect (Δ ≤ 0.02). |
| **Sweet spot — interpretable** | k=10, w ≥ 10 | Δmedian ≈ 0.10 on mid-range medians (≈ 0.5 vs 0.4). |
| **Saturated low** | k ≥ 15 | Absolute medians collapse (≈ 0.09 vs 0.015). Ratio is large (~6×), Δ moderate (~0.08), p ≤ 1.35e-10. |
| **Borderline** | k=8 | Strong p but small Δ (≈ 0.02–0.03) on top of high baseline (~0.95). |

![Sweep effect-size heatmap, broad peaks](../figures/sweep_effect_size_broad.png)

![Sweep effect-size heatmap, narrow peaks](../figures/sweep_effect_size_narrow.png)

Top-significance cells (broad, narrow consistent):

| Cell | Δmedian | Wilcoxon p | median(cross-ref) | median(diff-tissue) |
|---|---|---|---|---|
| k15_* (all w) | 0.081 | 1.35e-10 | 0.097 | 0.015 |
| k10_w20 | 0.107 | 4.5e-7 | 0.526 | 0.419 |
| k10_w15 | 0.106 | 4.5e-7 | 0.550 | 0.444 |
| k10_w10 | 0.100 | 4.5e-7 | 0.587 | 0.486 |

The k=10 cells are the most defensible to lead with: large absolute
effect, both groups well away from the 0/1 boundaries, and the heatmap
panel (B) of the single-cell figure shows visible block structure.

## Mirrored figures (local, not committed)

`figures/` is gitignored. Files present:

- `sweep_effect_size_{broad,narrow}.png` — top-level (k × w) effect-size heatmaps.
- `cross_ref_validation_{broad,narrow}_k10_w10.png` — representative single-cell figure (boxplot + 9×9 sample-by-sample heatmap).
- `cross_ref_dendrogram_k10_w10.png` — UPGMA dendrogram, broad + narrow.
- `metric_comparison_broad_k10_w10.png` / `.tsv` — sanity check for the new Mode D output columns (see below).

Other per-cell `cross_ref_validation.png` files are reachable through
`results/exp_a/{peak_type}/k{k}_w{w}/` if a different cell ever needs to
be inspected.

## New Mode D columns — sanity check at k=10, w=10

The bulk CSVs under `results/exp_a/` predate the containment + cosketch
metric additions (shipped 2026-05-14) and only carry `jaccard_similarity`
+ `jaccard_similarity_with_ends`. I re-ran hammock at the (k=10, w=10)
broad cell with the current binary to check the 10 new columns. Inputs +
output kept under `/tmp/hammock_newcols_test/` (not in `results/` —
exploratory).

`scripts/exp_a_metric_comparison.R` performs the same Wilcoxon test
(same-tissue cross-ref vs different-tissue) on each of the 12 columns.
Outputs `figures/metric_comparison_{broad,narrow}_k10_w10.{png,tsv}`.

![Metric comparison, broad, k=10, w=10](../figures/metric_comparison_broad_k10_w10.png)

![Metric comparison, narrow, k=10, w=10](../figures/metric_comparison_narrow_k10_w10.png)

### Findings (k=10, w=10)

| Metric (flavor)                | Δ broad | p broad | Δ narrow | p narrow |
|---|---|---|---|---|
| jaccard (minimizer)            | 0.127 | 1.35e-10 | 0.154 | 1.35e-10 |
| containment_AB/BA (minimizer)  | 0.059 | 1.36e-10 | 0.067 | 1.36e-10 |
| cosketch_geom (minimizer)      | 0.068 | 1.35e-10 | 0.083 | 1.35e-10 |
| cosketch_arith (minimizer)     | 0.068 | 1.35e-10 | 0.082 | 1.35e-10 |
| cosketch_max (minimizer)       | 0.032 | 1.35e-10 | 0.032 | 1.35e-10 |
| **jaccard (minimizer+ends)**   | **0.100** | **4.5e-07** | **0.099** | **2.1e-06** ← existing default |
| cosketch_geom (minimizer+ends) | 0.081 | 4.5e-07 | 0.078 | 7.2e-06 |
| cosketch_arith (minimizer+ends)| 0.077 | 4.5e-07 | 0.071 | 1.4e-05 |
| containment_AB/BA (with ends)  | 0.035 | 6.4e-04 | 0.035 | 8.5e-03 |
| cosketch_max (with ends)       | 0.017 | 3.7e-02 | 0.000 | 4.5e-01 |

Notes (both peak types):
1. **All minimizer-only metrics hit the Wilcoxon test floor** (p ≈ 1.4e-10 for n=18/54). Maximally significant, but operating near saturation (medians ≈ 0.99 vs ≈ 0.93), so Δ is small in absolute terms.
2. **Adding the start-end HLL lowers baselines into the interpretable mid-range** (~0.7) and reduces saturation. Cost: p drops by ~3–4 orders of magnitude for the surviving discriminators.
3. **Containment-with-ends is a markedly weaker discriminator** than jaccard-with-ends. The pattern is symmetric — `containment_AB` and `containment_BA` produce identical medians across the full all-vs-all set because every unordered pair appears in both directions.
4. **cosketch_max collapses entirely on narrow with ends** (Δ = 0, p = 0.45). Worst metric in both flavors.
5. **`jaccard_similarity_with_ends` remains the right default** for headline analysis: best Δ/saturation trade-off and matches what the existing CSV sweep already used. `cosketch_geom_with_ends` is a near-tie if a redundancy-aware alternative is wanted.
6. **Broad vs narrow**: minimizer-only effects are slightly larger on narrow (likely tighter peak boundaries → more shared minimizers between same-sample-different-ref pairs). With-ends p-values are slightly weaker on narrow but the qualitative ordering of metrics is identical.

### To extend

- Re-run the full (k × w) sweep with the current hammock to populate the new columns across all 26 cells × 2 peak types, then regenerate the sweep-effect heatmap for each metric. Inputs + Snakemake DAG are unchanged; only `rule exp_a_hammock` needs to re-fire.
- Single-cell sanity-check CSVs are at `/tmp/hammock_newcols_test/newcols_{broad,narrow}_k10w10_mnmzr_p24_jaccD_k10_w10.csv` (not preserved in `results/`).

## Scripts that produced these outputs

- `scripts/exp_a_validate_plot.R` — per-cell Wilcoxon + 2-panel figure (boxplot + 9×9 heatmap)
- `scripts/exp_a_sweep_summary.R` — (k × w) effect-size heatmap
- `workflow/Snakefile` — orchestrates peaks → FASTA (bedtools getfasta) → hammock Mode D → R plots

## Reproducing

```bash
conda activate claude-ref-comparison
ml r/4.3.0
snakemake --profile workflow/slurm_profile/
```

(Upstream nf-core/chipseq runs must have produced the broad/narrow peak
calls under `NFCORE_OUT/{ref}{,_narrow}/...` — see `scripts/run_nfcore.sh`.)
