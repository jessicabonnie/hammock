# Experiment A — Results

Cross-reference robustness validation: same-biological-sample H3K27ac
ChIP-seq peaks, aligned to GRCh37 / GRCh38 / CHM13, should sketch more
similarly *across references* than across tissues on any one reference.

3 tissues (heart, liver, lung) × 3 references = 9 sample × ref combinations.
2 peak callers (MACS3 broad + narrow). (k, w) sweep: 20 cells with w ≥ k.

Hammock Mode D with minimizer HLL, `--precision 24`. Similarity metric:
`jaccard_similarity_with_ends` (minimizer ∪ start-end HLL) for the stats
tables below; **the dendrogram was re-rendered on `jaccard_similarity` in
v0.6.0** (see that section). The `_with_ends` column no longer exists —
CLAUDE.md divergence #8, `docs/mode-d-ends-removal.md`.

> **Rerun 2026-05-14**: all hammock CSVs, stats, and figures below were
> regenerated after the Mode D minimizer-ingest bug fix (orig's slow-path
> `add_string(str(hash_val))` silently dropped most minimizers at k ≥ 12;
> `python/hammock/modes/sequence.py` now calls `add_hash64` directly —
> see CLAUDE.md "Intentional divergences" §6). Pre-fix the k=15 and k=20
> cells appeared to be a "saturated low" regime; post-fix they are the
> strongest discriminators in the sweep.

> **Corrected HLL rehash rerun 2026-08-25**: the complete 63-cell Figure 5
> union was rerun for broad and narrow peaks at p=24 with
> `--sequence-hll-hash rehash-selector64` and prespecified seeds 1, 42, and
> 99. All 378 matrices passed completeness, symmetry, diagonal, metadata, and
> boundedness validation. On `jaccard_similarity_ie`, k20_w20 remained the
> unique lead cell in both peak types at every seed. At the paper-facing seed
> 42, Δ is 0.5590 broad and 0.5954 narrow (legacy: 0.5632 and 0.5986), and the
> nine broad peak sets retain the same three tissue-first clades. Raw results
> and summaries are under `results/exp_a_rehash_selector64/`; the frozen
> manifest, array driver, and validator are in `rehash_rerun/`.

## Where to find the raw results

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
different-tissue Jaccard across every (k, w) cell with k ≥ 8 in both
peak types. The hypothesis ("reference choice contributes less than
tissue identity") holds across the entire usable parameter range, and
at k ∈ {15, 20} the two groups are **fully separated** — the minimum
same-tissue cross-ref similarity exceeds the maximum different-tissue
similarity.

### Dendrogram — direct visual proof

UPGMA clustering on `1 − jaccard_similarity`. Each tissue's
three references form a tight monophyletic clade in both peak types.
At k=15, w=15 the clades are deeply separated; at k=10, w=10 they
still partition cleanly but with smaller margin. Source:
`scripts/exp_a_dendrogram.R`.

> **Re-rendered 2026-08-03** on `jaccard_similarity`, after v0.6.0 removed
> `jaccard_similarity_with_ends` (CLAUDE.md divergence #8). The topology is
> **unchanged**: cutting either dendrogram at k=3 recovers the 3×3
> tissue partition exactly, under both columns and both peak types. Only the
> figure's provenance changed, not its claim. The stats tables further down
> this file were computed on the old column and have not been regenerated.

![Cross-reference dendrogram (broad + narrow, k=15, w=15)](../figures/cross_ref_dendrogram_k15_w15.png)

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
| **Saturated high** | k=5 (any w); k=8 (any w) | Both groups ≈ 0.95–0.99. Δ ≤ 0.02; k=5 mostly non-significant. |
| **Interpretable mid** | k=10, w ≥ 10 | Medians ≈ 0.55–0.65; Δ ≈ 0.09; p ≈ 1e-7. Groups overlap. |
| **Interpretable + fully separated** | k=15, k=20 (any valid w) | Medians cross-ref ≈ 0.55–0.78 vs diff-tissue ≈ 0.24–0.39. Δ ≈ 0.32–0.45, p ≈ 1.35e-10 (Wilcoxon floor). **min(xref) > max(diff-tissue)**. |

![Sweep effect-size heatmap, broad peaks](../figures/sweep_effect_size_broad.png)

![Sweep effect-size heatmap, narrow peaks](../figures/sweep_effect_size_narrow.png)

Top-effect cells (broad / narrow consistent ordering):

| Cell | peak | Δmedian | Wilcoxon p | median(cross-ref) | median(diff-tissue) |
|---|---|---|---|---|---|
| k15_w15 | broad  | 0.398 | 1.35e-10 | 0.783 | 0.385 |
| k15_w15 | narrow | 0.387 | 1.35e-10 | 0.729 | 0.342 |
| k20_w20 | broad  | 0.413 | 1.35e-10 | 0.725 | 0.313 |
| k20_w20 | narrow | 0.366 | 1.35e-10 | 0.643 | 0.277 |
| k15_w20 | broad  | 0.383 | 1.35e-10 | 0.743 | 0.361 |
| k15_w30 | broad  | 0.353 | 1.35e-10 | 0.676 | 0.323 |
| k20_w30 | broad  | 0.379 | 1.35e-10 | 0.657 | 0.277 |
| k10_w10 | broad  | 0.086 | 2.02e-07 | 0.640 | 0.554 |
| k10_w10 | narrow | 0.087 | 4.54e-07 | 0.640 | 0.553 |

**k=15, w=15 is the best lead cell**: largest Δ, maximally significant,
medians sit in the interpretable mid-range, and the groups are *fully
separated* (broad: min cross-ref 0.753 > max diff-tissue 0.443). k=10
cells are still useful when overlap (rather than clean separation) is
the story being told.

> **Recomputed on the surviving column, 2026-08-06.** Every number in the
> table above is on `jaccard_similarity_with_ends`, which v0.6.0 deleted.
> Recomputing the same statistic on `jaccard_similarity` from the same
> archived CSVs (18 same-tissue cross-ref pairs, 54 different-tissue,
> counting both orderings as the original did):
>
> | cell | peak | median(cross-ref) | median(diff-tissue) | Δ | min(xref) | max(diff) |
> |---|---|---:|---:|---:|---:|---:|
> | k15_w15 | broad  | 0.926 | 0.443 | 0.483 | 0.913 | 0.506 |
> | k15_w15 | narrow | 0.953 | 0.417 | 0.536 | 0.898 | 0.468 |
> | k10_w10 | broad  | 0.991 | 0.920 | 0.071 | 0.988 | 0.935 |
> | k10_w10 | narrow | 0.993 | 0.910 | 0.084 | 0.985 | 0.927 |
>
> **The headline survives the column change and gets stronger**: k=15, w=15
> remains the lead cell and Δ grows (0.398 → 0.483 broad). The k=10 cells
> shift into the saturated-high regime on this column (medians ≈ 0.99 vs
> 0.92), so the "interpretable mid-range" characterisation above is specific
> to `_with_ends`.
>
> **But full separation is no longer what distinguishes the lead cell.** On
> `jaccard_similarity` *all four* cells separate completely — k10_w10 broad
> manages 0.988 > 0.935 and narrow 0.985 > 0.927. On the deleted `_with_ends`
> column k=10 did **not** separate (0.579 < 0.610), so the column change
> quietly turned "k=10 fails the criterion" into "k=10 passes it with a thin
> margin". Rank the cells on Δ, not on whether they separate; separation is
> satisfied everywhere here and no longer carries information.
>
> Wilcoxon p-values were not recomputed. Two reasons to treat the archived
> p = 1.35e-10 as indicative only: the group ordering is unchanged, so it can
> only be at or below the same test floor — but the n=18/54 counts are
> **ordered** pairs, so every unordered comparison is entered twice. The
> independent units are 9 cross-reference and 27 different-tissue
> comparisons over 9 files, and even those share files. The medians in the
> table above are unaffected by the duplication; the significance test is not.
>
> Two standing caveats on all of these. `jaccard_similarity` is
> register-equality, not set Jaccard — it carries a chance-agreement floor
> and is not rank-faithful across pairs of differing set size
> (`CLAUDE.md` divergence #2). And every CSV here predates 2026-05-14 in
> schema terms only in that it lacks `jaccard_similarity_ie`; the calibrated
> value is exactly recoverable from the `containment_AB`/`containment_BA`
> columns these files do carry, via `J = 1/(1/C_AB + 1/C_BA − 1)`.

### Recomputed on `jaccard_similarity_ie` across all 20 cells, 2026-08-08

Generator: `experiments/ref-comparison/estimator_ie_crossref.py`, writing
`results/exp_a_estimator_delta.csv`. Same statistic, same 18/54 ordered pairs,
no sketching — `_ie` is recovered from the containments via the canonical
`runner._jaccard_ie_from_containments`. No p-values, for the reason given above.

**Two separate findings. Only the second is about the estimator.**

**1. "k=15, w=15 is the lead cell" does not survive the full grid — and this is
a `jaccard_similarity` result, not an `_ie` one.** The 2026-08-06 recompute
above covered only two cells (k15_w15, k10_w10), so it never tested whether
k15_w15 still won once every cell was scored. It does not. Ranked by Δ on
`jaccard_similarity`, broad:

| rank | cell | Δ (register-equality) | Δ (IE) |
|---:|---|---:|---:|
| 1 | k20_w30 | 0.5398 | 0.5609 |
| 2 | k20_w20 | 0.5326 | 0.5632 |
| 3 | k15_w30 | 0.4940 | 0.5081 |
| 4 | k15_w20 | 0.4883 | 0.5095 |
| 5 | **k15_w15** | **0.4827** | **0.5099** |

k15_w15 is fifth of twenty, on both columns; narrow peaks give the same
ordering. It led on the deleted `_with_ends` column (Δ 0.398 against k20_w30's
0.379) and lost that position when the column changed — the two-cell recompute
simply could not see it. **Treat "k=15, w=15" as "a cell in the k ≥ 15 plateau",
not as the argmax.** The top five are within 0.06 of each other and the top two
within 0.008, so no cell in that band is meaningfully best.

**2. The estimator change is close to a no-op for this ranking.** Spearman
ρ = 0.9925 between the two orderings on both peak types; 4 of 20 cells change
rank, all inside the top five near-tie band. Every cell from rank 6 down holds
its position exactly, and Δ agrees to four decimals there (k10_w10 broad:
0.0705 vs 0.0704). The claim Exp A actually rests on — that k ≥ 15 separates
from k ≤ 10 — holds under both, and the gap is slightly *wider* on `_ie`:

| peak | gap, register-equality | gap, IE |
|---|---:|---:|
| broad | +0.3891 | +0.4146 |
| narrow | +0.4287 | +0.4490 |

Full separation (`min(xref) > max(diff)`) holds under `_ie` at every cell where
it holds under `jaccard_similarity`, which is the point made above: it is
satisfied almost everywhere and carries no information. It is reported in the
CSV as a diagnostic only.

**Figure 5 is unaffected.** Rendering it on `_ie`
(`paper/figures/cross_reference_identity_ie.png`) leaves all three tissue clades
monophyletic — reference choice remains a within-tissue perturbation.

## New Mode D columns — at k=10, w=10

`scripts/exp_a_metric_comparison.R` runs the same Wilcoxon test on each
of the 12 metric columns the new hammock emits (6 minimizer + 6
minimizer+ends). Outputs: `figures/metric_comparison_{broad,narrow}_k10_w10.{png,tsv}`.

> **Stale as a recommendation, 2026-08-06.** Half of this section — every
> "(minimizer+ends)" / "with ends" row, and finding 5's "`jaccard_similarity_with_ends`
> remains the right default" — is about a column family v0.6.0 deleted
> (`CLAUDE.md` divergence #8). Mode D emits one similarity block now, so the
> live rows are the five minimizer-only ones plus `jaccard_similarity_ie`,
> which did not exist when this ran. The measurements are kept as the record
> of what the two families did on this corpus; do not read finding 5 as
> current guidance. Findings 1, 3 and 6 still apply to the minimizer-only
> family. Redoing this comparison at k=15, w=15 (the lead cell) on the
> post-v0.6.0 schema is the outstanding item, not just the k-shift noted
> under "To extend".

![Metric comparison, broad, k=10, w=10](../figures/metric_comparison_broad_k10_w10.png)

![Metric comparison, narrow, k=10, w=10](../figures/metric_comparison_narrow_k10_w10.png)

### Findings (k=10, w=10)

| Metric (flavor)                | Δ broad | p broad | Δ narrow | p narrow |
|---|---|---|---|---|
| jaccard (minimizer)            | 0.071 | 1.35e-10 | 0.084 | 1.35e-10 |
| containment_AB/BA (minimizer)  | 0.032 | 1.36e-10 | 0.036 | 1.36e-10 |
| cosketch_geom (minimizer)      | 0.037 | 1.35e-10 | 0.044 | 1.35e-10 |
| cosketch_arith (minimizer)     | 0.036 | 1.35e-10 | 0.043 | 1.35e-10 |
| cosketch_max (minimizer)       | 0.017 | 1.35e-10 | 0.016 | 1.35e-10 |
| **jaccard (minimizer+ends)**   | **0.086** | **2.02e-07** | **0.087** | **4.54e-07** ← existing default |
| cosketch_geom (minimizer+ends) | 0.065 | 1.54e-07 | 0.064 | 9.94e-07 |
| cosketch_arith (minimizer+ends)| 0.062 | 2.02e-07 | 0.059 | 4.45e-06 |
| containment_AB/BA (with ends)  | 0.029 | 8.81e-04 | 0.028 | 6.62e-03 |
| cosketch_max (with ends)       | 0.013 | 5.13e-02 | 0.000 | 4.92e-01 |

Notes (both peak types):
1. **All minimizer-only metrics hit the Wilcoxon test floor** (p ≈ 1.4e-10 for n=18/54). Maximally significant but operating near saturation (medians ≈ 0.99 vs ≈ 0.92), so absolute Δ is small.
2. **Adding the start-end HLL lowers baselines into the interpretable mid-range** (~0.55–0.78) and reduces saturation. Cost: p drops by ~3–4 orders of magnitude for the surviving discriminators.
3. **Containment-with-ends is a markedly weaker discriminator** than jaccard-with-ends. `containment_AB` and `containment_BA` produce identical medians across the all-vs-all set because every unordered pair appears in both directions.
4. **cosketch_max collapses on with-ends** (broad Δ 0.013 p ≈ 0.05, narrow Δ ≈ 0 p ≈ 0.49). Worst metric in both flavors.
5. **`jaccard_similarity_with_ends` remains the right default** for headline analysis at k=10: best Δ/saturation trade-off and matches what the existing CSV sweep uses. `cosketch_geom_with_ends` is a near-tie if a redundancy-aware alternative is wanted.
6. **Broad vs narrow**: minimizer-only effects are slightly larger on narrow (likely tighter peak boundaries → more shared minimizers between same-sample-different-ref pairs). With-ends p-values are slightly weaker on narrow but the qualitative ordering of metrics is identical.

### To extend

- Replicate this 12-metric comparison at k=15, w=15 (the new headline cell) — `metric_comparison.R` is parameterised on the input CSV.
- The same minimizer / minimizer+ends pattern likely holds, but with-ends may matter less when the minimizer-only signal is already fully separated.

## Scripts that produced these outputs

- `scripts/exp_a_validate_plot.R` — per-cell Wilcoxon + 2-panel figure (boxplot + 9×9 heatmap)
- `scripts/exp_a_sweep_summary.R` — (k × w) effect-size heatmap
- `scripts/exp_a_dendrogram.R` — UPGMA dendrogram (broad + narrow)
- `scripts/exp_a_metric_comparison.R` — metric Wilcoxon comparison at one (k, w).
  The committed figure and TSV cover 12 metrics. ~~hammock 0.5.0 emits 14 (adding
  `jaccard_similarity_ie` and its `_with_ends` twin)~~ — **superseded 2026-08-06**:
  v0.6.0 dropped the whole `_with_ends` family, so current hammock emits **7**
  Mode D similarity columns (`jaccard_similarity`, `jaccard_similarity_ie`,
  `containment_AB/BA`, `cosketch_{geom,arith,max}`). The script picks up
  whatever columns are present once the input CSVs are regenerated. Not yet re-run
- `workflow/Snakefile` — orchestrates peaks → FASTA (bedtools getfasta) → hammock Mode D → R plots

Two of the four are Snakemake-only (they read `snakemake@input`/`@params` and
cannot be invoked standalone): `exp_a_validate_plot.R` and
`exp_a_sweep_summary.R`. The other two take positional argv (verified
2026-08-06 against their `commandArgs` blocks):

```bash
ml r/4.3.0
Rscript scripts/exp_a_dendrogram.R <broad_csv> <narrow_csv> <metadata_tsv> <out_png> [kw_label]
Rscript scripts/exp_a_metric_comparison.R <csv> <metadata_tsv> <out_png> <out_tsv> [peak_type_label]
```

`<metadata_tsv>` is `config/exp_a_metadata.tsv`. `exp_a_metric_comparison.R`
already tolerates post-v0.6.0 CSVs — it drops the missing `_with_ends`
columns from the figure and prints which ones it dropped.

## Reproducing

```bash
conda activate claude-ref-comparison
ml r/4.3.0
snakemake --profile workflow/slurm_profile/
```

(Upstream nf-core/chipseq runs must have produced the broad/narrow peak
calls under `NFCORE_OUT/{ref}{,_narrow}/...` — see `scripts/run_nfcore.sh`.)
