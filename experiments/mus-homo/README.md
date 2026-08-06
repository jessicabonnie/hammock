# mus-homo — Tissue-over-Species Clustering (Human + Mouse, DNase-seq)

Tests whether minimizer-HLL sketching of DNase-seq peak FASTAs recovers the
literature's tissue-over-species clustering result at human↔mouse divergence
(~80 Mya).

Replaces the H3K27ac null result from the legacy `claude-ref-comparison`
Experiment B with **DNase-seq + brain** — open chromatin captures the full
regulatory repertoire (promoters + enhancers + insulators) with higher
cross-species peak overlap than H3K27ac, and brain is the most
transcriptomically distinct major tissue.

## Sample design

5 tissues × 2 species = **10 samples** (testis dropped — no adult mouse DNase-seq on ENCODE).

| tissue | human (GRCh38) | mouse (mm10, 2-month adult) |
|---|---|---|
| heart  | ENCSR792QQE | ENCSR671NPO |
| liver  | ENCSR909HFI | ENCSR000CNI |
| lung   | ENCSR757NJT | ENCSR702AVR |
| spleen | ENCSR567EEO | ENCSR087PCD |
| brain (frontal cortex) | ENCSR000EIY | ENCSR325RZJ |

Full sample manifest with peak BED URLs: `config/samples.tsv`.

## Pipeline shape

No nf-core/chipseq required. ENCODE publishes uniform-pipeline peak BEDs
as direct downloads, so the pipeline is much simpler than `ref-comparison`:

```
ENCODE peak BED (curl)  →  bedtools getfasta  →  hammock all-vs-all  →  cluster + plot
        ~1–5 MB / sample       seconds                10 min – 12 h         <2 min
```

> **Mode D timing note, 2026-08-06.** That "10 min – 12 h" was measured at
> `--threads 8`, as every archived Mode D run in this repo was. Mode D
> threading is a GIL convoy: v0.6.1 defaulted Mode D to `--threads 1` after
> measuring the 8-thread pool at **2× slower** (`CLAUDE.md`, Architecture).
> Treat the figure as inflated and re-baseline before sizing a new run.

(k, w) sweep: same as `ref-comparison` (`k ∈ {5, 8, 10, 15, 20}` × `w ∈ {5, 8, 10, 15, 20, 30}` filtered to `w ≥ k` = 20 cells).

## Running

```bash
conda activate claude-ref-comparison       # same env as ref-comparison; has hammock + snakemake
ml r/4.3.0                                  # clustering step is scripts/cluster_plot.R
snakemake --dryrun                          # check job plan
snakemake --profile workflow/slurm_profile/ # submit to SLURM (symlinked from ref-comparison)
# or, as one batch job:
sbatch sbatch_workflow.sh
```

Outputs land at `/vast/blangme2/jbonnie/hammock/mus-homo/results/`, one
`k{k}_w{w}/` directory per cell.

Scripts not wired into the Snakefile:

| script | what it does |
|---|---|
| `scripts/compute_column_comparison.{py,R}` | re-clusters every cell under each similarity column and tabulates tissue/species ARI (the 2026-05-13 addendum in `docs/experiment_design.md`) |
| `estimator_ie_tissue.py` | asks whether reading `jaccard_similarity_ie` (derived from the shipped containment columns — no re-sketching) rescues the cross-species tissue separation. Run as `python estimator_ie_tissue.py`. Conclusion not recorded here. |

## Pre-registered success criteria

| outcome | what it means |
|---|---|
| ≥ 5 of 20 (k, w) cells tissue-dominant; tissue ARI ≥ 0.3 in those cells | **Positive** — sketch + DNase-seq recovers tissue at human↔mouse divergence. The signal was there; previous H3K27ac null was a mark-choice problem. |
| 1–4 cells tissue-dominant, near-noise margins | **Marginal** — mark/data-type matters but doesn't fully close the gap. Frame as a divergence ceiling. |
| 0 cells tissue-dominant | **Null** — sequence-level sketching of regulatory peak FASTAs cannot recover cross-species tissue identity at ~80 Mya, full stop. Strong negative result. |

## Results (2026-05-14)

Sweep complete; 20 of 20 (k, w) cells produced dendrogram, PCA, and
cluster-assignment outputs. Clustering uses `complete`-linkage UPGMA on
`1 − jaccard_similarity_with_ends` cut into k = 5 clusters (matching the
number of true tissues), per `config/config.yaml`'s `primary_sim_col`.
All numbers below are post the Mode D `add_string` ingestion bug fix
(`CLAUDE.md` "Intentional divergences" §6).

> **Column note, 2026-08-06.** `jaccard_similarity_with_ends` was removed in
> hammock v0.6.0 (`CLAUDE.md` divergence #8), and `config/config.yaml`'s
> `primary_sim_col` has since been changed to `jaccard_similarity` — so the
> sentence above no longer describes what a re-run would do, only what the
> archived outputs were built from. The verdict does not move: the
> 2026-05-13 addendum in `docs/experiment_design.md` re-clustered all 20
> cells on `jaccard_similarity` and found the same max tissue ARI (+0.069)
> and the same 0/20 cells clearing the 0.3 threshold. Parse the archived
> CSVs by column *name*; the schema changed and field indices do not line up.
>
> Also note `jaccard_similarity` is register-equality, not set Jaccard — it
> carries a chance-agreement floor and is not rank-faithful across pairs of
> differing set size (divergence #2), which matters here because the
> human and mouse peak FASTAs differ in size. `estimator_ie_tissue.py` in
> this directory asks directly whether reading the calibrated
> `jaccard_similarity_ie` instead rescues the tissue signal; its answer is
> not recorded in these docs.

> **Cross-species Mode D caveat.** Human↔mouse Mode D Jaccard measures
> **shared k-mer content**, which at ~80 Mya is dominated by repeats and
> low-complexity sequence — it is not a homology measure (`CLAUDE.md`,
> "Cross-reference caveat"). hammock's default `k=8` is unsuitable for
> cross-species work; the sweep here runs k up to 20 for that reason.
> This is the same mechanism the Verdict section below identifies.

### Cluster structure across the (k, w) sweep

| regime | cells | clustering structure | tissues paired across species |
|---|---|---|---|
| **k ≥ 15** | 5 (k15_w15/20/30, k20_w20/30) | Perfect species partition — all 5 humans in clusters 1–3, all 5 mice in clusters 4–5 | **0 of 5** |
| **k = 10** | 4 (k10_w10/15/20/30) | Humans cluster in 1–2 groups; mice mostly singletons. Cluster 2 happens to include both lung samples (h + m) alongside 2 other samples | 1 of 5 (lung) |
| **k = 8** | 5 (k8_w8/10/15/20/30) | Same shape as k = 10; humans grouped, mice scattered, one accidental tissue coincidence per cell (typically lung) | 1 of 5 |
| **k = 5** | 6 (k5_w5/8/10/15/20/30) | Noisier; at high w one cluster collapses 5–6 mixed samples and accidentally pairs heart and liver | 0–2 of 5 |

Where a tissue does "pair" at small k, it sits inside a 4–6-sample
mixed cluster — the pairing is incidental, not a tissue-driven structure.
No (k, w) cell satisfies the pre-registered "tissue ARI ≥ 0.3" threshold;
at k ≥ 15 the tissue ARI is essentially zero while the species ARI is
1.0.

### Verdict — Null

Sequence-level minimizer-HLL sketching of DNase-seq peak FASTAs does
not recover cross-species tissue identity at ~80 Mya divergence. The
result lands in the pre-registered "0 cells tissue-dominant → strong
negative result" bucket. It is consistent with the cross-species
observation in `experiments/ref-comparison/` (within-human cross-reference
worked at high effect size; cross-species H3K27ac failed on the same
Mode D backbone); the DNase-seq + brain refinement does not lift it.

This is a property of k-mer sketching at long evolutionary distance
rather than a hammock-specific limitation — at minimizer-set turnover
this large any k-mer-set comparator (Mash, Sourmash, MinHash on the
same inputs) faces the same saturation. The practical conclusion is
that coordinate-space liftover (or some other orthology-aware
preprocessing) is the appropriate primitive for cross-species peak
comparison at this divergence, not direct sequence sketching.

Per-cell outputs (`dendrogram.png`, `pca.png`, `cluster_assignments.tsv`,
and the input Mode D CSV) are at
`/vast/blangme2/jbonnie/hammock/mus-homo/results/k{k}_w{w}/`.

### Paper-outline status

This experiment is intentionally excluded from `docs/paper_outline.md`;
the README is kept here as the internal record of the negative result.
