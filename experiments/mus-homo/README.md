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

(k, w) sweep: same as `ref-comparison` (`k ∈ {5, 8, 10, 15, 20}` × `w ∈ {5, 8, 10, 15, 20, 30}` filtered to `w ≥ k` = 20 cells).

## Running

```bash
conda activate claude-ref-comparison       # same env as ref-comparison; has hammock + snakemake
snakemake --dryrun                          # check job plan
snakemake --profile workflow/slurm_profile/ # submit to SLURM (symlinked from ref-comparison)
```

Outputs land at `/vast/blangme2/jbonnie/hammock/mus-homo/results/`.

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
