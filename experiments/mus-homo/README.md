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

## Success criteria

| outcome | what it means |
|---|---|
| ≥ 5 of 20 (k, w) cells tissue-dominant; tissue ARI ≥ 0.3 in those cells | **Positive** — sketch + DNase-seq recovers tissue at human↔mouse divergence. The signal was there; previous H3K27ac null was a mark-choice problem. |
| 1–4 cells tissue-dominant, near-noise margins | **Marginal** — mark/data-type matters but doesn't fully close the gap. Frame as a divergence ceiling. |
| 0 cells tissue-dominant | **Null** — sequence-level sketching of regulatory peak FASTAs cannot recover cross-species tissue identity at ~80 Mya, full stop. Strong negative result. |
