# mus-homo: Tissue-over-Species Clustering, Human + Mouse

**Status:** ✅ Run completed 2026-05-13. **Result: partial improvement over
H3K27ac but still not a clean tissue-over-species recovery.** See "Results
(2026-05-13)" section below.

**Created:** 2026-05-12 (split from `claude-ref-comparison/exp_b`).
**Locked:** 2026-05-12 — DNase-seq, 5 tissues (testis dropped), 10 samples.
**First run:** 2026-05-13.

---

## Motivation

The original `claude-ref-comparison` Experiment B asked whether hammock could
recover the well-known tissue-over-species clustering result of Yue 2014 /
Lin 2014 / Roadmap 2015 using H3K27ac peak FASTAs. After two rounds of
sample expansion (6 → 12 samples) and two orthology-aware preprocessing
redesigns (B2 conservation filter; B1 liftover, not completed), the answer
was **null**: 1 of 80 (k, w) cells was tissue-dominant, in the noise floor.

The most likely explanation, supported by Villar 2015 (Cell), is the
**choice of histone mark**. H3K27ac marks active enhancers — the
fastest-evolving class of regulatory element in vertebrate genomes — with
~10–20 % peak overlap between human and mouse liver. H3K4me3, by contrast,
marks active promoters anchored to orthologous TSSs, with ~60–70 % peak
overlap human↔mouse. A sketch over H3K4me3 peak FASTAs should have a
genuine chance to recover tissue identity across species.

## Hypothesis

> If DNase-seq open-chromatin peaks are sketched at small minimizer
> parameters (k = 5–10) and clustered, the clustering will recover
> **tissue** identity over **species** identity at one or more (k, w) cells
> across a 5-tissue × 2-species design with 10 samples.

This is the same null hypothesis that failed for H3K27ac, but with a data
type whose underlying biology — sequence overlap of orthologous peaks
across species — actually supports it. DNase-seq captures the full
regulatory repertoire (promoters + enhancers + insulators), giving the
sketch much richer input than any single ChIP-seq mark.

## Why DNase-seq (and not H3K4me3, ATAC-seq, etc.)

Three data types were considered and ranked:

1. **DNase-seq (chosen)** — open chromatin, all regulatory elements. ENCODE
   has matched human + mouse for 5/6 candidate tissues; mouse testis is
   unavailable but every other tissue lines up cleanly.
2. **ATAC-seq** — also open chromatin but mouse adult ATAC-seq for these
   tissues is patchier than DNase-seq on ENCODE; mouse testis ATAC-seq is
   also unavailable.
3. **H3K4me1** — enhancer-flavored ChIP-seq, more conserved across species
   than H3K27ac. Available for all 6 tissues, but a less-rich signal than
   open chromatin.
4. **H3K4me3** — promoter mark, most cross-species conserved, but
   tissue-discrimination is weak (most TSSs active across tissues).

## Design

**Samples:** 10 matched-tissue ENCODE DNase-seq experiments (5 tissues × 2 species).

### Human (Stamatoyannopoulos lab, single-replicate adult, GRCh38)

| Tissue | ENCSR | Donor |
|---|---|---|
| heart  | ENCSR792QQE | 54M heart RV |
| liver  | ENCSR909HFI | 53F right lobe |
| lung   | ENCSR757NJT | 60M lower left lobe |
| spleen | ENCSR567EEO | 61F spleen |
| brain  | ENCSR000EIY | 67F frontal cortex (Crawford, Duke) |

### Mouse (Stamatoyannopoulos lab, 2-month adult male, mm10)

| Tissue | ENCSR | Strain |
|---|---|---|
| heart  | ENCSR671NPO | B6CASTF1/J |
| liver  | ENCSR000CNI | C57BL/6J |
| lung   | ENCSR702AVR | B6CASTF1/J |
| spleen | ENCSR087PCD | B6CASTF1/J |
| brain  | ENCSR325RZJ | B6CASTF1/J frontal cortex |

### Caveats to flag in methods

- **Testis dropped.** ENCODE has 0 adult mouse testis DNase-seq (only 2
  postnatal samples at 10–14 days). Testis is highly tissue-distinct and
  its absence costs the design some discriminative power.
- **Mouse strain mix.** 4 of the 5 mouse samples are B6CASTF1/J; liver is
  C57BL/6J because it's the only clean (untreated) adult option. Minor
  strain heterogeneity vs species-level signal.
- **Human brain lab outlier.** All human samples are Stamatoyannopoulos
  lab *except* brain, which is Greg Crawford (Duke). Necessary because
  Crawford's lab is the source of human frontal-cortex DNase-seq matched
  to the mouse anatomy choice.
- **Single replicate per sample.** DNase-seq experiments are typically
  single-rep; this is the norm and isn't a design weakness.

## Pipeline shape

Crucial simplification over `ref-comparison`: **no nf-core/chipseq needed.**
ENCODE publishes uniform-pipeline-called peak BEDs as direct downloads.
For each ENCSR experiment, one `*.bed.gz` file (~1–5 MB) replaces the
multi-TB FASTQ + alignment workflow.

```
ENCODE peak BED (curl) → bedtools getfasta → hammock all-vs-all → cluster
```

Snakemake rules (built in `workflow/Snakefile`):
1. `download_peak_bed` — `curl` each peak ENCFF to `peaks_cache/`
2. `peak_bed_to_fasta` — `bedtools getfasta` from matching reference (stripped to BED3 for safety)
3. `fasta_list` — write the 10-sample FASTA list
4. `mus_homo_hammock` — hammock Mode D, (k, w) sweep
5. `cluster_and_plot` — `scripts/cluster_plot.R`: hierarchical clustering, dendrogram, PCA, cluster_assignments.tsv

(k, w) sweep: `k ∈ {5, 8, 10, 15, 20}` × `w ∈ {5, 8, 10, 15, 20, 30}` filtered to `w ≥ k` (20 cells).

## Success criteria

The experiment can land in one of three places:

| outcome | what it means |
|---|---|
| ≥ 5 of 40 (k, w × peak_type) cells tissue-dominant; tissue ARI ≥ 0.3 in those cells | **Positive** — sketch + H3K4me3 recovers tissue signal at human↔mouse divergence. Strong story. |
| 1–4 cells tissue-dominant, near-noise margins | **Marginal** — mark choice matters but H3K4me3 alone doesn't close the gap. Frame as a divergence-ceiling result. |
| 0 cells tissue-dominant | **Null** (parallel to H3K27ac result) — strong evidence that *no* peak-FASTA sketching recovers cross-species tissue identity at ~80 Mya, regardless of mark. Useful negative result. |

## Results (2026-05-13)

Full (k, w) sweep over 10 samples completed. 62/62 snakemake jobs done; 20
hammock CSVs, 20 dendrograms, 20 PCAs, 20 cluster_assignments TSVs.

### Headline numbers

| metric | mus-homo DNase-seq (N=10) | comparison: claude-ref-comparison/exp_b H3K27ac (N=12) |
|---|---|---|
| (k, w) cells in sweep | 20 (broad+narrow are one assay here) | 40 (broad + narrow) |
| tissue-dominant cells | **3 of 20 (15%)** | 1 of 40 (2.5%) |
| max tissue ARI in any cell | **+0.069** | +0.057 |
| species ARI saturation | **0.374–0.476** | 0.555 |

### Tissue-dominant cells

| k | w | ARI_t | ARI_s | margin | interpretation |
|---|---|---|---|---|---|
| 5 | 8  | +0.069 | +0.032 | +0.037 | both near zero — noise |
| 5 | 20 | +0.040 | −0.062 | +0.102 | species ARI **negative** — sketch fails to find species, tissue barely above zero |
| 5 | 30 | +0.040 | −0.062 | +0.102 | same as k=5,w=20 |

### Interpretation

The result is **marginally better than H3K27ac but still not a meaningful
tissue-over-species recovery.** Two effects from DNase-seq:

1. **Species ceiling drops.** Where H3K27ac saturates at ARI species ≈ 0.555,
   DNase-seq saturates at 0.374–0.476. The richer regulatory repertoire
   (open chromatin captures promoters + enhancers + insulators) makes the
   species axis less starkly separable.
2. **Tissue signal does NOT meaningfully emerge.** The 3 "wins" are
   technical victories — cases where the species ARI drops to zero or
   slightly negative while tissue ARI stays in the +0.04 to +0.07 range.
   By the convention `ARI ≥ 0.3 = real clustering`, none of the 3 wins
   qualify.

The k=5 corner remains the only region where any tissue signal appears,
just as for H3K27ac. The qualitative pattern is identical; only the
absolute numbers shift slightly.

### Verdict

Lands in the **"marginal"** outcome category (1–4 cells tissue-dominant,
near-noise margins) rather than positive. Data-type choice (DNase-seq vs
H3K27ac) makes a measurable but not decisive difference. The remaining
hypothesis is that sequence-level peak-FASTA sketching at human↔mouse
divergence fundamentally cannot recover tissue clustering — regardless of
mark/data-type — and the only path forward is either:
- coordinate-space liftover (B1-style)
- or shorter divergence (primate, `man-monkey` experiment)

### Output paths

- ARI summary table (regenerable): `python3 -c "..."` snippet from this session
- Per-cell figures: `results/k{k}_w{w}/{dendrogram.png, pca.png, cluster_assignments.tsv}`
  (where `results/` symlinks to `/vast/blangme2/jbonnie/hammock/mus-homo/results/`)

### Addendum (2026-05-13): column choice does not change the verdict

Re-clustered all 20 (k, w) cells with `jaccard_similarity` (minimizers only)
in addition to `jaccard_similarity_with_ends` (the default). Max tissue ARI
is identical (+0.069), 0/20 cells clear the 0.3 "real clustering" threshold
under either column, and the tissue-dominant count only nudges 3 → 4 (a
k=5,w=5 cell squeaks in at ARI_t = +0.040, still noise). The only material
difference is that `_with_ends` carries slightly more species signal at high
k (species ARI 0.374 vs 0.308 at k=15, w≥15). Per-cell comparison TSV at
`results/column_comparison.tsv`.

---

## Open questions / future work

- **H3K4me1 parallel arm.** Could compare "open chromatin" (DNase-seq) vs
  "enhancer mark" (H3K4me1) head-to-head on the same tissue panel. Likely
  to land in the same marginal/null regime but useful as a head-to-head.
  Deferred.
- **Coordinate-space liftover (B1).** Apply the lifted-coordinate approach
  to DNase-seq here, since it might give a stronger result than the
  abandoned-mid-run H3K27ac B1 experiment. Deferred.
- **Move to primate distance.** Run `man-monkey` to test whether shorter
  divergence (~6–25 Mya) actually recovers the tissue signal.
