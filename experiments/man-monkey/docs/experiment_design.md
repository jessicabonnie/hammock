# man-monkey: Tissue-over-Species Clustering at Primate Distance

**Status:** Data scouting. No samples committed. Question: does the
tissue-over-species result that failed at ~80 Mya human↔mouse divergence
work at the much shorter primate divergence (~6 Mya chimp, ~25 Mya
macaque)?

**Created:** 2026-05-12.

---

## Motivation

The `mus-homo` experiment tests whether mark choice (H3K4me3 instead of
H3K27ac) closes the cross-species tissue-clustering gap at long divergence
(~80 Mya). This experiment tests whether **shorter divergence** closes the
gap at the original mark (H3K27ac) or any mark.

There are two reasons to think primate distance might recover the result:

1. **Sequence divergence is small enough** that orthologous peaks share
   enough k-mers to overlap at small (k, w). Human-chimp peak sequences
   are ~99 % identical; human-macaque ~93 %; even at H3K27ac's fast
   enhancer turnover, primate-primate H3K27ac peak overlap is much higher
   than primate-rodent (Villar 2015, Reilly 2015).
2. **The literature claim is grounded at primate distance.** Reilly 2015,
   Vermunt 2016, and Prescott 2015 all report tissue-specific epigenomic
   signatures shared across primates.

## Data-shape constraint

ENCODE has **zero** non-human-primate H3K27ac / H3K4me3 ChIP-seq. This
experiment cannot reuse the `mus-homo` ENCODE-download pipeline. All data
must come from published primate ChIP-seq, deposited at GEO / ArrayExpress
/ supplementary material.

**The pattern in the literature:** primate ChIP-seq is overwhelmingly
**single-tissue, multi-species** — researchers source primate tissue with
difficulty and tend to focus on one organ per study.

### Candidate datasets

| dataset | accession | species | tissue(s) | mark(s) | n_samples |
|---|---|---|---|---|---|
| Vermunt et al. 2016 (*Nat Neurosci*) | GSE85354 | human, chimp, macaque | prefrontal cortex | H3K27ac, H3K4me1 | ~3–5 per species |
| Reilly et al. 2015 (*Science*) | GSE63648 | human, chimp, macaque | fetal cortex | H3K27ac, H3K4me2, H3K27me3 | ~3 per species |
| Prescott et al. 2015 (*Cell*) | GSE70751 | human, chimp, bonobo | cranial neural crest cells | H3K27ac, H3K4me1, H3K27me3 | ~2–3 per species |
| Roller et al. 2021 (*Nat Comm*) | ?, see Odom lab GitHub | mammals incl. human, chimp, macaque, marmoset | liver | H3K27ac | varies |
| Berthelot et al. 2018 (*Nat Comm*) | E-MTAB-2633 (Villar) extension | mammals incl. primates | liver | H3K27ac, H3K4me3 | varies |
| Garcia-Perez et al. 2018 / Marchetto 2013 | various | human, chimp, bonobo | iPSC, cranial NCC | various | varies |

(Accession numbers above are approximate — verify before use.)

## Design options

Three productive shapes depending on which dataset we commit to:

### Option A: Vermunt 2016 PFC

- 3 species × ~3 individuals × 1 tissue (PFC)
- **Question:** does H3K27ac sketch similarity within PFC pick up
  species identity in a way that mirrors the known phylogeny? Or do
  individual differences dominate at primate distance?
- **Strengths:** clean dataset, well-annotated, peaks in supplementary.
- **Weaknesses:** single tissue, so no "tissue vs species" question.

### Option B: Reilly 2015 fetal cortex + Vermunt 2016 PFC + Prescott 2015 CNCC

- 3 species × 3 tissues (~9–18 samples total)
- **Question:** does sketch cluster by tissue across species at primate
  distance, where the H3K27ac null at mouse distance was the alternative?
- **Strengths:** matched-species design across 3 tissues. Real
  tissue-vs-species at primate distance.
- **Weaknesses:** tissues differ in developmental stage (Reilly = fetal,
  Vermunt = adult, Prescott = stem-cell-derived). Strong methods caveat.

### Option C: Villar 2015 + Berthelot 2018 primates only, liver

- 4–6 primates × 1 tissue (liver)
- **Question:** species recovery within primates in liver — does sketch
  reconstruct the primate clade structure?
- **Strengths:** all peaks deposited and ready, same mark, same
  developmental stage. Cleanest dataset.
- **Weaknesses:** **single tissue**, so this is really a phylogeny
  question — duplicates `primate-phylogeny` on a subset of species.
  Maybe not its own experiment.

## Recommended direction

**Option B** is the only one that actually answers "tissue over species at
primate distance" — the original question that motivated the project. The
developmental-stage heterogeneity is a real concern that needs to be
addressed in design (e.g., either restrict to adult Vermunt + adult chimp /
macaque PFC equivalents from other studies, or accept the heterogeneity
and weight findings accordingly).

**Option A** is a useful sanity check that's complete on its own — even if
B isn't viable, A produces a publishable "does sketch cluster individuals
within / between species in a single tissue?" result.

## Pipeline shape (proposed)

Same as `mus-homo`: download peak BEDs from supplementary materials,
bedtools getfasta from the matching reference, hammock all-vs-all, cluster.

Additional complexity vs. `mus-homo`:
- Non-human primate reference genomes needed (panTro6, rheMac10, etc.)
- Peak BEDs may be in older assembly coords — liftOver as needed
- Peak files come from supplementary materials, not a single API endpoint

## Next steps

1. **Inventory data:** for each candidate dataset, pull the actual peak BED files (or links to them), count samples per species/tissue, note assembly versions.
2. **Decide on Option A, B, or C** based on what the inventory shows.
3. **Acquire primate reference FASTAs** — UCSC has them all.
4. **Build samples manifest.**
5. **Write Snakefile.**

## Open questions

- Is "developmental stage heterogeneity" across studies a blocker, or
  acceptable with caveats?
- Should we include human samples sourced from ENCODE to match each
  primate tissue, or only from the same studies as the primate samples
  (which would constrain to the same protocol)?
- Marker choice: H3K27ac or H3K4me3? Probably H3K27ac for tissue
  specificity within primates (where divergence is small enough that
  enhancer turnover is less of a problem).
