# Target-matched BED-file subsets

## Goal

Define biologically coherent groups of BED files that a researcher might reasonably compare in one analysis. The purpose is not merely to report repository size, but to give a concrete example of the number of coordinate-based datasets involved in a realistic similarity search or clustering task.

## Primary example: CTCF ChIP-seq across all references

CTCF is the preferred first example because:

- the target is a single, well-defined DNA-binding protein;
- peak BED files have a consistent biological interpretation across experiments;
- CTCF has been profiled across many cell types, treatments, studies, laboratories, species, and assemblies;
- multiple assembly generations are represented, exposing the reference-coordinate problem directly;
- comparing CTCF occupancy profiles is a plausible clustering, quality-control, retrieval, and annotation task.

The primary subset should include records meeting all of the following criteria:

1. assay is ChIP-seq or a repository-equivalent protein-DNA binding assay;
2. normalized target is exactly `CTCF`;
3. the record has an experiment- or sample-level peak BED file;
4. the assembly is known;
5. the BED represents called peaks rather than a repository-wide aggregate catalog;
6. controls, input DNA, blacklist files, and signal tracks are excluded.

The default analysis retains **all species and all non-empty assembly labels**. Human-only, species-specific, or single-assembly analyses are derived views created with explicit filters. This avoids hiding the reference fragmentation that the analysis is intended to quantify.

Counts should be reported at three levels:

- physical peak BED files;
- distinct experiments or processed samples represented;
- total genomic intervals across the included files, when available.

The results should be stratified by:

- repository;
- species;
- reference assembly;
- cell-type class;
- source study;
- original versus lifted coordinates, when the repository provides lifted files.

## Controlled derived views

The full all-reference inventory should be accompanied by narrower subsets that answer different questions:

- all human CTCF experiments across every human assembly;
- hg19 and hg38 separately for within-reference benchmarking;
- all mouse CTCF experiments across mm9 and mm10;
- one organism and one assembly for a strictly coordinate-compatible benchmark;
- cross-assembly or cross-species target-matched sets for demonstrating reference-independent comparison.

These subsets should be generated from the full manifest rather than collected independently.

## Secondary example: human H3K27ac ChIP-seq

H3K27ac is a strong secondary example because it is widely used to identify active enhancers and promoters. A researcher may plausibly compare all H3K27ac peak sets for tissue clustering, enhancer conservation, dataset retrieval, or quality control.

It is less tightly controlled than CTCF because H3K27ac profiles are highly responsive to cell state, treatment, peak-calling choices, and broad-versus-narrow peak representations. For that reason, H3K27ac should be reported separately and should not be mixed with CTCF in a single similarity benchmark.

Selection criteria should mirror the CTCF criteria, with the normalized target set to `H3K27ac` and the output restricted to peak BEDs.

## Additional candidate subsets

- `H3K4me3`: promoter-associated mark with comparatively sharp peaks and broad use across cell types.
- `H3K27me3`: broad repressive mark; useful for testing methods on long and diffuse intervals but more sensitive to file representation.
- ATAC-seq within a single tissue or cell-type class: same assay objective rather than same molecular target.
- CTCF within one cell-type class: a stricter subset for estimating the number of near-replicate or closely related BED files.

## Target normalization

Target labels must be normalized before counting. The manifest should preserve both the source label and normalized label.

For CTCF, acceptable source labels may include case and punctuation variants that unambiguously identify CTCF. Records involving fusion proteins, perturbation labels used as targets, multi-target mixtures, or ambiguous antibody descriptions should be excluded or reviewed manually.

For histone marks, normalization should standardize capitalization and punctuation while preserving biologically meaningful distinctions such as acetylation versus methylation state.

## Recommended manuscript table

| Target-matched comparison set | Repository | Species | Assembly | Peak BED files | Distinct experiments/samples | Cell types | Studies | Intervals |
|---|---|---|---|---:|---:|---:|---:|---:|
| CTCF ChIP-seq | ChIP-Atlas | Human | hg38 | — | — | — | — | — |
| CTCF ChIP-seq | ChIP-Atlas | Human | hg19 | — | — | — | — | — |
| CTCF ChIP-seq | ChIP-Atlas | Mouse | mm10 | — | — | — | — | — |
| CTCF ChIP-seq | ChIP-Atlas | Mouse | mm9 | — | — | — | — | — |
| CTCF ChIP-seq | ChIP-Atlas | Other | Other assemblies | — | — | — | — | — |
| Human CTCF ChIP-seq | Cistrome DB | Human | hg38 | — | — | — | — | — |
| Human CTCF ChIP-seq | ReMap | Human | hg38 | — | — | — | — | — |
| Human H3K27ac ChIP-seq | ChIP-Atlas | Human | hg38 | — | — | — | — | — |

## Pairwise-comparison interpretation

For a subset containing `n` BED files, an exhaustive all-pairs analysis requires `n(n-1)/2` file comparisons. The paper may report both the number of files and the corresponding number of pairwise comparisons, but the comparison count must be presented as a derived value rather than a repository statistic.

The summary should distinguish:

- comparisons possible directly within the same assembly;
- cross-assembly pairs that conventional coordinate overlap cannot compare directly;
- all possible pairs if a reference-independent representation is used.

## ChIP-Atlas implementation

The script `scripts/count_chip_atlas_target_subset.py` downloads the complete ChIP-Atlas experiment list and retains all non-empty assembly labels by default.

Run the full CTCF inventory:

```bash
python paper/database_counts/scripts/count_chip_atlas_target_subset.py
```

Create a human-only derived view:

```bash
python paper/database_counts/scripts/count_chip_atlas_target_subset.py \
  --species "Homo sapiens"
```

Create a conventional hg19/hg38 view:

```bash
python paper/database_counts/scripts/count_chip_atlas_target_subset.py \
  --assemblies hg19 hg38
```

Run the same workflow for H3K27ac:

```bash
python paper/database_counts/scripts/count_chip_atlas_target_subset.py \
  --target H3K27ac \
  --antigen-class Histone
```

The script preserves unknown assembly labels rather than silently dropping them. Species are inferred from a documented assembly-to-species mapping; unrecognized assemblies remain in the manifest with species labeled `unknown` for later review.

## Writing notes

Candidate motivating language developed from this analysis is collected in `motivating_statements.md`. Statements should not be moved into the manuscript until the relevant counts or published sources are available to support them.
