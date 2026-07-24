# Target-matched BED-file subsets

## Goal

Define biologically coherent groups of BED files that a researcher might reasonably compare in one analysis. The purpose is not merely to report repository size, but to give a concrete example of the number of coordinate-based datasets involved in a realistic similarity search or clustering task.

## Primary example: human CTCF ChIP-seq

CTCF is the preferred first example because:

- the target is a single, well-defined DNA-binding protein;
- peak BED files have a consistent biological interpretation across experiments;
- CTCF has been profiled across many cell types, treatments, studies, and laboratories;
- both hg19 and hg38 records are common, exposing the reference-coordinate problem directly;
- comparing CTCF occupancy profiles across cell types is a plausible clustering, quality-control, retrieval, and annotation task.

The primary subset should include records meeting all of the following criteria:

1. organism is *Homo sapiens*;
2. assay is ChIP-seq or a repository-equivalent protein-DNA binding assay;
3. normalized target is exactly `CTCF`;
4. the record has an experiment- or sample-level peak BED file;
5. the assembly is known;
6. the BED represents called peaks rather than a repository-wide aggregate catalog;
7. controls, input DNA, blacklist files, and signal tracks are excluded.

Counts should be reported at three levels:

- physical peak BED files;
- distinct experiments or processed samples represented;
- total genomic intervals across the included files, when available.

The results should be stratified by:

- repository;
- reference assembly;
- cell-type class;
- source study;
- original versus lifted coordinates, when the repository provides lifted files.

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

| Target-matched comparison set | Repository | Assembly | Peak BED files | Distinct experiments/samples | Cell types | Studies | Intervals |
|---|---|---|---:|---:|---:|---:|---:|
| Human CTCF ChIP-seq | ChIP-Atlas | hg38 | — | — | — | — | — |
| Human CTCF ChIP-seq | ChIP-Atlas | hg19 | — | — | — | — | — |
| Human CTCF ChIP-seq | Cistrome DB | hg38 | — | — | — | — | — |
| Human CTCF ChIP-seq | ReMap | hg38 | — | — | — | — | — |
| Human H3K27ac ChIP-seq | ChIP-Atlas | hg38 | — | — | — | — | — |

## Pairwise-comparison interpretation

For a subset containing `n` BED files, an exhaustive all-pairs analysis requires `n(n-1)/2` file comparisons. The paper may report both the number of files and the corresponding number of pairwise comparisons, but the comparison count must be presented as a derived value rather than a repository statistic.

The summary also distinguishes:

- comparisons possible within hg19;
- comparisons possible within hg38;
- cross-assembly pairs that cannot be compared directly without conversion or a reference-independent method;
- the hypothetical total if all coordinate systems were directly comparable.

## ChIP-Atlas implementation

The script `scripts/count_chip_atlas_target_subset.py` downloads the public ChIP-Atlas experiment list and selects records using exact, case-insensitive matching on:

- assembly: `hg19` or `hg38`;
- antigen class: `TFs and others`;
- antigen subclass: `CTCF`.

It treats one selected ChIP-Atlas experiment as one experiment-level peak BED at one specified q-value threshold. This avoids counting the same biological experiment three times merely because ChIP-Atlas publishes three thresholded BED files.

Run from the repository root:

```bash
python paper/database_counts/scripts/count_chip_atlas_target_subset.py
```

The script writes:

- `results/chip_atlas_ctcf_human_manifest.tsv`: every included experiment and its BED URL;
- `results/chip_atlas_ctcf_human_summary.tsv`: counts by assembly and derived comparison workloads;
- `results/chip_atlas_ctcf_human_provenance.json`: retrieval date, source URL, filters, and counting unit.

The same script can produce the secondary histone-mark examples by changing the target and antigen class, for example:

```bash
python paper/database_counts/scripts/count_chip_atlas_target_subset.py \
  --target H3K27ac \
  --antigen-class Histone
```

## Validation

The extracted experiment totals should be checked against the ChIP-Atlas antigen-count endpoint for the same assembly, antigen class, cell-type scope, and retrieval date. The manifest remains the authoritative audit trail because it records the included experiment identifiers rather than only an aggregate count.
