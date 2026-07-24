# Cross-database reference mismatch: Roadmap versus BLUEPRINT H3K27ac

## Why this example is more striking

This case study compares two distinct human epigenomics consortia that profiled the same histone modification but distributed their BED peak files on different reference assemblies:

- NIH Roadmap Epigenomics H3K27ac peaks on hg19/GRCh37
- BLUEPRINT H3K27ac peaks on hg38/GRCh38

Unlike alternate-reference copies of the same experiment, these collections contain different biological samples, donors, tissues, diseases, and cell types. A researcher studying active regulatory regions may reasonably want to search or cluster them together, but direct BED overlap is not defined across their original coordinate systems.

## Analysis question

How many Roadmap-by-BLUEPRINT H3K27ac BED-file pairs represent biologically related datasets that cannot be compared directly because one file is on hg19 and the other is on hg38?

## Counting policy

Use FILER's harmonized metadata tables as the source inventory because they expose:

- data source;
- genome build;
- assay;
- antibody/target;
- file format;
- processed download URL;
- cell type and tissue category;
- number of intervals.

Retain records satisfying all of the following:

1. `assay` is ChIP-seq;
2. `antibody` normalizes exactly to `H3K27ac`;
3. `output_type` is peaks;
4. the processed file is BED-like;
5. Roadmap records are native hg19 records;
6. BLUEPRINT records are native hg38 records;
7. lifted Roadmap copies are excluded from the primary count;
8. duplicate processed URLs are removed.

## Quantities to report

For each resource:

- BED-file count;
- total interval count;
- median intervals per file;
- cell-type count;
- tissue-category count;
- directly comparable within-resource pair count.

Across resources:

- Roadmap BED count multiplied by BLUEPRINT BED count;
- total Roadmap-by-BLUEPRINT pairs blocked by reference mismatch;
- optional counts stratified by tissue category or normalized cell type.

## Interpretation

The key result is not that hg19 and hg38 both exist. It is that two complementary collections of H3K27ac experiments cannot be included in one conventional coordinate-overlap analysis without first transforming one collection.

Candidate motivating statements:

> Distinct epigenomic consortia can profile the same regulatory mark in complementary biological samples while publishing their peak sets on incompatible reference assemblies.

> A researcher attempting to compare Roadmap and BLUEPRINT H3K27ac profiles cannot directly evaluate any cross-resource BED pair in the native files because Roadmap peaks are represented on hg19 whereas BLUEPRINT peaks are represented on hg38.

> Reference mismatch therefore excludes an entire block of biologically relevant comparisons, rather than merely duplicating a small number of files.

> Lift-over can create a common coordinate system, but it changes the analysis from direct comparison of the published interval sets to comparison after a reference-dependent coordinate transformation.

## Source notes

FILER publishes downloadable metadata templates separately for native hg19 and native hg38 tracks, and documents that its tracks are represented as sorted, indexed BED files. FILER also provides a distinct hg38-lifted collection; those lifted records must not be mixed with native hg38 records in the primary mismatch count.
