# BED-file counting methodology

## Goal

Build a reproducible companion to published repository-comparison tables that quantifies coordinate-dependent interval resources directly.

The analysis will report three related but non-equivalent quantities:

1. **Physical BED files**: distinct downloadable files with BED-compatible interval content.
2. **Biological datasets represented**: experiments, samples, or processed datasets represented by those files.
3. **Genomic intervals**: total records contained in the selected BED files, when practical to compute.

These quantities must remain separate. A repository may provide several BED files for one experiment, one aggregate BED file representing thousands of experiments, or multiple significance thresholds for the same peak calls.

## Primary output table

The manuscript-facing table will contain:

| Repository | Assay | Species | Assembly | Physical BED files | Represented datasets | Intervals | Download granularity | Snapshot date |
|---|---|---|---|---:|---:|---:|---|---|

Counts will be accompanied by a row-level manifest that records every included file or catalog entry.

## Counting rules

### Physical file

Count one physical file when all of the following hold:

- it is independently downloadable;
- it contains genomic intervals in BED, BED-like, narrowPeak, broadPeak, BEDPE, bigBed, or a documented BED-compatible format;
- its coordinate assembly is known or recoverable from repository metadata;
- it belongs to the selected repository snapshot.

Compressed and uncompressed versions of the same file are one file. Mirrors of the same file are one file when identity can be established.

### Dataset represented

Use the repository's native biological unit:

- ChIP-Atlas: SRX experiment accession;
- Cistrome DB: processed sample or database sample identifier;
- ReMap: source dataset identifier, not the aggregate catalog file;
- GTRD: experiment or peak set identifier as documented by the release;
- MethBank: sample or experiment identifier as documented by the release.

Aggregate BEDs must list, or be linked to metadata listing, the datasets represented. If no defensible mapping exists, report the physical file count but leave represented datasets unavailable.

### Interval count

For plain or compressed BED files, count non-comment, non-track, non-browser records. For bigBed files, use a documented conversion or summary utility. Aggregate interval counts should not be added to per-experiment interval counts when both represent the same underlying peaks.

## File classes

Each manifest record receives one `file_granularity` value:

- `experiment`: one experiment or processed sample;
- `replicate`: one biological or technical replicate;
- `threshold`: one experiment at a specific peak-call threshold;
- `factor_aggregate`: multiple datasets grouped by regulator or antigen;
- `biotype_aggregate`: multiple datasets grouped by cell or tissue class;
- `catalog_aggregate`: an entire repository or release catalog;
- `other_aggregate`;
- `unknown`.

The main manuscript count should emphasize experiment-level or sample-level BED files. Aggregate downloads are reported separately because they measure packaging rather than dataset abundance.

## Peak types

Record the output class without collapsing distinct biological objects:

- narrow peaks;
- broad peaks;
- gapped peaks;
- accessibility peaks;
- methylated regions;
- chromatin-interaction pairs;
- merged or consensus peaks;
- other genomic intervals.

## Assembly handling

Preserve the repository label and add a normalized label. Examples:

- `GRCh38`, `hg38` -> `GRCh38/hg38`;
- `GRCh37`, `hg19` -> `GRCh37/hg19`;
- `GRCm38`, `mm10` -> `GRCm38/mm10`;
- `GRCm39`, `mm39` -> `GRCm39/mm39`.

Lifted files are counted separately from native files and marked with `coordinate_origin = lifted`. A file must not be counted as an independent biological dataset merely because an additional lifted copy exists.

## Initial repositories

### ChIP-Atlas

ChIP-Atlas provides individual experiment peak calls and assembled BED files. The analysis will first count individual experiment-level peak outputs and then report assembled files separately. Multiple significance thresholds for the same experiment will be distinguishable rather than interpreted as independent experiments.

### Cistrome Data Browser

Cistrome exposes downloadable BED peak files at the processed-sample level and supports bulk peak-file downloads. Counts should use sample identifiers and preserve assay class and species/build metadata.

### ReMap

ReMap distributes BED files at several aggregation levels, including regulator, biotype, and whole-catalog files. Physical file counts alone would greatly understate the represented datasets, so source dataset counts and catalog-level files must be reported separately.

### GTRD and MethBank

These will be added after confirming stable release manifests and exact BED or BED-compatible download structure.

## Required provenance

Every repository snapshot must record:

- repository name and release/version;
- retrieval date;
- manifest or API URL;
- query parameters or commands;
- source-paper citation where applicable;
- checksums for downloaded metadata files when feasible;
- inclusion and exclusion counts with reasons.

## Planned analysis sequence

1. Inventory stable metadata and bulk-download endpoints.
2. Write one repository-specific manifest extractor per source.
3. Normalize manifest columns without altering source labels.
4. Validate assembly, file granularity, and dataset identifiers.
5. Generate counts by repository, assay, species, assembly, and granularity.
6. Optionally stream files to count intervals without retaining all BED data locally.
7. Generate manuscript tables from the normalized manifest.
