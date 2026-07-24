# NCBI GEO repository history

## Source

NCBI Gene Expression Omnibus (GEO), repository history table:

- https://www.ncbi.nlm.nih.gov/geo/summary/summary.cgi?type=history

Related publication:

- Clough E, Barrett T. NCBI GEO: archive for gene expression and epigenomics data sets: 23-year update. *Nucleic Acids Research*. 2024;52(D1):D138-D144. doi:10.1093/nar/gkad965.

## Retrieval and transcription

- Retrieved: 2026-07-24.
- The source page reports cumulative counts by calendar quarter for GEO Series, Platforms, and Samples.
- `source_data/geo_annual_q4_counts.tsv` contains the fourth-quarter snapshot for each complete year from 2001 through 2025.
- The incomplete 2026 calendar year is intentionally excluded.
- Values are transcribed directly from the official NCBI history table. No interpolation or digitization from a figure was used.

## Counting units

- **Series**: GEO Series records (GSE accessions), which group related Samples into a study-level record.
- **Platforms**: GEO Platform records (GPL accessions).
- **Samples**: GEO Sample records (GSM accessions).

These counts describe repository records, not unique biological experiments, donors, processed interval files, or reference assemblies. A GEO Series may contain multiple assays or supplementary files and may link to raw sequencing data in SRA.

## Why this source is useful

The official history table provides a long, first-party, reproducible time series for a major functional-genomics repository. The associated GEO update paper uses the same history resource to describe repository growth and reports that studies were increasing at roughly 15% per year at the time of publication.

## Validation checks

The plotting script verifies that:

1. years are unique and strictly increasing;
2. all three cumulative count columns are nondecreasing;
3. the first and last complete-year values match the transcribed source endpoints;
4. no partial 2026 observation is included.
