# Database counts for interval-sketching motivation

This directory contains reproducible analyses used to quantify the scale and reference-build fragmentation of public genomic resources discussed in the Hammock paper.

## Target tables

1. Sequence or assembly resources grouped by reference assembly or coordinate system.
2. ChIP-seq BED files grouped by reference assembly.
3. ATAC-seq BED files grouped by reference assembly.

## Reproducibility principles

- Every reported count must be generated from a versioned script or taken directly from a citable source.
- Repository-derived counts must record the data source, exact query or API request, retrieval date, inclusion and exclusion criteria, and counting unit.
- Raw metadata manifests should be retained so that every aggregate count can be traced to the records included in it.
- Results should distinguish file counts from experiment-level or dataset-level counts.
- Reference-build labels should be preserved as reported by the source and normalized separately for aggregation.
- Generated tables and figures should not be edited by hand.

## Structure

```text
paper/database_counts/
├── README.md
├── published_growth_figures.md
├── queries/
├── scripts/
├── source_data/
├── sources/
├── manifests/
└── results/
```

## Published ChIP-Atlas growth series

The first completed analysis reproduces the cumulative number of ChIP-Atlas experiments reported for 2015–2025 in Table 1 of the ChIP-Atlas 2025 update.

Files:

- `source_data/chip_atlas_annual_experiments.tsv`: direct transcription of the published annual counts and selected repository milestones.
- `sources/chip_atlas_2025.md`: citation, retrieval date, transcription policy, counting-unit limitations, and verification notes.
- `scripts/plot_chip_atlas_growth.py`: validates the source data and regenerates the figure.
- `results/chip_atlas_growth.svg`: version-controlled vector figure.
- `results/chip_atlas_growth.png`: generated locally by the plotting script; the PNG may be omitted from version control when the SVG is sufficient.

Run:

```bash
python paper/database_counts/scripts/plot_chip_atlas_growth.py
```

## NCBI GEO growth series

The second analysis uses the official GEO repository history table rather than digitizing a published plot. It takes the fourth-quarter snapshot for every complete calendar year from 2001 through 2025 and graphs cumulative GEO Series and Sample records.

Files:

- `source_data/geo_annual_q4_counts.tsv`: official year-end counts for GEO Series, Platforms, and Samples.
- `sources/geo_history.md`: first-party source, related publication, retrieval date, counting units, and transcription policy.
- `scripts/plot_geo_growth.py`: validates the cumulative series and generates the figure.
- `results/geo_growth.svg` and `results/geo_growth.png`: generated outputs.

Run:

```bash
python paper/database_counts/scripts/plot_geo_growth.py
```

## Shared-axis ChIP-Atlas and GEO comparison

The harmonized comparison places ChIP-Atlas cumulative experiments and GEO cumulative Series on the same linear y-axis for the overlapping years 2015–2025. GEO Samples are excluded because they use a different counting unit and are much larger in scale.

Files:

- `source_data/chip_atlas_geo_comparison_2015_2025.tsv`: joined annual series.
- `sources/chip_atlas_geo_comparison.md`: comparability and construction notes.
- `scripts/plot_chip_atlas_geo_comparison.py`: validates and generates the shared-axis figure.
- `results/chip_atlas_geo_comparison.svg` and `results/chip_atlas_geo_comparison.png`: generated outputs.

Run:

```bash
python paper/database_counts/scripts/plot_chip_atlas_geo_comparison.py
```

## GenBank and WGS growth series

The third analysis follows the convention used in the GenBank 2024 Update: one August release per year beginning with release 173 in 2009. It retains traditional GenBank and WGS counts separately and derives explicit combined totals for sequence records and base pairs.

Files:

- `source_data/genbank_august_growth_2009_2025.tsv`: official August release counts and derived totals.
- `sources/genbank_statistics.md`: source, publication, counting-unit definitions, and construction notes.
- `scripts/plot_genbank_growth.py`: validates component sums and generates log-scale growth plots.
- `results/genbank_growth.svg` and `results/genbank_growth.png`: generated outputs.

Run:

```bash
python paper/database_counts/scripts/plot_genbank_growth.py
```

Python dependencies for these analyses are `pandas` and `matplotlib`.

## Initial repository-counting scope

The first original counting implementation will focus on ENCODE ChIP-seq and ATAC-seq BED files because ENCODE exposes file-level metadata for assay, output type, file format, status, and assembly. Sequence/reference counts will be defined separately after choosing the exact counting unit, since a FASTA file is not itself normally described as being aligned to a reference.

## Open methodological decisions

- Whether the sequence table counts genome assemblies, reference FASTA files, or sequencing datasets aligned to each assembly.
- Whether the primary interval tables report all released BED files, one representative final peak set per experiment, or both.
- How to treat legacy assemblies, alternate loci, unplaced contigs, and assembly aliases.
- Whether later versions combine ENCODE with GEO, ArrayExpress/BioStudies, or other repositories, and how cross-repository duplicates will be identified.
