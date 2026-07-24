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

From the repository root, run:

```bash
python paper/database_counts/scripts/plot_chip_atlas_growth.py
```

Python dependencies are `pandas` and `matplotlib`.

## Initial repository-counting scope

The first original counting implementation will focus on ENCODE ChIP-seq and ATAC-seq BED files because ENCODE exposes file-level metadata for assay, output type, file format, status, and assembly. Sequence/reference counts will be defined separately after choosing the exact counting unit, since a FASTA file is not itself normally described as being aligned to a reference.

## Open methodological decisions

- Whether the sequence table counts genome assemblies, reference FASTA files, or sequencing datasets aligned to each assembly.
- Whether the primary interval tables report all released BED files, one representative final peak set per experiment, or both.
- How to treat legacy assemblies, alternate loci, unplaced contigs, and assembly aliases.
- Whether later versions combine ENCODE with GEO, ArrayExpress/BioStudies, or other repositories, and how cross-repository duplicates will be identified.
