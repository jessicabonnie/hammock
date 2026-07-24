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
- Generated tables should not be edited by hand.

## Planned structure

```text
paper/database_counts/
├── README.md
├── queries/
├── scripts/
├── manifests/
└── results/
```

## Initial scope

The first implementation will focus on ENCODE ChIP-seq and ATAC-seq BED files because ENCODE exposes file-level metadata for assay, output type, file format, status, and assembly. Sequence/reference counts will be defined separately after choosing the exact counting unit, since a FASTA file is not itself normally described as being aligned to a reference.

## Open methodological decisions

- Whether the sequence table counts genome assemblies, reference FASTA files, or sequencing datasets aligned to each assembly.
- Whether the primary interval tables report all released BED files, one representative final peak set per experiment, or both.
- How to treat legacy assemblies, alternate loci, unplaced contigs, and assembly aliases.
- Whether later versions combine ENCODE with GEO, ArrayExpress/BioStudies, or other repositories, and how cross-repository duplicates will be identified.
