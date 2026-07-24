# Motivating statements for the interval-sketching paper

This file collects candidate statements developed during the database-count analysis. They are working language, not final manuscript prose. Each statement should eventually be paired with the specific table, figure, or analysis that supports it.

## Repository growth and scale

> Public epigenomic repositories are growing rapidly, and the experiments they contain produce billions of genomic intervals.

> The growing volume of public interval data changes the computational problem from comparing a few hand-selected files to searching, clustering, and organizing collections containing hundreds or thousands of related BED files.

> Repository-scale interval analysis is therefore both a storage problem and a comparison problem: the number of datasets is increasing, while exhaustive pairwise comparison grows quadratically with the number of files.

## Reference fragmentation

> Public functional-genomics datasets are distributed across multiple reference assemblies. BED files aligned to different assemblies do not share a common coordinate system and therefore cannot be compared directly with conventional coordinate-overlap operations.

> Even among experiments targeting the same protein, the resulting BED files are partitioned across species-specific and version-specific coordinate systems.

> Large collections of genomic intervals are not only distributed among repositories; they are also partitioned among numerous incompatible reference coordinate systems.

> As genome assembly diversity grows through updated references, strain- and population-specific assemblies, and pangenomes, reference-dependent comparison becomes an increasingly important limitation.

## Target-matched example

> A researcher interested in one well-defined biological target, such as CTCF, may reasonably want to compare every available peak set across cell types, studies, laboratories, species, and reference assemblies.

> Target-matched BED files have a shared biological purpose but may differ in cell context, experimental design, peak-calling behavior, and coordinate system. This makes them a realistic collection for similarity search, clustering, quality control, and dataset retrieval.

> A human-only CTCF analysis is a useful controlled example, but restricting the primary count to hg19 and hg38 would hide the broader reference fragmentation that motivates reference-independent interval comparison.

> The primary inventory should therefore retain every reported assembly, with species- and assembly-restricted analyses treated as derived views rather than as the default dataset.

## BED files, experiments, and intervals

> File counts, experiment counts, and interval counts describe different aspects of repository scale and should be reported separately.

> The number of BED files reflects the operational number of coordinate-dependent objects a researcher may need to manage and compare; the number of experiments reflects biological dataset count; and the number of intervals reflects computational workload within those files.

> A raw BED-file count can be misleading when one repository distributes one file per experiment while another distributes target-level or catalog-level files that aggregate thousands of experiments. Reproducible comparisons must record the granularity represented by each file.

> A repository record that lists several supported assemblies does not by itself prove that a downloadable BED file exists for every experiment–assembly combination. Assembly labels define candidate coordinate-specific outputs; physical files must be verified before they are counted.

> Experiment counts and verified BED-file counts may therefore differ in either direction: one experiment may have multiple verified reference-specific BED files, while another may have no downloadable peak BED at the selected threshold.

## Possible narrative sequence

1. Public repositories contain a rapidly growing number of functional-genomics experiments.
2. These experiments produce billions of genomic intervals and large collections of BED files.
3. A biologically coherent question can require comparison of hundreds or thousands of target-matched BED files.
4. Those files are fragmented across incompatible reference assemblies.
5. Existing exact coordinate-overlap workflows therefore require remapping, lift-over, or separate within-reference analyses.
6. Compact reference-independent sketches offer a way to search and compare these collections without requiring every dataset to share one coordinate system.
