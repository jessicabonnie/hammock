# Published figures and tables on repository growth

Research date: 2026-07-24

This document records published figures and tables that may help motivate the scale and growth of genomic repositories relevant to interval sketching. The emphasis is on sources that either provide explicit time series or can support a reproducible redraw.

## Highest-priority candidates

### ChIP-Atlas 2025 update

**Citation**

Oki S. et al. *ChIP-Atlas 2025 update: 10-year anniversary of a data-mining platform for exploring epigenomic landscape.* Nucleic Acids Research. 2026;54(W1):W77-W84. DOI: 10.1093/nar/gkag378.

**Useful item**

- **Table 1: Ten-year progress of ChIP-Atlas.**
- Reports annual cumulative experiment counts from 2015 through 2025:
  - 2015: 37,720
  - 2016: 52,249
  - 2017: 72,199
  - 2018: 90,322
  - 2019: 131,903
  - 2020: 158,863
  - 2021: 196,136
  - 2022: 361,856
  - 2023: 408,697
  - 2024: 459,231
  - 2025: 464,655
- The same table records when additional assay types and genome assemblies were introduced, including ATAC-seq and newer assemblies.

**Why it is especially useful**

This is the closest published match to our motivation. It simultaneously demonstrates growth in public epigenomic experiments and expansion across assay types and coordinate systems.

**Potential use**

Replot the published annual cumulative counts with a citation and clearly state that values were transcribed from Table 1. Add annotations for assay and assembly expansion.

**Source**

- https://doi.org/10.1093/nar/gkag378
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13355075/

### ChIP-Atlas 2021 update

**Citation**

Oki S. et al. *ChIP-Atlas 2021 update: a data-mining suite for exploring epigenomic landscapes by fully integrating ChIP-seq, ATAC-seq and Bisulfite-seq data.* Nucleic Acids Research. 2022;50(W1):W175-W182. DOI: 10.1093/nar/gkac199.

**Useful items**

- **Figure 1B:** cumulative number of SRX-based experiments recorded in ChIP-Atlas.
- **Table 1:** counts by experiment type and database version.
- The paper states that more than 100,000 ChIP-seq experiments were added after 2018 and describes the addition of all ATAC-seq datasets archived in SRA.

**Potential use**

Useful as an earlier independent source for the growth trajectory and for assay-specific counts. The 2026 update is preferable for the longest time series.

**Source**

- https://doi.org/10.1093/nar/gkac199
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9252733/

### GEO 23-year update

**Citation**

Clough E. and Barrett T. *NCBI GEO: archive for gene expression and epigenomics data sets: 23-year update.* Nucleic Acids Research. 2024;52(D1):D138-D144. DOI: 10.1093/nar/gkad965.

**Useful items**

- **Figure 1:** GEO growth and data-type trends from 2013 through 2022. Stacked bars show the fraction of broad data-type categories, and a line shows total studies.
- **Figure 3:** cumulative supplementary-data storage and number of supplementary files from 2013 through 2022.
- The text reports more than 200,000 studies and 6.5 million samples, with studies increasing at approximately 15% per year and doubling about every five years.

**Why useful**

This establishes that public functional-genomics archives are growing rapidly and that next-generation sequencing now dominates submissions. It is broader than BED/interval data and should therefore be used as general context rather than as a direct count of comparable interval files.

**Reproducibility note**

The paper points to the GEO history summary, which may permit a current, independently reproducible extension of the time series.

**Source**

- https://doi.org/10.1093/nar/gkad965
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10767811/
- https://www.ncbi.nlm.nih.gov/geo/summary/?type=history

### Sequence Read Archive decade update

**Citation**

Katz K. et al. *The Sequence Read Archive: a decade more of explosive growth.* Nucleic Acids Research. 2022;50(D1):D387-D390. DOI: 10.1093/nar/gkab1053.

**Useful item**

- **Figure 1:** growth of SRA over the preceding decade.
- Through September 2021, the figure reports approximately 25.6 petabase pairs from more than 14.8 million public runs.

**Why useful**

This is a strong general-scale figure for raw sequencing repositories, but it does not directly quantify reference-aligned BED files. It may be useful in an introductory panel paired with a more targeted ChIP-Atlas or ENCODE panel.

**Source**

- https://doi.org/10.1093/nar/gkab1053
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8728234/

## Secondary candidates

### Cistrome Data Browser v3.0

**Citation**

Zheng R. et al. *Cistrome Data Browser: integrated search, analysis and visualization of chromatin data.* Nucleic Acids Research. 2024;52(D1):D61-D68. DOI: 10.1093/nar/gkad949.

**Relevant result**

- The update added approximately 30,000 samples.
- The paper notes especially strong growth in chromatin-accessibility samples due to adoption of ATAC-seq.

**Limitations**

The paper is more useful for before/after database-version comparisons than for a continuous annual time series.

**Source**

- https://pmc.ncbi.nlm.nih.gov/articles/PMC10767960/

### Cistrome Data Browser 2019 update

**Citation**

Mei S. et al. *Cistrome Data Browser: expanded datasets and new tools for gene regulatory analysis.* Nucleic Acids Research. 2019;47(D1):D729-D735. DOI: 10.1093/nar/gky1094.

**Relevant item**

- Figure 1 compares the earlier and expanded collections by assay/sample categories.

**Potential use**

Together with the 2017 and 2024 Cistrome papers, this could provide a sparse version-to-version growth series, but ChIP-Atlas provides a cleaner published annual series.

**Source**

- https://pmc.ncbi.nlm.nih.gov/articles/PMC6324081/

### ATACdb 2.0

**Citation**

Fang Q-L. et al. *ATACdb 2.0: a comprehensive chromatin accessibility database of human and mouse.* Nucleic Acids Research. 2026.

**Relevant result**

- The paper reports more than a 3.5-fold increase in samples and more than a 7.5-fold increase in accessible chromatin regions relative to ATACdb 1.0.

**Potential use**

This directly illustrates ATAC-seq resource growth, although it provides a two-version comparison rather than an annual trajectory.

**Source**

- https://pmc.ncbi.nlm.nih.gov/articles/PMC12807738/

### ChIP-Atlas 3.0

**Citation**

Oki S. et al. *ChIP-Atlas 3.0: a data-mining suite to explore chromosome architecture together with large-scale regulome data.* Nucleic Acids Research. 2024;52(W1):W45-W53. DOI: 10.1093/nar/gkae358.

**Useful items**

- **Table 2** reports experiment and interval counts by experiment type and database version.
- The unified pipeline identified more than 11 billion genomic intervals, including approximately 2.09 billion ChIP-seq binding intervals and 1.30 billion ATAC-seq/DNase-seq accessibility intervals.
- The paper reports about a 30% increase in total intervals relative to ChIP-Atlas 2.0.
- **Table 3** compares ChIP-Atlas with Cistrome DB, ReMap, GTRD, and MethBank across data sources, supported assays, preprocessing, experiment counts, organisms, genome assemblies, analysis tools, and access requirements.

**Why Table 3 is especially useful**

Table 3 demonstrates that the relevant interval ecosystem is distributed across multiple large services rather than contained in one archive. More importantly for Hammock, its genome-assembly row exposes the number and diversity of coordinate systems supported by these resources. It can therefore motivate both repository scale and cross-reference fragmentation.

The comparison includes, among others:

- ChIP-Atlas: hg19, hg38, mm9, mm10, rn6, dm3, dm6, ce10, ce11, and sacCer3.
- Cistrome DB: hg38 and mm10.
- ReMap: hg38, mm10, dm6, and TAIR10, with some legacy human and mouse data available through lift-over.
- GTRD: hg38, mm10, TAIR10, ce11, danRer11, dm6, rn6, sacCer3, and spo2.
- MethBank: multiple animal and plant assemblies.

**Potential use**

Rather than reproduce the full service-feature matrix in the main manuscript, derive a compact table with one row per service and columns for:

1. service;
2. assay classes;
3. reported number of experiments;
4. number of organisms;
5. number of explicitly listed assemblies;
6. assemblies or assembly families;
7. whether lift-over is used.

The full transcription should remain in source data, while the manuscript table should emphasize scale and coordinate-system diversity. Counts in the original table use different inclusion rules across services and should not be interpreted as directly harmonized measurements.

**Why useful**

This is highly relevant to Hammock because it counts genomic intervals themselves, not merely studies or files, and Table 3 connects those data to a fragmented multi-repository, multi-assembly ecosystem.

**Source**

- https://doi.org/10.1093/nar/gkae358
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11223792/

## Preliminary recommendation for the paper

A compact motivation figure could combine:

1. **Repository growth:** redraw the annual cumulative ChIP-Atlas experiment counts from Table 1 of the 2026 update.
2. **Interval scale:** report the billions of ChIP-seq and accessibility intervals from ChIP-Atlas 3.0.
3. **Reference fragmentation:** use Table 3 of ChIP-Atlas 3.0 as the published starting point, then add our own reproducible current counts grouped by genome assembly.

This would create a clear progression:

> the number of experiments is growing, those experiments produce billions of intervals, and the intervals are distributed across multiple repositories and incompatible reference coordinate systems.

## Caveats

- Published repository totals count different units: studies, experiments, samples, runs, files, bytes, bases, or intervals. They should not be plotted on a shared numeric axis without clear labeling.
- ChIP-Atlas counts SRX-based experiments and may include multiple assay classes. Assay-specific values should be taken from the appropriate table or recomputed from a released metadata snapshot.
- Counts reported for different services may use different inclusion, filtering, and deduplication rules and are not necessarily directly comparable.
- Published figures may be redrawn from reported numerical data with citation; copying the original image may require publisher permission depending on license.
- Any transcription of values from a published table should be stored as a small source-data file with the citation, table identifier, and transcription date.
