# ChIP-Atlas annual experiment counts

## Source

Zou Z, et al. **ChIP-Atlas 2025 update: 10-year anniversary of a data-mining platform for exploring epigenomic landscape.** *Nucleic Acids Research*. 2026;54(W1):W77–W84. doi:10.1093/nar/gkag378.

- Published online: 2026-04-29
- Source location: Table 1, “Ten-year progress of ChIP-Atlas”
- Values transcribed: row “Cumulative number of experiments,” 2015–2025
- Retrieved: 2026-07-24
- Primary article: https://doi.org/10.1093/nar/gkag378
- Open full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC13355075/

## Transcription policy

The values in `../source_data/chip_atlas_annual_experiments.tsv` are direct transcriptions from Table 1. Thousands separators were removed so that the count column is machine-readable. No interpolation or smoothing was performed.

The event columns preserve milestones from the same table:

- initial ChIP-seq and DNase-seq support in 2015;
- additional assemblies in 2017;
- hg38, mm10, dm6, and ce11 in 2020;
- ATAC-seq and Bisulfite-seq in 2021;
- annotation tracks in 2023.

## Interpretation limits

These are cumulative **experiment** counts reported by ChIP-Atlas, not counts of BED files, unique biological samples, publications, or genomic intervals. ChIP-Atlas uses SRX-based experiments, and one experiment may yield multiple processed files or interval sets. The series is therefore appropriate for illustrating repository growth but should not replace the file- and interval-level counts planned elsewhere in this directory.

## Verification

The plotting script checks that years are unique and increasing, counts are non-negative and cumulative, and the expected 2015 and 2025 endpoint values are unchanged. If the source data are revised, those assertions must be updated deliberately.
