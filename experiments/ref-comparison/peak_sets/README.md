# H3K27ac peak sets for the cross-reference analysis

MACS3 peak calls from nf-core/chipseq 2.1.0 (`../scripts/run_nfcore.sh`,
`../nextflow.config`, `../config/`) for three ENCODE H3K27ac experiments
(heart ENCSR175ABH, liver ENCSR864OOO, lung ENCSR954JMZ), each processed
against GRCh37, GRCh38, and CHM13 v2.0. `broad/` holds the
`*_peaks.broadPeak` files and `narrow/` the `*_peaks.narrowPeak` files,
gzipped and otherwise unmodified.

These are the exact inputs to `bedtools getfasta` (rule `peaks_to_fasta` in
`../workflow/Snakefile`); each file's line count equals the record count of
the corresponding FASTA used in the paper. Regenerate the FASTAs from these
peaks and the reference sequences listed in the paper's Data Availability
statement. Checksums of the gzipped files are in `SHA256SUMS`.
