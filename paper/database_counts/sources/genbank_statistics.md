# GenBank and WGS growth series

## Source

NCBI GenBank and WGS Statistics:

- https://www.ncbi.nlm.nih.gov/genbank/statistics/

Related publication:

- Sayers EW, Cavanaugh M, Clark K, et al. GenBank 2024 Update. *Nucleic Acids Research*. 2024;52(D1):D134–D137. doi:10.1093/nar/gkad903.

## Counting units

The official NCBI table reports, for each GenBank release:

- bases and sequence records in traditional non-set-based GenBank divisions;
- bases and sequence records in the Whole Genome Shotgun (WGS) collection.

The source data retain these components separately. `total_bases` and `total_sequences` are derived as the row-wise sums of the GenBank and WGS columns.

## Annual sampling rule

The GenBank 2024 Update states that its growth figure uses the August release from each year beginning with release 173 in August 2009. This analysis follows that convention through release 268 in August 2025.

## Reproducibility notes

- Values were transcribed directly from the official NCBI release table; no figure digitization was used.
- The plotting script checks that component sums equal the stored totals.
- The plotting script does not require a network connection.
- The NCBI table notes that traditional GenBank counts exclude CON-division records and set-based WGS/TSA/TLS records; WGS is reported separately and added here only in explicitly labeled derived totals.
