# ChIP-Atlas and GEO comparison series

This file documents the harmonized comparison series used to plot ChIP-Atlas and GEO on a shared axis.

## Purpose

The comparison figure aligns annual ChIP-Atlas cumulative experiment counts with annual GEO cumulative Series counts for the overlapping years 2015–2025.

## Why these units

- **ChIP-Atlas** reports cumulative **experiments** in Table 1 of the 2025 update paper.
- **GEO** reports cumulative **Series** in the official GEO repository history table.
- These are not identical counting units, but they are more comparable than pairing ChIP-Atlas experiments with GEO Samples.
- GEO Samples remain useful as a separate measure, but should not share this linear axis because both their scale and unit differ.

## Sources

- ChIP-Atlas 2025 update, Table 1.
- NCBI GEO repository history table.

## Construction rule

`source_data/chip_atlas_geo_comparison_2015_2025.tsv` is an inner join by calendar year over the overlapping interval 2015–2025. No interpolation or digitization was used.
