# SRA experiment growth series

## Counting unit

This analysis counts **SRA Experiment records**, not SRA Runs.

NCBI describes an SRA Experiment as a unique sequencing result for a specific sample and the main publishable unit in the SRA database. An Experiment is defined by the relevant combination of sample, library, sequencing strategy, layout, instrument model, and related experimental metadata. One Experiment may be associated with one or more Runs.

This makes SRA Experiments substantially more comparable to ChIP-Atlas experiments than SRA Runs would be. They are still not identical to GEO Series, which are study-level groupings.

## Primary sources

- NCBI SRA metadata overview: https://www.ncbi.nlm.nih.gov/sra/docs/submitmeta/
- NCBI SRA Entrez search documentation: https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/
- NCBI E-utilities documentation: https://www.ncbi.nlm.nih.gov/books/NBK25501/

## Query method

The script `scripts/fetch_sra_experiment_growth.py` queries the NCBI Entrez `sra` database using ESearch. NCBI states that SRA Entrez query results are Experiment records.

For each year, the script requests the number of publicly indexed SRA Experiment records with publication dates from 1900-01-01 through December 31 of that year:

```text
1900/01/01:YYYY/12/31[PDAT]
```

The returned count is therefore a cumulative year-end snapshot, not the number of newly released experiments during that year.

## Reproducibility and limitations

- The output records the exact query term, generated ESearch URL, and retrieval date for every row.
- Counts are live repository queries and may change if NCBI revises release dates, suppresses records, restores records, or changes indexing behavior.
- Re-running the script may therefore produce slightly different historical counts. The generated TSV should be committed when values are used in a manuscript.
- The published SRA growth paper primarily reports Runs and base pairs. Those quantities are intentionally not substituted for Experiment counts here.
- SRA Experiments and ChIP-Atlas experiments are related but not independent universes: ChIP-Atlas processes selected SRA experiments, so the databases overlap and their lines must not be summed.
