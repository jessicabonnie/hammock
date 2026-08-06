# man-monkey — Tissue-over-Species Clustering at Primate Distance

> **NEVER RUN. This directory is a scouting record, not an experiment.**
> No samples were committed, no data was acquired, and no code was written.
> `workflow/`, `config/` and `scripts/` exist but are **empty**; there is no
> Snakefile, no samples manifest, and no `results/` symlink. Nothing here has
> produced a number. Steps 1–5 of "Next steps" in
> `docs/experiment_design.md` are all outstanding.

## What it asks

The tissue-over-species clustering result failed at ~80 Mya human↔mouse
divergence (see `experiments/mus-homo/`, verdict: null). Does it work at
primate distance — ~6 Mya to chimp, ~25 Mya to macaque — where orthologous
peak sequences are 93–99 % identical?

## Contents

| file | what it is |
|---|---|
| `docs/experiment_design.md` | design options A/B/C over three candidate published datasets, with a recommendation (Option B) and open questions |
| `docs/primate_data_scouting_2026-05-12.md` | the raw scouting pass over GEO/ArrayExpress that produced the candidate table |

Accession numbers in the design doc are marked approximate and were **not**
verified. ENCODE has no non-human-primate H3K27ac / H3K4me3, so this cannot
reuse the `mus-homo` download pipeline — every peak set would have to be
pulled from a paper's supplementary material and, in most cases, lifted
between assemblies.

## Caveat that would apply if it were run

Cross-species Mode D similarity measures **shared k-mer content** — driven by
repeats and low-complexity sequence — and not homology (`CLAUDE.md`,
"Cross-reference caveat"). hammock's default `k=8` is unsuitable for
cross-species work; the sibling experiments sweep k up to 20 and still find
the small-k cells saturated. Sequence-level sketching across species also has
to contend with the finding that both `mus-homo` (~80 Mya) and
`primate-phylogeny` (deep mammalian) reached, respectively, a null and a
partial result — the design doc here was written before either closed, so its
framing of the question as still open should be read against those outcomes.

## Overlap with sibling experiments

Option C in the design doc (Villar/Berthelot primates, liver only) is largely
subsumed by `experiments/primate-phylogeny/`, which ran a 7-species Villar
panel and recovered the primate clade `((hsa, rheMac), calJac)` cleanly while
failing on deeper topology. If this experiment is picked up, Option A or B is
the only non-duplicative direction.
