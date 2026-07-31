# Multi-tissue primate data scouting — 2026-05-12

Findings from a background research agent. Treat accessions as approximate
pointers; the agent could not live-verify against GEO/ArrayExpress (web access
denied in its sandbox).

## Headline

**The single best public study that cleanly meets the ≥2 tissues per species
bar appears to be the NHP-PsychENCODE multi-brain-region ChIP-seq dataset**,
*if* we accept brain regions (DLPFC, V1, cerebellum, hippocampus, striatum)
as "tissues" — analogous to how Yue 2014 and Lin 2014 treat cortex and
cerebellum separately. For true cross-organ multi-tissue with multiple
primates, **no dataset confidently meets all bars in publicly-released form**
as of 2026-05; the "primate multi-organ epigenome" paper is still pending.

## Ranked candidates

### 1. NHP-PsychENCODE (Sousa, Sestan et al.)
- Multi-brain-region ChIP-seq (H3K27ac, H3K4me3) in rhesus + chimpanzee +
  matched human in same study.
- Regions: DLPFC, V1, cerebellum, hippocampus, striatum (5+ tissues if we
  count brain regions).
- Approx accession: GSE127898 and related. Processed peak BEDs available
  via Synapse (registration required, but public).
- **Verdict:** strongest candidate. Verify accession before committing.

### 2. NIH NHP Reference Epigenome Project / PRIMATE-FREEZE (Carbone, Ahituv, Pollen labs)
- Multi-tissue ATAC-seq + H3K27ac CUT&RUN in rhesus macaque (brain, liver,
  heart, kidney, testis).
- GSE2398xx / GSE243xxx range from Carbone lab (OHSU/ONPRC), 2023-2024.
- Status: partial release; some embargoed.

### 3. Garcia-Perez / Marques-Bonet great-ape multi-tissue ATAC-seq
- Garcia-Perez et al., Cell Genomics 2024 (approx): chimp, gorilla, orang
  + human iPSC-derived AND biobanked tissue panel. Brain + testis + liver
  in chimp/gorilla.
- Tissue coverage uneven across species — chimp has the most.

### 4. Marmoset Brain Atlas (RIKEN / Okano lab, Brain/MINDS)
- Multi-region ChIP-seq + ATAC-seq in *Callithrix jacchus* brain.
- Single species, brain-only — qualifies only if we accept brain regions
  as tissues.
- Data on DDBJ (DRA accessions) and Brain/MINDS portal.

## Dead ends (<2 tissues per primate)

- ENCODE rhesus pilot: single-tissue per assay mostly
- 4D Nucleome NHP: too sparse, mostly cell lines
- "Ape vs human" papers from Gallego-Romero, Kanton, Pavlovic: LCLs / iPSCs, not tissues
- Zoonomia regulatory companion: NHP tissue coverage is mostly liver (overlap with Villar/Berthelot)

## Pan-tissue NHP atlases underway (not yet released)

- NIH SenNet NHP arm (rhesus, marmoset) — multi-tissue ATAC + snATAC; preprints expected 2026
- "PriMAP" — primate multi-tissue molecular atlas (marmoset + macaque); RFA-RM-23-xxx
- Allen Institute marmoset whole-brain epigenome (multi-region snATAC); partial release on portal.allen.brain

## Honest summary from the agent

> The Yue/Lin-style analysis at primate divergence is genuinely under-served
> by current public data, which is why so many groups have just used liver
> as the lowest-common-denominator. This is itself a useful finding.

## Recommended next moves (for when we return to man-monkey)

1. **Verify accessions** in live GEO: `(Macaca mulatta[Organism] OR Pan troglodytes[Organism]) AND (ATAC-seq OR ChIP-seq) AND tissue`, sort by sample count.
2. Check Synapse for `syn4921369` (PsychENCODE) sub-projects with NHP tag.
3. Watch bioRxiv for "PriMAP" or "primate multi-tissue regulatory atlas".
4. **Tentative plan:** lead with NHP-PsychENCODE brain-region tissue-vs-species (Option A from the design doc, scoped to brain regions), and explicitly call it brain-region clustering rather than whole-organ tissue clustering.

## Caveats on this report

- Agent had no web access; findings are from training-data recall as of January 2026
- Accession numbers are approximations — every one needs live verification
- The 2026-05-2026 landscape may have moved further forward than the agent could capture
