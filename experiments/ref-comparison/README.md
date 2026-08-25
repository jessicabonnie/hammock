# ref-comparison — Cross-Reference Robustness

Validates minimizer-based sequence sketching (hammock Mode D) on ENCODE H3K27ac
ChIP-seq by asking whether **the same biological sample aligned to different
human reference genomes** (GRCh37, GRCh38, CHM13) sketches more similarly across
references than across tissues on any one reference.

**Outcome: positive.** The hypothesis holds across the whole usable
parameter range. On the corrected `rehash-selector64` sequence-HLL path and
the current inclusion--exclusion estimator, the lead cell is **k=20, w=20**
(seed-42 broad Δmedian 0.559). It remains the lead cell for broad and narrow
peaks at all three prespecified seeds (1, 42, and 99). Numbers, figures and
caveats: `docs/exp_a_results.md`.

> Two things a reader must carry. (1) The headline stats were computed on
> `jaccard_similarity_with_ends`, a column **hammock v0.6.0 removed**;
> `docs/exp_a_results.md` carries a 2026-08-06 recomputation on the surviving
> `jaccard_similarity` showing the result holds and strengthens (Δ 0.483,
> still fully separated). (2) A Mode D ingest bug was fixed on 2026-05-14 and
> everything was re-run; the pre-fix reading in which k=15/20 looked like a
> "saturated low" regime is wrong — post-fix they are the *strongest*
> discriminators. Any pre-2026-05-14 note claiming otherwise is stale.

> **Corrected HLL rehash rerun, 2026-08-25.** The 63-cell Figure 5 union was
> rerun at p=24 with `--sequence-hll-hash rehash-selector64` for seeds 1, 42,
> and 99. This separates minimizer selection from HLL register assignment.
> The tissue-first topology and k20_w20 optimum are unchanged; the seed-42
> broad Δ changed from 0.5632 to 0.5590. Reproducible drivers and the frozen
> cell manifest are in `rehash_rerun/`.

This experiment was originally part of `claude-ref-comparison` together with a
tissue-over-species experiment ("Exp B"). The tissue-over-species work has been
split into separate experiments (`mus-homo`, `primate-phylogeny`, `man-monkey`)
because it doesn't actually need the multi-reference alignment pipeline that
justifies this directory's existence.

---

## Repository structure

```
ref-comparison/
├── docs/
│   ├── experiment_design.md   ← design rationale, accessions, success criteria
│   ├── exp_a_results.md       ← results summary + figure index (start here)
│   ├── paper_outline.md       ← paper outline (Exp A only)
│   ├── data_selection_slides_prompt.md  ← slide-deck prompt for the sample-choice story
│   ├── workflow_diagram_prompt.md       ← diagram prompt for the pipeline figure
│   └── references.bib         ← BibTeX bibliography
├── workflow/
│   ├── Snakefile              ← downstream pipeline (peaks → FASTA → hammock → plots)
│   └── slurm_profile/
│       └── config.yaml        ← SLURM cluster profile for Snakemake
├── config/
│   ├── config.yaml            ← samples, refs, sketch parameters
│   ├── exp_a_metadata.tsv     ← sample × ref metadata for plotting
│   ├── confirmed_accessions.tsv
│   └── nfcore_samplesheet_human.csv
├── scripts/
│   ├── fetch_accessions.py       ← resolve GEO/SRA accession numbers
│   ├── run_nfcore.sh             ← launches nf-core/chipseq per reference
│   ├── exp_a_validate_plot.R     ← Wilcoxon + 2-panel boxplot + heatmap
│   ├── exp_a_sweep_summary.R     ← (k × w) effect-size heatmap across the sweep
│   ├── exp_a_dendrogram.R        ← UPGMA dendrogram (broad + narrow)
│   └── exp_a_metric_comparison.R ← per-metric Wilcoxon comparison at one (k, w)
├── setup_experiment.sh        ← one-time scaffolding (dirs + symlinks into /vast)
├── environment.yaml           ← conda env spec for `claude-ref-comparison`
├── nextflow.config            ← cluster profile consumed by scripts/run_nfcore.sh
├── results/                   ← fine-grained symlinks into legacy storage (Exp A subset)
│   ├── exp_a/ → .../claude-ref-comparison/results/exp_a
│   ├── fastas/ → .../claude-ref-comparison/results/fastas
│   └── logs/{exp_a,fastas}/   ← only the Exp A-relevant log subtrees
├── figures/                   ← local mirror of headline + representative cell figures (gitignored)
└── README.md
```

Under `results/exp_a/{broad,narrow}/` there are 20 `k{k}_w{w}/` cells (the
full `w ≥ k` sweep), each holding the all-vs-all CSV, `cross_ref_stats.tsv`,
and `cross_ref_validation.png`, plus a per-peak-type `sweep_effect_size.png`.

Hammock + nf-core outputs live at `/vast/blangme2/jbonnie/hammock/claude-ref-comparison/`
(kept on the legacy path so existing CSVs / plots remain valid without
re-running). The `results/` directory exposes only the Exp A-relevant
subtree; the Exp B artifacts under that legacy storage are not surfaced
here. See `docs/exp_a_results.md` for the headline numbers.

---

## Two-stage pipeline

1. **nf-core/chipseq** (upstream) — download FASTQs, align to each reference, peak-call. Multi-reference re-alignment is the part of the workflow that **needs** a full pipeline; the other experiments in `hammock_claude/experiments/` can use pre-called peak BEDs from ENCODE directly.
2. **Snakemake** (downstream, this directory) — peaks-to-FASTA → hammock all-vs-all → (k, w) sweep → plots.

---

## Quickstart

```bash
# Step 1: Upstream — run nf-core/chipseq across references
bash scripts/run_nfcore.sh                       # all references (GRCh38, GRCh37, CHM13)
bash scripts/run_nfcore.sh GRCh38                # one reference only

# Step 2: Downstream — Snakemake (peaks → FASTA → hammock → plots)
conda activate claude-ref-comparison
snakemake --dryrun                               # check job plan
snakemake --profile workflow/slurm_profile/      # submit to SLURM
```

---

## Sample design

| Tissue | ENCODE accession | References |
|--------|------------------|------------|
| Heart | ENCSR175ABH | GRCh37, GRCh38, CHM13 |
| Liver | ENCSR864OOO | GRCh37, GRCh38, CHM13 |
| Lung  | ENCSR954JMZ | GRCh37, GRCh38, CHM13 |

3 samples × 3 references = 9 sample × ref combinations. Two peak types
(MACS3 broad + narrow). The archived sweep used `k ∈ {5, 8, 10, 15, 20}`
× `w ∈ {5, 8, 10, 15, 20, 30}` filtered to `w ≥ k` (20 cells). The current
focused extension uses `k ∈ {20, 25, 30, 35, 40, 50}` ×
`w ∈ {20, 25, 30, 35, 40, 50, 60, 75}` with the same filter (33 cells),
retaining k=20 as an overlap check because the prior maximum occurred at the
boundary of the archived grid. Targeted follow-up runs added `w ∈ {25, 50,
75}` at `k ∈ {5, 8, 10, 15}` to make the low-k and k=15 trajectories directly
comparable. The merged table contains 63 unique cells per peak type and is
staged at `docs/data/exp_a_estimator_delta_expanded.csv`.

---

## Compute requirements

| Step | CPUs | RAM | Est. walltime |
|---|---|---|---|
| nf-core download + alignment (per sample × ref) | 16 | 32 GB | ~4 h |
| hammock per (k, w) | 1 | 16 GB | re-baseline on the extended grid |
| R plotting (per cell, sweep summary) | 1 | 8 GB | <2 min |

> **Mode D timing note, 2026-08-06.** The hammock row was measured at
> `--threads 8`. Mode D threading is a GIL convoy, not parallelism — v0.6.1
> changed the Mode D default to `--threads 1` after measuring the 8-thread
> pool at **2× slower** than single-threaded (`CLAUDE.md`, Architecture).
> So the walltime above is *inflated*, not a floor; re-baseline before
> sizing a new run off it.
