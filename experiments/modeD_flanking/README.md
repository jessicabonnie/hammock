# modeD_flanking — When (if ever) do boundary k-mers help Mode D?

Mode D emits two Jaccard estimates per pair:

| column | what it counts |
|---|---|
| `jaccard_similarity`           | minimizers only |
| `jaccard_similarity_with_ends` | minimizers **plus** the first / last `(k-1)` k-mers of every sequence ("ends" / "flanking k-mers") |

The `_with_ends` form was originally introduced because short sequences may
contain *no* minimizer — every interior k-mer's minimizer happens to fall
outside the sequence — so the sketch is empty without boundary content. The
question this experiment asks is the natural follow-up:

**Under what conditions do the flanking k-mers actually improve the
estimate, and under what conditions do they corrupt it?**

The Maurano DHS sweep (`experiments/maurano_dhs_validation/`) gave the
first hint: at small `w`, `_with_ends` is wildly inflated (mean pred ~0.97
vs bedtools truth ~0.34); at modest `w` the no-ends column tracks bedtools
much more cleanly (k=10, w=20, p=24 → r=0.966 vs `_with_ends` r=0.93). But
that's one corpus, one ground-truth notion (bedtools interval Jaccard),
and one regime of sequence lengths. This experiment widens that picture.

---

## Hypothesis

Flanking k-mers help when interior minimizers are *sparse* relative to the
information they need to convey, and hurt when they are dense enough to
already characterise the sequence.

Concretely:

1. **Help regime:** short sequences (length ≲ `w`), few intervals, large
   `w`. Every sequence contributes only its flanks; without them, the
   sketch is empty or near-empty and the Jaccard is unstable.
2. **Hurt regime:** many short intervals concatenated into one FASTA
   (e.g. ChIP-seq / DHS peak FASTAs — thousands of ~150 bp regions).
   Boundary k-mers are over-represented relative to information content,
   and any two corpora with similar length distributions will share
   their flanking-k-mer fingerprints by coincidence — inflating Jaccard.
3. **Indifferent regime:** long, single-contig sequences (e.g. assembled
   mitochondrial genomes). Boundaries are 2(k-1) k-mers in a sequence of
   ~16 kb, so they're noise-level either way.

If these are right, the deciding parameter is the ratio of "boundary
k-mers per FASTA" to "minimizers per FASTA" — call this the **flanking
fraction** φ:

> φ ≈ ( 2 · (k−1) · n_intervals ) / ( total_length / w )

When φ is small, the two columns should agree. When φ is large,
`_with_ends` inflates and decouples from interior content.

---

## Design

Two parts, decoupled so part 1 needs no new compute:

### Part 1 — re-analyse `maurano_dhs_validation`

The Maurano sweep has already produced ≥ 200 Mode D CSVs, each with both
similarity columns. Read those, join against the bedtools reference, and:

- compute `Δr = r(no_ends) − r(with_ends)` and `Δmae = mae(no_ends) − mae(with_ends)`
  per `(k, w, p)`,
- compute φ from the input FASTAs (a single scalar per pair, since the
  Maurano corpus has a single (n_intervals, length distribution) per file
  — but φ does vary across pairs because both files contribute),
- plot `Δr` and `Δmae` against `w`, `k`, and φ; the prediction is a sign
  flip around some threshold of φ × w.

This is the "free" part: no hammock re-runs, just a new R/Python analysis
pass over existing CSVs.

### Part 2 — synthetic corpora with exact k-mer Jaccard ground truth

bedtools Jaccard is the wrong yardstick for the *flanking-k-mer* question:
boundary k-mers come from sequence content, not interval coordinates, so
the right ground truth is **exact k-mer set Jaccard** (the same notion
the `mode-d-optimality` experiment used in the original hammock repo).

Generate synthetic FASTA corpora that systematically vary:

| axis | values |
|---|---|
| `n_intervals` (per FASTA) | 10, 100, 1 000, 10 000 |
| mean interval length     | 50, 150, 500, 5 000 bp |
| length distribution      | constant, lognormal (σ=0.5), heavy-tailed (σ=1.5) |
| mutation rate between pair | 0.01, 0.05, 0.15, 0.30 |

For each (n_intervals, length, dist, mutation), produce a pair of FASTAs
with known divergence. Compute:

- **exact k-mer Jaccard** at the test k (ground truth);
- Mode D Jaccard, both columns, at the same (k, w, p);
- φ (flanking fraction) for the pair.

Then ask: at what φ does `(with_ends − no_ends)` flip sign, and how does
that boundary move with mutation rate?

### Optional Part 3 — real long-contig negative control

The old hammock `fish-mito` experiment has 25 mitochondrial FASTAs (~16 kb
single contigs, almost zero boundary k-mers relative to interior). It's
the cleanest "ends shouldn't matter" datapoint. Symlink those FASTAs in,
run a small (k, w, p) grid, and confirm the two columns agree on this
corpus.

---

## Reuse from `maurano_dhs_validation`

| asset | how used |
|---|---|
| `results/raw_d/*.csv`        | source for Part 1; both similarity columns already present |
| `data/maurano_bedtools_ref.tsv` | ground truth for Part 1 |
| `data/fastas/*.fa`           | inputs for measuring φ (length distribution) per file |
| `workflow/Snakefile` style   | template for this experiment's Snakefile |

Nothing in this experiment overwrites or modifies anything in
`maurano_dhs_validation/`; we only read from it.

---

## Outputs

- `results/part1_flanking_vs_bedtools.csv` — per-(k,w,p) Δr / Δmae and φ stats
- `results/part2_synthetic_truth.tsv`     — exact k-mer Jaccards for synthetic pairs
- `results/part2_sweep.csv`               — Mode D output on synthetic pairs
- `results/part2_summary.csv`             — Δr / Δmae binned by φ, mutation, length
- `figures/maurano_delta_r_vs_w.png`      — Part 1 headline
- `figures/synthetic_delta_vs_phi.png`    — Part 2 headline
- `figures/synthetic_phase_diagram.png`   — sign(Δerror) over (φ, mutation)

---

## Layout

```
experiments/modeD_flanking/
├── README.md                  (this file)
├── prepare_data.sh            symlink maurano_dhs_validation outputs into data/maurano_link/
├── generate_synthetic.py      [Part 2] seeded synthetic FASTA pair generator
├── exact_kmer_jaccard.py      [Part 2] exact k-mer set Jaccard reference
├── run_sweep_synthetic.py     [Part 2] Mode D sweep driver on synthetic corpora
├── analyze_part1_maurano.R    [Part 1] re-analyse existing Maurano CSVs
├── analyze_part2_synthetic.R  [Part 2] flanking gain/loss on synthetic
├── workflow/
│   ├── Snakefile              optional — wires Part 2 + analyses together
│   └── slurm_profile/
├── data/   → /vast/.../modeD_flanking/data
│   ├── maurano_link/          symlink back to maurano_dhs_validation results/data
│   └── synthetic/             [Part 2] generated FASTA pairs
├── results/ → /vast/.../modeD_flanking/results
├── figures/                   in-repo (small PNGs)
└── logs/   → /vast/.../modeD_flanking/logs
```

---

## Status

**Design only.** Scripts are stubs/skeletons; nothing has been run. Start
with Part 1 — it's purely an R analysis on existing data and will tell us
whether the hypothesis even survives contact with one real corpus before
we spend cluster hours on the synthetic grid.
