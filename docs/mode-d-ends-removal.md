# Why the Mode D `_with_ends` columns were removed

**Status: done, not proposed.** The seven `*_with_ends` columns were deleted in
v0.6.0 and Mode D has emitted a single similarity block since. This file is the
record of why, including the counter-evidence and the honest limits of both.
Short version in `CLAUDE.md` divergence #8.

- **Question.** `jaccard_similarity_with_ends` merged the minimizer HLL with a
  per-record start/end-k-mer HLL. Was it measuring anything the minimizer column
  did not?
- **Measured.** Four structural properties of the construction (an uncontrolled
  size-weighted blend; 77–79% of the added elements chimeric; a 1 bp boundary
  cliff; a fallback that made the two columns bit-identical exactly where the
  column was supposed to help), plus a win-rate scoring of all 235 archived
  Maurano Mode D configs against `bedtools jaccard`, plus a v0.5.0-vs-v0.6.0
  wheel-to-wheel timing.
- **Answer.** Remove it. `no_ends` wins 185/235 outright and 93.2% under an
  outcome-independent filter, and removing the column bought 1.5–2.5×
  single-threaded on peak-shaped FASTA at byte-identical output on every
  surviving column.
- **Open.** Whether a *properly built* flank statistic would be useful on
  corpora with genuinely shared exact termini — see "What this does not settle".
  No corpus in this repo can test it.

## What the column was

Mode D built two HLLs per FASTA. `minimizer_hll` took every (k, w) window
minimizer. `startend_hll` took, per record,
`payload = canonicalize_kmer(s[:k] + s[-k:])` — the first k and last k bases
concatenated into one 2k string, lex-min'd against its reverse complement —
and then hashed **every** sliding k-mer of that payload, k+1 of them.
`jaccard_similarity_with_ends` compared `minimizer_hll ∪ startend_hll` on each
side; six containment/cosketch/`_ie` twins came along with it.

The construction was inherited from orig `hammock/lib/minimizer.py`, so orig
emits `jaccard_similarity_with_ends` too. Dropping it is a deliberate divergence.

## Four structural problems

### 1. It was a blend at a weight the user could not set

The two hash sets are disjoint, so the merged Jaccard is exactly the
size-weighted mediant of its parts:

    J_merged = (I_M + I_E) / (U_M + U_E) = (1−φ)·J_minimizer + φ·J_flank,
    φ = U_E / (U_M + U_E)

Measured against the shipped register-equality column (p=24, 3000 records of
400 bp, half shared):

| k | w | J_minimizer | J_flank | φ | predicted | actual |
|---|---|---|---|---|---|---|
| 10 | 20 | 0.6722 | 0.3442 | 0.418 | 0.5350 | 0.5352 |
| 10 | 100 | 0.6444 | 0.3442 | 0.782 | 0.4095 | 0.4097 |
| 10 | 200 | 0.6209 | 0.3442 | 0.886 | 0.3758 | 0.3759 |
| 20 | 100 | 0.4153 | 0.3338 | 0.794 | 0.3506 | 0.3508 |

**Caveat, measured:** the identity is tight only at high precision. Residual by
p at k=10, w=100: +0.0034 (p=8), −0.0136 (p=10), −0.0144 (p=16), −0.0012 (p=20),
−0.0002 (p=24). Below p≈18 the chance-agreement floor bends it. The
inclusion-exclusion flavor is cleaner throughout.

φ is a function of (k, w, record count, total length) — never a user parameter.
On one Maurano DHS file (198,621 records, 91.4 Mbp) the flank sketch is 12% of
the merged sketch at w=10, 55% at w=100, and 86% at w=500.

Note this identity is *algebraic*, not empirical — for disjoint sets it cannot
fail. It establishes that the weight is uncontrolled. It does **not** by itself
establish that `J_flank` is uninformative; that is the next section.

### 2. Most of the added elements are chimeric

Of the k+1 hashes inserted per record, only indices 0 and k are real sequence
k-mers. The other k−1 span the artificial `start||end` splice and occur in no
genome. On real Exp A H3K27ac peak FASTAs at k=10 they are **77–79%** of
distinct end elements.

Because canonicalization is applied to the whole 2k concatenation rather than
per k-mer, a record's start-k-mer orientation depends on its *other* end — so a
mutation at one end can change how the untouched end is hashed.

### 3. The flank term measures exact boundary identity, as a cliff

Same underlying genome, record boundaries jittered by δ bp (k=10, w=50, p=24):

| δ (bp) | 0 | 1 | 3 | 10 | 50 | 200 |
|---|---|---|---|---|---|---|
| J_minimizer | 1.000 | 0.997 | 0.995 | 0.986 | 0.932 | 0.783 |
| J_flank | 1.000 | **0.130** | 0.047 | 0.024 | 0.017 | 0.016 |

One base destroys it, and further jitter barely moves it — a step, not a
gradient. Under pure content divergence with boundaries fixed, J_flank instead
just tracks J_minimizer at lower resolution.

On real data the bit is essentially always off: two Maurano fetal-brain DHS BEDs
share **15 exactly-matching intervals out of 168,322 (0.01%)**; across all 190
Maurano pairs, 0.0029%–0.018%.

### 4. The stated rationale was false

The column was justified as rescuing records too short to produce minimizers.
But the no-minimizer fallback fed the whole record into **both** HLLs and
returned early, so the flank block was unreachable there. On a constant-length
sub-threshold corpus the two columns were **bit-identical** — verified at
(k=8, w=100, L=80), (k=10, w=200, L=150), (k=8, w=40, L=30) — and both collapse
to 0.0 under one SNP per record, because the fallback element is an exact-match
indicator. Dropout begins at `L < k + w − 1`.

## The empirical case

Scoring all 235 existing Maurano Mode D configs by Spearman ρ against the
`bedtools jaccard` reference (190 unordered pairs each):

| filter | family | `no_ends` win rate |
|---|---|---|
| all cells | register-equality | 185/235 (78.7%) |
| all cells | inclusion-exclusion | 165/235 (70.2%) |
| dropout < 10%, unsaturated | inclusion-exclusion | 109/117 (**93.2%**) |
| dropout < 5%, unsaturated | inclusion-exclusion | 83/87 (**95.4%**) |

The filter is deliberately defined on quantities independent of the outcome —
record-dropout fraction from the input length distribution and `k+w−1`, and
`mean(J_no_ends) < 0.9` — so it is not selection on the dependent variable. All
8 exceptions under it are cells where *both* columns have ρ < 0.39.

Global optimum: **k=20, w=20, p=24 → ρ = 0.9997** on `jaccard_similarity_ie`,
at 1.9% dropout.

### Testing the "large window" hypothesis directly

The original hope was that flanks would help where minimizers are sparse. Built
explicitly — records just above the dropout threshold, so minimizers per record
approaches 1 with *zero* dropout (k=10, p=22, A/B sharing 50% of records):

| w | L | minimizers/record | J_minimizer | J_flank | J_with_ends |
|---|---|---|---|---|---|
| 50 | 61 | 1.00 | 0.368 | 0.340 | 0.342 |
| 100 | 114 | 0.98 | 0.415 | 0.343 | 0.348 |
| 200 | 219 | 0.89 | 0.458 | 0.341 | 0.349 |
| 400 | 429 | 0.76 | 0.556 | 0.340 | 0.352 |

The flank component does not rescue the estimate in the sparse-minimizer
regime — it *takes over*. At large w the column stops tracking sequence content
and converges on record identity.

## Counter-evidence, stated honestly

**`with_ends` does win 12 cells on the inclusion-exclusion family**, at ρ up to
0.9498 vs 0.8950 for `no_ends` (k ∈ {15,20,25}, w=500, p ∈ {20,22,23,24}). This
was found by an adversarial reviewer and independently re-derived. Two things
qualify it: every one of those cells has **72.5–73.1% record dropout**, i.e. the
regime where the fallback makes the column a record-identity detector; and they
sit well below the 0.9997 the minimizer column reaches at 1.9% dropout.

**The register-equality and inclusion-exclusion families disagree** on the
winner in roughly 28 of 234 cells — a reminder that `jaccard_similarity` carries
a cardinality-dependent chance-agreement floor and merging changes the
cardinality (CLAUDE.md divergence #2).

**The bedtools yardstick is partly circular.** Mode D FASTAs are
`bedtools getfasta` extractions from a single reference, so two records have
identical sequence roughly iff they have identical coordinates:
ρ(`jaccard_similarity_ie`, bedtools) = 0.9997. Any column adding a term
bedtools does not measure can only lose ρ. This is why the removal rests on the
structural findings above rather than on the win rate.

**The ref-comparison Exp A result that originally justified the column does not
survive scrutiny — in either direction.** The "Metric choice at k=10, w=10"
block in `docs/paper_outline.md` §4.4 (Robustness to reference genome) ranked
metrics by Δmedian, which rewards de-saturation; on scale-free AUC,
`jaccard_similarity` gets 1.0000 (broad and narrow) vs 0.9012/0.8889 for
`_with_ends`. But that AUC is one biological sample deep — every
leave-one-sample-out refit restores 1.0000 for `_with_ends` — and a
**content-free** size statistic, `−|log(C_BA/C_AB)|`, scores AUC 0.9918/1.0000
on the same task. Exp A cannot support a comparative claim about either column.

## Performance effect

Measured by building v0.5.0 (`e84be5c^`) and v0.6.0 as separate wheels in
isolated venvs and timing both. All 17 surviving CSV columns are byte-identical
between the two builds over a 64-row run, so the speedup is free.

| cell | in-process sketch + parse (28 Mbp) | real CLI, 8 files, `--threads 1` | real CLI, 8 files, `--threads 8` |
|---|---|---|---|
| k=8, w=40 | 1.66× | 1.54× (p=20) | 2.07× (p=20) |
| k=20, w=100 | 2.70× | 2.45× (p=24) | 4.54× (p=24) |

Peak RSS at N=20, p=24: **3.45 GB → 2.14 GB**.

Post-sketch saving (the removed `merged()` construction plus the second
`pairwise_metrics_hll` pass) scales as 4^p: 0.02 s at p=18, 0.14 s at p=20,
0.78 s at p=22, 2.86 s at p=24 — under 1% of a sweep cell at every precision.
Sketching is ~99% of the work.

**Mechanism.** The removed per-record cost is `canonicalize_kmer` over 2k chars
plus k+1 Python-level `xxh64` calls: 10.3 µs/record at k=8, 22.6 µs/record at
k=20, with 84% of it in the hashes. It is flat in `w` and **roughly linear in
k**. The ratio therefore grows on both axes — numerator with k, denominator
(minimizer work) shrinking with w — so the large k=20, w=100 figure is mostly a
*k* effect, not a *w* effect.

**Caveat: the speedup is a function of mean record length, and does not
transfer to genome-scale FASTA.** Holding total bases at 28 Mbp and re-chunking
the same sequence:

| mean record length | k=8, w=40 | k=20, w=100 |
|---|---|---|
| 4,660 bp | 1.10× | 1.29× |
| 466 bp (peak-shaped) | 1.81× | 3.29× |
| 130 bp | 3.59× | 8.89× |

So this is a peak-FASTA result. On contig, transcript, or whole-genome input the
removal is worth ~1.1–1.3×. The table above also uses fBrain-DS14718 (mean 460 bp),
the second-shortest of the 20 Maurano files; corpus-wide (mean 574 bp) the
single-threaded figures are ~5–10% lower.

**Read the `--threads 8` column with care.** It is larger only because Mode D
threading is a GIL convoy that the removal partly relieves — not because the
removal parallelizes. Mode D is *faster at `--threads 1` than at `--threads 8`*
in both versions; see the Architecture note in CLAUDE.md. Fixing that is worth
more than this removal was.

## What this does not settle

Whether a *properly built* flank statistic would be useful on corpora whose
records genuinely share exact termini — restriction fragments, amplicon panels,
annotation versions, assembly contig ends. That would need its own column
(per-k-mer canonicalization, no junction k-mers, exposed separately so the user
sets the weight) and a corpus with a δ=1 bp jitter control to prove it is not
just a coordinate-equality detector. No such corpus exists in this repo.

## Consequences

Archived CSVs keep the columns; nothing on disk is invalidated. Re-pointed to
`jaccard_similarity`: `experiments/ref-comparison/workflow/Snakefile`,
`experiments/ref-comparison/scripts/exp_a_dendrogram.R`, all three
`config.yaml` `primary_sim_col` keys
(two of which were dead — ref-comparison hardcodes it, primate reads it
nowhere), and `experiments/primate-phylogeny/` (which dropped its three `_we`
trees; the `_mz` equivalents already exist on `/vast`).

**Not regenerated.** Figures and quoted numbers computed on the old column are
stale, including Fig 7 in `docs/paper_outline.md` (rendered by
`experiments/ref-comparison/scripts/exp_a_dendrogram.R`) and the §4.4
metric-choice claim, which `docs/paper_outline.md` has since marked WITHDRAWN
(2026-07-31) rather than regenerated. Also stale:
`experiments/mus-homo/` results. `experiments/mus-homo/scripts/compute_column_comparison.R`
compares the two columns and is now moot — it was never run (no
`column_comparison.tsv` on disk).
