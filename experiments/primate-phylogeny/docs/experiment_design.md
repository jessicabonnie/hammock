# primate-phylogeny: Sketch-Distance Recovery of Mammalian Phylogeny

**Status:** Run and closed as **partial** (accepted 2026-05-14, H3K27ac
refreshed 2026-05-16) on 7 of the planned 20 species. Primate clade recovered;
deep topology fails. See "Verdict: partial" below. The header below is the
original planning framing, kept as written.

**Question:** does pairwise minimizer-sketch similarity over H3K4me3 (or
H3K27ac) peak FASTAs in a single tissue recover the known mammalian phylogeny?

**Created:** 2026-05-12. **Status line corrected 2026-08-06** — it still said
"Planning. No code yet." long after the pipeline ran.

> **Schema note added 2026-08-06.** Everything below was measured before
> hammock v0.6.0, which removed the `jaccard_similarity_with_ends` column
> family (`CLAUDE.md` divergence #8). Passages that quote `_with_ends`
> values, and the `_we` trees / distance matrices in the output paths, refer
> to a column that no longer exists and cannot be regenerated. The archived
> files are unaffected; only re-runs are. Parse them by column *name* —
> positional field indices differ between the pre- and post-v0.6.0 schemas.
> Separately, `jaccard_similarity` is register-equality and is **not** set
> Jaccard: it has a chance-agreement floor and is not rank-faithful across
> pairs of differing FASTA size (divergence #2), which is directly relevant
> here because the per-species peak FASTAs differ in size by 2× and that
> asymmetry is already implicated below in the `cosketch_max` failure.
> `estimator_ie_topology.py` in the parent directory re-scores the topology
> on the calibrated `jaccard_similarity_ie` instead.

---

## Motivation

The cross-species tissue-clustering experiments (`mus-homo`, `man-monkey`)
ask whether the sketch can pick up tissue identity *across* species at a
given divergence. A complementary and arguably more tractable question is
whether the sketch picks up *species* identity in a way that tracks
phylogeny — i.e., does sketch similarity decay with evolutionary distance
in a recoverable, ordered way?

This question is well-served by **Villar et al. 2015** (*Cell* 160:554),
which generated H3K4me1 / H3K4me3 / H3K27ac ChIP-seq on **adult liver** from
20 mammalian species, deposited at ArrayExpress `E-MTAB-2633`. Liver is a
single, well-defined tissue across all 20 species, so the species axis is
isolated from tissue confounds.

## Hypothesis

> Pairwise hammock similarity over H3K4me3 peak FASTAs in adult liver
> across the Villar mammals will reconstruct a UPGMA / neighbor-joining
> tree whose Robinson-Foulds distance to the published mammalian phylogeny
> is significantly lower than randomized expectation.

A weaker but still meaningful version: pairwise sketch similarity will
correlate with phylogenetic divergence time (Spearman ρ).

## Why H3K4me3 (not H3K27ac)

H3K4me3 peaks are anchored to orthologous TSSs and conserved much more
strongly across mammals than H3K27ac enhancers. For a *phylogeny-recovery*
experiment, conservation is a feature — it means there's actually shared
signal to recover. Villar 2015 reports approximately:
- H3K27ac: ~10–20 % peak overlap across mammals
- H3K4me3: ~60–70 % peak overlap across mammals

Both marks are available in the Villar dataset, so H3K27ac can serve as a
**comparison** showing how the sketch handles a faster-evolving mark.

## Dataset

**Villar et al. 2015, ArrayExpress E-MTAB-2633.**

| species | common name | reference assembly |
|---|---|---|
| *Homo sapiens* | human | GRCh38 |
| *Pan troglodytes* | chimpanzee | panTro6 |
| *Gorilla gorilla* | gorilla | gorGor6 |
| *Pongo pygmaeus* | orangutan | ponAbe3 |
| *Macaca mulatta* | rhesus macaque | rheMac10 |
| *Callithrix jacchus* | marmoset | calJac4 |
| *Mus musculus* | mouse | mm10 |
| *Rattus norvegicus* | rat | rn7 |
| *Oryctolagus cuniculus* | rabbit | oryCun2 |
| *Cavia porcellus* | guinea pig | cavPor3 |
| *Canis familiaris* | dog | canFam5 |
| *Felis catus* | cat | felCat9 |
| *Mustela putorius furo* | ferret | musFur1 |
| *Sus scrofa* | pig | susScr11 |
| *Equus caballus* | horse | equCab3 |
| *Ovis aries* | sheep | oviAri4 |
| *Bos taurus* | cow | bosTau9 |
| *Monodelphis domestica* | opossum | monDom5 |
| *Macropus eugenii* | tammar wallaby | macEug2 |
| *Ornithorhynchus anatinus* | platypus | ornAna1 |

Peaks are deposited as per-species BED files in the ArrayExpress submission.
Some species' peaks are called against an older assembly than what's listed
here; we'll need to either liftOver or download the matching assembly.

> **What was actually run (note added 2026-08-06).** The assemblies in the
> table above are the *planning* choices and were not used. Phase 1 ran 7
> species against Ensembl release-73 (Sept 2013) assemblies, matching
> Villar's peak-call coordinates so no liftOver was needed:
> hg19, mm10, rheMac2, calJac3, canFam3, bosTau7, monDom5. The authoritative
> list is `config/config.yaml`'s `references:` block — the Snakefile only runs
> species whose `ucsc_assembly` key appears there. The remaining 13 species in
> `config/samples.tsv` were never staged.

> **Cross-species caveat (note added 2026-08-06).** At mammalian divergence,
> Mode D similarity between peak FASTAs measures **shared k-mer content** —
> repeat elements, low-complexity sequence, conserved promoter motifs — not
> homology (`CLAUDE.md`, "Cross-reference caveat"). This is the same
> mechanism the "Verdict: partial" section below diagnoses as the bottleneck,
> stated as a property of the method rather than of this dataset. hammock's
> default `k=8` is unsuitable for cross-species work; the sweep here runs k
> up to 20 for that reason, and the k=5 cells are saturated.

## Pipeline shape (proposed)

```
Villar peak BEDs (per species)  →  bedtools getfasta from species reference
                                ↓
hammock all-vs-all (single (k, w) cell, possibly small sweep)
                                ↓
  • distance matrix  →  neighbor-joining tree
  • compare to TimeTree reference  →  Robinson-Foulds distance
  • correlate pairwise distance to divergence time  →  Spearman ρ
```

Comparison versions:
1. H3K4me3 only (primary)
2. H3K27ac only (comparison — faster-evolving mark)
3. Combined or other mark (if useful)

## Success criteria

| outcome | what it means |
|---|---|
| Recovered tree has RF distance significantly below random expectation; primate clade + rodent clade visible | **Positive** — sketch recovers phylogeny. Quantify how much of the tree topology survives. |
| Sketch similarity correlates with divergence time (Spearman ρ > 0.5) but full tree topology not recovered | **Partial** — distance signal works, topology requires more sophisticated tree-building. Still useful. |
| No correlation with phylogeny | **Null** — peak-FASTA sketches don't carry phylogenetic signal in liver. Important negative result given Villar 2015 says the underlying data does. |

## Results — Phase 0 + Phase 1 (2026-05-13 / -14)

Pipeline runs end-to-end. Phase 0 (hsa + Mmus only) validated the
Villar-format normalization (header strip / CRLF strip / chr-prefix
adjust) and produced a biologically reasonable hsa↔Mmus H3K4me3 Jaccard
of 0.698 at k=10, w=10 — matching Villar 2015's reported ~60–70 % peak
overlap.

Phase 1 added rheMac, calJac, canFam, bosTau, and monDom from Ensembl
v73 archives (5 × 2.3–3.5 GB; chromosome names match Villar's bare-name
peak coordinates without liftOver). Seven-species all-vs-all completed
across the full (k, w) sweep.

### Both similarity columns are used (corrected 2026-05-14)

An earlier note in this file said `jaccard_similarity_with_ends` was
uniformly 0.000 on Villar peaks and had been dropped. That was a
header-lookup bug in my shell analysis script, not a hammock issue —
see `memory/project_villar_peaks_zero_with_ends.md`. Real values are
0.30–0.34 for cross-species pairs at k=10, w=10; both columns carry
information and the pipeline now emits a tree per column.

### Phase 1 topology readout

| (k, w) regime | behavior | usable? |
|---|---|---|
| k=5 (any w) | sketch saturated, topology near-random | no |
| **k=8–10, w=10–20** | **stable topology, primate clade recovered** | **yes** |
| k≥15 | cross-species Jaccard saturates at 0.000 → polytomy collapse | no |

> **The k≥15 row is a pre-fix artifact — checked 2026-08-06.** It was
> measured before the Mode D minimizer-ingest fix (`CLAUDE.md` divergence #6),
> which is exactly the bug that silently dropped most minimizers at k ≥ 12.
> In the post-fix `metric_spreads.tsv` files on disk, H3K4me3
> `jaccard_similarity` across species is **0.074–0.236 at k=15/w=15** and
> **0.025–0.146 at k=20/w=20** — not zero, and with *wider* spread (0.163 and
> 0.121) than the k=10/w=10 cell's 0.055. So "saturates at 0.000" is false
> post-fix. What was **not** re-derived is whether the usable-window verdict
> changes: the "Containment / cosketch addendum" below re-ran the full sweep
> post-fix and reports the Jaccard tree topology unchanged, but the
> usable-window table was never rewritten against those numbers. Treat the
> k=8–10 window as the pre-fix recommendation, not a re-confirmed one.

At the usable (k, w) cells, **the primate sub-tree ((hsa, rheMac),
calJac) is recovered cleanly** at every cell — hsa-rheMac sister
(~25 Mya, highest cross-species Jaccard 0.731 at k=10, w=10) and calJac
correctly placed as the New World outgroup.

**Long-branch artifact outside the primates.** Mmus pairs with monDom at
every usable cell, despite ~180 Myr divergence. Both species are
"equally unlike" the placental panel, so a narrow-dynamic-range distance
matrix groups them together by their shared dissimilarity from the rest.
canFam and bosTau also fail to form a Laurasiatheria clade and instead
sit as separate branches off the placental backbone.

**Dynamic range is the bottleneck — and HLL precision is not the cause.**
Cross-species similarities span 0.64–0.73 at p=24 (k=10, w=10) — a 0.09
window to encode ~180 Myr. A precision probe (2026-05-14) confirmed that
p=26 and p=28 produce *identical* spreads (0.0936 at k=10/w=10, 0.0360
at k=8/w=10), so the narrow band is the true Jaccard signal, not sketch
sampling noise. **The bottleneck is the distance metric / sketch content,
not its precision.**

`jaccard_similarity_with_ends` lives at 0.30–0.34 cross-species (k=10,
w=10) — a 0.04-wide window, even tighter than minimizers-only. Both
columns carry the same topology (hsa-rheMac highest, Mmus-monDom long-
branch artifact) but `_with_ends` compresses further.

### H3K27ac swap (2026-05-14): didn't widen the spread

Re-ran with `mark: H3K27ac` (Villar reports ~10–20 % cross-species peak
overlap vs ~60–70 % for H3K4me3, so I expected lower Jaccards / wider
spread). Result was the opposite: H3K27ac peak counts are 3–4× higher
per species (30–46k vs 11–17k), which pushed the sketch toward k-mer
saturation. At k=10/w=10 cross-species Jaccards moved to 0.774–0.810
(spread 0.037, **narrower** than H3K4me3's 0.093). Primate clade still
recovered, Mmus-monDom artifact persists. Mark choice doesn't unstick
this.

### Containment / cosketch addendum (2026-05-14, post-hammock-fix)

After hammock shipped containment_AB/BA + cosketch_{geom,arith,max} for
Mode D, re-ran k∈{5..20}/w∈{5..30} on H3K4me3 with the 12-column CSVs.
A separate Mode D arithmetic bug was discovered and fixed during this
work; pre-fix Jaccard values were ~10–15 % too low. After the fix, the
**Jaccard tree topology is unchanged** (primate clade ✓, Mmus-monDom
LBA ✗), absolute values shifted from 0.64–0.73 up to 0.80–0.86.

Containment did not fix the bottleneck either:

- `cosketch_max` gave the largest hsa↔rheMac vs Mmus↔monDom margin
  (+0.024 vs Jaccard's +0.013, ~1.9×) **but** broke the primate clade
  by pairing calJac with Mmus. The reason: calJac's FASTA is the
  smallest (20 MB) and Mmus's is among the largest (40 MB), so
  containment(calJac → Mmus) = 0.945 is the highest cosmax involving
  calJac — driven by FASTA size, not phylogeny.
- `cosketch_geom` (geometric mean — penalizes size asymmetry) recovers
  the same primate clade as Jaccard but with a *narrower* spread (0.032
  vs Jaccard's 0.055).

So containment-based metrics either don't help (geom) or actively
hurt (max) on this dataset because the size of the per-species peak
FASTAs varies more with peak count than with divergence.

### H3K27ac post-fix readout (2026-05-16)

Re-ran H3K27ac against the post-fix hammock with the 12-column CSVs and
6 trees per cell (jacc / cosmax / cosgeom × minimizer / with_ends).
Topology lands in the same place as H3K4me3:

- `jaccard_similarity` k=10/w=10: `(((hsa,rheMac), calJac), ((Mmus,monDom), canFam), bosTau)` — primate clade ✓, Mmus-monDom LBA ✗
- `cosketch_geom` matches `jaccard_similarity` topology
- `cosketch_max` breaks the primate clade exactly as it does on H3K4me3 (calJac size-asymmetry again)

But the **signal-to-noise on the critical pair is dramatically better
in H3K27ac**. The cross-species `jaccard_similarity` range at k=10/w=10
is only 0.021 wide (vs H3K4me3's 0.055 — H3K27ac saturates more
because there are 3–4× more peaks per species), but the hsa↔rheMac
vs Mmus↔monDom margin sits at +0.020, so the *primate-vs-LBA decision
spans 96 % of the cross-species range* (vs 23 % in H3K4me3). On
cosketch_geom the same ratio holds (96 %). So while H3K27ac is the
*smaller* sketch signal, it's a *cleaner* one — the right pair stands
above the noise with very little ambiguity.

| metric | H3K4me3 spread / margin / ratio | H3K27ac spread / margin / ratio |
|---|---:|---:|
| jaccard | 0.055 / +0.013 / 0.23 | 0.021 / +0.020 / **0.96** |
| cosketch_geom | 0.032 / +0.008 / 0.23 | 0.012 / +0.012 / **0.96** |
| containment_BA | 0.099 / +0.024 / 0.24 | 0.045 / +0.033 / 0.72 |
| cosketch_max | 0.040 / +0.024 / 0.59 | 0.029 / +0.027 / 0.92 |

This doesn't change the verdict — deep topology still fails the same
way — but for writing this up, the H3K27ac jaccard headline can be
stated cleanly: *the closest-primate pair is at the top of the
cross-species range and a long-branch noise pair is at the bottom,
with the right pair winning by margin covering 96 % of the range*.

### Verdict: partial (accepted, 2026-05-14)

Lands in the **partial** outcome category. The signal that exists —
the primate sub-tree ((hsa, rheMac), calJac) — is recovered cleanly and
*robustly* across every variation tried:

- (k, w) sweep — stable in the k=8–10, w=10–20 window
- HLL precision sweep — identical at p=24, 26, 28
- mark swap — same topology under H3K4me3 and H3K27ac
- metric swap — same topology under Jaccard, cosketch_geom, `_with_ends`
  variants

What does *not* recover, under any of those variations, is the deeper
topology: Mmus pairs with monDom (long-branch attraction), canFam and
bosTau don't form a Laurasiatheria clade, monDom isn't the marsupial
outgroup. The bottleneck is structural to peak-FASTA sketching at deep
mammalian divergence — the peak content is dominated by shared
mammalian k-mer background (conserved promoter motifs, repeat
elements, AT-rich segments) that crowds out species-specific signal.
No precision, no mark, no metric fixes that.

## Next direction (deferred): substrate engineering

The remaining unblocked hypothesis is that the *peak-FASTA substrate
itself* is the bottleneck — too much of the shared k-mer content is
mammalian-conserved background (repeats, common promoter motifs)
rather than peak-defining sequence. The follow-on experiment:

1. **Repeat-mask the peak FASTAs.** Strip RepeatMasker hits before
   sketching. The Ensembl v73 reference FASTAs we downloaded include
   soft-masked versions (`.dna_sm.toplevel.fa.gz`); switching `curl`
   target + adding a `bedtools getfasta -fold` (or hard-mask via
   `tr a-z N`) would convert this to a config switch.
2. **Restrict to non-promoter peaks.** Use a TSS-coordinate file to
   bedtools-subtract proximal-promoter peaks; only intergenic /
   intronic peaks go into the FASTA. These are the faster-evolving
   parts of the regulome, so any phylogeny signal should survive while
   the conserved-background floor drops.
3. **k-mer weighting (TF-IDF).** Down-weight k-mers shared across all
   7 species, up-weight species-specific ones. This requires a
   different sketch backend or a post-hoc reweighting step, not just
   a hammock CLI flag.

See `memory/project_primate_phylogeny_substrate_engineering.md` for the
pickup brief.

## Mash distance (excluded from the main paper; candidate here)

Mash distance turns a k-mer Jaccard `J` into an estimated per-base
mutation rate via `D = -(1/k)·ln(2J/(1+J))` (≈ 1 − ANI). It was
**deliberately left out of the main hammock paper** (`docs/paper_outline.md`):
every claim there is either rank-based (Pearson/Spearman/Wilcoxon) or
clustering-based (ARI/NMI, UPGMA topology on a single distance matrix),
and mash is a strictly monotone transform of `J`, so it moves none of
those numbers. It would only earn a place if the paper made a *calibrated
evolutionary-distance* claim — which is this experiment, not that one.

Why it's a genuine candidate **here** and not a no-op: neighbor-joining
is **not** invariant under a nonlinear monotone transform of the input
distances (its Q-criterion sums actual distance values, so the recovered
topology can change), unlike the rank/clustering metrics in the main
paper. And the `-ln` stretches small-`J` differences hardest — precisely
the deep-divergence regime where our signal is compressed (cross-species
`J` spans only 0.64–0.73 at k=10/w=10). So mash is worth a try as
*distance-metric* engineering, complementary to the *substrate*
engineering above.

Caveat / expectation: the diagnosed bottleneck is **dynamic range of the
true Jaccard signal**, not its scaling — the precision probe showed p=24/26/28
give identical spreads, i.e. the narrow band is real, not noise. A
monotone re-scaling can re-space a narrow band but cannot manufacture
separation that isn't in `J`, so mash is unlikely to fix the Mmus–monDom
long-branch artifact on its own. Most promising as a transform applied
*after* substrate engineering widens the underlying `J` spread. Revisit
if/when we decide to report a phylogenetic-distance (rather than
tissue-clustering) result.

## Output paths

- H3K4me3 (post-fix, 12-column CSVs + 6 trees per cell):
  `results/H3K4me3/k{k}_w{w}/{phylo_mnmzr_p24_jaccD_k{k}_w{w}.csv,
   metric_spreads.tsv, nj_tree_{jacc,cosmax,cosgeom}_{mz,we}.{newick,png},
   dist_matrix_{jacc,cosmax,cosgeom}_{mz,we}.tsv}`
- H3K27ac (post-fix, 12-column CSVs + 6 trees per cell — fully refreshed
  2026-05-16): `results/H3K27ac/k{k}_w{w}/...`

## Related work

- Roller et al. 2021 (*Nature Comm*, "LECIF") — similar data shape but with cCRE-based comparison; not a sketch method.
- Berthelot et al. 2018 (*Nat Comm*) — H3K27ac across mammals, extension of Villar.
- Schmidt et al. 2010 (*Science*) — TFBS conservation across vertebrates; supports the H3K4me3-anchored-to-TSS argument.
