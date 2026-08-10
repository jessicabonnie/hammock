# Submittability concerns (2026-08-09)

Written after a benchmarking bug hunt turned up four harness defects and then,
when the scope widened to "review the whole approach," several concerns about
the paper's scientific claims. Recorded here for review, not acted on.
Figure 3 itself is now correct and separate from this list — see
`docs/bedtools-baseline-retraction.md` for that.

## 1. No control isolates what the speedup is evidence for

hammock's cost is Θ(N) sketching + Θ(N²) cheap comparisons; the BEDTools
workflow is Θ(N²) full file re-reads (`bedtools jaccard` has no batch mode). A
tool that just read each BED once and did N² exact sorted-merge intersections
in one process — no sketching, no approximation — would recover most of
Figure 3's advantage. Until that control is run, the speedup is evidence for
*not re-parsing files N² times*, not for HyperLogLog specifically. Any
fixed-size sketch would show the same curve.

**Cost to build:** small — the merge algorithm is well understood, and the
gate (its output must equal `bedtools jaccard` exactly before any timing is
trusted) is cheap to write.

**Cost to run: not cheap, and should be estimated honestly rather than
guessed.** Tonight's Panel A rerun (N up to 512, 3 reps, both existing tools)
took 1h04m end to end (21:38–22:43, per `sacct` — an earlier draft of this
note said "23:43," an unverified guess, not a checked timestamp; corrected).
Most of that time is the tool arms themselves, not overhead: the N=512
bedtools cell alone (3 reps) is 35.7 minutes, N=256 is 9.3, N=512 hammock is
7.4 — 55 of the 64.5-minute total from the three largest cells. A batch-exact
control swept across the same N grid, as a genuinely new third arm, is
realistically comparable — tens of minutes to an hour or more, not a quick
check. Budget it like the Panel A/B jobs already run today, not like a unit
test.

## 2. Speed and accuracy are never measured at the same point

Figure 3's speed number is p=18 on synthetic data, which has no accuracy
measurement against BEDTools anywhere (though BEDTools runs in that harness
already — the exact answer is a free byproduct of the benchmark and is
currently discarded). Figure 4's accuracy number is p=21 on Maurano, where
hammock is measured *slower* than BEDTools without subsampling. There is no
single (N, p) cell in the repo where a reader can read off both numbers at
once.

**Cost to run:** does *not* require repeating the expensive N-scaling sweep —
accuracy only needs one representative N (say N=20–32), run once, with
BEDTools' jaccard output actually parsed and kept instead of discarded. Small
next to tonight's multi-hour Panel A job, but "small" here means minutes to
tens of minutes for corpus generation plus one BEDTools/hammock pass, not
zero.

## 3. The timed configuration isn't the recommended one

Figure 3's headline runs `--no-metrics`, i.e. the register-equality
`jaccard_similarity` column. The paper's own guidance is to read
`jaccard_similarity_ie` for anything that matters — and that column costs
8.5% more wall time at N=512, rising with N. The speed number shown is for a
column the paper tells readers not to trust.

**Cost to fix:** wording only — the `+IE` curve is already plotted in Figure
3A; state plainly which curve is "the recommended configuration."

## 4. `jaccard_similarity` is shipped as the default column, and isn't set Jaccard

It's column 1, unrenamed, no runtime warning unless `--verbose` is passed. It
carries a chance-agreement floor (~0.17, not constant — varies with size
ratio) and isn't rank-faithful across differently-sized pairs. Documented
internally as frozen for backward compatibility with the original tool's
values, not as a scientific choice. Figure 6's headline clustering result
(ARI=0.91) uses this column, against the paper's own stated rule to use the
calibrated one for anything that matters.

**Cost to fix:** the frozen constraint is on the *values*, not the *label* —
renaming the CLI flag/column is cheap. Recomputing Figure 6 on `_ie` and
checking the result holds is a rerun, not a redesign (a scoped check already
shows the two columns agree at that one cell).

## 5. HyperLogLog vs. MinHash is never argued

One descriptive sentence in the whole repo states the choice; none justifies
it against the standard alternative for Jaccard estimation. There's a real
argument available and unmade: HLL's inclusion-exclusion route estimates J as
a difference of two independently-noisy cardinalities, so relative error
grows as ~1/J — worse than MinHash's ~1/√J for low-similarity pairs. The
repo's own accuracy data (flat absolute MAE across J) shows exactly this
signature.

**Cost to fix:** writing, not experiments — the data to make the argument
already exists in the repo.

## 6. Mixed-stride's periodicity risk is asserted, not tested

The one genuinely novel algorithmic contribution here is a **deterministic**
per-chromosome stride sampler — same stride and phase for every file, every
run. At the recommended aggressive setting (`subB=0.01`), stride = 100 bp,
which exactly divides the ~200 bp nucleosome repeat and the ~300 bp Alu
element. The outline's dismissal of periodicity concerns addresses a
different question (cross-chromosome alignment) than the one it's cited
against (within-chromosome resonance), and no test in the repo constructs
periodic input to check.

**Cost to run:** one focused experiment — synthetic BEDs with controlled
periodic spacing, sweep period vs. stride through exact-divisor cases.

## 7. Accuracy claims rest on one biological corpus

Every accuracy/clustering number traces to the same 20-file Maurano set (one
species, one tissue type, one build). The synthetic corpus can't substitute —
its pairs all have nearly the same true J by construction, so it's unusable
for an accuracy curve. Figure 6's k=10, w=30 optimum is the argmax of a
39-cell sweep on that same 20-file corpus, an isolated spike over an
otherwise-flat landscape, with no held-out check or permutation test.

**Cost to fix:** a second real corpus is the expensive option; at minimum,
disclose the argmax-selection issue explicitly rather than reporting the
sweep optimum's own score as if it were an independent result.

## 8. An internal contradiction in checked-in guidance

CLAUDE.md states Mode D's hash-width bias as an established "−0.5% to −8.3%."
The document it cites for that number, `docs/seed-mode-d-hash-width.md`,
says explicitly it "has never been reproduced and its sign disagrees with the
mechanism." This is a direct contradiction inside the file every session
loads, independent of anything else on this list.

**Cost to fix:** minutes — resolve which document is right and correct the
other.

## What doesn't need rework

The calibration numbers that can be independently checked hold up: Figure 4's
MAE (4.3×10⁻⁴), Pearson r (0.999987), Kendall τ (0.9947) at p=21 all
reproduce from the CSVs, on a real dynamic range (true J 0.135–0.627), not a
narrow band. The internal documentation is unusually willing to correct
itself — it caught and retracted its own errors repeatedly over the life of
this project, including several during this review. The problem is
consistently that the rigor lags what's exposed on the CLI and in the
manuscript draft, not that the rigor is absent.

## Where this leaves the paper

None of the above says the tool doesn't work. It says the single speed claim
needs a scope statement (file count, and now also: which column, and
possibly interval-width regime — see the note below), the headline
similarity column needs renaming or replacing in the headline figure, and the
one novel algorithm needs a validation experiment it hasn't had. Items 3, 4,
and 8 are cheap. Items 1, 2, 6, and 7 are real experiments, in roughly that
priority order.

**One more finding, flagged separately because it wasn't independently
re-verified as carefully as the rest of this list**: a reviewer measured that
hammock's Mode B cost scales with total interval *width* (it hashes every
base position) while BEDTools' scales with interval *count*, and found the
speed ordering reverses for wide intervals (~270 kb mean width: BEDTools 4.8×
faster). Real Maurano data (mean width 460 bp) doesn't hit this regime, and
the synthetic corpus (100–10,000 bp) doesn't either, so it doesn't affect
Figure 3 as published — but broad ChIP peaks, CNVs, and TADs routinely have
mean widths in the range where this reviewer's numbers say hammock loses, and
nothing discloses that. Worth confirming with a direct measurement before
either citing or dismissing it.
