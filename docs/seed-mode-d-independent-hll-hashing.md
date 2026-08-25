# Seed: independent 64-bit HLL hashing for sequence-mode minimizers

## Objective

Determine whether feeding `digest`'s 32-bit minimizer selector hashes directly
into hammock's 64-bit HyperLogLog biases sequence-mode cardinality-derived
metrics, measure the effect across relevant precisions and datasets, and—if a
change is justified—implement and validate a principled independent 64-bit HLL
hashing path.

Do not begin by assuming that the existing sequence results are wrong or that
rehashing must improve them. The width mismatch and high-precision register
spike are confirmed; the sign and practical magnitude of any resulting
cardinality bias remain unresolved.

The primary scientific acceptance criterion is:

> Sequence-mode HLL inputs satisfy the estimator's uniform-hash assumption,
> cardinality and inclusion--exclusion estimates improve against exact
> feature-set truth, and every affected manuscript result is either reproduced
> or regenerated and explicitly reviewed.

## Repositories and baseline

- Code repository: `~/interval_sketch/hammock_claude`
- Manuscript repository: `~/interval_sketch/hammock_paper`
- Code baseline at seed creation:
  `325cbe7 Add mixed-stride v2 implementation seed`
- Manuscript baseline at seed creation:
  `786076e Correct sequence methods and document stride audit`
- Seed date: 2026-08-21

Read `docs/seed-mode-d-hash-width.md` completely before acting. It contains
the original measurements, corrected mechanism, and historical caveats. Treat
its measured facts as starting evidence, not as a decision.

The primary code worktree contains unrelated modified and untracked experiment
artifacts. They belong to the user. Do not alter, clean, stage, or delete them.

## Autonomous work only in disposable clones

Do not implement this task in either primary checkout. Create independent
local clones beneath a fresh directory in `/tmp`:

```bash
mktemp -d /tmp/hammock-mode-d-hash.XXXXXX
git clone --local ~/interval_sketch/hammock_claude WORK_ROOT/hammock
git clone --local ~/interval_sketch/hammock_paper WORK_ROOT/hammock_paper
git -C WORK_ROOT/hammock remote set-url origin \
  git@github.com:jessicabonnie/hammock.git
git -C WORK_ROOT/hammock_paper remote set-url origin \
  git@github.com:jessicabonnie/hammock_paper.git
```

Replace `WORK_ROOT` with the exact path returned by `mktemp`. Create a
dedicated local branch in each clone. A local clone initially points its
`origin` at the source filesystem checkout, so restore the canonical remote
URLs as shown.

Before editing the manuscript copy, run `git pull --ff-only` there to obtain
the current Overleaf state. Do not use `git worktree add`; it would mutate
administrative state in the primary repository.

The user authorizes autonomous editing, testing, local benchmarking, artifact
generation, manuscript drafting, and local commits inside these disposable
copies. Do not ask for confirmation for ordinary work confined to the copies.
This does not authorize:

- modifying, cleaning, staging, or committing either primary checkout;
- pushing any branch or commit;
- changing the canonical Overleaf manuscript;
- deleting or overwriting source datasets;
- submitting scheduler jobs or consuming substantial external compute without
  separate authorization;
- bypassing platform-enforced network or privilege approvals.

Leave the copies in place at handoff and report their exact paths, branches,
commits, tests, benchmarks, result comparisons, and proposed manuscript diff.

## Mandatory manuscript safeguards

The user has a standing instruction never to change the canonical manuscript
without explicit permission and always to pull it immediately before an
authorized change.

For this seed, manuscript edits are authorized only as drafts inside the
disposable manuscript copy. Do not push them. Before transferring any draft to
the primary manuscript, obtain explicit approval for the exact text and run a
fresh `git pull --ff-only` in the primary manuscript repository.

## Confirmed current behavior

The relevant paths are:

- `python/hammock/modes/sequence.py`
- `cpp/include/hammock/hll_sketch.hpp`
- `cpp/src/hll_sketch.cpp`
- `python/hammock/cli.py`
- `docs/seed-mode-d-hash-width.md`

Current sequence mode:

1. Calls `digest.window_minimizer(..., include_hash=True)`.
2. Receives pairs of minimizer position and canonical selector hash.
3. The returned selector hashes are at most 32 bits.
4. Passes each selector hash directly to
   `HLLSketch.add_hash64`, where it is widened to `uint64_t` but not
   independently rehashed.
5. Constructs the HLL with the default `hash_size=64`.

The selector hash is an order statistic: it is selected because it is the
minimum hash in a window. Its numerical distribution is therefore skewed
toward small values rather than uniform over the hash range.

The 2026-07-31 investigation established:

- no observed selector hash exceeded `2^32 - 1`;
- at HLL precision 24, approximately 3--6% of observed minimizer hashes were
  below `2^24`, depending on `k,w`;
- those hashes produce `rest == 0` in the current HLL update and enter the
  maximum rho bucket for a nominal 64-bit hash;
- the corresponding fraction was at most approximately 0.04% at precision 18
  and below.

The mechanism is minimizer small-value bias, not simply that a 32-bit value
leaves too few bits for rho.

## Important unresolved questions

Do not repeat unverified claims from the earlier investigation. Establish:

1. The sign and magnitude of cardinality bias at each precision.
2. Whether the `rest == 0` spike is the only practically important
   consequence of non-uniform selector hashes.
3. How union and inferred intersection estimates are affected.
4. How much `jaccard_similarity_ie`, containment, and cosketch outputs move.
5. Whether register equality changes materially under an independent hash.
6. Whether any numerical change improves exact feature-set estimation,
   agreement with coordinate Jaccard, cross-reference identity, or tissue
   recovery.
7. Whether the current `--seed` documentation is inaccurate in sequence
   mode. In the normal minimizer path, raw `digest` hashes currently bypass
   hammock's xxh64 seed; only the no-minimizer whole-record fallback uses it.

## Metrics affected

The current non-uniform inputs can affect HLL register states and therefore
potentially every sequence-mode output. Distinguish two mechanisms:

- Cardinality-derived outputs—`jaccard_similarity_ie`,
  `containment_AB`, `containment_BA`, and the cosketch summaries—depend on
  HLL cardinality estimates and are the primary concern.
- `reg_eq_similarity` does not call the cardinality estimator, but changing
  the ingested hashes changes register assignments and values. Its numerical
  output can therefore change after a fix even though it is not affected by
  cardinality-estimator bias in the same way.

Do not describe exact BEDTools coordinate Jaccard as ground truth for sequence
similarity. Use two distinct validation targets:

1. exact Jaccard/cardinality of the selected sequence feature sets, which is
   estimator ground truth;
2. BEDTools coordinate Jaccard and biological labels, which are external
   application objectives.

## Candidate designs

Keep the legacy behavior available during evaluation. Compare at least these
three candidates behind an explicit experimental option:

### A. Legacy raw selector hash

Feed the 32-bit `digest` selector hash directly to the nominal 64-bit HLL.
This is the baseline and compatibility mode.

### B. Independently rehashed selector identity

Serialize the 32-bit selector hash as fixed-width binary data with an explicit
domain tag, then hash it with seeded xxh64 before HLL ingestion.

Advantages:

- inexpensive;
- restores a full-width, approximately uniform HLL input;
- preserves the existing selected feature identity.

Limitation:

- collisions already created by the 32-bit selector hash cannot be recovered.

Never hash the decimal text representation of the selector hash. Use a
documented fixed byte order and domain separation.

### C. Independently hashed canonical minimizer sequence

Use the minimizer position returned by `digest`, recover the selected
`k`-mer from the FASTA record, canonicalize it against its reverse
complement, and hash that canonical sequence independently with seeded xxh64
for HLL ingestion.

Advantages:

- cleanly separates the hash used for minimizer selection from the hash used
  by HLL;
- supplies uniform 64-bit HLL inputs;
- avoids retaining 32-bit selector-hash collisions as feature identities;
- gives `--seed` meaningful sequence-mode semantics.

Risks and costs:

- additional substring extraction, reverse complementation, and hashing;
- case and ambiguous-base behavior must match `digest`;
- positions and tie behavior must be interpreted correctly;
- every sequence-mode sketch and result changes.

Candidate C is the principled preferred design if its correctness and runtime
are acceptable. Candidate B is a useful lower-cost comparator.

### D. A 32-bit HLL

Constructing the HLL with `hash_size=32` is a diagnostic, not the presumptive
fix. Selector hashes are still minima and therefore not uniform HLL inputs.
Evaluate it to isolate width effects, but do not select it merely because the
types match.

## Hash-domain design

If independent hashing is adopted:

- use a documented fixed-width binary encoding;
- apply a distinct domain tag for minimizer features;
- apply a separate domain tag for whole-record fallback features;
- normalize sequence case explicitly;
- preserve strand invariance;
- define behavior for ambiguous bases;
- ensure identical features hash identically across files and references;
- use `args.seed` consistently;
- test that changing the seed changes the sketch but not expected estimates;
- record hash mode and seed in provenance.

Domain separation prevents an ordinary minimizer and a whole-record fallback
feature from being treated as the same element merely because their numeric
hashes collide across feature types.

## Multistep execution plan

### 1. Freeze the current sequence-mode baseline

Record:

- code and manuscript commits;
- Python, compiler, `digest`, xxhash, and BioPython versions;
- exact commands and input checksums;
- all current sequence-mode CSVs and figure-source tables;
- checksums of rendered figures;
- current runtime measurements where available.

Identify every manuscript dependency using repository search and workflow
inspection. At minimum include:

- cross-reference identity analyses;
- Maurano tissue clustering;
- sequence precision sweeps;
- parameter-objective tradeoff analyses;
- supplementary sequence results;
- all captions and prose reporting sequence-mode metrics.

Gate:

- Reproduce existing baseline CSVs exactly before modifying ingestion.
- If baseline reproduction fails, stop and resolve environment or provenance
  drift first.

### 2. Reproduce and extend the hash-distribution diagnosis

Build a versioned diagnostic script that records, for each dataset and
`k,w,p`:

- total minimizer occurrences;
- number of distinct selector hashes;
- selector-hash bit-width distribution;
- mean and quantiles relative to `2^32`;
- fraction below `2^p`;
- HLL register histogram;
- count in the maximum rho bucket;
- HLL load `n / 2^p`.

Cover:

- synthetic random DNA;
- repetitive DNA;
- reverse-complement pairs;
- sequences containing ambiguous bases;
- the Maurano FASTAs;
- cross-reference FASTAs;
- default `k=8,w=40`;
- biological optimum `k=10,w=30`;
- numerical optimum `k=20,w=20`;
- all precisions used in the manuscript, especially 18, 21, 23, and 24.

Gate:

- Reproduce the confirmed 32-bit bound and high-precision spike.
- Save machine-readable output and plotting code.
- Do not infer cardinality bias from the register spike alone.

### 3. Establish exact feature-set truth

For manageable datasets, materialize:

- the set of distinct `digest` selector hashes;
- the set of distinct selected canonical minimizer sequences;
- exact individual-set, union, intersection, Jaccard, and containment values
  for both feature definitions.

Deduplicate exactly as HLL set semantics require. Process FASTA records
independently so no `k`-mer crosses a record boundary.

Gate:

- Verify reverse-complement invariance.
- Verify case normalization.
- Document ambiguous-base and no-minimizer behavior.
- Clearly label which exact feature definition is the comparator for each
  candidate.

### 4. Measure the legacy estimator directly

Feed known-distinct legacy selector hashes into HLL sketches across a precision
sweep and compare:

- estimated versus exact set cardinality;
- estimated versus exact union cardinality;
- inferred versus exact intersection size;
- inclusion--exclusion Jaccard versus exact feature-set Jaccard;
- containment estimates versus exact containment;
- register histograms and maximum-rho occupancy.

Use multiple synthetic inputs and, where the candidate supports it, multiple
independent seeds. Report signed bias, absolute error, relative error, and
variance.

Gate:

- Resolve the sign and magnitude of the legacy bias.
- Determine whether the effect is practically relevant at precision 18 and at
  higher precisions.
- If no meaningful error is found, preserve that result; do not manufacture a
  rationale for changing production behavior.

### 5. Implement candidate modes behind an experimental switch

Add explicit modes such as:

- `legacy-selector32`;
- `rehash-selector64`;
- `canonical-kmer64`;
- optionally `hll32-diagnostic`.

Do not silently replace the default. Keep output filenames or provenance
sufficiently distinct to prevent accidental comparison of incompatible
sketches.

Gate:

- Existing legacy-mode outputs remain byte-identical.
- Hash mode and seed are visible in verbose output or structured provenance.
- Incompatible sketches cannot be compared or merged without detection.

### 6. Add unit and property tests

Test:

- selector hashes remain the same across candidate modes;
- candidate B matches a fixed set of binary-encoding golden vectors;
- candidate C selects the same minimizer positions as legacy mode;
- canonical `k)-mers are strand invariant;
- uppercase and lowercase inputs agree;
- ambiguous-base behavior matches the stated contract;
- repeated minimizers have set semantics;
- record boundaries are respected;
- short-record fallback is deterministic and domain-separated;
- changing `--seed` changes independently hashed sketches;
- serial and threaded executions agree;
- saved golden HLL register arrays are stable for fixed mode and seed;
- union and pairwise comparisons reject incompatible hash modes.

Include randomized tests against exact feature sets.

Gate:

- Legacy tests remain green.
- New candidates pass exact feature-identity tests.
- No silent fallback occurs when `digest` is unavailable.

### 7. Compare estimator correctness across candidates

For every candidate, sweep relevant `k,w,p` and compare HLL estimates to the
corresponding exact feature-set values.

Primary outcomes:

- cardinality relative bias and MAE;
- union and intersection error;
- inclusion--exclusion Jaccard MAE;
- containment MAE;
- precision dependence;
- maximum-rho artifacts.

Register equality should be reported descriptively, not judged as calibrated
set Jaccard.

Gate:

- A production candidate must materially improve or eliminate the
  high-precision artifact.
- It must not worsen estimator accuracy at the command-line default precision.
- Candidate C should beat candidate B enough to justify its added complexity;
  otherwise prefer the simpler independently rehashed selector identity.

### 8. Benchmark runtime and memory

Measure separately:

- FASTA parsing;
- minimizer selection;
- canonical sequence recovery;
- independent hashing;
- HLL updates;
- pairwise comparison;
- total wall time and peak memory.

Use matched inputs, hardware, thread settings, and repeated runs. Candidate C
may need a batched C++ helper if Python substring/reverse-complement work is a
material bottleneck. Optimize only after correctness is established.

Gate:

- Predefine an acceptable runtime regression.
- Do not hide a substantial sketch-construction slowdown behind unchanged
  pairwise time.
- Reject any optimization that changes feature identity without a new
  validation cycle.

### 9. Re-run application-level sequence analyses

Run legacy and viable candidates on:

- exact feature-set validation;
- BEDTools coordinate-Jaccard comparison;
- cross-reference identity matching;
- Maurano tissue clustering;
- precision stability;
- parameter-objective tradeoffs;
- any supplementary sequence analysis.

Compare:

- every pairwise output;
- Pearson, Spearman, Kendall, and MAE where used;
- ARI and NMI;
- dendrogram topology, cluster cuts, and distance matrices;
- selected optimal `k,w,p`;
- figure-source data and rendered figures.

Gate:

- Improved HLL estimator accuracy is necessary but not sufficient.
- A production change must not silently damage the biological analyses.
- If the biological optimum moves, investigate whether the change reflects
  improved feature estimation or merely a new hash realization.

### 10. Perform seed-sensitivity analysis

Independent hashing makes `--seed` active in ordinary sequence mode. Run
multiple seeds for candidates B and C.

Measure:

- estimator variance;
- stability of pairwise rankings;
- stability of tissue clustering and cross-reference matches;
- stability of selected parameter optima.

Gate:

- Manuscript conclusions must not depend on a lucky single seed.
- If meaningful seed sensitivity remains, report uncertainty or aggregate
  across seeds rather than selecting the most favorable run.

### 11. Choose behavior with an explicit decision record

Possible outcomes:

1. **Legacy error is negligible at all supported precisions.**
   Keep legacy behavior, correct documentation, and consider restricting or
   warning on problematic high precision.
2. **Rehashed selector identity fixes the issue adequately.**
   Prefer candidate B if candidate C adds cost without meaningful scientific
   benefit.
3. **Canonical minimizer hashing clearly improves correctness.**
   Select candidate C and document the intentional feature-identity change.
4. **No candidate meets accuracy, stability, and runtime gates.**
   Do not change the default; continue investigation.

Write a short decision document containing evidence, rejected alternatives,
and consequences for backward compatibility.

### 12. Apply a manuscript result gate

Every independent-hash candidate will generally change HLL registers and
sequence-mode outputs, even at precision 18. Do not assume existing manuscript
numbers remain valid.

Classify each reported result:

1. **Numerically identical:** retain the value and record the exact comparison.
2. **Changed below displayed precision:** update archived source data and
   verify that rounded prose remains true.
3. **Meaningfully improved:** regenerate dependent figures and tables and
   propose exact text updates.
4. **Meaningfully worse or conclusion-changing:** do not update the manuscript
   until the method choice is reconsidered.

Draft manuscript changes only in the disposable copy. Include:

- separation of minimizer selection hash from HLL ingestion hash;
- canonical and strand-invariant feature definition;
- hash algorithm, width, seed, byte encoding, and domain tags;
- short-record fallback behavior;
- rerun statistics and figure provenance;
- any change to default parameters or biological conclusions.

### 13. Final verification

Before recommending merge:

- run the complete Python and C++ test suites from a clean build;
- run sanitizer and property tests;
- reproduce chosen-candidate outputs from a clean clone;
- build the complete manuscript draft;
- run `git diff --check` in both copies;
- produce an artifact manifest with input hashes, commands, environments,
  outputs, and figure checksums;
- make separate local commits for code, experiments/results, and manuscript
  drafts;
- report all pre-existing or unrelated failures explicitly.

## Recommended acceptance thresholds

Set final numerical thresholds before viewing candidate results. At minimum:

- no detectable systematic cardinality bias against exact feature sets at
  manuscript precisions;
- lower or equal inclusion--exclusion Jaccard MAE than legacy at precision 18;
- elimination of the anomalous maximum-rho spike at high precision;
- no material degradation in cross-reference identity classification;
- no material degradation in Maurano ARI/NMI or unexplained topology changes;
- stable conclusions across multiple ingestion seeds;
- deterministic results for fixed mode and seed;
- acceptable runtime and memory overhead.

If candidate C changes the exact feature universe relative to legacy, compare
both candidates against their own exact feature sets and state that distinction
clearly.

## Stop conditions

Stop and report rather than silently adapting if:

- the baseline cannot be reproduced;
- `digest` behavior differs from the recorded version;
- exact canonical feature recovery cannot be made consistent with selector
  positions;
- a fix improves synthetic cardinality but worsens application results
  materially;
- results depend strongly on ingestion seed;
- candidate sketches can be accidentally compared across hash modes;
- required manuscript changes lack explicit authorization for the canonical
  repository;
- completing validation would require unapproved scheduler or external
  compute.

## Required handoff

Report:

- disposable clone paths and branches;
- local commits;
- confirmed and rejected hypotheses;
- exact feature-set definitions;
- candidate comparison tables;
- test and benchmark results;
- all changed application-level metrics;
- proposed default and compatibility mode;
- proposed manuscript diff;
- remaining risks and unexecuted large-scale runs.
