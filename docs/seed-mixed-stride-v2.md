# Seed: chromosome-wide mixed-stride v2

## Objective

Redesign hammock's mixed-stride interval subsampling so that non-integral
reciprocal rates mix two gap lengths within each chromosome while preserving
the existing coordinate grid exactly whenever `1 / subB` is integral.

The principal acceptance criterion is:

> Improve arbitrary-rate behavior without changing any sampled coordinate,
> HLL register, pairwise value, or other numerical manuscript result for the
> published `subB=0.1` and `subB=0.01` configurations.

Do not assume that a changed result is acceptable merely because it appears
better. Freeze and reproduce the baseline, predefine comparison criteria, and
trace every changed manuscript value to regenerated source data.

## Terminology and methodological framing

The phrase `mixed-stride subsampling` does not appear to be an established,
standardized name in the sampling literature. Treat it as hammock-specific
terminology and define it explicitly rather than claiming that it is a
generally recognized algorithm.

The proposed v2 conforms better than v1 to the literal meaning of
mixed-stride sampling:

- v1 chooses either `S0` or `S1` for an entire chromosome, so an individual
  chromosome uses only one stride;
- v2 intersperses `S0`- and `S1`-length gaps within the same chromosome;
- the mean gap in v2 is `1/p`, so every sufficiently long chromosome
  approaches the requested density;
- a deterministic chromosome-specific phase preserves coordinate consistency
  across files and interval boundaries.

The most defensible general description is:

> a deterministic fractional-interval systematic sampling scheme that mixes
> the two adjacent integer strides required to realize the requested sampling
> rate.

Fractional-interval and unequal-interval designs are recognized extensions of
ordinary systematic sampling, although the name mixed-stride is not itself a
standard term. Useful background sources are:

- Hankin, Mohr, and Newman, *Sampling Theory: For the Ecological and Natural
  Resource Sciences*, systematic sampling chapter:
  https://academic.oup.com/book/43702/chapter-abstract/367153227
- Singh and Tailor, `Linear Systematic Sampling with Unequal Sampling
  Intervals in the Presence of Linear Trend`:
  https://ese-journals.unisalento.it/index.php/ejasa/article/view/18988

Like other systematic sampling designs, v2 can interact with periodic
structure in the sampled coordinates. A chromosome-specific phase reduces
fixed alignment but does not eliminate that general risk. Retain the planned
comparisons with hash-threshold sampling and explicitly test synthetic
periodic interval patterns.

## Repositories and current handoff state

- Code repository: `~/interval_sketch/hammock_claude`
- Manuscript repository: `~/interval_sketch/hammock_paper`
- Code baseline at handoff: `d685a65 Strengthen mixed-stride sampling tests`
- Manuscript baseline at handoff:
  `786076e Correct sequence methods and document stride audit`
- Date of this handoff: 2026-08-21

The code worktree contained unrelated modified and untracked experiment
artifacts at handoff. They belong to the user. Do not stage, alter, clean, or
delete them. Scope commits by explicit path.

The manuscript has a mixed-stride implementation-audit
`\begin{comment}` block under the mixed-stride Methods subsection. It records
the findings summarized below. The rendered mixed-stride prose has not yet
been corrected.

## Autonomous operation in disposable copies

Do not perform implementation work in either primary repository listed above.
At the beginning of the session, create independent disposable clones beneath
a newly allocated directory in `/tmp`. For example:

```bash
mktemp -d /tmp/hammock-mixed-stride-v2.XXXXXX
git clone --local ~/interval_sketch/hammock_claude WORK_ROOT/hammock
git clone --local ~/interval_sketch/hammock_paper WORK_ROOT/hammock_paper
git -C WORK_ROOT/hammock remote set-url origin \
  git@github.com:jessicabonnie/hammock.git
git -C WORK_ROOT/hammock_paper remote set-url origin \
  git@github.com:jessicabonnie/hammock_paper.git
```

Replace `WORK_ROOT` with the exact directory returned by `mktemp`; do not
use an unresolved environment variable in commands that modify or remove
files. Create a dedicated local branch in each clone. Do not use `git
worktree add`, because that would mutate administrative state in the primary
repository and share its object database. Each local clone is an independent
Git working tree and is the only place this task may make changes.

After cloning the manuscript, run `git pull --ff-only` inside the manuscript
copy against its restored GitHub remote before editing it so that the copy
includes the current Overleaf state. A local clone initially points `origin`
at the primary filesystem checkout, which is why the commands above restore
the canonical remote URL. If the pull conflicts or cannot fast-forward, stop
manuscript work and report the problem. The primary manuscript checkout must
remain untouched.

The user authorizes autonomous work inside these disposable copies. Without
asking for further confirmation, the session may:

- create branches and edit files;
- build and run tests;
- run local correctness and performance experiments;
- generate temporary fixtures and derived artifacts;
- draft and compile manuscript changes in the manuscript copy;
- make local commits in either copy;
- iterate on failures until the acceptance gates are met.

This authorization is limited to the disposable copies and local computations.
It does not authorize:

- editing, staging, committing, cleaning, or resetting either primary
  repository;
- pushing branches or commits to any remote;
- changing the canonical Overleaf manuscript;
- deleting or overwriting source datasets;
- submitting scheduler jobs or using substantial external compute unless
  separately authorized;
- bypassing platform-enforced approval prompts for network or privileged
  operations.

Leave the disposable copies in place at the end of the session unless the user
explicitly requests cleanup. Report their exact paths, local branches, commits,
tests, benchmarks, and proposed manuscript diff so another session or the user
can inspect them.

## Mandatory manuscript safeguards

The user has an explicit standing instruction:

1. Never change manuscript text without the user's explicit instruction.
2. Immediately before any authorized manuscript change, run
   `git pull --ff-only` in `~/interval_sketch/hammock_paper`.
3. Preserve unrelated changes and stop if the pull conflicts.

For this seed, the autonomy grant above is explicit permission to draft and
test manuscript edits only inside the disposable manuscript copy. It is not
permission to modify or push the primary manuscript. Before any later transfer
of those drafts into the primary manuscript, pull the primary manuscript again
and obtain the user's explicit approval for the exact text changes.

## Verified current implementation behavior

The production implementation is in:

- `cpp/src/stride.cpp`
- `cpp/include/hammock/stride.hpp`

For `0 < p < 1`, the normal mixed-stride path computes:

```
S0 = floor(1 / p)
S1 = ceil(1 / p)
```

When `S0 != S1`, it computes:

```
q = (p - 1/S0) / (1/S1 - 1/S0)
```

The code then hashes the chromosome name and chooses `S1` with probability
`q`; otherwise it chooses `S0`. It chooses only one stride for the entire
chromosome. It does not mix gap lengths within a chromosome.

A separate hash of `chr + "|res"` selects a residue modulo the chosen
stride. The retained coordinates therefore satisfy one chromosome-wide
congruence relation:

```
x mod S = residue
```

Every interval on the same chromosome recomputes the same stride and residue,
so normal-path sampling is independent of interval boundaries.

The current term "mixed" describes a mixture of chromosome-level assignments
or assignments across gate seeds, not a mixture of gap lengths within one
chromosome.

At the manuscript's evaluated rates:

- `p=0.1`: `S0=S1=10`
- `p=0.01`: `S0=S1=100`

Consequently, a correctly designed v2 can and should retain exactly the same
coordinate grids for both reported rates.

## Confirmed prose and implementation issues

### Active manuscript probability sentence

The rendered Methods prose says `S0` is selected with probability `q` and
`S1` with probability `1-q`. This is reversed. The equation and code use
`q` as the probability of `S1`.

### Rate semantics

The current method matches `p` only in expectation over the pseudorandom
chromosome-to-stride assignment. A fixed chromosome is sampled at
approximately `1/S0` or `1/S1`, not at `p`. A fixed finite collection of
chromosomes can therefore deviate from the requested rate, particularly when
datasets have different chromosome composition.

### Long chromosome-name fallback

When `chr + separator + coordinate` cannot fit the fixed stack buffer
(roughly chromosome names longer than 43 characters with the default
separator), the fallback path:

- uses `round(1/p)` rather than mixed-stride interpolation;
- starts stepping at each interval's start;
- does not use a chromosome-wide residue.

This violates interval-boundary invariance and the manuscript description.

### Zero sampling rate

Both CLIs accept `subB=0`, while mixed-stride evaluates `1/p`. Decide and
implement explicit semantics. Returning an empty sample is recommended, but
rejecting zero specifically for mixed-stride is also defensible.

### Similarity invariance

The manuscript's qualitative explanation is sound: when the same
coordinate-anchored sampling grid is applied to both sets, individual-set,
intersection, and union cardinalities tend to shrink together, leaving
Jaccard and containment comparatively stable. This is approximate because
systematic sampling can interact with short intervals or periodic structure.

## Tests already added

`cpp/tests/test_stride.cpp` now includes:

- a `p=0.3` test showing that chromosome `1` selects stride 3 and
  chromosome `chrX` selects stride 4 with the default gate seed;
- a test showing that partitioning a chromosome span into adjacent intervals
  produces the same HLL as sketching the span as one interval.

At handoff, the focused tests passed:

```bash
/tmp/hammock-codex-tests/cpp/tests/hammock_tests --test-case="*mixed-stride*"
```

Result:

```
3 test cases passed
7 assertions passed
```

The complete C++ suite had one unrelated pre-existing failure:

```
Mode A: chr/non-chr prefixes normalize identically
```

The full result was 23 of 24 test cases passing. Do not attribute that failure
to mixed-stride without new evidence.

The temporary test build was configured with:

```bash
cmake -S . -B /tmp/hammock-codex-tests \
  -DHAMMOCK_BUILD_TESTS=ON \
  -DHAMMOCK_BUILD_CLI=OFF \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build /tmp/hammock-codex-tests --target hammock_tests -j2
```

CMake fetches pinned doctest v2.4.11 and may require network approval on a
fresh machine.

## Proposed v2 mathematical contract

For `0 < p <= 1`, define:

```
S0 = floor(1 / p)
S1 = ceil(1 / p)
```

If `S0 == S1`, preserve the current chromosome-wide residue grid exactly.
This branch is required for bitwise parity with the published rates.

If `S0 != S1`, mix gaps within each chromosome. The fraction `r` of gaps
having length `S1` must satisfy:

```
(1 - r) * S0 + r * S1 = 1 / p
r = (1/p - S0) / (S1 - S0)
```

Do not reuse the current `q` formula. The current formula mixes reciprocal
densities across whole-chromosome assignments; the new formula mixes gap
lengths within one coordinate grid.

Implement the grid as a chromosome-anchored Beatty/Bresenham-style sequence
using fixed-point integer arithmetic. The implementation must:

- use both `S0` and `S1` when required;
- derive a deterministic phase from chromosome and gate seed;
- be anchored to genomic coordinates rather than interval starts;
- calculate the first sampled coordinate in an interval directly;
- jump only between retained positions, preserving the speed advantage over
  per-position hash-threshold sampling;
- avoid platform-dependent floating-point accumulation;
- use the same path for chromosome names of every length;
- remain deterministic across input order and thread count.

## Execution plan and gates

### 1. Freeze and reproduce the baseline

Before changing the algorithm:

- record exact code and manuscript commits;
- locate every mixed-stride-dependent workflow, result table, figure, caption,
  and prose value;
- archive source CSVs, figures, checksums, commands, environment information,
  and scheduler resources;
- reproduce the current reported analysis once using the current algorithm.

Gate:

- Numerical outputs must reproduce exactly.
- Timing must reproduce within normal run-to-run variation.
- If the baseline cannot be reproduced, stop and diagnose it before coding.

### 2. Write an executable specification

Specify behavior for:

- `p=1`;
- integral reciprocal rates;
- non-integral reciprocal rates;
- `p=0`;
- short and empty intervals;
- overlapping intervals;
- long contig names;
- negative coordinates if the parser permits them;
- different seeds;
- thread and input-order independence.

Manually work examples for `p=0.3`, `0.1`, `0.01`, and rates near zero
and one.

Gate:

- State a density-discrepancy bound.
- Demonstrate why interval partitioning cannot change the coordinate grid.

### 3. Create a slow reference oracle

Implement a transparent coordinate generator for tests, preferably with exact
or fixed-point arithmetic.

Gate:

- Density converges to `p` for non-integral rates.
- Integral reciprocal rates match the legacy coordinate grid exactly.
- Partitioned and unpartitioned intervals yield identical coordinates.

### 4. Expand tests before switching production behavior

Cover:

- `p=1`, `0.1`, `0.01`, `0.3`, and irregular values;
- both gap lengths within one chromosome;
- exact legacy parity for integral reciprocal rates;
- split, merged, overlapping, and reordered intervals;
- chromosome and seed changes;
- long chromosome names;
- `subB=0`;
- short intervals and coordinate boundaries;
- serial and multithreaded execution;
- randomized property comparisons against the reference oracle.

Gate:

- Existing tests remain green except for independently documented
  pre-existing failures.
- New non-integral behavior tests fail against v1 for the intended reason.

### 5. Implement behind a temporary strategy name

Initially expose the implementation as `mixed-stride-v2` while retaining v1.
Add method/version provenance to verbose output or structured output metadata.

Gate:

- Unit, property, integration, and sanitizer tests pass.
- Integral-rate sampled coordinates and HLL register arrays are byte-identical
  between v1 and v2.
- Results are identical across thread counts.

### 6. Run compact end-to-end comparisons

On small synthetic BED collections, compare v1 and v2 for:

- retained coordinate sets;
- realized rate by chromosome and overall;
- HLL registers;
- set, intersection, and union estimates;
- Jaccard and containment;
- runtime.

Gate:

- Require exact equality at `p=0.1` and `p=0.01`.
- Require lower rate error and interval-boundary invariance at non-integral
  rates.
- Stop if accuracy or runtime degrades materially without explanation.

### 7. Re-run all reported manuscript configurations

At minimum, reproduce the Maurano analysis using:

- `subB=1.0`, `0.1`, and `0.01`;
- HLL precision 18;
- eight threads;
- the same input files, estimator, environment, and output processing.

Compare every pairwise value, MAE, summary statistic, figure source, and
rendered figure. For runtime, use the same hardware and allocation and compare
replicate medians and ranges rather than requiring exact equality.

Gate:

- Pairwise numerical outputs at 0.1 and 0.01 must be identical.
- If they differ, do not update manuscript values; fix the parity violation.
- Treat runtime as changed only when the difference exceeds ordinary
  replicate variation.

### 8. Evaluate where v2 should improve

Run rates such as `0.3`, `0.15`, `0.07`, and `0.025`. Compare v1, v2,
full sampling, and hash-threshold sampling on:

- per-chromosome and global rate error;
- Jaccard and containment error;
- chromosome-composition bias;
- runtime;
- determinism and interval-boundary invariance.

Predefine improvement as:

- smaller deviation from requested density;
- no worse, preferably lower, similarity error;
- no material runtime regression;
- exact partition invariance;
- deterministic results;
- no pathological degradation on synthetic periodic interval patterns.

### 9. Apply a manuscript decision gate

There are four possible outcomes:

1. Reported numerical results are identical: update only the authorized
   Methods description.
2. Raw values differ below displayed precision: update archived source data,
   verify all rounded claims, and document the comparison internally.
3. Results meaningfully improve: regenerate every dependent figure and table
   and propose exact prose/caption updates to the user.
4. Results regress or behave inconsistently: do not make v2 the default.

Before any manuscript modification, obey the mandatory pull and explicit-text-
authorization safeguards above.

### 10. Finalize

- Make v2 the default only after all gates pass.
- Decide whether v1 remains available under a legacy name.
- Update CLI documentation and method provenance.
- Resolve or remove the manuscript audit comment only with user authorization.
- Run a clean full build and test suite.
- Build the complete manuscript after authorized edits.
- Produce a manifest linking reported statistics and figures to commands and
  source artifacts.
- Commit code and manuscript changes separately, staging only explicit paths.

## Recommended stop conditions

Stop and report rather than silently adapting if:

- integral reciprocal rates are not byte-identical to v1;
- the baseline manuscript results cannot be reproduced;
- v2 reduces rate error but worsens similarity accuracy materially;
- timing comparisons were run on unmatched hardware or scheduler conditions;
- implementing the grid requires changing estimator semantics;
- manuscript changes are needed but have not been explicitly authorized.
