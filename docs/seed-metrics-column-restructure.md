# Seed: restructure default/`--register-equality`/`--metrics` output columns and filenames

**Status: spec accepted 2026-08-11, implementation not started.** This is an
implementation plan, not a research seed — the "evidence gathered" is the
current-state inventory below; "what still needs establishing" is whether the
concrete design in this doc survives contact with the ~36 dependent files it
touches. Work happens on a worktree (see "Where to work," below) so the two
currently-running SLURM jobs against `main` are undisturbed. Follow the Parts
in order; each is scoped to be picked up cold in a fresh session with cleared
context — read this whole doc first, then jump to the Part you're resuming.
**Part 1 is a read-only payoff analysis — read it before starting Part 2's
implementation work, since it scopes what Part 9 (the benchmark-retargeting
follow-up) is actually for.**

## The spec (user request, 2026-08-11)

Three output shapes, selected by (new) mutually-exclusive flags, replacing
the current two-shape `--metrics`/`--no-metrics` contract:

1. **Default (no flag): `jaccard_similarity_ie` only.** Today's default is
   the full 7-metric-column block; that block is no longer the default for
   *either* front-end.
2. **`--register-equality` / `--re`** (replaces `--no-metrics`, which is
   removed, not aliased — see "Decisions," below): two columns,
   `jaccard_similarity` and `register_equality_similarity`, the latter a
   literal duplicate of the former. Cheap: this is (and must remain) the
   register-equality-only timing arm, i.e. it must not compute the
   union/containment path.
3. **`--metrics`**: everything — the current 7-column block
   (`jaccard_similarity`, `jaccard_similarity_ie`, `containment_AB`,
   `containment_BA`, `cosketch_geom`, `cosketch_arith`, `cosketch_max`) plus
   `register_equality_similarity` appended, duplicating `jaccard_similarity`.
   8 metric columns total.

Output filenames must reflect the shape (see "Filename tags," below).

## Decisions made in this seed (treat as settled; revisit only if a Part
hits a concrete problem with one)

Resolved via `AskUserQuestion` before any code was touched:

1. **Filename tags: all three shapes get an explicit tag, none stays
   bare.** Today only the reduced shape (`--no-metrics`) is tagged (`_j3`);
   the default is untagged. Reusing "untagged" for the *new* default would
   collide with the meaning "untagged" has carried in every archived,
   pre-existing full-block CSV — a script doing positional column access
   (`row[2]`) against an old untagged file expects `jaccard_similarity`
   there; against a new untagged file it would silently get
   `jaccard_similarity_ie` instead. Tagging all three removes that failure
   mode entirely. This is exactly the class of bug this repo's memory
   already documents once (the Villar `_with_ends` awk bug,
   `memory/project_villar_peaks_zero_with_ends.md` — not in `CLAUDE.md`
   itself) — parse by name, and here, never let two semantically different
   shapes share a bare filename.
2. **`--no-metrics` is removed outright, not kept as a deprecated alias.**
   Matches standing project practice (`memory/project_no_external_api_consumers.md`:
   "no API consumers outside this repo — delete rather than shim; removals
   are not breaking changes"). The ~36 dependent files that pass
   `--no-metrics`/`--metrics` today are inventoried below and updated in
   Part 6, in this same worktree, before merge — so nothing on `main` ever
   observes a window where the flag is gone but callers still expect it.
3. **Both front-ends move together, kept in bit-for-bit sync**, matching
   existing repo convention (`CLAUDE.md` divergence #9; the cross-tool `==`
   gate in `tests/test_hammock_cpp_metrics.py`). The Python CLI currently has
   *no* `--metrics`/`--no-metrics` flag at all (see "Current state," below)
   — this seed adds the full three-shape surface to it for the first time,
   not just a rename.

### Column layout, concretely

| Shape | Flag | `similarity_measures` (A/B/C; Mode D is the same list, `ref1`/`ref2` still appended after) | Cost |
|---|---|---|---|
| default | *(none)* | `["jaccard_similarity_ie"]` | needs the fused union pass + per-sketch cardinalities — same cost family as `--metrics`, NOT the cheap arm (see note below) |
| register-equality | `--register-equality`/`--re` | `["jaccard_similarity", "register_equality_similarity"]` (2nd is `= 1st`) | cheap — skips the union pass entirely, must stay that way |
| full | `--metrics` | `["jaccard_similarity", "jaccard_similarity_ie"] + _CONTAINMENT_COLS + ["register_equality_similarity"]` (8 cols; last `= jaccard_similarity`) | current full cost, plus one more formatted field for the duplicate — not separately measured, see Part 1's caveats |

**Cost note for whoever implements Part 4 (C++):** `HLLSketch::jaccard_and_union_cardinality`
computes register-equality Jaccard *and* the union cardinality in one fused
pass (`cpp/include/hammock/hll_sketch.hpp:67`; written up in `CLAUDE.md`'s
"Open seeds" section, under the `docs/seed-mode-d-hash-width.md` bullet's
"pairwise metrics loop is one fused register pass" paragraph — that's where
the fused-pass mechanism is actually documented, not in the seed file
itself) — so the IE-only default gets `jaccard_similarity` for free as a
side effect of computing the union it needs anyway; it simply isn't written
to the CSV. The register-equality arm is the only one of the three that can
skip the union pass, which is precisely why it must be its own branch and
not "full minus some columns."

### Filename tags, concretely

Appended last in the suffix chain, same position `_j3` occupies today.
`cpp/app/hammock_cli.cpp`'s `outprefix_with_suffix` builds
mode/expA-or-subA/subB, in that order, then (today) `_j3` — it has **no**
k/w component, since `hammock-cpp` only supports Modes A/B/C, never Mode D.
The new tag goes last in that chain too. The Python side additionally has a
k/w component (Mode D only, `outprefix.py::get_new_prefix`), and the new tag
goes after that on the Python side.

| Shape | Tag |
|---|---|
| default | `_ie` |
| `--register-equality`/`--re` | `_re` |
| `--metrics` | `_full` |

E.g. `myrun_hll_p18_jaccB_ie.csv`, `myrun_hll_p18_jaccB_re.csv`,
`myrun_hll_p18_jaccB_full.csv`. This applies to **both** front-ends —
currently only `hammock-cpp` tags at all; the Python CLI (`outprefix.py`)
gains a metrics-mode tag for the first time.

## Current state (inventory, established 2026-08-11 by direct code read —
re-check line numbers before editing, this doc will drift)

**Python CLI (`hammock`) has no `--metrics`/`--no-metrics` flag today.**
`python/hammock/cli.py` defines no such argument. `runner.py` unconditionally
writes the full 7-column block in both the A/B/C path (`runner.py:579`
region, `similarity_measures = ["jaccard_similarity", "jaccard_similarity_ie"] + _CONTAINMENT_COLS`)
and the Mode D path (`_write_mode_d_csv`, `runner.py:397` region, identical
list). `python/hammock/outprefix.py`'s `get_new_prefix` has no metrics
parameter and never tags for it.

**`hammock-cpp` (`cpp/app/hammock_cli.cpp`) has the current two-shape
contract:** `Args.metrics` (bool, default `true`). `--metrics`/`--no-metrics`
parsed at line ~202-205. Header + stride branch at line ~360-370 (`stride = args.metrics ? 7 : 1`).
Per-pair branch at line ~405-430 (`if (!args.metrics) { cell[0] = ...jaccard_similarity...; continue; }`
else the fused/containment path filling `cell[0..6]`). Filename tag at line
~285 (`if (!a.metrics) out += "_j3";`). Help text at line ~90-155 documents
the current two-shape contract and must be rewritten.

**Test inventory** (files referencing `jaccard_similarity`/`containment_AB`/
`cosketch_*`/`jaccard_similarity_ie`, i.e. touched by column-shape changes):
`tests/test_hammock_cpp_metrics.py` (the cross-tool bit-for-bit gate — full
rewrite, Part 4), `tests/test_bed2fasta_cli.py`, `tests/test_containment_estimator.py`,
`tests/test_jaccard_ie.py`, `tests/test_mode_d_parity.py`, `tests/test_mode_d.py`,
`tests/test_parity_against_original.py` (Part 3 — each needs auditing for
whether it invokes the CLI with no explicit flag and then expects
now-non-default columns).

**Dependent scripts/docs referencing `--no-metrics` or `--metrics`**
(`grep -rln -- "--no-metrics"` / `"--metrics"`, 2026-08-11; union of both
lists, repo-relative):
```
README.md
CLAUDE.md
cpp/app/hammock_cli.cpp
experiments/bedtools_benchmark/RESULTS.md
experiments/bedtools_benchmark/README.md
experiments/bedtools_benchmark/sbatch_cli_overhead.sh
experiments/bedtools_benchmark/sbatch_subB_nsweep.sh
experiments/bedtools_benchmark/sbatch_fig3_panelA.sh
experiments/bedtools_benchmark/sbatch_fig3_panelA_v2.sh
experiments/bedtools_benchmark/sbatch_fig3_panelB.sh
experiments/bedtools_benchmark/sbatch_fig3_largeN.sh
experiments/bedtools_benchmark/sbatch_fusion_ab.sh
experiments/bedtools_benchmark/sbatch_files_t16.sh
experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py
experiments/bedtools_benchmark/measure_cli_overhead.py
experiments/bedtools_benchmark/measure_threads_isolation.py
experiments/bedtools_benchmark/sweep.py
experiments/bedtools_benchmark/fusion_ab.py
experiments/bedtools_benchmark/pairwise_cost_by_precision.py
experiments/subB_mixed_stride/sbatch_synthetic_ie.sh
experiments/subB_mixed_stride/sbatch_maurano.sh
experiments/subB_mixed_stride/run_sweep.py
experiments/subB_mixed_stride/summarize_synthetic_ie.py
experiments/subB_mixed_stride/RESULTS.md
experiments/subB_mixed_stride/README.md
paper/pairwise_scaling/plot_pairwise_scaling.R
paper/pairwise_scaling/plot_pairwise_scaling_supplement.R
paper/pairwise_scaling/plot_subB_nscaling_supplement.R
docs/metrics-by-default.md
docs/seed-hammock-cpp-file-dispatch.md
docs/seed-benchmark-methodology.md
docs/seed-subsampling-synthetic-supplement.md
docs/figure3-panel-a-rebuild.md
docs/paper_outline.md
docs/submittability-concerns.md
docs/bedtools-baseline-retraction.md
```
`docs/seed-*.md` files are historical records of past measurements — do not
edit their prose to match the new flag name (that would misrepresent what
was actually run); only `docs/submittability-concerns.md` and
`docs/metrics-by-default.md` are live design docs that need updates (Part 7).
Everything under `experiments/` and `paper/` is live tooling and needs a real
audit, not a find-replace (Part 6) — see below for why.

**Not mechanical:** every script above passing `--no-metrics` today wants
"the cheap register-equality timing arm," which under the new contract is
exactly `--register-equality`/`--re` — rename only. But every script passing
`--metrics` today wants it because that was the *only* way to get
`jaccard_similarity_ie` at all (the Python CLI had no opt-out, so
"`--metrics` on hammock-cpp" was how the two front-ends were made
comparable). Now that bare (no-flag) is the cheap-ish IE-only shape, some of
those `--metrics` call sites may want to switch to the new default instead
of continuing to pay for the full 8-column block — but that's a timing
characteristic change requiring new benchmark numbers, not a safe mechanical
edit. Flag it per-script in Part 6; don't decide it there — Part 1, below,
scopes whether chasing it (Part 9) is worth doing at all.

## Where to work

Worktree already created, 2026-08-11, off `main` at `256931c`:
```
/home/jbonnie1/interval_sketch/hammock_claude_wt_metrics   (branch: metrics-restructure)
```
Two SLURM jobs were running against `main` at seed-creation time —
`29761450` (`subB_nsw`) and `29763124` (`cli_over`) — both launched from
already-built binaries/checkouts of `main`; nothing in this worktree can
affect them, and there is no need to check `squeue` before working here. Do
not run `pip install -e .` or rebuild `_core`/`hammock-cpp` in the **main**
checkout while those (or any other) jobs may be reading its `build/` — do
all builds inside the worktree.

Baseline to compare against: last known full-suite result before this seed,
`244 passed / 8 skipped` (`docs/seed-hammock-cpp-file-dispatch.md`, v0.7.1,
2026-08-11). Re-run `pytest tests/` in the worktree at the start of Part 2
(the first Part that touches code) to confirm that baseline still holds
there before making changes — if it doesn't, that's pre-existing drift, not
something this seed's changes caused.

## Parts

Each Part below states what it touches, what "done" looks like, and the test
command that proves it. Commit at the end of each Part (small, reviewable
commits on `metrics-restructure`) so a context-cleared session can `git log`
to see exactly how far execution got.

### Part 1 — Payoff analysis: anticipated time improvement, and is this worth doing

**Read-only. No code changes.** Written up front, before Part 2, because the
answer scopes Part 9 rather than being a nice-to-have retrospective — no
sense implementing Parts 2-8 without being honest about what the eventual
benchmark-retargeting payoff actually is.

**Two separable questions, deliberately not conflated:**

1. *Is the core three-shape restructure (Parts 2-8) worth doing?* Yes, and
   this isn't really a cost/benefit call — it's the user's explicit product
   spec (see "The spec," above), independently motivated by a real,
   already-documented correctness hazard: `CLAUDE.md`'s divergence #2 and
   `docs/jaccard-definitional-gap.md` exist because `jaccard_similarity` and
   `jaccard_similarity_ie` are routinely confused, on different scales, not
   substitutable, and the repo has needed multiple corrective passes to keep
   its own docs and figures from mixing them up (`CLAUDE.md`'s "counter-claim
   ... retracted as stated" paragraph is one such correction). A default that
   returns exactly the column this repo already recommends
   (`jaccard_similarity_ie`) removes an entire class of "which column did I
   actually read" mistakes at the source. The engineering cost is bounded and
   enumerated (Parts 2-8, ~2 CLIs + 7 test files + ~36 dependent files, plus
   the no-flag Python-CLI call sites Part 6 item 4 adds) — proceed regardless
   of what question 2 below concludes.
2. *Is Part 9 (retargeting benchmark call sites to the cheaper shape and
   re-measuring) worth doing?* This is genuinely discretionary — it is a
   benchmarking-harness efficiency question, not a correctness one, and it
   is the part this analysis actually needs to earn its keep on.

**Anticipated time improvement for Part 9: direction only — a specific
percentage was attempted here and withdrawn.** An earlier draft of this
section combined two `CLAUDE.md` numbers into a per-precision "modeled
default/metrics ratio" table claiming ~25-32% savings. A three-reviewer pass
over this seed caught that the combination was invalid, not merely
imprecise — recorded here so the mistake isn't repeated by whoever next
touches this section:

- The **block-cost table** (job 29628907, N=64/side, p=12-24) is built from
  `hammock-cpp`'s `pair_time` timer — `cpp/app/hammock_cli.cpp`'s `Pairwise:`
  checkpoint, taken *before* the serial write loop starts (`t_w0` is a
  separate, later checkpoint). By construction, `pair_time` **excludes**
  write/`fprintf` time entirely.
- The **72%/28% write/compute split** was computed from a *different*
  column, `comparison_time` (`Pairwise+write:`, i.e. pairwise compute *and*
  write combined) — the very thing `pair_time` was built to exclude.

Apportioning 72% of a `pair_time` delta to "write cost" therefore
apportions a quantity that delta structurally does not contain — not a
cross-experiment combination with some slack in it, an operation on the
wrong column. Independently, the 72/28 figure traces to a single sweep
fixed at one precision (p=14, per `docs/seed-benchmark-methodology.md`'s
"08-08 files sweep"), and this repo's own `docs/metrics-by-default.md`
pre-fusion table shows that split moving from roughly 80%+ write-share at
p=12 down to under 1% at p=24 (write time ~flat, pairwise compute scaling
~Θ(2^p)) — so treating 72/28 as one constant usable across p=12/18/24, as
the withdrawn table did, contradicts data already in this repo.

**What can honestly be said without a new measurement:**
- **Direction**: a call site that only needs `jaccard_similarity_ie` and
  switches from `--metrics` (8 columns under the new contract) to the bare
  default (1 column) should be faster — it stops paying to `fprintf` 7
  unused columns, a real Θ(N²) serial cost. (That the write phase can
  dominate at some precision/N combination is itself established by the
  72%/28% experiment; only the specific split value, and its transfer to
  other precisions, is what's unsupported.)
- **Magnitude**: unknown. It is both N- and precision-dependent by the
  argument above, and cannot responsibly be stated as a percentage without a
  fresh, purpose-built measurement: same paired/interleaved design as job
  29628907, timing the bare default directly against `--metrics` (not
  reconstructed from other flags' deltas), swept across precision rather
  than read off one point — the metrics-by-default.md table above shows
  precision alone moves the answer by roughly two orders of magnitude.
  `register_equality_similarity`'s own formatting cost (one more `%.17g`
  field on `--metrics`/`--re`, absent from any measurement cited above)
  should be captured by that same fresh run rather than estimated.
- This is exactly the class of mistake `CLAUDE.md`'s own retraction section
  warns against ("a 0.7s floor was accepted as the explanation for a 1262s
  gap... check a fix against the size of the thing it claims to explain") —
  the difference here is the mismatched quantities were caught before
  landing in `CLAUDE.md`, not after.

**Verdict on Part 9: still worth doing eventually, still low priority —
defer, don't schedule dedicated compute for it. This verdict does not rest
on the withdrawn numbers above**, and holds on narrower grounds:
- It's a **harness runtime** improvement, not a **published-number**
  improvement — Figure 3 already reports only `jaccard_similarity_ie`, so
  nothing in the paper's claimed values changes; only how long it takes to
  regenerate them does. That's real but modest value (faster iteration on
  reruns), not urgent value.
- Regenerating headline-adjacent figures/CSVs carries its own risk
  independent of correctness — it invites another round of "did anything
  actually change" scrutiny right as `docs/submittability-concerns.md`
  tracks submission readiness. Check that doc's current status before
  scheduling Part 9; if the paper is close to frozen, defer past submission
  entirely.
- **Do it opportunistically instead**: `docs/seed-hammock-cpp-file-dispatch.md`
  Part 2 already needs the `--no-metrics`/N-sweep benchmarks re-run for an
  unrelated reason (the sketch-phase threading fix invalidated their prior
  numbers). Piggyback Part 9 on that re-run rather than spending a separate
  SLURM allocation solely on this seed — and use the opportunity to actually
  measure the write/compute split at the precisions that re-run covers,
  rather than guessing at it again.

### Part 2 — Python CLI: add the three-shape flag/column/filename contract

Touches: `python/hammock/cli.py` (new mutually-exclusive
`--metrics`/`--register-equality`/`--re` group, default neither),
`python/hammock/runner.py` (parameterize `similarity_measures` and the
per-row value list in both the A/B/C write path and `_write_mode_d_csv` by
the new flag; `register_equality_similarity` is a plain duplicate — no new
math, just write `jaccard`'s value twice under two names), `python/hammock/outprefix.py`
(`get_new_prefix` gains a tag parameter for the three suffixes above; thread
it through both call sites in `runner.py`).

Keep `--register-equality`/`--re` cheap: its code path must not call
`_core.pairwise_metrics_hll`'s containment machinery if that's avoidable, or
must at minimum not depend on it being run — check whether `pairwise_metrics_hll`
already always computes containments as a fused side effect on the Python
binding side (it does, per `cpp/include/hammock/hll_sketch.hpp`'s fused pass) —
if so, the *Python* CLI's register-equality arm may not be able to skip the
union cheaply without a binding change, unlike the C++ standalone binary
(Part 4), which computes it per-branch. **Establish this before assuming
parity in cost characteristics between the two front-ends** — the *column
contract* must match bit-for-bit; the *cost profile* of `--re` on the two
front-ends is a separate, secondary question and may legitimately differ (the
Python CLI already always pays the fused-pass cost for every shape today, so
`--re` not being cheaper than `--metrics` on the Python side wouldn't be a
regression — only C++ standalone's `--re` must stay cheap, since that's the
one all the timing benchmarks in Part 6/9 use).

Done when: a new `tests/test_metrics_flags.py` (write it here) exercises all
three shapes on Mode B and Mode D, asserting exact header + exact filename
suffix for each, and that `jaccard_similarity` is byte-identical across the
`--re` and `--metrics` shapes (same value, computed the same way) and
`jaccard_similarity_ie` is byte-identical across default and `--metrics`.
Pre-existing tests are *expected* to start failing at this checkpoint (they
assume the old contract) — that's Part 3's job, not this one's. Don't try to
make the full suite green yet.

### Part 3 — Fix pre-existing Python-side tests for the new default

Touches: `tests/test_parity_against_original.py`, `tests/test_containment_estimator.py`,
`tests/test_jaccard_ie.py`, `tests/test_mode_d.py`, `tests/test_mode_d_parity.py`,
`tests/test_bed2fasta_cli.py` — audit each CLI invocation with no explicit
metrics flag; if the test's actual intent is the full block (e.g. it asserts
on `containment_AB`), add `--metrics` explicitly rather than rewriting the
assertions, so the test keeps testing what it always tested. If a test's
intent was genuinely about `jaccard_similarity_ie` alone, the new bare
default may let you simplify it — judgment call, not mechanical.

Done when: `pytest tests/ --ignore=tests/test_hammock_cpp_metrics.py` is
green at (at least) the Part-2 baseline pass count, modulo tests
intentionally rewritten.

### Part 4 — `hammock-cpp`: mirror the same contract in C++

Touches: `cpp/app/hammock_cli.cpp` — replace `Args::metrics` (bool) with a
3-way `enum class MetricsMode { Ie, RegisterEquality, Full }` (default `Ie`);
rewrite arg parsing (`--metrics` → `Full`, `--register-equality`/`--re` →
`RegisterEquality`; no `--no-metrics` case — falls through to the existing
unknown-argument path, same pattern as the removed `--peak-height`, see
`test_peak_height_is_gone`); rewrite the header/stride/per-pair branch
3-way, not 2-way: `Ie` (stride 1, `cell[0] = jaccard_similarity_ie`, needs
the fused union pass + hoisted cardinalities, same as today's `--metrics`
branch minus the containment/cosketch arithmetic), `RegisterEquality`
(stride 2, `cell[0] = cell[1] = jaccard_similarity`, skips the union pass —
today's `if (!args.metrics)` branch, just duplicated into a second cell),
`Full` (stride 8, today's fused/containment path filling `cell[0..6]` plus
`cell[7] = cell[0]` for the duplicate). Rewrite `outprefix_with_suffix`'s
tag logic (3-way, see filename table above); rewrite the `--help` text
block (lines ~90-155).

Then rewrite `tests/test_hammock_cpp_metrics.py` to the new contract —
`_METRIC_COLS`, `_FULL_HEADER`/`_SLIM_HEADER` (now three headers), the
flag strings passed in `_run_cpp_path`, and every filename assertion
(`_j3` → `_ie`/`_re`/`_full`). This file is the cross-tool bit-for-bit gate
(`test_metrics_block_matches_python_bit_for_bit` and friends) — keep that
property for all three shapes, not just `--metrics`.

Build inside the worktree per `CLAUDE.md`'s compiler caveat (conda env's own
compiler, not a spack gcc module):
```bash
cd /home/jbonnie1/interval_sketch/hammock_claude_wt_metrics
CR=$CONDA_PREFIX; rm -rf build/
CC=$CR/bin/x86_64-conda-linux-gnu-gcc CXX=$CR/bin/x86_64-conda-linux-gnu-g++ \
  pip install -e . --no-build-isolation
```
Done when: `pytest tests/test_hammock_cpp_metrics.py -q` passes.

### Part 5 — Full suite green

Run `pytest tests/` end-to-end in the worktree; fix stragglers Parts 2-4
missed. Done when: full green at ≥ the baseline pass count noted in "Where
to work" above (some skips are expected/environment-gated, per `CLAUDE.md`'s
Parity environments section — don't chase those).

### Part 6 — Update dependent experiment/benchmark scripts

Per-file judgment call (see "Not mechanical," above), not a find-replace.
For each file in the inventory under `experiments/` and `paper/`:
1. If it passes `--no-metrics` to get the cheap register-equality timing
   arm — rename to `--register-equality` (or `--re`). Mechanical, safe.
2. If it passes `--metrics` — read *why*, and record it with a single
   mandated mechanism: an inline `# PART9:` comment at the call site (shell
   or Python, whichever fits) for whichever consumers only ever read
   `jaccard_similarity_ie` downstream. **Use the comment, not a commit-message
   note** — Part 9 will run cold in a context-cleared session and its first
   move is `grep -rn "PART9:"`; anything recorded only in a commit message is
   invisible to that grep and will be silently missed. **Don't act on the
   flag here** — leave `--metrics` as-is at every such call site; Part 9
   decides whether and when to retarget it, per Part 1's verdict above (low
   priority, piggyback on other reruns rather than a dedicated pass).
3. `.sh`/`.py` files with hardcoded output-filename globs (e.g. anything
   matching `*_j3.csv`, or constructing the old suffix by hand rather than
   reading it off a `--verbose` "Wrote ..." line the way
   `tests/test_hammock_cpp_metrics.py`'s `_run_cpp_path` already does) need
   the glob/construction updated to `_re`/`_full`/`_ie`. Grep each file for
   `_j3` specifically, not just the flag names.
4. **Separately, and easy to miss: audit for invocations of either
   front-end that pass *no* metrics flag at all** — not just Python-CLI
   callers. The `--no-metrics`/`--metrics` grep this Part's inventory was
   built from is structurally blind to these: the Python CLI has never had
   the flag (see "Current state," above), so every existing Python-CLI
   caller relies implicitly on today's always-full-block default and
   contains neither string; and any pre-v0.7.0 or otherwise-unflagged
   `hammock-cpp` call site relies the same way, since `--metrics` is
   currently the (already-full) default there too — Part 2 flips the
   Python CLI's implicit default and Part 4 flips `hammock-cpp`'s explicit
   one, so a bare invocation of *either* binary is equally exposed once both
   land.

   **Default policy for every no-flag call site found this way: add
   `--metrics` explicitly, to preserve exactly the columns it already
   depends on — unless the script's purpose is specifically
   benchmarking/timing, in which case apply item 1/2's logic instead** (pick
   the flag that matches what it's actually timing, which may be the new
   bare default or `--register-equality`, not necessarily `--metrics`).
   Most call sites are the former case: something downstream reads named
   columns and has no interest in which arm produced them, so the
   column-preserving default is `--metrics`, not a case-by-case guess at
   which of the 8 columns is actually used.

   **Judge this by what each script actually logs/reads downstream, not by
   its directory** — `experiments/bedtools_benchmark/` and
   `experiments/subB_mixed_stride/` each contain a mix of genuine timing
   scripts (`measure_cli_overhead.py`, `measure_threads_isolation.py`) *and*
   accuracy/analysis scripts that happen to live in the same directory:
   `experiments/bedtools_benchmark/estimator_compare.py` compares
   `jaccard_similarity` against `jaccard_similarity_ie` against `bedtools`
   truth (reads `containment_AB`/`containment_BA` as data, never logs wall
   time — needs `--metrics`, item-2 treatment, despite the directory) and
   `experiments/subB_mixed_stride/run_ie_subb.py` asks whether subsampling
   degrades `jaccard_similarity_ie` (an accuracy question, also needs
   `--metrics`, also in a directory that sounds timing-related). Directory
   is a weak prior at best; open the script and see what it does with its
   output before deciding.

   At minimum, audit and fix: `experiments/maurano_dhs_validation/run_sweep_abc.py`
   and `run_sweep_d.py` (invoked by `sbatch_sweep.sh`/`sbatch_fill_highk_w.sh`;
   their output feeds `experiments/maurano_dhs_validation/analyze.R`, which
   reads `jaccard_similarity`, `containment_AB`, `containment_BA`,
   `cosketch_geom`/`arith`/`max` — all of which vanish from the bare default
   in Part 2/4, silently breaking the R script's column lookups on the next
   sweep — none of this is a timing script, so `--metrics` per the policy
   above), `experiments/primate-phylogeny/scripts/precision_probe.sh`, and
   `experiments/modeD_flanking/run_sweep_synthetic.py` (archival, lower
   priority). These three are a starting list, not a complete one — this
   Part's original grep-based inventory can't see no-flag call sites by
   construction, so re-grep for bare `hammock`/`hammock-cpp` invocations
   across `experiments/` (and any other directory that shells out to either
   binary) rather than trusting that list as exhaustive.

Suggested split if a fresh session wants a smaller unit: 6a
`experiments/bedtools_benchmark/*`, 6b `experiments/subB_mixed_stride/*`,
6c `paper/pairwise_scaling/*.R`, 6d the no-flag Python-CLI audit (item 4).
Don't re-run any SLURM jobs as part of this Part — updating a script to
launch correctly is in scope; re-generating archived results is not (that's
Part 9's job, deliberately deferred).

Done when: `grep -rn -- '--no-metrics'` returns nothing outside
`docs/seed-*.md` (historical, left alone per "Current state" above); every
live `experiments/`/`paper/` script that reads named similarity/containment
columns (not just timings) passes `--metrics` explicitly, on whichever
front-end it calls, per item 4's default policy; every script that's
genuinely about timing carries the flag matching what it actually measures;
`grep -rn 'PART9:'` shows one marker per `--metrics` call site flagged for
Part 9.

### Part 7 — Docs: README, CLAUDE.md, live design docs

- `README.md`: CLI flags section and "Output columns" section need the new
  three-shape contract (currently describe none of this, since the Python
  CLI had no flag at all).
- `CLAUDE.md`: the "Since v0.7.0 it emits the full 9-column block by
  default" paragraph (in the top-level "Build / test" section) is now wrong
  in every particular (default, flag names, filename tags) — rewrite it,
  keep the surrounding cost-measurement paragraphs (block-cost table,
  fusion-benefit table) since those numbers are unaffected by this
  restructure — only which flag reaches which shape changed, not what any
  shape costs relative to another. Add a dated changelog-style note (this
  repo's convention, e.g. divergence #9's "SHIPPED" framing) once merged.
- `docs/metrics-by-default.md`: superseded by this restructure (it documents
  the two-shape, always-full-by-default contract as current) — add a
  retraction/superseded banner at the top pointing here, matching the
  pattern in `docs/bedtools-baseline-retraction.md`. Don't delete it; it's
  cited from `CLAUDE.md`'s Architecture section as historical per-precision
  data.
- `docs/submittability-concerns.md`: line ~54 ("Figure 3's headline runs
  `--no-metrics`, i.e. the register-equality...") needs the flag name
  updated to `--register-equality`/`--re`.
- `pyproject.toml`: bump `version` from `0.7.1` to **`0.8.0`** — a minor
  bump, not a patch. This repo's own version history draws that line
  consistently: minor bumps mark CLI-contract/default changes (v0.4.0 added
  a column, v0.5.0 added `jaccard_similarity_ie`, v0.6.0 removed the
  `_with_ends` column family, v0.7.0 changed `hammock-cpp`'s defaults —
  `git log --oneline` on each shows exactly one minor-version commit per
  behavior change), while patch bumps are reserved for fixes to already-shipped
  behavior with no contract change (v0.6.1 changed Mode D's *default thread
  count*, not its output; v0.7.1 bounded `--threads`, not its output). This
  seed changes what every invocation without special flags outputs and
  removes a flag outright — squarely the v0.7.0-class change, not the
  v0.6.1/v0.7.1-class one. Note the bump in a `CLAUDE.md` changelog-style
  entry alongside the rest of Part 7's doc updates. `--version` on both
  front-ends is load-bearing, not cosmetic: the benchmark harnesses probe it
  and refuse anything older than the version they expect, so a future
  harness gating on "≥ 0.8.0" is meaningful only if this bump actually
  happens as part of the merge, not deferred to some later commit.

Done when: `grep -rln -- '9-column\|--no-metrics'` across `README.md`,
`CLAUDE.md`, `docs/metrics-by-default.md`, `docs/submittability-concerns.md`
returns nothing describing the old contract as current; `docs/metrics-by-default.md`
carries a superseded banner; `pyproject.toml`'s version reads `0.8.0`.

### Part 8 — Merge

Once Parts 2-7 are done and `pytest tests/` is green in the worktree: rebase
`metrics-restructure` onto current `main` (it will have moved — check
`git log main` for what landed on it while this seed was in progress, e.g.
the CLI-vs-hammock-cpp seed's Part 2 re-measurement), re-run the full suite
once more post-rebase, then merge to `main` per `memory/feedback_work_on_main.md`.
Confirm the two SLURM jobs noted above (or whatever is running by then) are
still on `main` and untouched — `git worktree remove` the worktree only
after the merge lands and nothing references its path.

Done when: `metrics-restructure` is merged into `main`, `pytest tests/` is
green on `main` post-merge, and `git worktree list` no longer shows
`hammock_claude_wt_metrics`.

### Part 9 — Re-run benchmarks with the new targeted flags, update results

**Per Part 1's verdict: worth doing, low priority — piggyback on another
scheduled rerun rather than spending dedicated compute.** Part 6
deliberately left every `--metrics` call site alone and just recorded which
ones only need `jaccard_similarity_ie`. This Part is that follow-up: run
`grep -rn 'PART9:'` to find Part 6's markers, and for each flagged call site
that only ever reads `jaccard_similarity_ie` downstream, switch it to the
new bare default and **re-run** the benchmark rather than just editing the
invocation.

This changes timing numbers, so it is not safe to do silently:
- Regenerate the affected `docs/data/*.csv` files and any
  `experiments/*/results/*.csv` they were staged from, using the same
  methodology already established (paired/interleaved runs, `--exclusive`
  where the harness used it, same node-selection discipline) per
  `docs/seed-benchmark-methodology.md` — don't introduce a new confound
  while fixing this one.
- Update the corresponding `RESULTS.md` files
  (`experiments/bedtools_benchmark/RESULTS.md`,
  `experiments/subB_mixed_stride/RESULTS.md`) and any figure this feeds
  (`paper/pairwise_scaling/*.R` outputs, `paper/figures/*.png`,
  `paper/outline.md`'s numeric callouts) — grep each old number before
  overwriting it, the way `CLAUDE.md`'s retraction section models: name the
  old value, the new value, and why they differ (flag/shape change, not a
  methodology fix), so a reader doesn't mistake this for another baseline
  retraction.
- **Do not** describe the resulting speedup/timing deltas as a "fix" or
  improvement to hammock itself — the underlying computation is unchanged;
  only which columns get computed (and therefore what gets `fprintf`'d) is
  different. Frame it as "this measurement now uses the arm that matches
  what it actually needed," matching how divergence #9's `_j3` tag was
  introduced without claiming a speedup.
- If a call site's old number and new number differ by less than this
  repo's established run-to-run noise floor (±2-4% per the paired-design
  control in `docs/seed-benchmark-methodology.md`), it's fine to leave the
  archived figure as-is and just note in `RESULTS.md` that the invocation
  was retargeted with no material effect — don't force a regeneration pass
  that isn't going to change anything readers would notice. Per Part 1, the
  write phase this Part is chasing shrinks fast at higher precision (well
  under 1% of pairwise time by p=24 in the metrics-by-default.md table), so
  expect a near-zero effect at high-p call sites specifically, not just
  small-N ones.
- Part 1 deliberately withdrew its numeric estimate rather than leave a
  wrong one in place — there is no percentage to "trust" here. Run one
  small paired measurement per call site (bare default vs `--metrics`,
  mirroring job 29628907's design, at that site's actual precision) before
  deciding whether a regeneration pass is worth it.

Done when: every call site Part 6 flagged has either been retargeted +
re-measured with results updated in place, or explicitly noted as
"retargeted, no material change" per the noise-floor bullet above — no
flagged site should be silently left on `--metrics` with a stale
"TODO" comment.

## Open question for whoever picks up Part 2

Should `--register-equality`/`--re` and `--metrics` be genuinely mutually
exclusive (argparse error if both given) or should `--metrics` simply win if
both are passed? Mutually exclusive is the safer default (matches how
`--ref` vs `--ref1`/`--ref2` is already handled per `CLAUDE.md`'s bed2fasta
section — "combining X with Y is an error") and was assumed throughout this
seed, but wasn't explicitly asked. Confirm or flag before Part 2 lands if it
matters.
