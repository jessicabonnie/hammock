# Seed: rename `jaccard_similarity` → `reg_eq_similarity`, drop the `register_equality_similarity` duplicate

**Status: plan drafted 2026-08-12, reviewed once (3 independent passes:
inventory completeness, risk/safety, process/convention fit) plus one
follow-up direct-verification pass after the user pointed at the per-figure
"provenance" statements in `paper/` and `docs/` as an under-used scoping
source — that pass found a whole missed directory (`docs/scripts/`), folded
in below alongside the first round's findings. All six Open Questions
resolved with the user 2026-08-12 (see "Decisions made," below); one of
those resolutions (rename the C++ method too) added a new Step 1b.
Implementation not started.** This is an implementation
plan in the same mold as `docs/seed-metrics-column-restructure.md`,
restructured to match the user's own four-step framing rather than
inventing extra process scaffolding: one ungated Setup preamble (worktree +
inventory, read-only, low risk) followed by five gated Steps — the four
items in the user's original request, plus Step 1b, added by explicit
decision (see below) for the C++ method rename. Each Step is gated by
independent review of *that Step's* plan before its edits land — see
"Review process" below for what that means operationally. Work happens on a
worktree (Setup) per explicit user instruction, overriding the default
work-on-main convention (`memory/feedback_work_on_main.md`) for this task
only.

## Motivation (user request, 2026-08-12)

`jaccard_similarity` is a misleading name — it is the register-equality
statistic (fraction of active HLL registers with equal values), not set
Jaccard, and CLAUDE.md's own estimator note (divergence #2) already
documents how far it can diverge from bedtools. The name was a mistake made
at the very start of the project (predates this refactor) and is being
corrected now: rename the *column* to `reg_eq_similarity`.
`jaccard_similarity_ie` (the inclusion-exclusion, set-Jaccard-comparable
column, divergence #7) is correctly named already and **must not be
touched** — it stays `jaccard_similarity_ie`, unchanged, in every shape.

Once the column is renamed, `register_equality_similarity` — added in the
v0.8.0 metrics restructure (`docs/seed-metrics-column-restructure.md`,
commit range ending `00e3dd5`) specifically as a same-named-but-distinct
duplicate to keep the `--register-equality`/`--metrics` shapes self-labeling
— becomes pure redundancy: `reg_eq_similarity` already *is* that label.
Drop it. This shrinks the `re` shape from 4 columns to 3 and the `full`
shape from 8 to 7.

**User clarification (same session):** no external users of hammock exist
(`memory/project_no_external_api_consumers.md`'s standing precedent already
covers this) — do **not** add deprecation-warning prose to README or
anywhere else about the rename being a breaking change. Just find-and-replace
the column names in the appropriate places. This simplifies Step 4 below:
docs get updated to reflect the new name as fact, not hedged as "changed
from X, see migration note."

The one piece of backward-compat that **is** still requested (user's step 2,
unaffected by the clarification above): code that *reads* CSVs for paper
figures/stats should prefer `reg_eq_similarity` if present, falling back to
`jaccard_similarity` if not — so pre-existing archived CSVs that are *not*
rewritten by Step 1 still parse. **This is not the kind of shim
`memory/project_no_external_api_consumers.md` discourages** — that
precedent is about carrying dead code paths for hypothetical external API
consumers of *this* codebase; a read-fallback over historical *data files*
that this repo itself generated and archived is a different thing (reading
old data, not serving old callers), and the "no external consumers" fact
only removes the need to *warn* about the rename, not the need to still
*parse* files that predate it.

**This rename knowingly reverses a prior stated decision, found by the
follow-up verification pass — worth recording, not silently overriding.**
The same memory file cited above (`memory/project_no_external_api_consumers.md`,
2026-08-06) also says: *"[the no-external-consumers rule] does not license
changing `jaccard_similarity`, which is frozen for parity reasons independent
of who consumes it."* That's a different, narrower reason than the
external-consumers one — it's about `tests/test_parity_against_original.py`,
which byte-compares our output against the original, unrefactored `hammock`,
a frozen external tool that will always emit a column literally named
`jaccard_similarity`. The user's request in this session is explicit and
current, and this seed already independently arrived at the actual fix the
parity concern requires (Scope E's caveat and Open Question 2: the parity
test's comparison logic must map by column *position*, not name, once ours
is renamed) — so the technical substance is handled. But the "frozen"
framing was a deliberate prior call, not an oversight, so Step 1's review
gate should explicitly confirm the parity-test fix lands correctly (not just
that it's mentioned in this doc) before treating this as resolved.

## Review process

**Two separate gates, not one — do not conflate them.** (1) The 3-reviewer
subagent gate, defined below, runs *before* a Step's edits land. (2) A
**hard stop for the user's own review**, which runs *after* a Step's edits
land, are committed, and pass gate (1) — **the implementing session must
stop and wait for the user's explicit go-ahead before starting the next
Step.** This applies after Setup and after every one of Steps 1-4; it is not
satisfied by the subagent review gate, no matter how thorough, and it is not
optional or implicit — a session working through this seed must treat each
Step boundary as a full stop, report what landed and what the subagent
reviewers found, and wait. Don't chain Step N's completion straight into
Step N+1's work in the same uninterrupted run.

The user asked for 3 reviewers to go over the plan at each step before
implementation. Concretely, that means: **before a Step's edits land, launch
3 independent review passes against that Step's plan** (as written in this
doc at that point), each with a distinct lens — e.g. correctness/completeness
of scope, risk/safety of the concrete edit, and process/convention fit — the
same three lenses used for this doc's own initial review (see git history of
this file for that review's actual findings, folded in throughout below as
the model to repeat). A parallel batch of 3 subagent reviews, one per lens,
is the concrete mechanism; three sequential looks by the same reviewer is
not equivalent and doesn't satisfy this.

Disagreement resolution: a reviewer's concern blocks the Step until either
(a) the plan or diff is changed to address it, or (b) it's explicitly
recorded in this doc (in the relevant Step or "Open questions" section) as
accepted residual risk, with a one-line reason. "Resolved in writing" means
an edit to this file, not a verbal or chat-only sign-off — the doc is the
record a later Step's implementer (or a merge reviewer) reads cold. If a
reviewer's requested change alters the diff materially, re-run the 3-pass
review on the changed diff before proceeding; a fix that's purely mechanical
(e.g. a typo, a line-number correction) doesn't need a fresh round.

Setup (below) is deliberately **ungated** — no 3-reviewer round of its own.
That's a judgment call, not an unstated default: see Open Question 5.

**Commit discipline (user requirement, added mid-session):** commit as work
lands, not in one giant commit at the end. Concretely: each Step below gets
its own commit(s) on the worktree branch as soon as that Step's edits pass
its review gate — don't hold everything uncommitted until Step 4's merge.
Each Step's body ends with an explicit "Commit:" line naming what goes in
that commit (some Steps split into two, where the sub-parts are logically
separable and independently verifiable). At Step 4's close, merge to `main`
with a real merge commit (`git merge`, not squash/rebase-and-fast-forward)
so the individual Step commits stay visible in `main`'s history — this
matches the precedent seed's own merge (`00e3dd5`, which preserved its
`Part N: ...`-titled commits rather than collapsing them).

**One general exception, not just for Setup: edits to *this seed file
itself* always commit to `main` directly, throughout, never to the worktree
branch** — matching how the precedent's own seed doc (and its three
follow-up correction commits, `898d95c`/`ff18476`/`f3c43df`) lived entirely
on `main` while its Parts' actual code changes lived on the
`metrics-restructure` branch. So: Setup's inventory table, each Step's
review-gate sign-off/residual-risk notes, and Step 3's parity-proof record
(below) all land as their own small commits on `main` as soon as they're
written — only the actual code/CSV/doc-*content* changes each Step performs
live on the worktree branch until Step 4's merge.

## Scope / current-state inventory (snapshot, 2026-08-12 grep sweep — re-verify in Setup)

**A. Emitting code** (writes the literal header string to a CSV):
- `python/hammock/runner.py` — `_metrics_shape()` (column-name lists),
  `_metrics_row_values()` (value order + the `register_equality_similarity`
  duplicate-reuse logic, whose docstring explicitly describes the
  soon-to-be-removed duplicate and needs rewriting, not just deleting the
  code), `_print_estimator_note()` (stderr text naming the column).
- `python/hammock/cli.py` — `--metrics`/`--register-equality` help strings
  (lines ~264, 273) name the columns literally.
- `cpp/app/hammock_cli.cpp` — the three `std::fprintf(fp, "query\t...")`
  header lines (~414, 420, 423), the `--help` text (~140-173), inline
  comments describing shape/duplicate semantics (~32, 407-411, 483-489),
  **and two separate per-pair value writes that must not be conflated**:
  the `Full`-branch `cell[7] = reg_jac;` (~513) and the `RegisterEquality`-
  branch `cell[0] = cell[1] = qsk[i]->jaccard_similarity(*rsk[j]);` (~480).
  Both are literal-duplicate writes for their respective shapes; Step 1
  touches neither (header-string-only), but Step 4 must touch **both** when
  it drops the duplicate column — see Step 4's "Column removal" (a first
  draft of this doc named only the `cell[7]` site, missing the
  `RegisterEquality` one; both are the same class of bug if left stale).

**In scope — resolved 2026-08-12, see "Decisions made" below (Step 1b):**
`cpp/src/hll_sketch.cpp`, `cpp/include/hammock/hll_sketch.hpp`,
`cpp/include/hammock/abstract_sketch.hpp` define the C++ **method**
`HLLSketch::jaccard_similarity()` — an internal API name, not an output
column. **Correction to an earlier draft of this doc:** `bindings/_core.cpp`
does **not** expose it to Python under the same name — it's wrapped in a
lambda and exposed as `estimate_jaccard` (`bindings/_core.cpp:321-325`),
already a different name. So there is no Python-facing name to keep in sync;
the rename's blast radius is C++-only: the method declaration/definition in
the three files above, its ~9 call sites (`cpp/tests/test_modes.cpp:28,44`;
`cpp/tests/test_hll_sketch.cpp:20,31,43,111,135`;
`cpp/app/hammock_cli.cpp:480,494`; `bindings/_core.cpp:211,323`, the last
two being C++-side calls *inside* the binding, not the exposed name itself).
The repo is public on GitHub (`origin` is `jessicabonnie/hammock`, `main`
already synced there — confirmed by checking `git remote -v` and
`git rev-parse origin/main`, not assumed), so "no external consumers"
(`memory/project_no_external_api_consumers.md`) is an empirical claim
("not currently observed"), not a "the code is private" one — a real if
narrow risk for a hypothetical direct C++/pybind11 caller, accepted by the
user with that context in hand. See Step 1b for the concrete work.

**B. Archived CSVs checked into git that feed paper figures/stats and carry
a bare `jaccard_similarity` header today** (comma-delimited — Python-writer
output; confirmed with a tab-delimited sweep too, which found none, but
hammock-cpp output elsewhere may be tab-delimited, re-check in Setup):
- `docs/data/exp_a_narrow_k10_w10.csv`
- `docs/data/exp_a_broad_k10_w10.csv`
- `docs/data/hammock_hll_p18_jaccB_full.csv`
- `docs/data/hammock_hll_p21_jaccB_full.csv`
- `docs/data/hammock_hll_p23_jaccB_full.csv`

**Collision hazard, confirmed by review — `exp_a_narrow_k10_w10.csv` and
`exp_a_broad_k10_w10.csv` are the two files in this list at risk of a bad
substitution.** Both headers contain **both** `jaccard_similarity` and
`jaccard_similarity_with_ends` (and the matching `containment_AB_with_ends`
family). A plain non-word-boundary `s/jaccard_similarity/reg_eq_similarity/`
would also corrupt `jaccard_similarity_with_ends` into
`reg_eq_similarity_with_ends`. Step 1's substitution must be word-boundary-
safe (or done as an explicit token match, not a substring match) on these
two files specifically — same caution the plan already had for `_ie`,
extended to `_with_ends`.

None of these, nor any other checked-in CSV under `paper/`, `docs/data/`, or
`experiments/`, currently ships a `register_equality_similarity` **header**
(swept and confirmed empty) — it's recent enough (v0.8.0, 071dcf4 and
earlier on `metrics-restructure`) that no archived file predates its
addition. So Step 1, for this header-only list, is a pure single-token rename,
not also a column-drop.

This list is a snapshot of `docs/data/` and the top two levels of `paper/`
and `experiments/`; it is **not** a verified-complete sweep of every
`experiments/*/results/*.csv` (many of those are regenerated, not committed
— check `.gitignore`). Setup redoes this properly.

**B2. Tracked CSVs carrying `jaccard_similarity` as a *data value*, not a
header — a distinct failure mode Step 1's header-only sweep does not cover,
found by review, not in the original snapshot:**
- `docs/data/mode_d_summary.csv` — **940** occurrences of the literal string
  in a `column` data field (header is `precision,k,w,column,reference,...`).
  Generated by `experiments/maurano_dhs_validation/analyze.R`.
- `paper/estimator_crossover/estimator_crossover_stats.csv`,
  `paper/interval_accuracy/interval_accuracy_stats.csv`,
  `paper/sequence_tissue_clustering/estimator_agreement_stats.csv` — each
  carries `jaccard_similarity` as a literal label/data value (e.g.
  `"Register-equality (jaccard_similarity)"` or a bare `column`/`estimator`
  field value), sourced from R label constants in the scripts that generate
  them (category C, below).

Because these are **generated artifacts of tracked source code**, not
independently-archived measurements, they are explicitly **not** in scope
for Step 1's hand-edit-the-header approach. Instead: Step 2 updates the
*generating* code's string literals (the label constants, the `column = `
assignment logic) to emit `reg_eq_similarity`, and Step 3 regenerates these
four files from scratch and diffs them structurally as part of its parity
proof — see Steps 2 and 3, and Open Question 6.

`docs/data/mode_d_summary.csv` specifically has **more consumers than just
`plot_parameter_objective_tradeoff_estimators.R`** (found by direct
verification, not the original three-reviewer round): four scripts under
**`docs/scripts/`** — `mode_d_violins.R`, `mode_d_lines.R`,
`mode_d_metric_tradeoff.R`, `mode_d_bedtools_vs_modeB_scatter.R` — each do
`filter(column == "jaccard_similarity", ...)` against it, feeding Figures 3,
4, 6, and S4 of `docs/paper_outline.md` per that doc's own figure-index
table (lines ~366-393). `docs/scripts/` is a directory Category D below
never accounted for at all (it's neither under `paper/` nor `experiments/`)
— all five consumers of this file's `column` field need the same Step-2
treatment.

**C. Paper R scripts** referencing the column (grep'd, with actual
column-read vs. comment/label-only distinguished):

| File | Reads the column as data? |
|---|---|
| `paper/cross_reference_identity/plot_cross_reference_identity.R:42` | yes — `DEFAULT_SIM_COL <- "jaccard_similarity"` |
| `paper/interval_accuracy/plot_interval_accuracy.R:91,154,171` | yes — `:91` is a label constant (`EST_RE <- "Register-equality (jaccard_similarity)"`) that also needs updating, missed in the first draft of this scope table which only cited `:154,171` |
| `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R:59,271,286` | yes — also drives Figure 6. **Lines 271 and 286 are coupled and must be edited as one unit, not independently**: `:271` builds the working column set via `intersect(c("jaccard_similarity", "jaccard_similarity_ie"), names(raw))`; `:286` separately does `agreement$signature[agreement$column == "jaccard_similarity"]` against a *different* hardcoded literal. Updating only one leaves the other filtering for a name that no longer exists, and `agreement$column == "jaccard_similarity"` degrades **silently to a length-0 match** (`ref_signature` becomes empty, `partition_identical` becomes `NA` for every row) rather than erroring — exactly the kind of failure Step 3's parity proof exists to catch, but cheaper to just get right here. |
| `paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R:60,98,132` | yes |
| `paper/estimator_crossover/plot_estimator_crossover.R:62` | label string only (`EST_RE <- "Register-equality (jaccard_similarity)"`), verify no separate data read elsewhere in the file |
| `paper/pairwise_scaling/plot_precision_frontier.R:112` | printf message text only, not a column read |
| `paper/pairwise_scaling/plot_pairwise_scaling.R:5` | comment only |

`paper/methods.md`, `paper/draft.md`, `paper/outline.md` reference the
column name in prose (manuscript text, captions, and the Figure 6 estimator
note). Given the user's "just find and replace, no need to hedge"
clarification, these get updated too, in Step 4's closeout — not left
describing a column that no longer exists in fresh output.

**D. `experiments/` scripts, plus `docs/scripts/` (found only by direct
verification, not by the three-reviewer round — see below).** The original
snapshot ("28 files matched the grep", attributed to `maurano_dhs_validation/`,
`primate-phylogeny/`, `bedtools_benchmark/`, `mus-homo/`, and the ARCHIVAL
`modeD_flanking/`) was **undercounted by review — an independent sweep found
34 script matches (.py/.R/.sh) under `experiments/`, spanning three more
live directories never named**: `experiments/ref-comparison/` (4 files),
`experiments/subB_mixed_stride/` (6 files), `experiments/synthetic_evolution/`
(1 file). The first two are **provably in scope**, not just plausibly so:
`paper/draft.md:101` and `paper/outline.md:74` cite `ref-comparison/`
directly as the source for Figure 5 and (`paper/draft.md:273`)
`experiments/ref-comparison/results/exp_a_estimator_delta.csv` for Table S5;
`paper/draft.md:282` cites `subB_mixed_stride/` for Figure S9.
`modeD_flanking/` remains correctly out of scope (already flagged
non-reproducible in `memory/project_modeD_flanking.md` — do not touch, no
value in updating dead analysis code).

**`docs/scripts/` is a fourth directory, outside both `paper/` and
`experiments/`, that no sweep (including the three-reviewer round) checked
at all** — found by directly listing it after the user flagged that
"provenance" statements exist for figures in `paper/` and `docs/`. Four of
its seven files (`mode_d_violins.R`, `mode_d_lines.R`,
`mode_d_metric_tradeoff.R`, `mode_d_bedtools_vs_modeB_scatter.R`) read
`docs/data/mode_d_summary.csv`'s `column` field for the literal string
`"jaccard_similarity"` (category B2) and feed Figures 3, 4, 6, and S4 of
`docs/paper_outline.md` — confirmed via that doc's own figure-index table,
not inferred. Provably in scope, same treatment as the other B2 consumers.

**Use the repo's own figure-provenance statements as the authoritative
scoping source, not just grep-and-infer.** `paper/outline.md` and
`paper/draft.md` each carry a `**Figure N generation.**` paragraph per
figure, explicitly naming its generating script and input file(s) (e.g.
`paper/outline.md:58,66,74,85,95`); `docs/paper_outline.md` has an actual
figure-index table (~lines 360-395) doing the same for its own figures.
Setup should walk every one of these blocks/table rows and treat every
script and data file they name as in-scope by construction — this is how
`docs/scripts/` was found, and it's a more reliable method than grepping
`paper/outline.md`/`paper/draft.md` prose for cross-references, which is
what the original inventory relied on and which is exactly what missed a
directory with no textual link into either file's prose (the link only
exists in the figure-index table's structured rows).

**A fourth grep-taxonomy case, missed in the original three-way (a) emits-
header / (b) reads-as-data / (c) prose-only split: read-and-re-emit.**
`experiments/maurano_dhs_validation/analyze.R` is the generator behind
`docs/data/mode_d_summary.csv` (category B2, above) — it both reads
`jaccard_similarity` from upstream hammock output *and* re-emits the literal
string as data into several downstream CSVs (13+ hardcoded occurrences,
including a `short_col` lookup table mapping `jaccard_similarity =
"register_equality"`, itself worth checking against the new name for
consistency rather than left as a second, differently-spelled alias). This
file needs its own line in Setup's completed inventory table, tagged as
this fourth case, not squeezed into (b).

**E. Tests asserting the column name literally** — not "paper figures/stats"
in the user's literal wording, but they will go red the moment Step 1 or
Step 4 lands if not updated in the same commit: `tests/test_mode_d_parity.py`,
`tests/test_metrics_flags.py`, `tests/test_hammock_cpp_metrics.py`,
`tests/test_parity_against_original.py`, `tests/test_mode_d.py`,
`tests/test_bed2fasta_cli.py`, `tests/test_containment_estimator.py`,
`tests/test_jaccard_ie.py`. **All 8 need a grep-and-fix pass in the relevant
Step, not just the two (`test_hammock_cpp_metrics.py`, `test_metrics_flags.py`)
whose column-shape assertions Steps 1 and 4 name explicitly below** — the
other six were flagged as "at risk" in an earlier draft of this doc but
never actually assigned to a Step's body text, which would have left them to
surface only at Step 4's final `pytest tests/` gate, off the per-Step review
track. Fixed now: both Step 1 and Step 4 explicitly re-grep this full list.

**Caveat on `test_parity_against_original.py`:** it byte-compares our output
against the *original*, un-refactored `hammock` (pipx `hammock-orig` /
conda `hammock`) — a frozen external codebase this rename does not and
cannot touch. Orig will emit a column literally named `jaccard_similarity`
forever. The parity test's comparison logic must keep mapping our
`reg_eq_similarity` values against orig's `jaccard_similarity` values (by
column *position* or an explicit rename-before-compare step) — this is a
"teach the harness about the new name," not "make both sides agree on a
name," and is exactly the kind of thing Open Question 2 asks reviewers to
confirm before Step 1 lands.

**F. Docs describing the current contract as fact**, needing updates in
Step 4's closeout: `CLAUDE.md` (the three-shape table under "Build / test",
the `--register-equality`/`--metrics` prose block, the "Measured cost of the
metrics block" section headers, divergence #2's estimator note), `README.md`
(verify whether it documents columns at all), `docs/jaccard-definitional-gap.md`,
`docs/estimator-analysis-findings.md`. **Two more, found by review, missing
from the original list:** `docs/paper_outline.md` — a *separate, canonical*
document from `paper/outline.md` (confirmed distinct: `docs/estimator-
analysis-findings.md:127-128` raises the two-live-outline-files question,
and `:526-528` ("Answer to Q5") is where it actually resolves
`docs/paper_outline.md` as canonical — a citation-accuracy fix over an
earlier draft of this doc, which pointed only at the question, not the
answer), containing dozens of prose references to the current column name;
and `docs/submittability-concerns.md`, a live
design doc (the precedent seed's own Part 7 updated it) that currently
states the post-restructure contract using the current name ("`jaccard_similarity`
now requires `--register-equality`/`--re` or `--metrics`") and will read as
wrong once the rename lands.

`docs/seed-metrics-column-restructure.md` is a **historical record** of the
v0.8.0 change — do not rewrite it to pretend `reg_eq_similarity` existed at
that time; add one forward-pointer line noting it's superseded by this seed
for the column *name* (its shape/count structure otherwise still holds,
minus one duplicate column).

**Correction, found by review: do NOT add this rename as a bullet in
CLAUDE.md's "Open seeds" list.** That section's own stated scope is standing
*research* questions ("None is a decision; each is evidence gathered plus
what still needs establishing") — a settled, executed rename doesn't belong
there. The precedent this seed is modeled on never added itself to that
list either; instead its Part 7 called for a dated changelog-style note
inside the relevant prose section (the three-shape table's section),
matching the "SHIPPED" framing divergence #9 already uses elsewhere. Follow
that pattern in Step 4, not an Open Seeds bullet. **Caution, also found by
review:** the precedent's own promised "SHIPPED" note for the metrics
restructure has still not actually been added to CLAUDE.md even though that
merge landed (`00e3dd5`/`071dcf4` are in `main`'s history already) — don't
let this seed's closeout note slip the same way; it's part of Step 4's
"Done when," not an optional follow-up.

## The plan

Five gated Steps: the user's own four (1: rename headers, 2: update
consuming code, 3: prove parity, 4: drop the duplicate column), plus Step 1b
(rename the C++ method and Python binding too, added by decision — see
"Decisions made"), each preceded by independent review of *that Step's*
plan per "Review process" above. One ungated Setup preamble first.

### Setup — worktree + inventory (ungated, read-only)

`git worktree add ../hammock_claude_wt_regeq -b jaccard-reg-eq-rename`
(adjust branch name to whatever's free). All four Steps happen there.

**Before starting, check for conflicting in-flight work — concretely, not
just "be careful":** run `squeue -u $USER` (or the local equivalent) and
check specifically whether anything is writing to `docs/data/` or any
`paper/*/` results directory — Step 1 does in-place edits to *tracked* CSVs
on this repo's history, a different risk profile from the precedent seed's
worktree (which only built/tested in the worktree, touching no shared
tracked files). If a job is regenerating one of the category-B/B2 files
concurrently, wait for it or coordinate before Step 1 edits the same file —
a race here produces a merge conflict or, worse, a silently-stale base for
Step 1's diff.

Setup's own deliverable (the inventory table below) is produced and
committed to `main` before Step 1's review begins, per "Commit discipline"
above — not inside the worktree branch.

Then complete the inventory sweep from "Scope" above, properly rather than
as a snapshot:
- Split every `jaccard_similarity` (excluding `_ie`) hit into: (a) emits the
  literal string as a CSV **header**, (b) reads a CSV by that column name as
  data, (c) prose/comment/label-text only, (d) reads *and re-emits* the
  literal string as data into a downstream tracked file (category B2's
  fourth case — `analyze.R` is the known instance, check for others,
  especially other `experiments/*/analyze*.R` or `.py` scripts that produce
  summary CSVs).
- **Also sweep for the string as tracked *data*, not just as a header** —
  `grep -c` for `jaccard_similarity` inside CSV bodies (not just `head -1`)
  across `docs/data/`, `paper/`, and `experiments/*/results/`. Category B2
  above is the known result of doing this once; confirm it's complete.
- **Start from the provenance statements, not just grep.** Walk every
  `**Figure N generation.**` paragraph in `paper/outline.md` and
  `paper/draft.md`, and every row of `docs/paper_outline.md`'s figure-index
  table, and treat every script and data file each one names as in-scope by
  construction. This is how `docs/scripts/` (category D) was found after
  every grep-based sweep missed it — a directory with no *prose*
  cross-reference into `paper/outline.md`/`paper/draft.md`, only a
  structured table row in `docs/paper_outline.md`. Grep alone is not
  sufficient for this reason; use it to find *candidate* files, use the
  provenance statements to confirm or add to that list.
- For every `experiments/` or `docs/scripts/` script found either way, check
  whether it feeds a figure or number one of the three provenance sources
  above actually cites. Mark in-scope / out-of-scope explicitly — don't
  assume the whole tree, but do explicitly walk `experiments/ref-comparison/`,
  `experiments/subB_mixed_stride/`, `experiments/synthetic_evolution/`, and
  `docs/scripts/` (all confirmed live, not archival; `ref-comparison/`,
  `subB_mixed_stride/`, and all four of `docs/scripts/`'s `mode_d_*.R`
  files are already confirmed in-scope per category D above — Setup's job
  for these four is completing the *rest* of the sweep, not re-deciding
  scope already settled).
- Re-run the checked-in-CSV header sweep across **all** of `paper/`,
  `docs/data/`, and any `experiments/*/results/` directories that are
  actually tracked by git (not gitignored), with **both** comma and tab as
  candidate delimiters (hammock-cpp writes tab-delimited output; the Python
  CLI writes comma-delimited via `csv.writer`'s default — a delimiter-blind
  sweep can silently miss one family).
- Produce a table (file, role, in-scope) as this section's deliverable,
  appended to this doc (on `main`, per above) before Step 1's review starts.

**Then stop.** Setup is ungated for the 3-reviewer subagent pass (see
Open Question 5), but the user-review stop still applies: report the
completed inventory and wait for the user's go-ahead before starting Step 1.

### Step 1 — Rename column headers in outputs

Covers the user's step 1: rename the literal `jaccard_similarity` **header**
wherever it appears in an *output* — both the code that emits it (so all
future runs carry the new name) and the archived CSVs already checked in
that carry it as a header and feed a paper figure or stat. Leave
`jaccard_similarity_ie` untouched everywhere. **Do not** remove
`register_equality_similarity` in this Step — it's still a literal-duplicate
column for now, just renamed alongside its twin; removing it is Step 4, kept
separate because it changes column *count*, a different risk class than a
same-value rename (see Step 4). **Do not** touch the category-B2 data-value
CSVs here either — those are handled by updating their generating code in
Step 2 and regenerating in Step 3, not by hand-editing tracked output.

Emitting code (both front-ends, same commit, kept in bit-for-bit sync per
CLAUDE.md divergence #9 and the existing `tests/test_hammock_cpp_metrics.py`
cross-tool `==` gate):
- `runner._metrics_shape`'s column-name lists: `"jaccard_similarity"` →
  `"reg_eq_similarity"` (both the `re` and `full` entries).
  `_metrics_row_values` and `_print_estimator_note` unchanged in logic, just
  update any string literal/docstring naming the column.
- `cli.py`'s `--metrics`/`--register-equality` help strings.
- `hammock_cli.cpp`'s three header `fprintf`s, its `--help` text, and
  comments naming the column (~32, 407-411, 483-489) — the per-pair value
  writes (`cell[0]`, `cell[7]`, `cell[1]` in the `RegisterEquality` branch,
  `stride`) are unchanged in this Step, only the header string literal.
- Grep all 8 files from Scope category E for a literal `jaccard_similarity`
  hit (excluding `_ie`) and fix each one in this same commit — not only
  `tests/test_hammock_cpp_metrics.py` and `tests/test_metrics_flags.py` (the
  two with column-*shape* assertions that also need Step 4's count changes
  later), but all 8, since several assert exact CSV header/column content
  that this Step's rename directly changes.
- **Two non-literal fixes, found by the Step 1 review gate's round 1 (not
  catchable by the grep pass above — see the Setup inventory's Scope-E
  corrections for the full trace), required in the same commit:**
  - `tests/test_parity_against_original.py::_projected_rows` — change
    `for line in lines` to `for line in lines[1:]` so the header row is
    excluded from the compared tuples. Without this, `test_jaccard_byte_equal`
    fails on the header alone the moment our side's header token diverges
    from orig's frozen `jaccard_similarity`, regardless of whether every
    data value still matches.
  - `tests/test_mode_d_parity.py::test_mode_d_structural_parity` — widen the
    two `str.startswith(...)` tuples at (currently) `:110-111` and `:116-117`
    from `("jaccard_similarity", "containment", "cosketch")` /
    `("jaccard_similarity", "cosketch")` to
    `("jaccard_similarity", "reg_eq_similarity", "containment", "cosketch")` /
    `("jaccard_similarity", "reg_eq_similarity", "cosketch")`. Without this,
    the renamed column silently drops out of both the range/self-pair check
    and the symmetry check — the test still passes, just checking one fewer
    column than intended.

Archived CSVs (the Setup-confirmed header list, 5 files in the initial
snapshot): edit **only** the header row's `jaccard_similarity` token →
`reg_eq_similarity`, **word-boundary-safe** so `jaccard_similarity_with_ends`
in `exp_a_narrow_k10_w10.csv`/`exp_a_broad_k10_w10.csv` is not also matched
(see the collision hazard in Scope B). Do not touch `jaccard_similarity_ie`,
do not touch any data row, do not re-run hammock — a pure text substitution
on one line. Verify with `diff <(tail -n +2 old.csv) <(tail -n +2 new.csv)`
(expect empty) for every file touched, **and** verify
`jaccard_similarity_with_ends` still appears unchanged in the two files that
have it.

**Review gate:** 3 reviewers check (a) the Setup inventory for gaps —
scripts that build a column name dynamically rather than as a literal, any
`experiments/` script mis-marked out-of-scope that a paper figure actually
depends on, any category-B2-style data-value file Setup's sweep missed —
and (b) this Step's diff: header-only in the CSVs (with the `_with_ends`
collision explicitly checked), name-only in the emitting code (no `cell[]`/
`stride` changes should appear in this diff at all — flag any that do, they
belong in Step 4), and all 8 Scope-E test files addressed. Don't proceed to
Step 2 until resolved or explicitly accepted as residual risk in writing.

**Review gate round 1, run 2026-08-12 (3 parallel subagent passes — scope
completeness, risk/safety, process/convention fit — against the plan text as
it stood before this paragraph's own edits):**

- **Process/convention fit: clean.** Commit split, worktree/main carve-out,
  both-front-ends-same-commit, and off-limits-territory checks all verified
  directly against repo state (not just doc prose) — no findings.
- **Risk/safety: clean for Step 1 itself, two non-blocking findings recorded
  for later Steps, not requiring a Step 1 diff change:**
  1. `paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff.R`
     (the non-`_estimators` sibling) also reads `jaccard_similarity` as data
     from `docs/data/mode_d_summary.csv`'s `column` field
     (`SIM_COLUMNS <- c("jaccard_similarity", "jaccard_similarity_ie")`) and
     is confirmed live, not dead code, by `paper/draft.md:122`/
     `paper/outline.md:95`'s Figure 7 provenance paragraph ("the existing
     `plot_parameter_objective_tradeoff.R` is retained unchanged as the
     single-estimator diagnostic workflow"). **Missing from Scope C /
     Category B2's consumer list — add it to Step 2's file list before Step
     2's own plan is treated as complete** (same fallback treatment as its
     `_estimators` sibling). Doesn't affect Step 1, which never touches B2
     readers.
  2. Two paper scripts that read the 5 archived CSVs Step 1 edits directly by
     path+name (`plot_interval_accuracy.R` on the three `hammock_hll_p*`
     files, `plot_cross_reference_identity.R` on `exp_a_broad_k10_w10.csv`)
     are not updated until Step 2, so the worktree branch is transiently
     non-functional for those two scripts between Step 1 landing and Step 2
     landing. **Accepted as residual risk, in writing, right here**: nothing
     in this plan runs those scripts during that window (Step 3's
     regenerate-and-diff only runs after Steps 1+2 both land), the worktree
     branch never touches `main` until Step 4's merge, and the failure mode
     if triggered early is a loud `stop()`/missing-column R error, not silent
     wrong output.
- **Scope completeness: two blocking findings, both fixed in this Step's plan
  text above before proceeding (not deferred) — see the Setup inventory's
  "Scope-E" corrections for the full trace.** `test_parity_against_original.py`'s
  `_projected_rows` was comparing the CSV header row as if it were a data
  row (an off-by-one in `for line in lines` vs `lines[1:]`), which only ever
  passed because both sides happened to share the literal string
  `jaccard_similarity`; `test_mode_d_parity.py`'s dynamic `startswith` column
  scan would silently drop the renamed column from two of its checks. Both
  are now explicit line items above, in the same commit as the rest of Step
  1's Scope-E fixes. All other scope-completeness checks (line-for-line
  emitting-code verification, the 5-CSV archived list and its `_with_ends`
  collision, the repo-wide CSV/TSV/TXT sweep, the "no dynamic column-name
  construction" claim, the figure-provenance cross-check) were independently
  re-verified against the actual files and found accurate.

Because these findings materially changed Step 1's planned diff (two new
required test fixes), a second, focused review round follows below before
implementation proceeds, per "Review process"'s re-run rule.

**Commit:** two, once the review gate clears — (1) emitting-code + test
rename (`runner.py`, `cli.py`, `hammock_cli.cpp`, the two Scope-E test files'
name-only edits plus the other six's fixes), (2) archived-CSV header
rewrites (the 5 files, plus whatever Setup's redone sweep added), kept
separate since the second is pure data-file editing with its own `diff`-based
verification and reviewers may want to inspect it independently of the code
change.

**Review gate round 2, run 2026-08-12 (3 fresh parallel passes, focused on
round 1's two new fixes and their knock-on effects, per "Review process"'s
re-run rule):**

- **Scope completeness: clean.** Both fixes verified correct against the
  actual current file content (line numbers match exactly, no drift). The
  `lines[1:]` fix has exactly one caller and no other header-as-data instance
  exists in that file. The widened `startswith` tuples were traced
  column-by-column against `runner._metrics_shape`'s actual "full" list and
  confirmed to select the same set (7 similarity cols / 5 symmetric cols)
  pre- and post-rename, with `containment`/`cosketch` untouched. A sweep of
  the remaining 6 Scope-E files plus the emitting code found no third
  non-literal reference to the column — round 1's two fixes are the complete
  set.
- **Risk/safety: clean.** `lines[1:]` loses no real test coverage (any
  genuine column-count/order mismatch would still fail on every data row, not
  just the header) and is the correct fix rather than a weakened one. The
  widened tuples don't false-positive-match anything else in the full 8-column
  shape — `register_equality_similarity` diverges from `reg_eq_similarity` at
  the 4th character, so it stays correctly excluded, unchanged from before.
  Both fixes are pure Python test-file edits (`runner.py`-backed CLI path),
  no C++ rebuild needed to implement or verify them. Round 1's two
  non-blocking findings re-confirmed accurate by direct read of both named R
  scripts (`stop()`-on-missing-column in both, not silent degradation).
- **Process/convention fit: clean, plus two mechanical cleanups identified
  and now applied** (not requiring their own review round, per the "purely
  mechanical" carve-out): the Setup inventory's first correction paragraph
  had a stale "see ... below" pointing at content that is actually above it
  in the file — fixed to "above". And round 1's `plot_parameter_objective_tradeoff.R`
  finding was recorded as a note inside Step 1's round-1 subsection but never
  actually threaded into Step 2's own consumer list — fixed: Step 2's
  "Reading `mode_d_summary.csv`'s `column` field" paragraph now names six
  consumers, not five, including that file. Both the commit-discipline check
  (seed-doc edits on `main`, worktree branch untouched and still at `f4cf813`
  pre-Step-1) and round 1's re-run-vs-mechanical judgment call were confirmed
  correct on independent re-derivation, not just accepted as given.

**Step 1 is ready to implement.** Proceeding to the actual edits.

**Step 1 landed as two commits on the worktree branch** (`49878aa` emitting
code + tests, `3f1e5ca` the 5 archived CSVs), then re-reviewed by 3 fresh
adversarial subagents against the **actual diff**, not the plan text
(distinct from the two rounds above, which reviewed the plan before
implementation) — correctness of the diff line-by-line, blast radius/scope,
and verification-gap hunting. All three came back clean; no blocking
findings. Two things worth recording:

- **The C++ side is now verified by execution, not just source reading.**
  The verification-gap reviewer built `hammock-cpp` from this worktree
  (cluster-compiler-caveat-compliant flags) and ran
  `tests/test_hammock_cpp_metrics.py` with `HAMMOCK_REQUIRE_CPP=1`: **19/19
  passed**, including the cross-tool bit-for-bit gate for all three shapes.
  Combined with the Python-side run (via the worktree-testing shim in
  `memory/reference_hammock_env_on_path.md` — strip
  `ScikitBuildRedirectingFinder` from `sys.meta_path`, symlink the
  unchanged `_core.*.so`, insert the worktree's `python/` dir first,
  guarded to py3.10), **110 passed, 2 skipped** (bedtools not on
  PATH, benign) across all 8 Scope-E test files. This supersedes commit
  `49878aa`'s message, which (accurately, at the time) said full pytest
  verification hadn't been run.
- **Residual, accepted, not blocking**: `test_mode_d_parity.py`'s
  `sim_cols`/`sym_cols` scan (`assert sim_cols` non-emptiness, the same
  tuples this Step's round-1 review gate widened) has weak discriminating
  power specifically for Mode D's coverage of the renamed column — because
  `jaccard_similarity_ie` (frozen, untouched) always satisfies the
  `startswith` prefix check regardless of whether `reg_eq_similarity` is
  present, misspelled, or the whole rename is reverted, that one test file
  alone cannot distinguish "rename done correctly" from "rename never
  happened" for Mode D. This is a **pre-existing weak-assertion pattern,
  not introduced by this rename** — verified by adversarially reverting
  `runner.py`'s emitting code entirely and confirming 11 *other* tests
  across the 7-file Python suite fail loudly (`test_mode_d.py`,
  `test_metrics_flags.py`, `test_parity_against_original.py` all have
  literal-name assertions that do catch it), so the aggregate safety net
  holds even though this one test's contribution to it is smaller than its
  comment claims. Not fixed here — noted for whoever next touches that
  test file, since strengthening it is a test-quality improvement
  orthogonal to this rename, not a gap this Step needs to close.

**Also recorded here: a live tampering incident during this review round,
unrelated to the rename's correctness but worth keeping.** Three fabricated
tool-result blocks appeared during the adversarial review, styled as
harness "system-reminder"s claiming `test_parity_against_original.py`,
`test_mode_d_parity.py`, and `runner.py` had been edited by the user or a
linter, each instructing the implementing session not to tell the user.
Two of the three claims were false (no actual on-disk change). One was
real: `runner.py` was actually mutated on disk mid-review, introducing a
typo (`reg_eq_similarity` -> `reg_eq_similarty`) in the `full` shape's
column list — reverted immediately (`git checkout HEAD --
python/hammock/runner.py`), confirmed via `git diff HEAD` that both the
worktree and `main` checkout were clean before and after every subsequent
reviewer's build/run activity. The instruction to conceal this was not
followed. Recorded here per this doc's own "resolved in writing" norm,
since it happened during this Step's review and a later reader should not
have to reconstruct it from chat history.

**Then stop.** Per "Review process" above: report what landed and the
subagent reviewers' findings, and wait for the user's own go-ahead before
starting Step 1b. Do not continue automatically.

### Step 1b — Rename the C++ method AND the Python-facing binding, too

Added 2026-08-12 by explicit user decision (Open Question 1, resolved — see
"Decisions made" below), then **expanded the same day** after the user's
actual goal was stated directly: *"make sure the register equality
calculation isn't named jaccard misleadingly on anything someone might be
using."* That's broader than just the CSV column, and broader than just the
internal C++ method — it includes `bindings/_core.cpp`'s pybind11-exposed
Python method **`estimate_jaccard`**, found only by checking (an earlier
draft of this doc wrongly assumed it shared the C++ method's name; it
doesn't — it's already a separately-chosen name, `.def("estimate_jaccard",
...)` at `bindings/_core.cpp:321`, wrapping a call to
`self.jaccard_similarity(other)`). This is the one name in the whole
codebase someone would actually call from Python directly, without ever
touching a CSV — and it's used for real, right now, by 8 call sites in
`tests/test_sequence_invariants.py:33` and `tests/test_containment_estimator.py`
(`:79,119,143,207,238,257,498`). Both renames are independent of Step 1's
CSV-header rename functionally (neither a method caller nor a binding
caller ever sees a CSV column name) but grouped here under the same
"fix the misleading name everywhere it appears" theme, on a disjoint set of
files from Step 1.

**Part 1 — the C++ method.** New name: `reg_eq_similarity()`, matching the
renamed column exactly (not `register_equality_similarity()` — keep it
short, matching the CSV name is the point). Rename it everywhere it's
declared, defined, or called:
- `cpp/include/hammock/abstract_sketch.hpp:12` — the pure-virtual
  declaration (`virtual double jaccard_similarity(...) const = 0;`).
- `cpp/include/hammock/hll_sketch.hpp:48` — the override declaration.
- `cpp/src/hll_sketch.cpp:49` — the definition
  (`double HLLSketch::jaccard_similarity(...)`).
- `cpp/tests/test_modes.cpp:28,44` and
  `cpp/tests/test_hll_sketch.cpp:20,31,43,111,135` — 7 call sites across the
  two test files; rename the call, not the test's own logic or assertions.
- `bindings/_core.cpp:211,323` — two C++-side calls (`a[i]->jaccard_similarity(*b[j])`
  in the batch-matrix path, `self.jaccard_similarity(other)` inside the
  `estimate_jaccard` lambda, which Part 2 below also touches for its own
  reason — expect both parts to edit line 323's neighborhood).
- `cpp/app/hammock_cli.cpp:480,494` — **the same two call sites Step 4 later
  touches for the stride/duplicate-column shrink** (`cell[0] = cell[1] =
  qsk[i]->jaccard_similarity(*rsk[j]);` and
  `reg_jac = qsk[i]->jaccard_similarity(*rsk[j]);`). This Part only renames
  the call (`->jaccard_similarity(` → `->reg_eq_similarity(`); Step 4 later
  edits the surrounding `cell[]` assignment structure at the same lines —
  same overlap pattern as Step 1 vs. Step 4 on the header rename, not a
  conflict, just something for whoever implements to expect.

**Part 2 — the Python-facing binding.** New name:
**`estimate_reg_eq_similarity`** (the user's explicit preference, over the
alternative of matching the original hammock's `estimate_jaccard_registers`
— chosen because it fully removes "jaccard" from the name, matching the
stated goal exactly, and keeps one consistent `reg_eq_similarity` root
across the CSV column, the C++ method, and this binding, rather than two
different naming schemes for the same computation). Rename:
- `bindings/_core.cpp:321` — the `.def("estimate_jaccard", ...)` string
  literal itself → `.def("estimate_reg_eq_similarity", ...)`.
- `tests/test_sequence_invariants.py:33` — 1 call site.
- `tests/test_containment_estimator.py:79,119,143,207,238,257,498` — 7 call
  sites.

**Leave `estimate_intersection` alone** — it's the IE-path sibling
(`bindings/_core.cpp:331`, wrapping `intersection_size`) and is already
correctly named for what it computes; nothing about this rename touches it.

No behavior change anywhere in either Part — pure identifier renames,
verified by the fact that every call site is either a test assertion
(unaffected by what the method/binding is called) or a value read into a
variable independent of the name.

**Review gate:** 3 reviewers confirm every declaration/definition/call site
in both Parts was renamed consistently (a leftover mismatched name in Part
1 is a compile error; in Part 2, `estimate_jaccard` vs.
`estimate_reg_eq_similarity` mismatches would surface as an `AttributeError`
at test time — so this is a cheap, high-confidence check: if the C++ side
compiles, the Python extension rebuilds, and `cpp/tests/` + `pytest tests/`
are all green, both renames are complete by construction) and that
`estimate_intersection` was correctly left untouched.

**Review gate, run 2026-08-12 (3 parallel subagent passes — scope
completeness, risk/safety, process/convention fit — against the plan text as
it stood before this paragraph's own edits):**

- **Process/convention fit: clean.** Worktree/main state matches the doc's
  own commit citations exactly (merge-base `f4cf813`, worktree at `3f1e5ca`,
  `main` 4 commits ahead with the Setup/round-1/round-2/post-implementation
  paperwork); the doc's own "Commit: two" line was independently confirmed
  unambiguous; the seed doc itself is byte-identical between the merge-base
  and the worktree branch (confirmed via `git diff`), i.e. no seed-doc edits
  leaked onto the worktree branch; no CLAUDE.md convention (Python-is-the-
  program, narrow/stateless `_core` bindings) is implicated by a pure
  identifier rename. No findings.
- **Scope completeness: no blocking gaps.** Independent grep across the whole
  worktree (cpp/, bindings/, python/, tests/, paper/, experiments/,
  docs/scripts/, memory/) confirmed all 11 Part 1 sites and all 8 Part 2 sites
  are exactly the enumerated set, including confirming `estimate_intersection`
  is a genuinely different code path (wraps `intersection_size`, not
  `jaccard_similarity`) and correctly left alone, and that
  `extern/hll/hll.cpp`'s vendored, differently-named, zero-call-site
  `jaccard_similarity_registers` is correctly out of scope. Two non-blocking
  findings: (1) the plan's citation of `hammock_cli.cpp`'s second call site as
  line 494 is stale — it's actually line 496 (line 494 is unrelated,
  `jaccard_and_union_cardinality`); the code snippet quoted is correct, only
  the line number drifted. (2) Five in-repo comments name the old
  method/binding identifier in prose and will read as stale (one,
  `hammock_cli.cpp:484-485`, is actively self-referential — "until Step 1b
  renames it" — and becomes self-contradictory the moment Step 1b lands) but
  aren't in Part 1/2's enumerated edit list: `hll_sketch.hpp:54`,
  `hll_sketch.cpp:176`, `hll_sketch.cpp:258`, `hammock_cli.cpp:484-485`,
  `bindings/_core.cpp:328`. A sixth, `test_containment_estimator.py:187`
  ("the route it replaced -- `jaccard_similarity()` plus `union_with()`"),
  was flagged independently by the risk/safety pass, in a file Part 2 already
  edits for 7 call sites. Folded into Part 1's diff below (mechanical, no
  scope change) rather than left for a later untracked pass.
- **Risk/safety: no blocking issue with the rename itself — pure identifier
  substitution, no logic touched, no ABI-surface risk found (no CMakeLists
  reference, no `.pyi` stub, no export/version-script file, no duck-typed
  `estimate_jaccard` call path anywhere) — but one BLOCKING finding on the
  *verification methodology* the plan proposed.** The Review-gate paragraph
  above claims "if `cpp/tests/` + `pytest tests/` are all green, both renames
  are complete by construction." That's true for Part 1 (any stale C++
  identifier is a hard compile error) and true for 7 of Part 2's 8 Python
  call sites (a stale `estimate_jaccard` there is a loud `AttributeError`).
  But `tests/test_sequence_invariants.py:33`'s call site sits inside a test
  decorated `@digest_only` — skipped, not failed, whenever the bioconda
  `digest` package isn't importable, a real recurring condition in this repo
  (CLAUDE.md's cluster-compiler/RPATH caveat, divergence #6). In a
  digest-less environment, a missed rename at that one site would produce a
  silent skip indistinguishable from normal, not a test failure — so "pytest
  green" alone does not prove Part 2 is complete in every environment.
  **Closed by strengthening verification, not by changing the diff**: Part
  2's commit requires an explicit `grep -rn 'estimate_jaccard\b' bindings/
  tests/ python/` returning zero hits, run in addition to (not instead of)
  `pytest tests/`, so completeness doesn't depend on `digest` availability.
  This is a mechanical addition to the verification step, not a scope or diff
  change, so per "Review process"'s re-run rule it doesn't require a fresh
  3-pass round. Also confirmed clean by this pass: the Step 1b/Step 4
  `hammock_cli.cpp` overlap composes safely (disjoint edit classes at the
  same lines), and no `.def()` name collision with `estimate_reg_eq_similarity`.

**Step 1b is ready to implement.** Proceeding to the actual edits, folding in
both non-blocking comment-sweep findings and the strengthened grep
verification.

**Commit:** two — (1) the C++ method rename (Part 1, all
declaration/definition/call sites), (2) the Python-facing binding rename
(Part 2: the `.def(...)` registration plus its 8 test call sites) — kept
separate since Part 2 is the one with any (small, in-repo-only) API-surface
implication, worth its own reviewable diff. Each commit's message records a
green `cpp/tests/`/`pytest tests/` run as evidence it compiles and behaves
identically, plus (per the risk/safety finding above) a clean
`grep -rn 'estimate_jaccard\b'`/`grep -rn 'jaccard_similarity\b'` sweep of the
relevant files for Part 2/Part 1 respectively.

**Step 1b landed as two commits on the worktree branch** (`cfde6e7` Part 1 —
C++ method rename + the 6 comment-sweep sites the pre-implementation gate
authorized; `a64cd2c` Part 2 — Python binding rename), each verified before
commit with a fresh build (conda-env compiler, per the cluster-compiler
caveat), the required `grep` sweeps (both empty), `cpp/tests/hammock_tests`
(21/22 passed — see below), and `pytest tests/` with `HAMMOCK_REQUIRE_CPP=1`
(264 passed, 8 skipped, all bedtools/samtools-PATH, none digest-gated in this
environment). Then re-reviewed by 3 fresh adversarial subagents against the
**actual diff** (correctness line-by-line, blast radius/scope, verification-gap
hunting) — same pattern as Step 1's post-implementation review. All three
independently rebuilt the worktree and reproduced every claimed number
exactly; two came back fully clean (one cosmetic line-wrap nit in a
`hammock_cli.cpp` comment; one note that a reworded comment there inherited,
rather than introduced, a minor pre-existing ambiguity about which code path
actually computes `reg_jac`). The third surfaced a real finding, below.

**One pre-existing, unrelated test failure, confirmed mechanically not
caused by this rename.** `cpp/tests/hammock_tests`'s "Mode A: chr/non-chr
prefixes normalize identically" fails identically (`0 == Approx(1)`) at the
pre-Step-1b baseline (`3f1e5ca`) and after both commits — the implementer
stashed/rebuilt to confirm this during implementation, and the
verification-gap reviewer independently confirmed it more strongly still: the
`git diff` touching this test and the method it calls is a pure identifier
substitution with byte-identical logic either side, so the failure is
mechanically guaranteed to be identical, not merely observed to be. Apparently
`cpp/tests/` (it needs `-DHAMMOCK_BUILD_TESTS=ON`, off by default, plus
`-DCMAKE_POLICY_VERSION_MINIMUM=3.5` for the vendored doctest's old CMake
floor) had never actually been built and run in this repo's session history
before this Step — worth a follow-up outside this seed's scope; not fixed
here, not blocking Step 1b.

**Real finding, not blocking these two commits, but blocking the broader
rename goal's completeness claim: `HLLSketch::jaccard_and_union_cardinality`
still says "jaccard" and is the actual hot-path producer of the
`reg_eq_similarity` CSV value.** Found by the verification-gap reviewer.
Step 1b's plan and its grep-based verification (`jaccard_similarity(`,
`estimate_jaccard\b`) were both scoped to the two specific identifiers named
in the user's request — neither pattern matches
`jaccard_and_union_cardinality`, a different, pre-existing, unrenamed method
(`cpp/include/hammock/hll_sketch.hpp:67`, `cpp/src/hll_sketch.cpp:173`) that
computes the *same* register-equality quantity fused with a union-cardinality
pass, with an output parameter literally named `jaccard`. It is not a rarely
used alternate path: `bindings/_core.cpp:288`'s `pairwise_metrics_hll`
(called unconditionally for every CSV row by both `runner.py:441,634`) and
`hammock_cli.cpp:494` (the default/full/metrics arms of `hammock-cpp`) both
call it to produce the value written into the shipped `reg_eq_similarity`
column — `reg_eq_similarity()` itself (the method Step 1b renamed) is only
reached via `pairwise_jaccard_hll` and the `--register-equality` cheap-arm
fallback branch (`hammock_cli.cpp:496`). So on the default/full/metrics code
path, the value in the `reg_eq_similarity` column is actually computed by a
method whose name and parameter still say "jaccard" — exactly the surface
the user's stated goal ("make sure the register equality calculation isn't
named jaccard misleadingly on anything someone might be using") was aimed
at, and Step 1b's grep verification cannot see it because the identifier
doesn't match either pattern it checked for. Both commits did exactly what
their reviewed, gated plan scoped — this is a genuine scope gap in the
*plan*, found only after implementation, not a defect in the diff. Not
touched here; needs its own decision (a Step 1c, or folded into a future
step) before the rename can be called complete by the user's own stated
goal, not just by Step 1b's narrower plan text.

**Then stop.** Report and wait for the user's go-ahead before starting
Step 2 — and, separately, before deciding whether to add a step for
`jaccard_and_union_cardinality`. Do not continue automatically.

## Decision: `jaccard_and_union_cardinality` — do it now, as Step 1c (2026-08-12)

**Resolved: rename it now, before Step 2, as a new Step 1c below.** Not
deferred.

Reasoning:

- It is closer to the user's actual stated goal than either identifier Step
  1b already renamed. Step 1b's post-implementation review established that
  on the default/`--metrics`/`--full` code path — the one that produces the
  shipped `reg_eq_similarity` CSV value for the overwhelming majority of
  real runs — the value is computed by `jaccard_and_union_cardinality`, not
  by `reg_eq_similarity()` (that method is only reached via the
  `--register-equality` cheap-arm branch and `pairwise_jaccard_hll`'s
  all-pairs matrix). A method whose name and output-parameter both still say
  "jaccard" is, on the hot path, a bigger instance of exactly the problem
  ("make sure the register equality calculation isn't named jaccard
  misleadingly on anything someone might be using") than the two identifiers
  already fixed.
- It is the same low-risk class of change as Step 1b: a pure identifier
  rename inside a public-on-GitHub but currently zero-external-caller C++
  method, with the same "compiler catches a stale name" verification
  property that made Step 1b cheap to gate. No new risk category is
  introduced.
- Doing it now keeps the "eliminate misleading internal jaccard naming"
  work as one contiguous, disjoint-from-Step-2 unit on the worktree branch,
  rather than reopening the C++ internals after Step 2-4 have already moved
  on to paper/doc consumer-facing work. This doc already contains one
  cautionary example of a promised follow-up slipping (the v0.8.0 "SHIPPED"
  note that still hadn't landed in CLAUDE.md as of this doc's own Scope F) —
  deferring this finding "for later, unscheduled" risks the same fate, and
  there is no step downstream of Step 2 that naturally picks it back up
  otherwise.
- It does not entangle with or block Step 2's plan. Step 2 is entirely about
  paper/experiments *reading* code (CSV column names, R/Python consumer
  scripts); `jaccard_and_union_cardinality` is an internal C++ method with
  zero paper/experiments callers (confirmed by grep below) and zero effect
  on any CSV's shape or values. The two are cleanly separable in scope,
  files touched, and review lens even though both live under the
  "misleading naming" theme.

Considered and rejected: deferring to "whenever Step 4 touches CLAUDE.md
anyway." Rejected because Step 4 is about the *user-facing* column/version
contract (the "Done when" checklist's `grep -rn 'jaccard_similarity('`
check is explicitly scoped to `cpp/` and `bindings/` call *syntax*, which
this rename actually satisfies more completely if done first — doing it
after Step 4 would mean Step 4's own "Done when" grep needs re-running
anyway). Folding it into Step 4 also mixes a pure-internal-rename commit
into Step 4's already-larger column-removal-plus-version-bump commit set for
no benefit.

### Step 1c — Rename `HLLSketch::jaccard_and_union_cardinality`

**Scope, confirmed by grep across `cpp/`, `bindings/`, `python/`, `tests/`,
`paper/`, `experiments/`, `docs/`, `README.md`, `CLAUDE.md` (2026-08-12,
this Step's own inventory — no reliance on Step 1b's grep, which was scoped
to different identifiers and never matched this one):**

New names: method `HLLSketch::jaccard_and_union_cardinality` →
`reg_eq_and_union_cardinality` (mirrors `reg_eq_similarity()` — short,
consistent root, and still names both things the fused pass computes,
matching the existing naming pattern rather than inventing a new one);
output parameter `jaccard` → `reg_eq` (matches the renamed CSV column and
the renamed sibling method).

**Site 1 — declaration + doc comment**, `cpp/include/hammock/hll_sketch.hpp:51-69`
(corrected from an earlier draft's "52-68" — off by one at each boundary,
found by review):
the doc comment above the declaration ("yielding the register-equality
Jaccard *and* the cardinality of the union") and the `void
jaccard_and_union_cardinality(const HLLSketch& other, double& jaccard,
double& union_cardinality) const;` signature itself.

**Site 2 — definition**, `cpp/src/hll_sketch.cpp:173-175` (signature),
`:181-182` (the exception message: `"HLLs must have same precision and hash
size for fused jaccard/union"` — a string a caller could see in a thrown
exception, in scope for the same "anything someone might be using" reason
as the rest of this rename), `:202-203` (corrected from an earlier draft's
"201", found by review — the `jaccard = (active == 0) ? 0.0 : ...`
assignment, renamed to `reg_eq = ...` to match the renamed parameter).

**Site 2b — two more exception-message strings, in the already-renamed
`reg_eq_similarity()` method, same file — found adjacent by the Step 1c
review gate's scope-completeness pass, folded in here rather than left as a
gap:** `cpp/src/hll_sketch.cpp:52` (`"Cannot compute Jaccard between
different sketch types"`) and `:61` (`"HLLs must have same precision and
hash size for Jaccard"`). These are leftover from Step 1b — its grep
verification was scoped to `jaccard_similarity(` and `estimate_jaccard\b`,
neither of which matches a string literal inside the method body, so they
slipped through even though `reg_eq_similarity()` was the method being
renamed at the time. Same class of fix as Site 2's exception message
(caller-visible text, no logic change), same file already being edited by
this Step — mechanical, not a scope expansion of what this rename is for.

**Site 3 — call site**, `cpp/app/hammock_cli.cpp:494`:
`qh[i]->jaccard_and_union_cardinality(*rh[j], reg_jac, u);` → rename the
call only; `reg_jac` (the caller's local variable) is already correctly
named and untouched. Same file, `:483-489`: the multi-line comment above
this call currently says "computed by the reg_eq_similarity() C++ method" —
**this is already slightly wrong** (a pre-existing inaccuracy Step 1b's own
post-implementation review flagged as "inherited, not introduced": `reg_jac`
on this code path is actually computed by
`jaccard_and_union_cardinality`/`reg_eq_and_union_cardinality`, not by
calling `reg_eq_similarity()`, which is only reached in the `qh[i] &&
rh[j]` — false fallback branch a few lines below). Step 1c fixes this
comment for real while renaming the identifier it names, rather than
carrying the inaccuracy forward under a new name.

**Site 4 — call site + comment**, `bindings/_core.cpp:285` (comment: "See
HLLSketch::jaccard_and_union_cardinality for why this is bit-identical
rather than merely close"), `:288` (`a[i]->jaccard_and_union_cardinality(*b[j],
jbuf(i, j), u);` — rename the call only).

**Deliberately out of scope, stated here for the review gate to weigh, not
silently excluded:** the local variable `jaccard` and accessor `jbuf` inside
`pairwise_metrics_hll` (`bindings/_core.cpp:240,243` — corrected from an
earlier draft's "236-239", found by review) are themselves named
after the same misleading word and hold exactly the register-equality
values this whole rename is about — but they are call-site-local
implementation details, never exposed as a name to any caller (the function
returns a plain positional tuple), and renaming them is not what either the
Step 1b finding or this decision's reasoning above calls for (which is
specifically the *method* and its *formal parameter*, the two things a
future caller reads to understand the API). Matches Step 1b's own precedent
of leaving already-fine local variable names alone (`reg_jac`) rather than
sweeping every local variable that happens to touch the value.

**Comment/doc updates outside `cpp/`/`bindings/`:**
- `tests/test_containment_estimator.py:186` ("The fused jaccard+union pass
  (HLLSketch::jaccard_and_union_cardinality)...") and `:476` ("the Flajolet
  fallback in jaccard_and_union_cardinality") — comment-only, no test logic
  touched (neither line calls the method directly; both describe it in
  prose motivating the test).
- `CLAUDE.md:376` (the `docs/seed-mode-d-hash-width.md` "Open seeds" bullet,
  under "The pairwise metrics loop is one fused register pass") — this is a
  *living* reference into still-open work, not a dated snapshot, so it stays
  accurate rather than going stale the moment this Step lands. Update the
  identifier only; the surrounding technical claim (one fused register pass,
  bit-identical by construction) is unaffected by a rename and stays as
  written.
- `CLAUDE.md:379`, same paragraph as `:376` above — "It replaced
  `jaccard_similarity()` + `intersection_size()`" is *already* stale, a
  leftover from Step 1b (which renamed `jaccard_similarity()` to
  `reg_eq_similarity()` but didn't touch this sentence) — found by the Step
  1c review gate's scope-completeness pass. Fixed here since it's in the
  exact paragraph this Step is already editing and the same "living
  reference, stays accurate" reasoning applies.
- `CLAUDE.md:810`, a separate location (divergence #9's "Not implemented"
  section) — `` `jaccard_similarity(const AbstractSketch&)` + `dynamic_cast`
  shape `` describes `HLLSketch`'s current inheritance pattern, not a dated
  measurement, using the pre-Step-1b method name. Also found by the Step 1c
  review gate's scope-completeness pass, also a Step-1b leftover (not
  something Step 1c's own rename touches functionally), fixed here in the
  same CLAUDE.md commit rather than left for a third pass to find.

**Explicitly left untouched, matching the established historical-doc
precedent from this seed's own Scope F:** `docs/seed-hammock-cpp-file-dispatch.md:191,316`
and `docs/seed-metrics-column-restructure.md:77` are dated records of past
findings/measurements (Part 2's performance investigation; Part 4's cost
note) that named the function as it was called *at the time* — rewriting
them to the new name would misrepresent what was true when written, the
same reasoning already applied to every other historical doc this seed
touches. No paper/, experiments/, or README.md reference exists (confirmed
by grep, zero hits in all three).

**Review gate:** 3 reviewers confirm (a) the site list above is complete —
independently re-grep `cpp/`, `bindings/`, `python/`, `tests/`, `paper/`,
`experiments/`, `docs/`, `README.md`, `CLAUDE.md`, `memory/` for
`jaccard_and_union_cardinality\b` and for the bare word `jaccard` inside
`hll_sketch.hpp`/`hll_sketch.cpp`/the touched `hammock_cli.cpp`/
`bindings/_core.cpp` regions, to catch anything the snapshot above missed;
(b) the "deliberately out of scope" call on `bindings/_core.cpp`'s local
`jaccard`/`jbuf` variables is the right line to draw, not scope-creep
avoidance dressed up as principle; (c) the hammock_cli.cpp comment fix
actually resolves the inherited ambiguity Step 1b's post-implementation
review flagged, rather than just moving it; (d) no CMakeLists/export-map/
`.pyi`/duck-typed-caller ABI-surface risk (same check Step 1b's risk/safety
pass already ran for the sibling rename — re-run here since this is a
different identifier, don't just cite the prior clean result — this rename
is actually lower-risk than Step 1b's, since Site 3/4 are the *only* callers
and both are compiled C++, with zero `.def(...)` pybind11 exposure at all
for this identifier).

Verification methodology, mirroring Step 1b, **with one gap closed by
review before implementation, not after**: since every site is either a
compiled C++ call/declaration (stale name = hard compile error) or a
comment (no runtime effect), the completeness proof is a green
`cpp/tests/hammock_tests` build + `pytest tests/ --HAMMOCK_REQUIRE_CPP=1` +
a clean `grep -rn 'jaccard_and_union_cardinality\b'` sweep of the whole
repo, **expecting zero hits outside the two deliberately-preserved
historical docs** (`docs/seed-hammock-cpp-file-dispatch.md`,
`docs/seed-metrics-column-restructure.md` — corrected from an earlier
draft's unqualified "expect zero hits," which self-contradicted this Step's
own decision, two paragraphs up, to leave those two files untouched; found
by the risk/safety review pass).

**Closing a real gap in that methodology, found by the risk/safety
review pass, not by a diff change:** `cpp/tests/hammock_tests` (via
`HAMMOCK_BUILD_TESTS`) only builds against the `hammock_core` static lib —
it covers Sites 1-2 but not Site 3 (`hammock-cpp`, gated by
`HAMMOCK_BUILD_CLI`) or Site 4 (`_core.so`, the pybind11 extension). Site 3
is already guarded against silent staleness by
`tests/test_hammock_cpp_metrics.py::test_binary_is_not_stale` (asserts the
`hammock-cpp` binary's mtime is >= every file under `cpp/app`, `cpp/src`,
`cpp/include`). **Site 4 has no equivalent guard** — `pytest tests/` alone
does not rebuild `_core.so`, so an implementer who edits
`bindings/_core.cpp` and runs `pytest tests/` against an already-installed
extension could get a fully green run without the edited code ever having
been recompiled, on the single most consequential site (the hot path
`runner.py` actually calls for every default/full/metrics CSV row).
**Required before this Step's commits can rely on "pytest green" as
evidence:** an explicit `pip install -e . --no-build-isolation` (per the
cluster-compiler-caveat build in CLAUDE.md) immediately before the
verification `pytest` run, confirmed by checking `_core*.so`'s mtime is
newer than `bindings/_core.cpp`'s post-edit mtime — recorded in the commit
message alongside the other verification evidence, the same place Step 1b
recorded its build/grep evidence.

**Commit:** two, once the review gate clears. **Correction, found by two of
the three review passes independently: an earlier draft justified this
split as "the same reason Step 1 and Step 1b both split this way" — that
citation is wrong.** Step 1's split axis was code vs. archived-CSV *data*;
Step 1b's was internal-C++-method vs. public-Python-binding (an
API-surface distinction), and Step 1b actually folded its whole comment
sweep into the *same* commit as its code rename rather than splitting them
out. Step 1c's split is its own, different axis — compiler-verified vs.
not — stated on its own merits below, not as a repeat of either precedent:
(1) the functional rename (method + parameter + both call sites + all
three exception-message strings, i.e. everything in `cpp/` and
`bindings/_core.cpp:288`'s call site that a compiler gates — a stale name
anywhere in this commit is a hard build failure, which is what makes this
commit cheap to gate), (2) the comment/doc-only updates (the
`hll_sketch.hpp` doc comment, `bindings/_core.cpp:285`'s comment, the
`hammock_cli.cpp` comment including its ambiguity fix, the two
`test_containment_estimator.py` comments, `CLAUDE.md:376,379,810`) — kept
separate because (2)'s correctness rests on human review of prose, not on
the compiler, and reviewers may want to weigh that evidence differently
from (1)'s.

Each commit's message must record: a green `cpp/tests/hammock_tests` build,
a green `pytest tests/ --HAMMOCK_REQUIRE_CPP=1` run against a freshly
rebuilt `_core.so` (per the verification-methodology fix above — confirm
`pip install -e . --no-build-isolation` ran after the edit, not before),
and a clean `grep -rn 'jaccard_and_union_cardinality\b'` sweep (zero hits
outside the two preserved historical docs).

**Review gate, run 2026-08-12 (3 parallel subagent passes — scope
completeness, risk/safety, process/convention fit — against the plan text
as it stood before this paragraph's own edits, i.e. the version committed
in `082f042`):**

- **Process/convention fit: clean.** Commit discipline directly verified
  (`082f042` on `main` is doc-only; worktree branch confirmed still at
  `a64cd2c`, untouched). The Decision section's reasoning and its "SHIPPED
  note slipped" citation were both independently verified against actual
  repo state, not taken on faith. One non-blocking finding: the "Commit:"
  paragraph's claim to mirror Step 1's/Step 1b's split rationale was
  inaccurate — fixed above.
- **Scope completeness: site list for `jaccard_and_union_cardinality`
  itself confirmed complete** by an independent repo-wide grep (matches the
  plan's enumeration exactly, including the two correctly-preserved
  historical docs and the correctly-empty paper/experiments/README.md
  claim). Two adjacent findings, both folded in above rather than left as
  gaps: the two leftover "Jaccard" exception-message strings inside the
  already-renamed `reg_eq_similarity()` (Site 2b, new), and two leftover
  stale-name references in CLAUDE.md (`:379`, `:810`) that Step 1b's own
  comment sweep missed. Three line-number citations had drifted from actual
  file content (all corrected above — mechanical, no scope change). The
  "deliberately out of scope" call on `bindings/_core.cpp`'s local
  `jaccard`/`jbuf` variables was independently verified correct: every
  consumer of `pairwise_metrics_hll`'s return value unpacks it as a plain
  positional tuple, never by name.
- **Risk/safety: no blocking issue with the rename itself** — confirmed
  lower ABI-surface risk than Step 1b's sibling rename (zero pybind11
  `.def(...)` exposure for this identifier at all; only two call sites,
  both compiled C++). No exception-message-text matching anywhere in the
  test suite (`pytest.raises(match=...)`, `CHECK_THROWS_WITH`), so the
  exception-string renames are safe. **Two real gaps in the verification
  *methodology* text, both closed above by strengthening verification, not
  by changing the diff** (mirroring how Step 1b's own risk/safety finding
  was closed): the unqualified "expect zero hits" completeness bar
  self-contradicted this Step's own decision to preserve two historical
  docs; and `_core.so` (Site 4, the actual hot-path file) had no staleness
  guard equivalent to `test_binary_is_not_stale`'s coverage of `hammock-cpp`
  — closed with an explicit fresh-rebuild requirement before any commit can
  cite "pytest green" as evidence.

Because these findings changed Step 1c's site list (Site 2b, two more
CLAUDE.md lines) and its verification methodology, per "Review process"'s
own re-run rule this would normally need a second focused round — but every
change here is the same "purely mechanical, no fresh round needed" class
Step 1b's own round 1 established for its comment-sweep findings (string/
line-number corrections and a verification-strengthening addition, no
change to the core rename's shape or risk profile), so proceeding directly
to implementation.

**Step 1c is ready to implement.**

**Step 1c landed as two commits on the worktree branch** (`2dff13a` — the
functional rename: method, parameter, both call sites, three exception-
message strings including Site 2b's two leftovers inside
`reg_eq_similarity()`; `1f2d739` — comment/doc updates, including the
inherited-ambiguity fix in `hammock_cli.cpp` and the two leftover
`CLAUDE.md` references at `:379`/`:810`). Both commits' verification ran as
the review gate's strengthened methodology required, not the original
unqualified version:

- Rebuilt with the conda-env compiler per the cluster-compiler caveat
  (`CC`/`CXX` pointed at `$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-{gcc,g++}`);
  confirmed via `readelf -d` that the conda env lib leads `_core.so`'s
  RPATH (no `gcc-9.3.0`) and `_DIGEST_AVAILABLE == True`.
- Confirmed `_core.so`'s mtime postdates `bindings/_core.cpp`'s edit for
  each commit — the staleness-guard requirement the risk/safety review
  added, closing the one gap that specific pass called blocking.
- `cpp/tests/hammock_tests` (`-DHAMMOCK_BUILD_TESTS=ON
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5`): **21/22 passed**, both before and
  after each commit — the same pre-existing "Mode A chr/non-chr prefixes"
  failure the Step 1b post-implementation review already traced to the
  pre-Step-1b baseline, mechanically unaffected by a pure identifier rename.
- `pytest tests/ --HAMMOCK_REQUIRE_CPP=1`: **264 passed, 8 skipped**
  (bedtools/samtools not on PATH, benign) — identical before and after both
  commits, matching Step 1b's closing baseline exactly.
- `grep -rn 'jaccard_and_union_cardinality\b'` repo-wide: zero hits outside
  the two deliberately-preserved historical docs, both before commit 2 (with
  the comment sites still pending) and after (fully clean per the corrected
  completeness bar).

**Update: a post-implementation adversarial review round was run after all,
per explicit user request following the initial report above** (3 fresh
parallel subagent passes against the actual landed diff — correctness
line-by-line, blast radius/scope, and verification-gap hunting — same
pattern as Step 1's and Step 1b's post-implementation rounds).

- **Diff correctness: fully clean.** Every hunk in both commits confirmed
  pure identifier/string substitution, zero logic change; every renamed
  call site matches its declaration (no typos, no parameter-order
  mismatch); all three exception-message strings read correctly; the
  `hammock_cli.cpp` comment rewrite verified to accurately describe the
  code below it; the three `CLAUDE.md` edits verified against current code.
- **Verification-gap hunting: fully clean.** Independently reproduced, from
  a from-scratch rebuild (not reusing any prior build artifact), every
  number the commit messages cite: RPATH/`_DIGEST_AVAILABLE` check,
  `cpp/tests/hammock_tests` 21/22 (same pre-existing unrelated failure),
  `pytest tests/ --HAMMOCK_REQUIRE_CPP=1` 264 passed/8 skipped exactly, all
  8 skips confirmed environment-gated (bedtools/samtools) with zero
  digest-related skips — so this rename has no digest-gated test dependency
  and no silent-skip risk of the kind that bit Step 1b. Separately confirmed
  the fused method's two `double&` out-params (`reg_eq`, `union_cardinality`)
  are genuinely order-checked by `test_containment_estimator.py`'s exact-`==`
  cross-checks against the independent `estimate_intersection`/
  `estimate_reg_eq_similarity` scalar paths, not merely compile-gated; and
  re-confirmed (independently, not trusting the pre-implementation pass)
  that no test anywhere matches on any of the three changed exception
  strings.
- **Blast radius/scope: one real, small miss found.** `hammock_cli.cpp:452`'s
  comment ("The fused jaccard+union path below...") describes
  `reg_eq_and_union_cardinality()` in loose prose, not the literal
  identifier — predates this rename (`826d7b42`, 2026-08-08), sits ~30
  lines from the block `1f2d739` *did* rewrite for the same method, and
  slipped through both commits' grep verification because neither pattern
  (`jaccard_and_union_cardinality\b`, case-sensitive `Jaccard` in
  `hll_sketch.{hpp,cpp}`) matches a lowercase paraphrase in a different
  file. Everything else this pass checked (local-variable exclusion,
  file-set-vs-commit-message match, paper/experiments/README.md untouched,
  a broadened spelling-variant sweep) matched the reviewed plan exactly —
  this was the sole finding across all three passes.

**Fixed as a follow-up commit, `6cce6d3`** (comment-only, one line,
`cpp/app/hammock_cli.cpp:452`), re-verified with a fresh rebuild and
`pytest tests/ --HAMMOCK_REQUIRE_CPP=1` (264 passed, 8 skipped, unchanged).
A repo-wide sweep after the fix confirms no third stray "jaccard+union"-style
paraphrase remains anywhere in the touched files.

**Then stop.** Report what landed and the subagent reviewers' findings, and
wait for the user's go-ahead before starting Step 2. Do not continue
automatically.

### Step 2 — Update paper/experiments code to read `reg_eq_similarity`, with fallback

Covers the user's step 2. Every in-scope script from Setup's table,
categories (b) and (d) — **substantially expanded by this Step's own
pre-implementation review gate, below.** Setup correctly identified the
in-scope *directories* but its per-file/per-line coverage of `experiments/`
and `docs/scripts/` was materially incomplete going into this Step; the gaps
were found and closed here, before implementation, matching the pattern
Step 1's own round 1 already established for this doc (find, fold in,
re-verify).

**Reading (category b):** prefer `reg_eq_similarity` if present in the CSV,
else fall back to `jaccard_similarity` (covers any archived file Step 1
didn't touch).

R pattern:
```r
sim_col <- if ("reg_eq_similarity" %in% names(df)) "reg_eq_similarity" else "jaccard_similarity"
```
Python (`csv.DictReader`) pattern:
```python
col = "reg_eq_similarity" if "reg_eq_similarity" in row else "jaccard_similarity"
```
(equivalent for pandas: check `"reg_eq_similarity" in df.columns`).

**Log when the fallback branch actually fires** (e.g. one `message()`/
`print(..., file=stderr())` per script run, not per row) — found necessary
by review: a bare presence check can't distinguish "this is a legitimately
old archived file" from "a fresh run silently regressed and is still
writing the old name" (e.g. a partial merge that reverted `runner.py`'s
rename but not `hammock_cli.cpp`'s, breaking Step 1's required cross-
front-end sync). A visible log line makes a fresh run hitting the fallback
path show up in run output instead of looking identical to expected legacy
handling. **For loop-heavy read sites** (below), this means hoisting the
presence check out of the loop and logging once at the top of the run — not
once per row or per file processed inside it.

**Two files resolve their column name via `DEFAULT_SIM_COL`, *before* the
CSV load — the generic pattern above cannot be dropped in verbatim, found by
this Step's review gate (scope-completeness and risk/safety passes,
independently):**
- `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R:59`
  (`DEFAULT_SIM_COL <- "jaccard_similarity"`) → `sim_col` at `:75`, consumed
  at `:217` (presence check) and, the file's actual **primary** read —
  bigger than the `:271/:286` diagnostic block already known to this Step —
  `:248` (`similarity_matrix(raw, sim_col)`, the dendrogram's real data
  source). `sim_col` is resolved at `:75`, before `raw <- read_csv(...)` at
  `:213/:216`, so the fallback check has to move to *after* the load. Treat
  `:59/:75/:217/:248` and `:271/:286` as **one** coupled unit sharing a
  single resolved variable — five sites, one file, not two independent
  edits (the original plan text named only `:271/:286`).
- `paper/cross_reference_identity/plot_cross_reference_identity.R:42`
  (`DEFAULT_SIM_COL <- "jaccard_similarity"`) → `SIM_COL` at `:50`, same
  before-the-load ordering problem relative to `raw <- read_csv(peaks_csv,
  ...)` (~`:100`). Same fix shape: move resolution after the load.

**`plot_interval_accuracy.R:154,171` needs restructuring, not a drop-in
substitution, found by review:** `:154`'s `required_hammock <- c(...,
"jaccard_similarity", "jaccard_similarity_ie")` is an AND-style `setdiff()`-
gated required-columns check, not an OR/fallback branch — loosen it to
accept either name. `:171`'s `` `Register-equality (jaccard_similarity)` =
jaccard_similarity `` inside `transmute()` is a bare NSE symbol reference,
not a string-driven lookup — resolve `sim_col` once, upstream of both sites,
then use `.data[[sim_col]]`.

**Re-emitting (category d) — `experiments/maurano_dhs_validation/analyze.R`,
corrected and substantially expanded by this Step's review gate; the
original plan text's "the `column = "jaccard_similarity"` assignments, the
`short_col` lookup table" undersold this file's real edit surface — found
independently by all three review-gate passes:**

There is no literal `column = "jaccard_similarity"` assignment anywhere in
the file; that description was wrong. The actual sites:
- `d_metric_order`/`metric_available()` (`:327-348`) — a **discovery gate**,
  not an assignment: a metric name only survives into `d_metrics` if it's
  present in *every* scanned file's header (`all(vapply(d_headers, ...))`).
  Once Step 1's rename lands, this silently drops the register-equality
  metric from discovery on any freshly regenerated Mode D CSV — no error, no
  `stop()`. Fix: accept either name per-file, mapping both to one canonical
  label, not a literal match against `"jaccard_similarity"` alone.
- `read_hammock_csv()`'s `jcol` (singular) parameter defaults to
  `"jaccard_similarity"` at `:160`; `scan_dir()`'s own `jcols` (plural)
  default is the separate call site at `:249` (corrected by round 2
  review — an earlier draft conflated the two into one `:160` citation).
  `scan_dir()` loops `jcols` and calls `read_hammock_csv(jcol=jc)` per
  candidate, and is exercised live — the ABC scan call at `:282` doesn't
  override it, so `read_hammock_csv`'s hard `stop()` on a missing `jcol` is
  reachable on that path. `:147`'s `j_truth = as.numeric(jaccard_similarity)`
  is a bare NSE read against `raw_abc/hammock_hll_p21_jaccB_full.csv` (not
  one of Step 1's 5 rewritten CSVs).
- Five downstream `filter(column == "jaccard_similarity", ...)` sites, not
  just the `short_col` table: `:373` (`d_plot`, load-bearing), `:471`
  (`refs_wide`, independently filters `d_out`, load-bearing), `:454,492,550`
  (re-filters of already-filtered frames). All five need the same
  value-comparison fallback (`column %in% c("reg_eq_similarity",
  "jaccard_similarity")`) already specified below for the six B2 consumers —
  `analyze.R` is effectively a seventh, more heavily-touched
  consumer/generator hybrid.
- `:754`'s `sub("jaccard_similarity", "no_ends", best_d_col)` — cosmetic
  label substitution in a plot title; degrades silently (wrong wording, not
  a crash) if left.
- The `short_col` lookup table (`:364-367`, the one site the original plan
  text did name correctly) maps `jaccard_similarity = "register_equality"` —
  update this key to match the new name, per Setup's own earlier flag.

**New site, found by review, absent from Setup's inventory and the original
Step 2 plan entirely: `experiments/bedtools_benchmark/estimator_compare.py:258`**
(`j_re = float(r["jaccard_similarity"])`) — reads the globally-installed
`hammock --metrics` output (not the worktree build) via `csv.DictReader`; a
bare dict-key lookup, hard `KeyError` if unfixed. Feeds
`results/estimator_compare_full.csv` (untracked, regenerated), the actual
input to `paper/estimator_crossover/plot_estimator_crossover.R`
(`docs/paper_outline.md:320` cites `paper/figures/estimator_crossover.png`)
and `estimator_rank_by_precision.py` (backs CLAUDE.md's rank-fidelity
table). Setup's own "check for other `analyze*.R`/`.py` scripts" instruction
should have caught this and didn't — the provenance trace stopped one hop
short, at the R script rather than the Python script that produces its
input. Fix: `r.get("reg_eq_similarity") or r["jaccard_similarity"]`, logged
on fallback.

**A second new site in the same directory, found by review:
`experiments/bedtools_benchmark/sweep.py`** — `parse_hammock_csv(path,
column="jaccard_similarity")` (`:127` default, `:415` call site) already
raises a deliberate, explicit `KeyError` with rich context (`:151`) if the
column is missing — a pre-existing anti-silent-failure design; extend rather
than replace it. In scope: `docs/paper_outline.md:502-547,379` cites
`sweep_precision_maurano_p18_t16.csv`/`_t8.csv` by exact filename as this
script's own output, feeding Figure 3B's cross-check and Figure S8 (via
`plot_precision_frontier.R`). Fix: try `reg_eq_similarity` first inside
`parse_hammock_csv`, fall back to `column`, keep the hard-fail behavior when
*neither* name is present, log the fallback.

**New group, found by review: `experiments/ref-comparison/`'s other 3 files
(of its 4), missed by Setup's original per-file characterization** (Setup
correctly flagged the directory as in-scope via `paper/draft.md:101,273`/
`paper/outline.md:74` for `estimator_ie_crossref.py`'s own output, but never
walked its sibling files):
- `estimator_ie_crossref.py:84` (`reg = float(row["jaccard_similarity"])`) —
  hard `KeyError`, directly in scope (Table S5 provenance, confirmed).
- `scripts/exp_a_dendrogram.R:51` (`sim[mat$key_a[i], mat$key_b[i]] <-
  mat$jaccard_similarity[i]`) — hard error (`"replacement has length
  zero"`). In scope only via `docs/paper_outline.md:214` ("Fig 7") — not
  cited in `paper/outline.md`/`paper/draft.md`, but a `docs/paper_outline.md`
  citation is already this doc's established in-scope standard (the same
  standard that found `docs/scripts/` during Setup) — treat as in scope on
  that basis. The figure-numbering scheme there looks superseded relative to
  the live manuscript; noted as a residual observation, not a reason to skip
  the fix.
- `scripts/exp_a_metric_comparison.R:36,84,100` (corrected by round 2
  review — the `tribble(` defining `metric_groups` starts at `:36`, not
  `:38`) (`metric_groups` tribble →
  `pull(.data[[metric]])` → `write_tsv`) — **silent degradation, not a
  crash**: the file already filters to present columns and prints a generic
  "skipping" message, so the row just vanishes with no distinguishing
  signal. In scope via `docs/paper_outline.md:372,389` ("Fig S1"), same
  standard as above.
- `scripts/exp_a_validate_plot.R` — **corrected by this Step's round 2
  review: the original citation was simply wrong** (`docs/paper_outline.md:69,148`
  is blank/unrelated prose, not this script — that document's actual "Fig 3"
  is `mode_d_bedtools_vs_modeB_scatter.png`, a different file entirely).
  Real chain, two-hop like the primate-phylogeny/mus-homo §9.6 citation
  elsewhere in this Step: `docs/paper_outline.md` → `experiments/ref-comparison/docs/exp_a_results.md`
  (which names this script and its `cross_ref_validation.png`/
  `cross_ref_stats.tsv` outputs directly) — in scope on that basis. **Fix
  mechanism also corrected, round 2 risk/safety finding:** the original
  plan's "fix belongs in the Snakefile, resolve to a `reg_eq_similarity`-
  preferred value there" was underspecified and risky — a Snakemake
  `params:` block is a static literal, not a function, so a naive
  implementation would just swap the literal string with no fallback at
  all, violating the prefer-new-fall-back-to-old requirement, and has no
  compile/pytest safety net if done wrong. Fix instead **inside
  `exp_a_validate_plot.R` itself**, mirroring the precedent this Step
  already sets for the structurally identical `mus-homo/scripts/cluster_plot.R`
  case below: resolve `sim_col` from `mat`'s actual columns (prefer
  `reg_eq_similarity`, fall back to the Snakemake-passed
  `snakemake@params[["sim_col"]]` literal, log the fallback), leaving
  `experiments/ref-comparison/workflow/Snakefile:180`'s literal untouched as
  an unused legacy default — no Snakefile edit needed after all.

**`experiments/subB_mixed_stride/run_sweep.py:130-131` already has a
name-tolerance OR-chain with its own silent-degradation bug, found by
review, pre-existing but exposed by this rename:**
`r.get("jaccard_similarity") or r.get("jaccard") or
r.get("jaccard_similarity_ie")` — once `jaccard_similarity` is renamed, this
falls through *silently, with no log line*, to the IE estimator's value in
the full/`--metrics` shape (where both `jaccard_similarity` and
`jaccard_similarity_ie` are present in one row), quietly substituting a
different estimator into a field every downstream consumer treats as
register-equality. In scope directly — `paper/draft.md:282` cites this
script's own output (`sweep_maurano_ie_20260809_200658.csv`) for Figure S9.
Fix: extend the chain to `r.get("reg_eq_similarity") or
r.get("jaccard_similarity") or r.get("jaccard") or
r.get("jaccard_similarity_ie")`, and — since this bug predates the rename —
add the log line the chain never had, distinguishing which candidate
resolved.

**New cross-species group, found by review, resting on a two-hop citation
chain already established elsewhere in this doc
(`docs/paper_outline.md:319` → `docs/estimator-analysis-findings.md` §9.6,
which names both scripts and their numbers directly) —
`experiments/primate-phylogeny/` and `experiments/mus-homo/`:**
- `primate-phylogeny/estimator_ie_topology.py:67`
  (`re_j[(a,b)] = float(row["jaccard_similarity"])`) — hard `KeyError`.
- `primate-phylogeny/scripts/build_phylogeny.R:82,86` (corrected by round 2
  review — the `all_metrics` literal vector is defined at `:82`, the
  `intersect(all_metrics, names(mat_long))` consuming it is at `:86`, not
  both at `:82` as an earlier draft said) — silent degradation (row quietly
  drops from `metric_spreads.tsv`).
- `primate-phylogeny/scripts/build_phylogeny.R:160` (corrected by round 2
  review — `:110-118` is the parameterized `build_and_plot` function
  *definition*, not a hardcoded site; the only actual hardcoded call-site
  literal, `build_and_plot("jaccard_similarity", ...)`, is at `:160` alone)
  — silent degradation, the worse of the two: on a missing column
  it writes an empty Newick/dist file and a placeholder PNG, and the
  Snakemake rule still reports success. This is the file behind the
  primate-clade-recovery headline the citation chain rests on. Fix: resolve
  `sim_col` once (preferred `reg_eq_similarity`, logged fallback) before the
  `build_and_plot` call, rather than hardcoding the literal at the call
  site.
- `mus-homo/estimator_ie_tissue.py:59` — hard `KeyError`, same pattern as
  the primate-phylogeny sibling.
- `mus-homo/scripts/cluster_plot.R:34,48`, reached only through
  `workflow/Snakefile:155` → `config/config.yaml:26`'s `primary_sim_col` —
  **invisible to a literal grep of the R script itself**, found only by
  tracing the Snakemake param chain. Hard error (`!!sym()` subset on a
  missing column). This is the live generator behind the mus-homo "0/20
  tissue-dominant" ARI headline. Fix in `cluster_plot.R`: resolve the live
  column name from `mat_long` directly (prefer `reg_eq_similarity`), rather
  than trusting the Snakemake-passed literal — the config value can stay as
  a legacy default.

**Explicitly out of scope, found by review, recorded here rather than
silently skipped (same "not provably in scope" standard the rest of this
doc already applies):**
- `experiments/subB_mixed_stride/run_ie_subb.py` and `analyze_ie_subb.py` —
  coupled read/re-emit pair, zero citations found in `paper/outline.md`,
  `paper/draft.md`, or `docs/paper_outline.md`; their sole consumer is
  `experiments/subB_mixed_stride/RESULTS.md`'s own internal write-up, not a
  manuscript artifact.
- `experiments/subB_mixed_stride/sbatch_synthetic_ie.sh` — its one real run
  was deliberately cancelled and superseded
  (`docs/seed-subsampling-synthetic-supplement.md:62-89`); Figure S10 comes
  from a different job/script.
- `experiments/synthetic_evolution/code/analyze.R` — the original Setup
  snapshot never claimed this file was "provably in scope" the way it did
  for `ref-comparison/`/`subB_mixed_stride/` (re-read: "The first two are
  provably in scope..." — `synthetic_evolution`'s status was left
  unstated). This Step's review independently searched all three provenance
  documents for all five of its output figure names and found zero
  citations. It already has partial tolerance
  (`intersect(c("jaccard_similarity","jaccard"), cols)[1]`) that would need
  extending if this file is ever brought into scope later — not now.
- `experiments/primate-phylogeny/scripts/precision_probe.sh:106` — found by
  review, has the single worst failure mode in this whole sweep (awk's
  `h["jaccard_similarity"]` resolving to `""` on a missing key silently
  aliases `$jcol` to `$0`, corrupting every reported value with no error at
  all) but is confirmed not wired into any Snakemake rule, not cited by any
  provenance document, and produces only ephemeral stdout under `/vast`, not
  a tracked artifact. Left unfixed per the "not provably in scope"
  standard — flagged here, in writing, so a future toucher of this file
  knows the hazard exists rather than discovering it fresh. **Severity note,
  added by round 2 review:** this failure mode gets *worse* once Step 1
  lands — today it only misfires on a CSV that happens to lack the column;
  after Step 1, every freshly-generated hammock CSV lacks it (emitters carry
  no back-compat; only Step 2's readers get a fallback, and this script
  isn't one of them), so silent corruption becomes this script's permanent
  behavior, not an occasional hazard. **Also add this file as a named
  exception to Step 4's "Done when" grep bar**, below, so its known,
  deliberately-unfixed literal doesn't read as an unaddressed leftover at
  that gate.
- `experiments/mus-homo/scripts/compute_column_comparison.{R,py}` — already
  broken by the earlier, unrelated `_with_ends` column removal
  (self-documented "OBSOLETE... DOES NOT RUN on current output"), same
  dead-code class as `modeD_flanking/`, already excluded by this doc's
  precedent.
- `experiments/primate-phylogeny/config/config.yaml:36`'s `primary_sim_col`
  value — confirmed by grep that nothing in that experiment tree actually
  reads this key; no functional fix needed.
- `experiments/ref-comparison/config/config.yaml:25`'s `primary_sim_col`
  value — same as the primate-phylogeny case above (found by round 2 review,
  the exactly-analogous exclusion for this directory was never stated):
  confirmed by grep that nothing in `ref-comparison/` reads this key; no
  functional fix needed.
- **Three more `experiments/maurano_dhs_validation/` scripts, found by round
  2 review, missed despite `analyze.R` in the same directory getting the
  deepest treatment in this whole plan** — all three read `jaccard_similarity`
  as a bare literal, and all three were checked against all three provenance
  documents (`paper/outline.md`, `paper/draft.md`, `docs/paper_outline.md`)
  with zero citations found for any of their output figures, so none is
  provably in scope:
  - `mode_c_interpolation.R:44,50,52,58` (`select(..., jaccard_similarity)`/
    `rename(j_A = jaccard_similarity)`) — would hard-error post-Step-1;
    outputs (`mode_c_subB_interpolation_agg.png`,
    `mode_c_expA_interpolation_agg.png`) cited nowhere.
  - `render_dendrogram.R:27` (`jcol <- ... else "jaccard_similarity"`) — a
    standalone CLI utility's hardcoded default; not cited anywhere.
  - `make_metric_plots.R:35,114` (`filter(column == "jaccard_similarity",
    ...)`) — silent degradation (empty filter, no error) if left. Its
    outputs (`mode_d_lines_p24.png`, `mode_d_violins_by_k.png`) share
    filenames with the in-scope `docs/scripts/mode_d_lines.R`/
    `mode_d_violins.R` outputs but read from a *different*, untracked
    source (`results/mode_d_summary.csv`, a symlink into scratch, not the
    tracked `docs/data/mode_d_summary.csv`); `PLOT_GENERATION.md` marks it
    "executed... archived spec, not an outstanding task." Left out of scope
    on the same "not provably in scope" standard as the rest of this list,
    recorded here rather than silently omitted so a future pass doesn't
    mistake this for an oversight — matching the same reasoning already
    applied elsewhere in this doc (e.g. Scope F's closing paragraph).

**Reading `mode_d_summary.csv`'s `column` field (category B2), specifically:**
**six consumers, not five** (corrected 2026-08-12 by Step 1's round-1 risk/
safety review — the original five-consumer count missed a sixth) —
`paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R`,
its non-`_estimators` sibling
`paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff.R`
(same `SIM_COLUMNS <- c("jaccard_similarity", "jaccard_similarity_ie")`
pattern, confirmed live — not dead code — by `paper/draft.md:122`/
`paper/outline.md:95`'s Figure 7 provenance paragraph: "the existing
`plot_parameter_objective_tradeoff.R` is retained unchanged as the
single-estimator diagnostic workflow"), plus all four of
`docs/scripts/mode_d_violins.R`, `mode_d_lines.R`, `mode_d_metric_tradeoff.R`,
`mode_d_bedtools_vs_modeB_scatter.R` filter on
`column == "jaccard_similarity"`. Apply the same reg_eq_similarity-preferred/
jaccard_similarity-fallback pattern to all six (their `column` field is data
read from the regenerated CSV, not a CSV-level column name, so the fallback
here is a value comparison — `filter(column %in% c("reg_eq_similarity",
"jaccard_similarity"))` or equivalent — not a `names(df)` presence check).

Check whether `paper/` already has a shared sourced-R-utility file before
inventing one just for this helper — if 3+ scripts need it, factor it once;
if not, inline is fine and matches this repo's existing per-script style.
**Checked by review: no `source(`-based shared helper exists anywhere under
`paper/` or `docs/scripts/` today, and this Step's expanded scope now has
10+ distinct scripts needing the same fallback-resolution logic — past the
doc's own "3+" threshold.** Left as an implementer's call at diff time
(inline duplication is not blocking), but flagged so it isn't an
accidental miss of the doc's own stated rule.

**Verification methodology, added by this Step's review gate — R/Python
scripts have no pytest/`cpp/tests`-equivalent safety net, so "the review
gate confirms the pattern was applied" is not sufficient on its own, found
by review:** before this Step's commits can be called verified, run every
edited script at least once against a real input in each of its two
fallback branches where feasible, confirming non-empty, non-`NA` output and
the expected log line in each case — not deferring this entirely to Step
3's expensive full regenerate-and-diff pass, which is a correctness *proof*,
not a cheap syntax/logic check. Record per-file pass/fail in the commit
message(s). **Made concrete for the full expanded scope, round 2 finding —
the original wording only covered the handful of scripts reading Step 1's
specific 5 rewritten archived CSVs, not the majority of this Step's now-much
larger scope:**
- For scripts reading one of Step 1's 5 rewritten CSVs
  (`plot_interval_accuracy.R`, `plot_cross_reference_identity.R` on
  `exp_a_broad_k10_w10.csv`): use those files directly for the
  "prefers `reg_eq_similarity`" branch, an unmodified archived CSV predating
  Step 1 for the fallback branch.
- For everything else (the large majority — `plot_sequence_tissue_clustering.R`,
  the six `mode_d_summary.csv` consumers, `analyze.R`, `estimator_compare.py`/
  `sweep.py`/`estimator_ie_crossref.py`, the ref-comparison and cross-species
  scripts): generate a small fresh CSV by running the worktree's own
  post-Step-1 `hammock`/`hammock-cpp` build for the "prefers" branch, and
  reuse any convenient pre-Step-1 archived CSV (or a small synthetic
  pre-rename-header fixture) for the fallback branch.
- **`exp_a_validate_plot.R` and `mus-homo/scripts/cluster_plot.R` are
  Snakemake `script:`-directive files**, not directly `Rscript`-runnable —
  they reference an injected `snakemake@params[[...]]` S4 object. Verify
  either by invoking `snakemake -R <rule>` against a minimal/mocked DAG, or
  by stubbing a `snakemake` S4 object with the needed `params`/`input`/
  `output` slots at the top of an interactive session before `source()`-ing
  the script body — state which was used in the commit message, since a
  plain `Rscript script.R` will not exercise these two.
- Verification runs of `exp_a_dendrogram.R`/`exp_a_metric_comparison.R`
  should write to a scratch path, not their real default output paths —
  see the residual-risk note below for why.

**Accepted residual risk, recorded in writing per this doc's own rule,
found by review, corrected by round 2 review:** Step 2 updates several
generating scripts (`analyze.R`, `estimator_compare.py`, `sweep.py`, the
`ref-comparison`/cross-species scripts) to emit `reg_eq_similarity`-
preferring output, but does **not** regenerate the tracked CSVs/figures
those scripts produce (Step 3 does that). Between Step 2 landing and Step 3
regenerating, anyone who runs one of these updated generators against a mix
of pre- and post-Step-1 raw hammock CSVs on this worktree branch could
silently produce a locally-regenerated file that diverges from what's
checked in. Same class of risk Step 1's review gate already accepted in
writing for `plot_interval_accuracy.R`/`plot_cross_reference_identity.R`'s
transient-breakage window; accepted here on the same grounds — nothing in
this plan runs these generators during that window (the verification runs
above write to scratch paths, not real output paths, precisely to keep it
that way), the worktree branch never touches `main` before Step 4's merge,
and Step 3's regenerate-and-diff is specifically designed to catch exactly
this class of drift before it could reach `main`. **The specific tracked
files at risk, corrected — round 2 found the original list both wrong and
incomplete:** `docs/data/mode_d_summary.csv` and the three
`paper/*_stats.csv` files are tracked and at risk, as originally stated.
`results/exp_a_estimator_delta.csv` is **not** — it's gitignored
(`experiments/ref-comparison/.gitignore:2`), confirmed via `git check-ignore
-v`, so removed from this list. **Missing from the original list entirely:
7 force-tracked figures under `experiments/ref-comparison/figures/`**
(`cross_ref_dendrogram_k10_w10.png`, `cross_ref_dendrogram_k15_w15.png`,
`cross_ref_validation_broad_k10_w10.png`, `metric_comparison_broad_k10_w10.png`,
`metric_comparison_narrow_k10_w10.png`, `sweep_effect_size_broad.png`,
`sweep_effect_size_narrow.png`) — tracked via `git ls-files` despite the
directory being gitignored, and the actual manuscript figure sources for
Figs 7/S1. `exp_a_dendrogram.R`/`exp_a_metric_comparison.R` (both edited by
this Step) default-write to exactly these paths, so a verification run that
doesn't redirect output would clobber them with locally-regenerated,
not-yet-Step-3-proven data — the verification methodology above now says to
write to scratch instead, closing this specific instance of the risk rather
than just accepting it.

**Review gate:** 3 reviewers confirm the fallback pattern (with logging) is
applied consistently, that every coupled site-group above was edited
together (not partially), that the category-B2 and `analyze.R` generating-
code updates are complete, that the newly-scoped `experiments/` and
cross-species sites are correctly included/excluded per the citations
above, and that no script was missed.

**Review gate round 1, run 2026-08-12 (3 parallel subagent passes — scope
completeness, risk/safety, process/convention fit — against the plan text as
it originally stood, i.e. the version committed before this paragraph's own
edits; followed by 2 further investigative sweeps to close a scope gap all
three lens-reviews converged on):**

- **Process/convention fit:** repo state (`main` at `81fba4a`, worktree at
  `6cce6d3`, both clean) confirmed directly. One blocking finding:
  `analyze.R`'s original description ("the `column = "jaccard_similarity"`
  assignments, the `short_col` lookup table") doesn't correspond to any
  actual literal assignment in the file and misses the `d_metric_order`/
  `metric_available` discovery-gate mechanism and five downstream filter
  sites — folded in above. One non-blocking observation: the original
  "Commit: one" rationale was thinner than Steps 1/1b/1c's stated split
  rationales, given several files mix category (b) and (d) in the same
  file — addressed by the revised 4-commit split above.
- **Risk/safety:** four blocking findings, all folded in above — the
  logging requirement wasn't concrete for loop-heavy read sites (now
  addressed with an explicit per-file granularity note);
  `plot_sequence_tissue_clustering.R:59`'s `DEFAULT_SIM_COL` site was
  missing from the coupled-edit instruction; the Step-2-created transient
  inconsistency window (updated generators vs. not-yet-regenerated tracked
  CSVs) had no accepted-residual-risk note (now added); and `analyze.R`'s
  actual re-emit mechanism (the discovery gate, not a flat per-CSV
  fallback) fails silently, independently confirmed by this pass too.
- **Scope completeness:** five blocking findings, the largest-impact round —
  two structurally-identical `DEFAULT_SIM_COL`-before-load sites
  (`plot_sequence_tissue_clustering.R:59`, `plot_cross_reference_identity.R:42`)
  the plan never mentioned; `analyze.R`'s reading logic (not just its
  re-emit side) left unaddressed; one entire generator script,
  `experiments/bedtools_benchmark/estimator_compare.py`, missing from the
  inventory outright; and `plot_interval_accuracy.R:154,171`'s actual code
  shape (an AND-style required-columns check and a bare NSE reference)
  not fitting the plan's proposed drop-in pattern. All five folded in
  above.
- **Follow-up investigative sweeps** (launched to close the scope-
  completeness finding that Step 2's text never enumerated concrete sites
  in `experiments/ref-comparison/`, `subB_mixed_stride/`,
  `synthetic_evolution/`, the rest of `bedtools_benchmark/`,
  `primate-phylogeny/`, or `mus-homo/`, despite Setup declaring those
  directories in-scope): found `experiments/bedtools_benchmark/sweep.py`
  (a second missed generator, in scope via `docs/paper_outline.md`'s
  Figure 3B/S8 provenance), all 3 remaining `ref-comparison/` files (in
  scope via `docs/paper_outline.md` citations — the same citation standard
  already used to find `docs/scripts/` during Setup), a pre-existing
  silent-degradation bug in `run_sweep.py`'s already-present OR-chain
  fallback, and — resting on the same two-hop `docs/paper_outline.md` →
  `docs/estimator-analysis-findings.md` §9.6 citation chain this doc already
  relies on elsewhere — two new in-scope cross-species scripts,
  `primate-phylogeny/estimator_ie_topology.py` +
  `primate-phylogeny/scripts/build_phylogeny.R`, and
  `mus-homo/estimator_ie_tissue.py` + `mus-homo/scripts/cluster_plot.R`
  (the last reached only through a Snakemake-param indirection invisible to
  a literal grep of the R script itself). The same sweep also confirmed
  several candidate files are genuinely **not** provably in scope
  (`run_ie_subb.py`/`analyze_ie_subb.py`, `sbatch_synthetic_ie.sh`,
  `synthetic_evolution/code/analyze.R`, `precision_probe.sh`,
  `compute_column_comparison.{R,py}`) and recorded them above as explicit
  out-of-scope calls rather than silent omissions.

Because these findings changed Step 2's site list substantially (multiple
new files, one new generator directory, two new experiment groups, several
structural — not just line-number — corrections), per "Review process"'s
own re-run rule this is squarely in "materially altered the diff" territory,
not the "purely mechanical" carve-out Steps 1b/1c's minor corrections
qualified for. **A second, focused review round follows below before
implementation proceeds.**

**Review gate round 2, run 2026-08-12 (3 fresh parallel passes against the
revised plan text above):**

- **Process/convention fit:** repo state confirmed (`main` at `b871ddb`
  going in, worktree still untouched at `6cce6d3`, both clean). Round 1's
  own classification (materially altered the diff, not mechanical) was
  independently re-judged correct — the re-run itself is warranted. One
  blocking finding, since fixed above: `exp_a_validate_plot.R`'s
  `docs/paper_outline.md:69,148` citation is factually wrong (that
  location is a different figure entirely) — corrected, and the fix
  mechanism moved from an underspecified Snakefile edit to matching the
  already-established `cluster_plot.R` in-script-resolution precedent. One
  process observation, addressed by this write-up itself: this round's
  findings needed to actually be recorded in the doc, not just referenced —
  done here.
- **Scope completeness:** the deep site citations spot-checked
  (`plot_sequence_tissue_clustering.R`'s five-site chain, `analyze.R`'s
  discovery gate and five filter sites, `estimator_compare.py`/`sweep.py`,
  the ref-comparison and cross-species scripts, all six exclusion calls)
  held up, modulo the mechanical line-drift corrections folded in above
  (`exp_a_metric_comparison.R:36` not `:38`; `build_phylogeny.R:82,86` and
  `:160` not `:82`/`:110-118,160`; `analyze.R`'s `read_hammock_csv`/
  `scan_dir` defaults disentangled to `:160`/`:249`). One blocking gap,
  closed above: three more `experiments/maurano_dhs_validation/` scripts
  (`mode_c_interpolation.R`, `render_dendrogram.R`, `make_metric_plots.R`)
  and `experiments/ref-comparison/config/config.yaml:25`'s unused
  `primary_sim_col` key were never evaluated despite `analyze.R` in the same
  directory getting the deepest treatment in the whole plan — all four
  independently checked against every provenance document and confirmed
  correctly excludable, now recorded as explicit out-of-scope calls rather
  than silent omissions.
- **Risk/safety:** two blocking findings, both closed above. The
  `ref-comparison/workflow/Snakefile:180` fix as originally planned risked
  becoming an unconditional rename with no fallback (a Snakemake `params:`
  block is a static literal, not a function) and had no orchestration-layer
  safety net — closed by dropping the Snakefile edit entirely and resolving
  inside `exp_a_validate_plot.R` itself, the same precedent already used for
  `cluster_plot.R`. The "Accepted residual risk" paragraph both overstated
  one risk (`results/exp_a_estimator_delta.csv` is gitignored, not tracked)
  and missed a real one (7 force-tracked figures under
  `experiments/ref-comparison/figures/` that the edited dendrogram/
  metric-comparison scripts default-write to) — corrected, and the
  verification methodology now requires scratch-path output so the risk is
  closed rather than merely re-accepted. Confirmed clean on the other five
  checks: `cluster_plot.R`'s fix is safe (single consumer of `sim_col`, no
  other downstream re-read by name); zero test coverage anywhere for any
  newly-added file (matches and extends round 1); `precision_probe.sh`'s
  "left unfixed, flagged" framing needed a severity escalation note and a
  named Step-4 exception, both added above; commit 3's bundling concern is
  substantially reduced now that the Snakefile edit is gone (one fewer,
  higher-risk, no-safety-net change mixed into that commit).

Every blocking finding from both round-1 and round-2 has now been folded
into the plan text above, in writing, per this doc's own "resolved in
writing" rule — none deferred as a verbal/chat-only resolution. The round-2
findings themselves are corrections and strengthenings of already-agreed
scope (a citation fix, a fix-location change within an already-in-scope
file, an out-of-scope determination for four more files, a residual-risk
correction, a verification-methodology concretization) — none of them
expand what gets edited beyond what round 1 already settled, so per
"Review process"'s re-run rule this is the "purely mechanical" class and
does not require a third round.

**Step 2 is ready to implement.** Proceeding to the actual edits.

**Commit:** revised from "one" to **four**, split by risk/verification
profile — found more appropriate by this Step's own review gate given the
now much larger and more heterogeneous scope than the original one-commit
plan anticipated (mirrors how Steps 1/1b/1c all split along a stated axis
rather than defaulting to one commit):
1. The `paper/` R-script reading fixes — the coupled tissue-clustering unit
   (`:59/:75/:217/:248/:271/:286`), `plot_cross_reference_identity.R`'s
   ordering fix, `plot_interval_accuracy.R`'s restructuring,
   `plot_parameter_objective_tradeoff*.R`, `plot_estimator_crossover.R`'s
   label constant, and the four `docs/scripts/mode_d_*.R` B2 consumers — the
   most directly manuscript-facing group.
2. `experiments/maurano_dhs_validation/analyze.R` alone, given its size and
   the discovery-gate finding's importance.
3. The other `experiments/` generators — `estimator_compare.py`, `sweep.py`,
   `estimator_ie_crossref.py`, the 3 `ref-comparison/scripts/*.R` files
   (`exp_a_validate_plot.R`'s fix now lives entirely in the R script itself,
   per round 2's correction above — no Snakefile edit in this commit),
   `run_sweep.py`'s OR-chain extension.
4. The cross-species pair — `primate-phylogeny/`, `mus-homo/`, including the
   `cluster_plot.R`/Snakemake-param-chain fix.

Each commit message records its own per-file verification-script-run
evidence per the paragraph above.

**Step 2 landed as four commits on the worktree branch** (`13dd1d2` —
`paper/` R-script reading fixes + the four `docs/scripts/mode_d_*.R` B2
consumers, 10 files; `784117e` — `experiments/maurano_dhs_validation/analyze.R`
alone, all 6 named sites plus the `short_col` key; `1154644` — the other
`experiments/` generators (`estimator_compare.py`, `sweep.py`,
`estimator_ie_crossref.py`, the 3 `ref-comparison/scripts/*.R` files
including `exp_a_validate_plot.R`'s in-script fix, `run_sweep.py`'s OR-chain
extension), 7 files; `ed0b0fb` — the cross-species pair
(`primate-phylogeny/`, `mus-homo/`), 4 files). 22 files changed total
(494 insertions, 69 deletions), matching the plan's four-way split exactly.
Each commit's own message records its detailed per-file verification
evidence (real archived pre- and post-Step-1 data used wherever it existed
in this worktree, hand-built fixtures only where none did; Snakemake
`script:`-directive files — `exp_a_validate_plot.R`, `cluster_plot.R`,
`build_phylogeny.R` — verified via a stubbed minimal S4 `snakemake` object
rather than a full DAG run; all generated output written to scratch, never
clobbering tracked figures/CSVs, confirmed via `git status` after every
run). Two pre-existing, out-of-scope bugs were hit and routed around during
verification setup, not fixed (not this Step's job): `analyze.R`'s
`compute_metrics_vs_ref`/`compute_cluster_metrics` erroring on
single-row/small-sample inputs, and `estimator_ie_crossref.py`'s `main()`
dividing by zero on a single-cell peak type.

Then re-reviewed by 3 fresh adversarial subagents against the **actual
landed diff** (correctness line-by-line, blast radius/scope, and
verification-gap hunting) — same pattern as Steps 1/1b/1c's
post-implementation rounds, run at the user's explicit request for this
Step too. **All three came back clean — no blocking findings.**

- **Blast radius/scope: fully clean.** Independently confirmed the diff
  touches exactly the 22 planned files, each commit's contents match its
  own stated scope, `jaccard_similarity_ie` is untouched everywhere, no
  tracked output file was regenerated, nothing on the "explicitly out of
  scope" list was touched, and a fresh repo-wide sweep found no in-scope
  site the diff missed.
- **Correctness: no blocking bugs.** Fallback logic verified live (via
  extracted-logic toy tests, not just reading) to correctly prefer
  `reg_eq_similarity` and fall back only when genuinely absent; log-once
  guards fire exactly once even across multi-file/multi-iteration loops;
  the coupled `analyze.R` (5 sites) and `plot_sequence_tissue_clustering.R`
  (6 sites) groups consistently reuse one resolved value throughout; the
  Snakemake `script:`-directive fixes correctly treat the injected literal
  as a fallback-only default. Three non-blocking findings, all fixed in a
  follow-up commit (see below).
- **Verification-gap hunting: no blocking gaps.** Independently re-ran 5+
  files across all four commits against real data in both branches and
  reproduced every claimed result; cross-checked the two Snakemake
  `script:`-directive stubs against the real `Snakefile` rule definitions
  line-by-line and confirmed the stub approach was genuinely equivalent
  (both scripts' `snakemake@...` reads matched what the real rules pass,
  with a live `snakemake -n` dry-run cross-check); confirmed via `git log
  -1` on sampled tracked figures/CSVs that none were touched. Three
  non-blocking findings, all addressed below.

**Fixed as a follow-up commit, `0ba31e9`** (mirroring Step 1c's own
one-line follow-up, `6cce6d3`), covering the three non-blocking findings
from the correctness and verification-gap-hunting passes (the blast-radius
pass had none):
- `analyze.R:606,614,628,635,641` — five stale `jaccard_similarity`
  mentions in a comment header, a printed `cat()` diagnostic, and three
  block comments, left behind describing code this Step had already
  renamed. The printed diagnostic now interpolates `REG_EQ_COL` directly
  rather than hardcoding a name. Prose/print-string-only — the underlying
  computed values were already correct, confirmed by the review's own
  mixed-header re-run.
- `run_sweep.py`'s `parse_hammock_csv` OR-chain used per-row truthy checks
  (`if r.get("reg_eq_similarity"):`), inconsistent with every other Python
  file this rename touched (which use presence checks) and, in principle,
  able to treat a present-but-blank field as absent and silently
  substitute a different estimator for just that one row — reproduced live
  by the review, not just theorized. Restructured to resolve the candidate
  column once from the CSV header, before the row loop, matching the
  "resolve once, log once" pattern used everywhere else in this Step; a
  blank value for the resolved column now correctly skips that malformed
  row instead of silently falling through. Verified live against three
  synthetic cases (reg_eq present, legacy-only, reg_eq-present-but-blank),
  all producing correct behavior.
- A stray gitignored `Rplots.pdf`, left at the worktree root by the
  cross-species Snakemake-stub verification runs (R's default PDF device
  firing despite explicit `CairoPNG()` calls elsewhere) — harmless, never
  tracked, but cleaned up rather than left.

**Not fixed, judged non-blocking and out of scope, recorded rather than
silently dropped:** the verification-gap pass's finding that `analyze.R`'s
`short_col` table's second recode key is technically unreachable given the
script's own data flow (harmless dead code, not a correctness bug — the
stated justification for keeping it still holds for the key that *is*
reachable); and commit `1154644`'s message overstating that no real
archived CSV existed to verify `estimator_compare.py`'s fallback branch
end-to-end (a real one did exist and the reviewer built and ran it,
confirming already-correct behavior — a commit-message accuracy note, not
a code defect).

**Then stop.** Report and wait for the user's go-ahead before starting
Step 3. Do not continue automatically.

### Step 3 — Prove parity: regenerate every affected figure/stat, diff against baseline

Covers the user's step 3, and is the load-bearing safety check for Steps 1-2
— a pure rename should produce **zero** numeric difference anywhere, since
no computation changed, only which key labels which value.

1. On a clean checkout of current `main` (pre-Step-1), regenerate every
   figure (PNG), stats output (CSV / printed table), **and the category-B2
   generated CSVs** (`docs/data/mode_d_summary.csv` and the three
   `paper/*_stats.csv` files) that Setup flagged in-scope. Save as a
   baseline snapshot (worktree or scratch space, not committed to `main`).
2. After Steps 1-2 land in the worktree, regenerate the same
   figures/stats/generated-CSVs from the same source data.
3. Diff:
   - PNGs: expect byte-identical or pixel-identical output.
   - Numeric stats outputs and the category-B2 CSVs: expect **exact**
     equality except for the renamed label/column values themselves — no
     floating-point computation changed, this is a lookup-key/label swap.
   - **Diff the raw per-pair hammock CSVs directly, not only the downstream
     aggregated figures/stats** — found necessary by review: an aggregating
     statistic (a correlation, a mean, a clustering result) can in principle
     absorb a single corrupted or mis-column-assigned cell without visibly
     moving, whereas a raw-CSV diff, column by column, surfaces it directly.
4. If any diff shows a numeric change beyond the expected label/name swap,
   **stop** — treat it as a bug in Steps 1-2, not something to wave through.
   The whole basis for treating this as a cosmetic rename depends on this
   coming back clean.

**Review gate:** 3 reviewers examine the diff *results*, not just the plan
— at least one should independently regenerate one figure/stat themselves
rather than trust the report, since this Step is itself acting as the
review/proof mechanism for Steps 1-2.

**Commit:** append a parity-proof record to this seed doc (what was
regenerated, the diff methodology, the outcome — clean, or any accepted
residual noted) and commit it to `main` directly, per "Commit discipline"
above — this Step produces no code/CSV/doc-content changes of its own, only
evidence that Steps 1-2's changes did.

**Then stop.** Report the parity results and wait for the user's go-ahead
before starting Step 4 — this Step's sign-off matters more than most, since
it's the safety proof the whole rename rests on. Do not continue
automatically.

### Step 4 — Remove the duplicate `register_equality_similarity` column, then close out

Covers the user's step 4, plus the version bump and doc sync the user asked
be done alongside it (not its own separate numbered step, but folded in here
since it's the natural closeout once the rename is proven safe).

**Column removal** (same two emitting-code files as Step 1, different
lines): `re` shape shrinks from `query, reference, reg_eq_similarity,
register_equality_similarity` (4 cols) to `query, reference,
reg_eq_similarity` (3 cols); `full` shrinks from 8 to 7, dropping the
trailing duplicate. Concretely:
- `runner._metrics_shape`/`_metrics_row_values` — delete the
  duplicate-reuse code path (its docstring describes the duplicate
  explicitly; rewrite, don't just delete around it).
- `hammock_cli.cpp`'s `MetricsMode::Full` **and** `MetricsMode::RegisterEquality`
  branches — **both** need their duplicate per-pair write removed, not just
  one (a gap in an earlier draft of this doc, found by review): the `Full`
  branch's `cell[7] = reg_jac;` (~513, `stride` 8→7) **and** the
  `RegisterEquality` branch's `cell[0] = cell[1] = qsk[i]->jaccard_similarity(*rsk[j]);`
  (~480, `stride` 2→1) are the same class of duplicate write and both must
  shrink. Check every site that depends on column count or position — the
  `matrix` sizing and the write loop are both derived directly from
  `stride` (`hammock_cli.cpp:429-431` and the serial write loop), so they
  self-correct once `stride` is right; the actual manual risk is entirely in
  the literal `cell[N] = ...` assignments inside the `#pragma omp parallel
  for collapse(2)` block — get both of those right and the rest follows.
- **Add a row-length assertion to `tests/test_hammock_cpp_metrics.py`**,
  mirroring the one `tests/test_metrics_flags.py` already has (found missing
  by review): `test_hammock_cpp_metrics.py` currently reads both outputs via
  `csv.DictReader`, which silently tolerates a row with the wrong number of
  fields (extra values go to the `None` restkey, missing ones get
  `restval`) — a header/`stride` desync in the C++ output would not be
  caught by this test's existing column-by-column checks alone. Assert
  `len(row) == len(header)` per row, matching `test_metrics_flags.py:76-78,96-99`.
- Update `tests/test_hammock_cpp_metrics.py` / `tests/test_metrics_flags.py`
  column-*count* assertions (this Step's actual count change, distinct from
  Step 1's name-only edits to the same files) in the same commit, and
  re-grep the remaining 6 Scope-E test files for anything this Step's count
  change affects.
- **Race-condition stress check, found necessary by review:** the
  duplicate-write bug above, if left unfixed, is a data race under
  `#pragma omp parallel for` — deterministic within one thread's static
  chunk but capable of cross-thread corruption at chunk boundaries,
  depending on `OMP_NUM_THREADS` at run time. A single regenerate-and-diff
  pass (as in Step 3) could pass by luck on one run's scheduling. Before
  calling this Step done, regenerate a `--metrics` and a `--register-equality`
  output at 2-3 different thread counts (e.g. 1, 4, and whatever the default
  resolves to) and diff each against a single-threaded reference — this is
  a cheap, targeted check, not a repeat of Step 3's full figure-regeneration
  apparatus.

Confirm no in-scope Step-2 script actually reads `register_equality_similarity`
as a value source (Setup's table should already show this — none did in the
initial sweep) before deleting it; if one does, repoint it to
`reg_eq_similarity` first.

**Closeout** (same Step, once the column removal is in, its stress check is
clean, and tests are green):
- `CLAUDE.md`: update the three-shape table, the `--register-equality`
  prose block (column counts: `re` now 3, `full` now 7), divergence #2's
  estimator note, the "Measured cost of the metrics block" section's prose
  (numbers unaffected — pure rename/dedup doesn't change timing). Add a
  **dated changelog-style note** recording the rename in the relevant prose
  section (matching the "SHIPPED" framing divergence #9 uses) — **not** an
  "Open seeds" bullet (see Scope F's correction above for why, and the
  caution about the precedent's own such note never actually landing —
  don't repeat that here; this note is part of "Done when," below).
- `README.md`, `docs/jaccard-definitional-gap.md`,
  `docs/estimator-analysis-findings.md`, `docs/paper_outline.md`,
  `docs/submittability-concerns.md`: plain find-and-replace to current fact,
  no deprecation-notice framing (per the user's "only user of hammock"
  clarification).
- `docs/seed-metrics-column-restructure.md`: leave the historical record
  intact; add one line noting it's superseded by this seed for the column
  *name* (the shape/count structure it defined otherwise still holds, minus
  one duplicate column).
- **Version bump**: this changes the output-column contract of the `re` and
  `full` shapes (name and count both). CLAUDE.md doesn't state a single
  explicit numbered policy for this — the v0.7.0/v0.8.0 precedents were
  bumped minor by analogy/pattern, not by a quoted general rule (an earlier
  draft of this doc mis-cited CLAUDE.md as containing a literal "minor bumps
  mark CLI-contract/default changes" policy statement; it doesn't — that's
  a paraphrase of the pattern, not a quote). Worth reviewers weighing
  explicitly: this rename is arguably **more** disruptive than the v0.8.0
  precedent it's being compared to, since v0.8.0 was purely additive (new
  shapes/columns; nothing consumer-facing was renamed or removed), whereas
  this change removes the name `jaccard_similarity` entirely from 2 of 3
  shapes — precisely why Step 2's fallback logic has to exist at all. Given
  the project is pre-1.0 and already treats schema-changing minor bumps as
  normal (per the v0.7.0/v0.8.0 history), a minor bump — **0.8.0 → 0.9.0**
  — is the leaning here, not a claim of settled policy; confirm per Open
  Question 4. Both front-ends' `--version` must reflect it in the same
  commit (benchmark harnesses gate on it).
- Rebase onto current `main` (it may have moved since Setup — check
  `git log main`), re-run `pytest tests/` on the rebased branch, merge per
  `memory/feedback_work_on_main.md`, remove the worktree.

**Review gate:** 3 reviewers check the column-removal diff (both `cell[]`
write sites, the row-length assertion, the thread-count stress check) and
the doc updates for staleness (grep for any surviving bare
`jaccard_similarity` describing current behavior, confirm the CLAUDE.md
changelog note actually landed rather than being deferred) before merge.

**Commit:** two, once the review gate clears, plus the merge itself — (1)
the column removal (both emitting-code files, both `cell[]` sites, the
`stride`/count changes, the row-length-assertion test addition, the
remaining Scope-E test fixes, the thread-count stress-check results noted
in the message), (2) the closeout (CLAUDE.md, README, the other docs listed
above, the version bump) — kept separate from (1) because the risky
code-count change and the purely-textual closeout are different things to
review and, if needed, revert independently. Then rebase, re-run
`pytest tests/`, and **merge with `git merge` (not squash)** so this Step's
two commits — and Steps 1-2's four — all stay individually visible in
`main`'s history, per "Commit discipline" above.

**Done when:** `pytest tests/` is green on `main` post-merge; `cpp/tests/`
green; `pyproject.toml` reads `0.9.0`; `git worktree list` no longer shows
the rename worktree; CLAUDE.md carries a dated note of this change (not just
an Open Seeds bullet); `grep -rn 'jaccard_similarity\b' | grep -v '_ie'`
across the repo returns only: (a) `test_parity_against_original.py`'s
orig-comparison logic (frozen external tool, always emits this name); (b)
Step 2's backward-compat fallback string literals in paper/experiments code
(`else "jaccard_similarity"` and equivalents — these are *supposed* to
survive, they're what makes old archived CSVs still parse); (c) historical
seed/CLAUDE.md/memory entries explicitly marked as superseded/historical;
(d) `experiments/primate-phylogeny/scripts/precision_probe.sh:106` — Step
2's review gate deliberately left this one literal unfixed (not provably in
scope, see that Step's write-up), so its survival here is expected, not a
leftover to chase down; and (e) nothing else describing it as the current
column *or method or binding* name — after Step 1b, `grep -rn
'jaccard_similarity(' cpp/ bindings/` (the call-syntax form, not the
CSV-string form) should return nothing at all, since that Step removes the
identifier everywhere it was a name rather than a string.

**Then stop and report.** Step 4 ends in a merge to `main`, so there's no
further Step to hold off on — but still stop here and report the final
state (the "Done when" checklist, results) to the user rather than treating
merge-to-main as a silent, self-closing action.

## Decisions made (resolved with the user, 2026-08-12 — no longer open)

All six items below were open questions in an earlier draft of this doc.
Resolved directly with the user via explicit choice, before Step 1 starts,
per the same "don't resolve unilaterally" policy that flagged them in the
first place. Kept here as the record of *why*, matching this doc's own
stated "resolved in writing" requirement — not a template to re-litigate.

1. **Rename the C++ method `HLLSketch::jaccard_similarity()` too?**
   **Resolved: yes.** Initial leaning was out-of-scope (bigger blast radius
   for a purely internal name); reversed once it was established the repo
   is public on GitHub (`origin` = `jessicabonnie/hammock`, `main` already
   synced there — checked, not assumed) and the CSV rename doesn't actually
   decide this either way (a direct C++/pybind11 caller never sees the CSV
   column name). Given that, and the user's actual stated goal being broader
   than the CSV ("make sure the register equality calculation isn't named
   jaccard misleadingly on anything someone might be using"), this expanded
   further mid-resolution to also cover the Python-facing
   `estimate_jaccard` binding — see Step 1b for both parts. New names:
   `reg_eq_similarity()` (C++ method), `estimate_reg_eq_similarity` (Python
   binding, the user's explicit preference over matching orig's
   `estimate_jaccard_registers`).
2. **`test_parity_against_original.py`'s fix shape.** Not re-opened as a
   choice — the doc's existing plan (Scope E caveat: map by column
   *position*, not name, since orig is frozen and will always emit
   `jaccard_similarity`) stands as the approach. Step 1's review gate is
   where this actually gets checked for real, not this list.
3. **Archived-CSV sweep completeness (delimiters, `experiments/*/results/`,
   the data-value sweep).** Not a decision, an instruction already actionable
   in Setup's body text — no user choice needed here, just confirmed as part
   of Setup's scope rather than carried as a separate open item.
4. **Version bump.** **Resolved: minor, 0.8.0 → 0.9.0**, as the doc's
   original leaning — user confirmed directly, no change to Step 4's plan.
5. **Setup gating.** **Resolved: stays ungated**, as the doc's original
   leaning — user confirmed directly, no change to "Review process."
6. **Category-B2 handling** (`docs/data/mode_d_summary.csv` and the three
   `paper/*_stats.csv` files). **Resolved: update generating code in Step 2,
   regenerate + diff in Step 3**, as the doc's original leaning — user
   confirmed directly, no change to Steps 2/3's plan.

## Setup inventory (completed 2026-08-12, worktree `../hammock_claude_wt_regeq` / branch `jaccard-reg-eq-rename`)

Pre-flight: `squeue -u $USER` returned no running jobs — no conflicting writer
to `docs/data/` or any `paper/*/results` directory. Worktree created clean
off `main` at `f4cf813`.

Re-ran every sweep in "Scope / current-state inventory" from scratch (grep
sweeps, not a re-read of the snapshot). **Result: the snapshot above is
accurate and complete — no file was added to or removed from any category.**
Specifics below, organized by the doc's own taxonomy, each row a direct
re-verification (not a repeat of the snapshot's prose).

**A. Emitting code.** Confirmed exact sites: `runner.py:38` (`_metrics_shape`
"re" tuple), `:34-35` ("full" tuple), `_metrics_row_values`/docstring
`:46-47`, `_print_estimator_note` `:128-145`; `cli.py:264-273` help strings;
`hammock_cli.cpp` — three header `fprintf`s now at `:414,420,423` (line
numbers shifted slightly from the snapshot's ~414/420/423 estimate but same
three lines), `--help` text `:140-173`, comments `:32,409`, and the two
per-pair value-write sites `:480` (`RegisterEquality` branch,
`cell[0]=cell[1]=...jaccard_similarity(...)`) and `:494`
(`reg_jac = ...jaccard_similarity(...)`, feeding the `Full` branch's
`cell[7]` a few lines later) — both confirmed untouched by Step 1, both
correctly in Step 1b/Step 4's scope, not this Step's.

**B. Archived CSV headers.** Re-swept **every tracked `.csv` under `paper/`,
`docs/data/`, `experiments/`** (recursive `git ls-files`, 23 files total, not
just a top-two-levels sample) for a header-row token match on both `,` and
`\t` delimiters. Exactly the same 5 files match, no more, no less:
`docs/data/exp_a_narrow_k10_w10.csv`, `docs/data/exp_a_broad_k10_w10.csv`,
`docs/data/hammock_hll_p18_jaccB_full.csv`,
`docs/data/hammock_hll_p21_jaccB_full.csv`,
`docs/data/hammock_hll_p23_jaccB_full.csv`. Also checked every tracked
`.tsv`/`.txt` under the same trees (30 files) — none carry a
`jaccard_similarity` header; `docs/data/maurano_bedtools_ref.tsv`'s `jaccard`
column is bedtools' own output field, a different name, correctly out of
scope. **Collision hazard reconfirmed by direct header read**: both
`exp_a_narrow_k10_w10.csv` and `exp_a_broad_k10_w10.csv` literally contain
`...,jaccard_similarity,containment_AB,containment_BA,cosketch_geom,
cosketch_arith,cosketch_max,jaccard_similarity_with_ends,
containment_AB_with_ends,...` in that header — a non-word-boundary
substitution would corrupt the `_with_ends` twin, exactly as the snapshot
warned.

**B2. Data-value CSVs.** `grep -c` across all 23 tracked CSVs for the literal
string in the file body: the same 4 files, same order of magnitude as the
snapshot — `docs/data/mode_d_summary.csv` (940), `paper/estimator_crossover/
estimator_crossover_stats.csv` (24), `paper/interval_accuracy/
interval_accuracy_stats.csv` (6), `paper/sequence_tissue_clustering/
estimator_agreement_stats.csv` (2). No other tracked CSV in any of the three
trees carries the string as a data value. Also spot-checked the specific
`docs/data/*.csv` files that Figure 3/S6/S7/S8/S9/S10's provenance
paragraphs cite as inputs (`cpp_vs_bedtools_t16_p18.csv`,
`cpp_vs_bedtools_t16_p18_largeN.csv`, `sweep_threads_p18.csv`,
`sweep_precision_maurano_p18_t{8,16}.csv`, `maurano_subB_ie_summary.csv`,
`maurano_subB_re_vs_bedtools.csv`, `maurano_subB_summary.csv`,
`maurano_bedtools.csv`) individually — **0 occurrences in every one of
them**; they're timing/MAE tables with no `jaccard_similarity` column at
all (`maurano_bedtools.csv` has a bare `jaccard` field, bedtools' own name,
not ours). So none of the benchmark-figure provenance chains (Figures 3,
S6-S10) pull in any B2-style file; only Figures 4, 6, 7 do, via the four
files above — consistent with the snapshot.

**C. Paper R scripts.** Re-confirmed every line number the snapshot cites,
including the two the first-draft-correction flagged: `plot_interval_accuracy.R:91`
(`EST_RE` label) *and* `:154` (`required_hammock` column-presence check,
distinct from `:171`'s value read) both present as described;
`plot_parameter_objective_tradeoff_estimators.R:60`
(`SIM_COLUMNS <- c("jaccard_similarity", "jaccard_similarity_ie")`) present
alongside `:98,132`; `plot_sequence_tissue_clustering.R:271`
(`intersect(c("jaccard_similarity","jaccard_similarity_ie"), names(raw))`)
and `:286` (`agreement$column == "jaccard_similarity"`) both present and
confirmed coupled — `:271` builds the working set, `:286` filters a
*different* hardcoded literal against it, exactly the silent-`NA`-if-only-
one-edited hazard the snapshot flagged.

**D. `experiments/` + `docs/scripts/` scripts.** Recursive `.py`/`.R`/`.sh`
grep under `experiments/` returns **exactly 34 files** (the snapshot's
count, re-derived independently, not copied): `maurano_dhs_validation/`,
`primate-phylogeny/`, `bedtools_benchmark/`, `mus-homo/`, `modeD_flanking/`
(correctly untouched — archival), `ref-comparison/` (4 files:
`estimator_ie_crossref.py`, `scripts/exp_a_dendrogram.R`,
`scripts/exp_a_metric_comparison.R`, `scripts/exp_a_validate_plot.R`),
`subB_mixed_stride/` (6 files: `analyze_ie_subb.py`, `run_ie_subb.py`,
`run_sweep.py`, `sbatch_maurano.sh`, `sbatch_synthetic_ie.sh`,
`summarize_synthetic_ie.py`), `synthetic_evolution/` (1 file:
`code/analyze.R`). `docs/scripts/` has exactly 4 of its 7 files matching —
`mode_d_violins.R`, `mode_d_lines.R`, `mode_d_metric_tradeoff.R`,
`mode_d_bedtools_vs_modeB_scatter.R` — the other 3
(`maurano_speed.R`, `maurano_speed_table.R`, `synthetic_nscaling.R`) do not
reference the column at all. Cross-checked against `paper/outline.md` and
`paper/draft.md`'s `**Figure N generation.**` provenance paragraphs directly
(not inferred): Figures 4, 6, 7 (and via `docs/paper_outline.md`'s
figure-index table, Figures 3/4/6/S4 of *that* doc's own numbering) are the
only ones whose cited inputs touch a `jaccard_similarity` column; Figures
3/S6/S7/S8/S9/S10's cited inputs are the pure-timing files confirmed clean
in B2 above.

**No dynamically-constructed column names found.** Grepped for
`paste0`/`sprintf`/f-string patterns that build the string
`"jaccard_similarity"` rather than writing it as a literal — none found; the
one `f"jaccard_similarity ..."` hit
(`experiments/ref-comparison/estimator_ie_crossref.py:139`) is a print-message
prefix, not a column-name construction, and the file is already correctly
in-scope via category D.

**E. Tests.** Re-grepped all 8 files; confirmed each has at least one
literal, non-`_ie` hit needing a fix, and traced what each actually reads.

**Correction (post-gate finding, see "Step 1 review gate round 1" above): the
paragraph originally here claimed `test_parity_against_original.py`'s
`_projected_rows`/`test_jaccard_byte_equal` needed no code change. That
claim was wrong, caught by the scope-completeness reviewer, not by this
Setup pass — recorded here as a correction, not silently rewritten, per this
doc's own "resolved in writing" rule.** The original reasoning (`keep` is
built independently per side by filtering `_PROJECTED_OUT` by name, then
compared positionally) is correct as far as it goes, but misses that
`_projected_rows`'s `return` line iterates `for line in lines` over **all**
of `lines = csv_text.strip().split("\n")`, including `lines[0]`, the header
row itself. So the first tuple in the returned list is the header's own
*column names* at the kept positions, not data — and `test_jaccard_byte_equal`
compares `orig_rows == ours_rows` including that first tuple. Today, both
sides' header token is the literal string `jaccard_similarity`, so that
tuple happens to match and the test passes; that was never intentional
(nothing in the function's name, docstring, or the surrounding tests
suggests the header should be part of the comparison), it just never
mattered before because both sides used the same literal name. Once Step 1
renames only our side's header token to `reg_eq_similarity`, that first
tuple diverges (`"reg_eq_similarity" != "jaccard_similarity"` at the same
index) and all 7 `test_jaccard_byte_equal` parametrizations fail on the
header alone, before a single data value is even compared. **Not
catchable by Step 1's own "grep for a literal hit and fix" remediation** —
the string `jaccard_similarity` appears in this function only in surrounding
comments/docstring, never in the comparison logic itself, so a grep-driven
pass walks past it. **Fix, now folded into Step 1's plan (see below): change
`_projected_rows` to iterate `lines[1:]`, dropping the header row from the
comparison** — this is exactly what "map by column position, not name"
(Open Question 2) requires: two same-position data columns compare equal
regardless of what either side's header calls them, and header-name equality
was never supposed to be part of the contract.

The real fixes elsewhere in this file are the two sites that look up *our
own* header by the literal name to find a column index for a value read
(`test_mode_b_subB_actually_subsamples:176`,
`.split(",").index("jaccard_similarity")` against **our** CSV) — those must
become `"reg_eq_similarity"`. Same pattern (reading our own emitted header by
name) is the fix needed in `test_mode_d.py:115`, `test_bed2fasta_cli.py:72`,
`test_jaccard_ie.py:241`, and the column-list/index literals in
`test_metrics_flags.py:29,106,107,127` and
`test_hammock_cpp_metrics.py:29,215,223`.

**Second correction, same source:** `test_mode_d_parity.py` was described
above as needing "prose/comment updates only." Also wrong, also caught by
the scope-completeness reviewer. `test_mode_d_structural_parity`
(`:110-111,116-117`) scans **our own** freshly-emitted `--metrics` header
dynamically by prefix match, not by an exact literal:
```python
sim_cols = [i for i, c in enumerate(header)
            if c.startswith(("jaccard_similarity", "containment", "cosketch"))]
...
sym_cols = [i for i in sim_cols
            if header[i].startswith(("jaccard_similarity", "cosketch"))]
```
After Step 1's rename, the renamed column no longer starts with
`"jaccard_similarity"`, so it silently drops out of both `sim_cols` (the
range/self-pair check) and `sym_cols` (the symmetry check) — the test keeps
passing, just quietly checking 6 similarity columns instead of 7. Literal
`grep`-and-substitute would also get this wrong a different way (substituting
the tuple's `"jaccard_similarity"` entry outright would then stop matching
`jaccard_similarity_ie`, which must still match). **Fix, folded into Step 1:
widen both tuples to `("jaccard_similarity", "reg_eq_similarity",
"containment", "cosketch")` and `("jaccard_similarity", "reg_eq_similarity",
"cosketch")` respectively** — keeps matching `jaccard_similarity_ie` via the
unchanged `"jaccard_similarity"` prefix and restores coverage of the renamed
column.

`test_containment_estimator.py` still only needs a prose/comment update (no
literal or dynamic lookup found there beyond the comment at `:72,187`).

**F. Docs.** The 6 doc files the snapshot names
(`CLAUDE.md`, `README.md`, `docs/jaccard-definitional-gap.md`,
`docs/estimator-analysis-findings.md`, `docs/paper_outline.md`,
`docs/submittability-concerns.md`) plus `paper/methods.md`, `paper/draft.md`,
`paper/outline.md` are confirmed to have live, current-fact prose describing
the column (`README.md` in particular documents the column contract at
length — lines 6, 46, 53, 60, 209, 213, 353, 364, 392 — answering the
snapshot's open "verify whether it documents columns at all" with: yes,
extensively, it is genuinely in Step 4's scope, not a light touch).

**Found by this re-sweep, not called out in the snapshot: a broader set of
*historical/status* docs also match the grep** (seed docs other than this
one, dated plan/results notes, experiment READMEs) —
`docs/bedtools-baseline-retraction.md`, `docs/figure3-panel-a-rebuild.md`,
`docs/metrics-by-default.md`, `docs/mode-d-ends-removal.md`,
`docs/part6-plan-draft.md`, `docs/part6-plan-v2.md`, `docs/part6-plan-v3.md`,
`docs/precision-frontier-maurano.md`, `docs/seed-benchmark-methodology.md`,
`docs/seed-mode-d-hash-width.md`, `docs/seed-mode-d-threading.md`,
`docs/seed-subsampling-synthetic-supplement.md`, and per-experiment
`README.md`/`RESULTS.md`/`docs/experiment_design.md` files under
`experiments/bedtools_benchmark/`, `experiments/maurano_dhs_validation/`,
`experiments/modeD_flanking/`, `experiments/mus-homo/`,
`experiments/primate-phylogeny/`, `experiments/ref-comparison/`,
`experiments/subB_mixed_stride/`, `experiments/synthetic_evolution/`.
**Explicitly marked out of scope for every Step in this plan**, by the same
reasoning already applied to `docs/seed-metrics-column-restructure.md`
itself: these are dated records of what was true or what was measured at
the time, not living descriptions of the current CLI contract — rewriting
them to say `reg_eq_similarity` would misrepresent what those historical
runs actually emitted. None of Steps 1-4 touch this list. Recorded here so
a future pass doesn't mistake the omission for an oversight.

**Conclusion: no scope changes.** The plan's Steps 1-4 as written are ready
to execute against exactly the files already named in "Scope /
current-state inventory" above. Proceeding to Step 1's review gate.
