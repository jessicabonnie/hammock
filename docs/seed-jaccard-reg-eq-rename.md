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

**Commit:** two, once the review gate clears — (1) emitting-code + test
rename (`runner.py`, `cli.py`, `hammock_cli.cpp`, the two Scope-E test files'
name-only edits plus the other six's fixes), (2) archived-CSV header
rewrites (the 5 files, plus whatever Setup's redone sweep added), kept
separate since the second is pure data-file editing with its own `diff`-based
verification and reviewers may want to inspect it independently of the code
change.

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

**Commit:** two — (1) the C++ method rename (Part 1, all
declaration/definition/call sites), (2) the Python-facing binding rename
(Part 2: the `.def(...)` registration plus its 8 test call sites) — kept
separate since Part 2 is the one with any (small, in-repo-only) API-surface
implication, worth its own reviewable diff. Each commit's message records a
green `cpp/tests/`/`pytest tests/` run as evidence it compiles and behaves
identically.

**Then stop.** Report and wait for the user's go-ahead before starting
Step 2. Do not continue automatically.

### Step 2 — Update paper/experiments code to read `reg_eq_similarity`, with fallback

Covers the user's step 2. Every in-scope script from Setup's table,
categories (b) and (d):

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
handling.

`paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R:271,286`:
edit as **one unit** sharing a single resolved column-name variable, not two
independent literal replacements — see the Scope C caveat for exactly what
breaks (a silent `NA` degradation, not an error) if only one line is
updated.

**Re-emitting (category d, and the category-B2 CSVs it produces):**
`experiments/maurano_dhs_validation/analyze.R`'s literal `jaccard_similarity`
data-value writes (the `column = "jaccard_similarity"` assignments, the
`short_col` lookup table) and the equivalent label constants in
`paper/estimator_crossover/plot_estimator_crossover.R`,
`paper/interval_accuracy/plot_interval_accuracy.R:91`, and
`paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R` (its
own label text, separate from the coupled 271/286 read-logic above) — update
these to emit `reg_eq_similarity` as the value they write into
`docs/data/mode_d_summary.csv` and the three `paper/*_stats.csv` files.
These four generated files are **not** hand-edited (Step 1 explicitly
excluded them); Step 3 regenerates and diffs them instead.

**Reading `mode_d_summary.csv`'s `column` field (category B2), specifically:**
five consumers, not one — `paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R`
plus all four of `docs/scripts/mode_d_violins.R`, `mode_d_lines.R`,
`mode_d_metric_tradeoff.R`, `mode_d_bedtools_vs_modeB_scatter.R` filter on
`column == "jaccard_similarity"`. Apply the same reg_eq_similarity-preferred/
jaccard_similarity-fallback pattern to all five (their `column` field is data
read from the regenerated CSV, not a CSV-level column name, so the fallback
here is a value comparison — `filter(column %in% c("reg_eq_similarity",
"jaccard_similarity"))` or equivalent — not a `names(df)` presence check).

Check whether `paper/` already has a shared sourced-R-utility file before
inventing one just for this helper — if 3+ scripts need it, factor it once;
if not, inline is fine and matches this repo's existing per-script style.

**Review gate:** 3 reviewers confirm the fallback pattern (with logging) is
applied consistently, that the coupled tissue-clustering lines were edited
together, that the category-B2 generating-code updates are complete against
Setup's table, and that no script was missed.

**Commit:** one, once the review gate clears — all in-scope
`paper/`/`experiments/`/`docs/scripts/` consuming-code edits from this Step
together, since they're one coherent "teach the readers about the new name"
change and splitting it further wouldn't add independent verifiability the
way Step 1's code/data split did.

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
and (d) nothing else describing it as the current column *or method or
binding* name — after Step 1b, `grep -rn 'jaccard_similarity(' cpp/
bindings/` (the call-syntax form, not the CSV-string form) should return
nothing at all, since that Step removes the identifier everywhere it was a
name rather than a string.

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
Worth recording for Step 1's review gate (this is the parity-test check the
motivation section asked to be verified for real, not just asserted):
`test_parity_against_original.py`'s `_projected_rows`/`test_jaccard_byte_equal`
comparison is **already position-effective, not name-matching** — it builds
`keep` independently per side by filtering out only the names in
`_PROJECTED_OUT` (which contains neither `jaccard_similarity` nor
`reg_eq_similarity`), then compares the surviving *tuples* positionally. So
once ours emits `reg_eq_similarity` in the same column position orig emits
`jaccard_similarity`, the byte-equality check keeps working with **no code
change needed** to that function — Open Question 2's "map by position, not
name" is already satisfied by the existing implementation, not something
Step 1 needs to newly add. The real fixes needed in this file are the two
sites that look up *our own* header by the literal name to find a column
index for a value read (`test_mode_b_subB_actually_subsamples:176`,
`.split(",").index("jaccard_similarity")` against **our** CSV) — those must
become `"reg_eq_similarity"`. Same pattern (reading our own emitted header by
name) is the fix needed in `test_mode_d.py:115`, `test_bed2fasta_cli.py:72`,
`test_jaccard_ie.py:241`, and the column-list/index literals in
`test_metrics_flags.py:29,106,107,127` and
`test_hammock_cpp_metrics.py:29,215,223`. `test_mode_d_parity.py` and
`test_containment_estimator.py` only need prose/comment updates (no literal
`.index("jaccard_similarity")`-style lookups found in either).

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
