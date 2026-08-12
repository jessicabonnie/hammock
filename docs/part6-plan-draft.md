# Part 6 implementation plan (DRAFT — under adversarial review, not yet executed)

Scope: `docs/seed-metrics-column-restructure.md` Part 6 — update dependent
`experiments/`/`paper/` scripts for the three-shape metrics contract Parts
2-4 already landed (`28f53cd`/`c7bdba1`/`89febf5` on this worktree,
`metrics-restructure` branch). New contract: no flag → `jaccard_similarity_ie`
only (tag `_ie`); `--register-equality`/`--re` → `jaccard_similarity` +
`register_equality_similarity` duplicate (tag `_re`), replaces the removed
`--no-metrics`; `--metrics` → full 8-column block (tag `_full`), unchanged
flag name/meaning.

**User-confirmed scope decision (this session):** for benchmark harnesses in
`experiments/bedtools_benchmark/` and `experiments/subB_mixed_stride/` whose
`--metrics-arm`/`--metrics-all`/second-pass machinery exists solely to obtain
`jaccard_similarity_ie` (previously only reachable via the full block): **mark
with `# PART9:` comments only, do not retarget flags or re-run anything in
Part 6.** Part 9 owns the actual switch to the bare default + fresh paired
measurement, per the seed doc's existing constraint.

Inventory was built by four independent read-only agents (bedtools_benchmark,
subB_mixed_stride, paper/pairwise_scaling R, repo-wide bare-invocation sweep)
and cross-checked against live `grep -n` output before this plan was written.
Line numbers below are verified against the current worktree state as of this
session, not just agent-reported.

## Category 1 — Mechanical rename `--no-metrics` → `--register-equality` (literal CLI flag strings passed to hammock-cpp)

1. `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py`
   - Line 383 (`harness_floor_ms()`): `"--no-metrics"` → `"--register-equality"`.
   - Line 506 (`run_hammock()`): `cmd += ["--metrics"] if metrics else ["--no-metrics"]` → `[..., "--register-equality"]`. This is the shared function `fusion_ab.py`, `sweep.py`, `pairwise_cost_by_precision.py` all import — fixing it here propagates to all four call sites without touching those files' own logic.
2. `experiments/subB_mixed_stride/run_sweep.py`
   - Line 173: `"--metrics" if metrics else "--no-metrics"` → `[..., "--register-equality"]`.

Verified via direct `grep -n` that no other file in `experiments/`/`paper/`
constructs a literal `"--no-metrics"` argument to pass to a subprocess — every
other `--no-metrics` grep hit in these two directories is prose/comment
(Category 2) or an internal Python variable/dict-key label unrelated to the
actual CLI flag (`fusion_ab.py`'s `"no_metrics"` arm-identifier strings,
`pairwise_cost_by_precision.py`'s same pattern — see "Explicitly not
touching," below).

## Category 2 — Prose/docstring corrections (describes CURRENT behavior wrongly, not historical record)

Scoped to comments/docstrings that assert something false about the *live*
contract going forward — not to comments recording what a past, already-run
job actually used (those are left alone, matching `RESULTS.md`/`docs/seed-*.md`
treatment). Two flavors: (a) literal `--no-metrics` string that will confuse a
reader since the flag no longer exists, (b) description of `--metrics`
*relative to the default* that's now wrong because the default changed
meaning (flagged explicitly by the user this session — the flag name is
unchanged but its relationship to "no flag" flipped).

1. `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py`
   - Line 69: version-gate error string `"This harness passes --no-metrics and parses microsecond timings, "` → update flag name.
   - Lines 479-484 (`run_hammock()` docstring): *"metrics selects the output shape explicitly: True passes --metrics (the binary's default since 0.7.0), False passes --no-metrics."* — both clauses are wrong now (full block is no longer the default; False now means `--register-equality`). Rewrite.
   - Lines 1131-1144 (`--metrics-arm`/`--metrics-all` argparse help text): "WITHOUT --no-metrics" / "instead of --no-metrics" / "replaces the default no-metrics arms" — all reference the removed flag and the stale "default = full block" assumption. Rewrite together with the Category 5 PART9 marker at the same location (single edit, not two passes).
2. `experiments/bedtools_benchmark/sweep.py`
   - Lines 124-131 (`parse_hammock_csv()` docstring) and line 143 (error message `"Re-run hammock-cpp without --no-metrics."`): rewrite column-count/flag-name claims against the actual current header shapes (verify exact column counts from `runner.py`/`hammock_cli.cpp` at edit time, not from memory).
   - Lines 491-492, 853, 859-866 (`--metrics-arm` docstring/help): same treatment as item 1's third bullet, combined with Category 5's marker C.
3. `experiments/bedtools_benchmark/measure_cli_overhead.py`
   - Lines 13-16: *"The hammock CLI has no --metrics/--no-metrics opt-out... it always emits the full 9-column block... so hammock-cpp here is run WITH --metrics for an apples-to-apples comparison."* Premise is false post-restructure (CLI now has the flag; bare CLI no longer emits the full block). Rewrite to describe why `--metrics` is now added explicitly to *both* arms (see Category 4).
4. `experiments/bedtools_benchmark/fusion_ab.py`
   - Lines 7, 31-38, 228: docstring/printed-output text literally names `--no-metrics` as a real flag ("the `--no-metrics` arm runs byte-identical code," the four-row ASCII arm diagram, `"CONTROL -- post/pre on --no-metrics."` print statement). Rename these specific prose references to `--register-equality`.
   - **Explicitly not touching:** the internal Python arm-identifier strings (`"no_metrics"`/`"metrics"` used as dict keys / printed column labels at lines 173, 192, 233-234, 243, 247-248) — these are this script's own internal vocabulary for labeling the two arms in its own output, not a hammock CLI flag, and renaming them would change this script's own CSV/printed-table schema for no behavioral reason. Flagged for reviewers to confirm this line is drawn in the right place.
5. `experiments/bedtools_benchmark/pairwise_cost_by_precision.py`
   - Line 5: docstring `` `--no-metrics` is the opt-out used `` → update to `--register-equality`. Same "leave internal `no_metrics`/`metrics` arm labels alone" reasoning as item 4.
6. `experiments/subB_mixed_stride/run_sweep.py`
   - Lines 59, 66-69, 113-118, 169-172, 229, 264-269, 291: multiple comments/help-text/print-statements describing "the 3-column --no-metrics arm." Rewrite flag name and column count together (verify exact new register-equality column count at edit time — spec says 2 metric columns, `jaccard_similarity` + `register_equality_similarity`, so likely 4 total incl. `query`/`reference`, not 3; confirm against `outprefix.py`/`hammock_cli.cpp` rather than asserting a number here).
7. `experiments/bedtools_benchmark/sbatch_fig3_largeN.sh` (~lines 39-41) and `sbatch_fig3_panelA.sh` (~lines 33-38): *"Both arms are kept (--no-metrics and the default 9-column block)"* — "the default" is no longer the 9-column block; rewrite to name the two arms correctly (`--register-equality` baseline + explicit `--metrics-arm`).
8. `experiments/subB_mixed_stride/sbatch_maurano.sh`: rename the line-78 `--metrics` forward is unaffected (flag unchanged), but the stale comment at lines 28-29 (*"Both legs are --no-metrics throughout"*) already contradicts the script's own body (line 78 runs `--metrics`) even before this restructure — fix the rename **and** correct this pre-existing inaccuracy while in the file.

**Explicitly deferred / left alone (historical record, not live guidance):**
`RESULTS.md` in both `bedtools_benchmark/` and `subB_mixed_stride/` (dated,
archival, describes what flag *was* passed to already-completed jobs);
`sbatch_synthetic_ie.sh` lines 16-17 (describes a pre-existing CSV predating
the `--metrics`/IE feature entirely — a historical fact, not current
guidance); all of `paper/pairwise_scaling/*.R` (see Category 6 — the
dedicated read-only agent found essentially nothing needing a change there).

## Category 3 — Add `--metrics` to bare (no-flag) call sites that read named similarity/containment columns as scientific data (not benchmark timing)

Per the seed doc's Part 6 item 4 default policy: preserve exactly the columns
already depended on. These are accuracy/science pipelines, not the
"benchmarks primarily concerned with register-equality/default" the user
scoped Category 5 to — they read `containment_AB`/`containment_BA`/
`cosketch_*` columns that exist in no shape but `--metrics`, so `--register-equality`
would not suffice even if it were cheaper.

1. `experiments/bedtools_benchmark/estimator_compare.py` — line 113-114 `cmd = [...]`, add `"--metrics"`. Reads `containment_AB`, `containment_BA`, `jaccard_similarity` (lines ~249-251).
2. `experiments/subB_mixed_stride/run_ie_subb.py` — line 45-47 `cmd = [...]`, add `"--metrics"`. Hard-checks for `jaccard_similarity`, `jaccard_similarity_ie`, `containment_AB`, `containment_BA` in the header (lines ~58-61), `raise SystemExit` if any missing.
3. `experiments/maurano_dhs_validation/run_sweep_abc.py` — line 154-161 `cmd = [...]`, add `"--metrics"`. Feeds `analyze.R`, which reads the full 7-metric block by name (`d_metric_order`, `analyze.R:325-327`).
4. `experiments/maurano_dhs_validation/run_sweep_d.py` — line 140-148 `cmd = [...]`, add `"--metrics"`. Same downstream (`analyze.R`).
5. `experiments/maurano_dhs_validation/workflow/Snakefile` — four rules, each missing the flag on its `{HAMMOCK} ...` shell line: `run_mode_ab` (~132), `run_mode_c_expA` (~153), `run_mode_c_subB` (~175), `run_mode_d` (~198). Add `--metrics` to each. Same `analyze.R` downstream, reached via a fully independent Snakemake path from items 3/4 (not visible to a `.py`-only audit).
6. `experiments/primate-phylogeny/scripts/precision_probe.sh` — line 47, the `sbatch --wrap=` invocation. Add `--metrics`. **Currently a silent-wrong-column bug, not a loud failure**: the script's own `awk` does `h["jaccard_similarity"]` (lines ~73-94, ~100-118); under the bare default that key is absent, `jcol=""`, and `$jcol`/`$0`-style awk semantics silently coerce to `0` rather than erroring — highest-priority fix in this category.
7. `experiments/primate-phylogeny/workflow/Snakefile` — `rule phylo_hammock`, line ~172. Add `--metrics`. Strongest dependency found in the whole audit: `scripts/build_phylogeny.R` needs `cosketch_max`/`cosketch_geom` (lines ~160-162), which exist in no shape but `--metrics` — not a judgment call.
8. `experiments/mus-homo/workflow/Snakefile` — `rule mus_homo_hammock`, line ~129. Add `--metrics`. `config/config.yaml`'s `primary_sim_col: "jaccard_similarity"` is read via `cluster_plot.R`.
9. `experiments/ref-comparison/workflow/Snakefile` — `rule exp_a_hammock`, line ~153. Add `--metrics`. `sim_col = "jaccard_similarity"` hardcoded in the Snakefile, read via `exp_a_validate_plot.R`.
10. `experiments/synthetic_evolution/Makefile` — four recipes (Mode A ~172, Mode B ~186, Mode C no-expA ~202, Mode C expA ~219). Add `--metrics` to each `$(HAMMOCK) ...` line. `code/analyze.R` requires `jaccard_similarity` or legacy `jaccard` (hard `stop()` if neither present).

**Note on uniform `--metrics` (not `--register-equality`) even where only
`jaccard_similarity` is read** (items 6, 8, 9): the seed doc explicitly
discourages "a case-by-case guess at which of the 8 columns is actually
used" for bare accuracy-reading call sites — `--metrics` is the
column-preserving default across this whole category, deliberately not
optimized per-script. This is called out for reviewers since it's a real
choice, not the only defensible one.

**Explicitly NOT touched (archival, already broken for an unrelated,
pre-existing reason, per the seed doc's own "archival, lower priority"
framing):** `experiments/modeD_flanking/run_sweep_synthetic.py` and
`value_demo.py` — both hard-require `jaccard_similarity_with_ends`, a column
removed in v0.6.0 (divergence #8), so both already `raise SystemExit`/fail
against any hammock build newer than that, independent of this restructure.
Adding `--metrics` would not make either script functional again, so it's
not a productive use of a Part 6 edit.

## Category 4 — Bare CLI calls in genuine timing scripts: add matching `--metrics` to preserve the historical measurement's apples-to-apples premise

1. `experiments/bedtools_benchmark/measure_cli_overhead.py` — line 102-104, `cli_cmd` (Python `hammock` CLI, currently bare). `cpp_cmd` (line 99-101) already hardcodes `--metrics`. Add `"--metrics"` to `cli_cmd` too — this restores exactly the historical comparison (the CLI used to always emit the full block with no flag; now it needs the flag to do the same), not a benchmark redesign. Combine with the Category 2 docstring fix and a Category 5 PART9 marker noting this pairing could eventually test register-equality/default arms instead.
2. `experiments/bedtools_benchmark/measure_threads_isolation.py` — line 96-98, same pattern, same fix, same PART9 marker.

## Category 5 — `# PART9:` markers (leave `--metrics` as-is; comment only, per user's confirmed decision this session)

Every one of these exists specifically to obtain `jaccard_similarity_ie` via
the (now unnecessarily expensive) full block, in a **benchmark/timing**
script — the class of call site the user's framing this session scoped to
"primarily concerned with `--register-equality` and default IE."

1. `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py` — `--metrics-arm`/`--metrics-all` (argparse defs ~1131-1144, `arms_for()`/`ie_tool_name_for_subb()` ~599-636). One marker at the argparse definitions.
2. `experiments/bedtools_benchmark/sweep.py` — the untimed second `--metrics` pass on the precision axis (~409-440, `if args.ie:` block) — the cleanest case: reads **only** `jaccard_similarity_ie` from its CSV and nothing else.
3. `experiments/bedtools_benchmark/sweep.py` — `--metrics-arm` for the threads/intervals axes (argparse ~859-866, `IE_ARM_TOOL` usage ~535/581).
4. `experiments/bedtools_benchmark/measure_cli_overhead.py` and `measure_threads_isolation.py` — at the (newly-explicit, per Category 4) `--metrics` on both arms; note these two arms could eventually be redesigned around register-equality/default instead of the full block.
5. `experiments/subB_mixed_stride/run_sweep.py` — at line 173 / the `--metrics` argparse definition (~264-269). **Caveat noted in the marker itself**: this arm's CSV also captures `containment_AB`/`containment_BA` (not just `jaccard_ie`) per the parse function at lines ~104-130 — Part 9 should verify those are actually unused downstream before assuming a clean retarget; if some consumer does use them, this one may need to stay on `--metrics` rather than drop to the bare default.

No SLURM re-runs, no flag value changes, no filename/output changes in this
category — comment-only.

## Category 6 — Confirmed no action needed

- `paper/pairwise_scaling/plot_largeN_supplement.R`, `plot_precision_frontier.R`, `plot_threading_supplement.R` — zero references to any flag/tag terms; read pre-existing hardcoded-path CSVs only.
- `paper/pairwise_scaling/plot_pairwise_scaling.R`, `plot_pairwise_scaling_supplement.R` — every `--no-metrics`/`--metrics-arm` reference is either (a) a dated historical description of what a past run used, or (b) a reference to `benchmark_cpp_vs_bedtools.py`'s own `--metrics-arm`/`--metrics-all` wrapper flags, which are not renamed by this restructure. No live filename/glob dependency on the old `_j3` tag exists in either file (both use hardcoded `file.path(...)`, no `Sys.glob`).
- `paper/pairwise_scaling/plot_subB_nscaling_supplement.R` — exists on `main`, not yet in this worktree (added after the branch point). Confirmed via `git show main:...` that its references are the same two categories as above (historical + unrelated wrapper-flag name); no action needed now, will be picked up automatically at the Part 8 rebase.
- `experiments/bedtools_benchmark/RESULTS.md`, `experiments/subB_mixed_stride/RESULTS.md` — dated, archival, describe already-completed jobs by job ID; matches the `docs/seed-*.md` "don't edit historical prose" treatment even though not literally under `docs/`.
- `experiments/bedtools_benchmark/estimator_rank_by_precision.py` — never invokes hammock; reads only `estimator_compare.py`'s own derived-column CSV (`j_register_equality`/`j_incl_excl`), insulated from the raw hammock header rename.
- No `_j3` hardcoded filename/glob reference exists anywhere in `experiments/bedtools_benchmark/`, `experiments/subB_mixed_stride/`, or `paper/pairwise_scaling/` (verified by direct grep, zero hits) — Part 6 item 3 (filename-glob fixes) is a non-issue for every file this plan touches. Output-CSV paths are captured dynamically at runtime (parsed from stderr/stdout `Wrote ...` lines or `run_hammock()`'s own return dict), never reconstructed from the old suffix.

## README.md updates (live docs, in scope per Part 6 — not Part 7, which is CLAUDE.md/root README/design docs)

- `experiments/bedtools_benchmark/README.md` — "Conventions"-style section describing the old two-shape contract needs the new three-shape contract.
- `experiments/subB_mixed_stride/README.md` — same.

## Done-criteria checklist (from the seed doc, mapped to this plan)

- [ ] `grep -rn -- '--no-metrics'` under `experiments/`/`paper/` returns nothing outside `RESULTS.md`/dated-historical prose explicitly listed above as deferred.
- [ ] Every live script reading named similarity/containment columns (Category 3) passes `--metrics` explicitly.
- [ ] Every genuine timing script (Category 4) carries the flag matching what it measures.
- [ ] `grep -rn 'PART9:'` shows exactly one marker per flagged call site (Category 5) — 5 files, ~7 markers as scoped above.
- [ ] Full `pytest tests/` still green after these edits (none of them touch `tests/`, `python/`, or `cpp/`, so this is a regression check, not an expected-change area).

## Open questions for reviewers

1. Is the Category 2 line drawn correctly between "internal arm-label strings" (left alone: `fusion_ab.py`/`pairwise_cost_by_precision.py`'s `"no_metrics"` dict keys) and "prose asserting a real CLI flag" (fixed)? Could a downstream reader be misled by the untouched internal labels?
2. Category 3 applies `--metrics` uniformly even to scripts that only read `jaccard_similarity` (items 6, 8, 9) rather than `--register-equality`, which would be cheaper and sufficient today. Right call, or should those three specifically get `--register-equality` instead, accepting the "case-by-case column guess" the seed doc otherwise discourages?
3. Any bare call site missed? The repo-wide sweep covered `.py`/`.sh`/Makefiles/Snakefiles specifically because the seed doc's own methodology couldn't see those formats — worth one more pass over other non-`.py`/`.sh` formats (e.g. any `.nf`/Nextflow, `.smk` includes, CI config) before considering Category 3/4 complete?
4. Category 5's scope (mark, don't retarget) was confirmed with the user this session for the *benchmark* directories specifically — confirm no Category-3 item (accuracy pipelines) was miscategorized into needing a PART9 marker instead of an outright `--metrics` fix, or vice versa.
