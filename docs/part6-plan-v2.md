# Part 6 implementation plan v2 (revised after round-1 adversarial review; under round-2 review, not yet executed)

Supersedes `docs/part6-plan-draft.md` (v1). Scope unchanged: `docs/seed-metrics-column-restructure.md`
Part 6 — update dependent `experiments/`/`paper/` scripts for the three-shape
metrics contract Parts 2-4 already landed (`28f53cd`/`c7bdba1`/`89febf5`,
`metrics-restructure` branch). Contract: no flag → `jaccard_similarity_ie`
only (tag `_ie`); `--register-equality`/`--re` → `jaccard_similarity` +
`register_equality_similarity` (tag `_re`); `--metrics` → full 8-column
block, i.e. `jaccard_similarity`, `jaccard_similarity_ie`, `containment_AB`,
`containment_BA`, `cosketch_geom`, `cosketch_arith`, `cosketch_max`,
`register_equality_similarity` (tag `_full`).

## What changed from v1, and why (round-1 review + user decisions)

1. **Critical correctness bug found and fixed**: `python/hammock/outprefix.py`
   now unconditionally tags every Python-CLI output filename (`_ie`/`_re`/`_full`)
   — new behavior from Part 2, the Python CLI never tagged before. v1's "just
   add `--metrics`" fix for 8 bare call sites would have broken them (Snakemake
   `MissingOutputException`, Make silently never re-running, R regex `stop()`
   or silent skip). v2 fixes the filename expectation everywhere it's
   hardcoded, not just the flag.
2. **`fusion_ab.py` fix corrected**: v1's shared-`--no-metrics`-rename in
   `run_hammock()` would have broken `fusion_ab.py`'s primary use case (A/B
   testing archival pre-0.8.0 `hammock-cpp` binaries that don't understand
   `--register-equality`). v2 version-gates the flag choice per binary
   (user-confirmed direction: "Version-gate the flag in run_hammock()").
3. **Pre-existing bugs unrelated to the rename, found during review, now
   in scope** (user-confirmed: "Fix them in Part 6"): `estimator_compare.py`'s
   already-broken hardcoded output path; four scripts silently selecting a
   stale file via `sorted`/unsorted `glob(...)[0]` once old- and new-format
   CSVs coexist (a live hazard, not hypothetical — this worktree already has
   both).
4. **Category 4 folded into Category 5** (spec-fidelity finding + user
   clarification): the seed doc's own item 4 says bare timing/benchmark call
   sites should NOT default to `--metrics` — v1 did this for
   `measure_cli_overhead.py`/`measure_threads_isolation.py`, contradicting the
   spec. These are now PART9-marked instead, with the stale "apples-to-apples"
   docstring corrected to honestly describe the now-mismatched arms rather
   than claim a fix that isn't Part-6-safe.
5. **PART9 marker content made concrete** (user-confirmed, stated three times
   this session): every benchmark call site currently paying for the full
   `--metrics` block solely to obtain `jaccard_similarity_ie` gets a marker
   stating the exact target — **drop the flag entirely; bare/no-flag is the
   IE arm now** — not a vague "reconsider this," and not `--register-equality`
   (that's for the OTHER arm, the mechanically-renamed one).
6. **Stale pre-restructure result files (348 total) get renamed, not left
   orphaned or regenerated** (user-confirmed, and independently verified safe
   by reading actual downstream column usage — see Category 3a below).

## Category 1 — Mechanical rename `--no-metrics` → `--register-equality`, version-gated

1. `experiments/subB_mixed_stride/run_sweep.py` line 173: `"--metrics" if metrics else "--no-metrics"` → `"--register-equality"`. This binary is always the current build (no archival-binary use case here, verified — `run_sweep.py` takes one `--binary` argument, never two), so a straight rename is safe.
2. `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py` — **not a blind rename**. Full design (verified against actual `run_hammock()`/`check_binary_version()` source, all call sites checked):
   - Add `REGISTER_EQUALITY_MIN_VERSION = (0, 8, 0)` (this restructure's version per Part 7's planned bump, `pyproject.toml` currently `0.7.1`).
   - Add `_probe_binary_version(binary)`, `lru_cache`-wrapped, refactored out of the existing `check_binary_version()`'s inline parsing (so the one-time startup gate and every later `run_hammock()` call against the same path share one cache entry/one exec).
   - Add `_metrics_off_flag(binary) -> str`: returns `"--register-equality"` if the probed version `>= (0, 8, 0)`, else `"--no-metrics"`.
   - `run_hammock()` line 506: `cmd += ["--metrics"] if metrics else [_metrics_off_flag(binary)]`. No signature change — `binary` was already `run_hammock`'s first positional parameter, so no call-site changes needed in `fusion_ab.py`, `sweep.py`, or `pairwise_cost_by_precision.py`.
   - This correctly handles `fusion_ab.py`'s dual-binary loop (`builds[build]` alternates `pre`/`post` per call) — each gets the flag vocabulary its own vintage understands.
   - Full diff verified ready-to-implement (docstrings, `check_binary_version` refactor, caching rationale) — apply as designed.

## Category 2 — Prose/docstring corrections (unchanged from v1, still needed)

Same scope as v1: fix comments/docstrings that assert something false about
the *live* default (not historical-record prose, which stays untouched).
Full list: `benchmark_cpp_vs_bedtools.py` (version-gate error string, `run_hammock()`
docstring, `--metrics-arm`/`--metrics-all` help text — this last one now also
gets its Category 5 PART9 marker in the same edit), `sweep.py` (`parse_hammock_csv()`
docstring + error message, `--metrics-arm` docstring/help), `fusion_ab.py`
(prose naming `--no-metrics` literally, not the internal `"no_metrics"` arm-label
strings — leave those, they're this script's own vocabulary, not a CLI flag),
`pairwise_cost_by_precision.py` (docstring), `run_sweep.py` (subB_mixed_stride,
multiple comments — verified column counts this time: `--register-equality`
is 4 total columns not 3, bare default is 3 total not 9, `--metrics` is 10
total not 9 — confirmed via `outprefix.py`+`runner.py` read, not guessed),
`sbatch_fig3_largeN.sh`/`sbatch_fig3_panelA.sh` ("the default 9-column block"
claim), `sbatch_maurano.sh` (rename + a pre-existing, unrelated stale-comment
fix noted as incidental).

Also new in v2: `measure_cli_overhead.py`'s docstring (lines 13-16) gets
rewritten to honestly state the current mismatch (cpp arm explicit `--metrics`,
CLI arm bare/IE-only) rather than describe a fix that isn't happening in
Part 6 — paired with its Category 5 marker.

## Category 3 — Bare Python-CLI call sites reading named columns as data: add `--metrics` AND fix the filename-tag propagation

Every item below needs **both** parts, verified via direct source read (not
estimated):

### 3a. Existing stale-data files: rename to add `_full` tag (do this first, no code change)

Verified safe by reading actual downstream column usage, not assumed:
- `analyze.R` (maurano)'s ABC scan (`scan_dir(abc_dir, parse_abc_name, refs=..., do_clustering=FALSE)`, line 279) never overrides `jcols`, so it uses the default `jcols = "jaccard_similarity"` **only** — confirmed by reading `scan_dir`'s signature (line 246) and the call site. The old `raw_abc/*.csv` files (header: `...,jaccard_similarity,containment` — pre-2026-05-14 placeholder format) already satisfy this.
- `analyze.R`'s Mode D path has a **pre-existing fallback** (lines 161-169) that reconstructs `jaccard_similarity_ie` from `containment_AB`/`containment_BA` when the `_ie` column itself is absent — built for exactly this situation. The old `raw_d/*.csv` files have `containment_AB`/`containment_BA`/`cosketch_*` (confirmed via header read) but predate `jaccard_similarity_ie` and still carry the removed `_with_ends` family (harmless extra columns, ignored by name-based readers).
- `experiments/synthetic_evolution/code/analyze.R` only reads `jaccard_similarity`/legacy `jaccard` (line 128) — confirmed via grep, no containment/cosketch/IE dependency at all.

So: rename (not regenerate, not leave orphaned) —
- `experiments/maurano_dhs_validation/results/raw_abc/*.csv` → append `_full` before `.csv` (29 files)
- `experiments/maurano_dhs_validation/results/raw_d/*.csv` → append `_full` before `.csv` (235 files)
- `experiments/synthetic_evolution/results/hammock_g100000000/*.csv` → append `_full` before `.csv` (84 files)

(`experiments/ref-comparison/`'s one "untagged .csv" hit is `config/nfcore_samplesheet_human.csv` — a sample sheet, not a hammock output. Confirmed, not touched. `primate-phylogeny`/`mus-homo` have zero existing result CSVs — pipelines haven't been run in this worktree yet, nothing to rename.)

This makes the "require `_full`" regex fixes below safe immediately, with no
regeneration needed and no orphaned files.

### 3b. Code fixes, per file (flag + filename-tag propagation together)

1. **`experiments/maurano_dhs_validation/workflow/Snakefile`** — five pattern-string locations must gain `_full` before `.csv`, kept in sync (Snakemake requires wildcarded `output:` to literally match): the `expand()`/list-comprehension pattern strings (`AB_CSVS` line 54, `C_EXPA_CSVS` line 57, `C_SUBB_CSVS` line 60, `D_CSVS_BASE` line 64, `D_CSVS_EXT` line 66) AND the four rules' own `output:` declarations (`run_mode_ab` line 123, `run_mode_c_expA` line 144, `run_mode_c_subB` line 166, `run_mode_d` line 188). `rule analyze`'s `input:` (line 208-211) references `ABC_CSVS`/`D_CSVS` by variable, so it inherits the fix automatically — no separate edit. Add `--metrics` to all four rules' `shell:` blocks (after `--threads {threads}`, before the outprefix arg, matching the existing arg order style).
2. **`experiments/primate-phylogeny/workflow/Snakefile`** — `rule phylo_hammock`'s `output.csv` (line 161) and `rule build_phylogeny`'s `input.matrix` (line 190) are two separate hardcoded strings with the identical pattern (`{OUTDIR}/{MARK}/k{{k}}_w{{w}}/phylo_mnmzr_p{PRECISION}_jaccD_k{{k}}_w{{w}}.csv`) — both need `_full` added, kept in sync. Add `--metrics` to the `shell:` block (line 172-179).
3. **`experiments/mus-homo/workflow/Snakefile`** — same pattern: `rule mus_homo_hammock`'s `output.csv` (line 118) and `rule cluster_and_plot`'s `input.matrix` (line 144), both need `_full`. Add `--metrics` to the `shell:` block (line 129-136).
4. **`experiments/ref-comparison/workflow/Snakefile`** — same pattern: `rule exp_a_hammock`'s `output.csv` (line 141) and `rule exp_a_validate_and_plot`'s `input.matrix` (line 167), both need `_full`. Add `--metrics` to the `shell:` block (line 153-161). (`sim_col = "jaccard_similarity"` at line 176 needs no change — present in `_full`.)
5. **`experiments/synthetic_evolution/Makefile`** — `HAMMOCK_AB` (lines 82-83, feeds the `%`-pattern rules at lines 165 and 179 — **both the variable AND the two pattern-rule target declarations need `_full` added, kept in sync**, since Make's `%` pattern matching requires the literal target string to match); `HAMMOCK_C_DEF` (line 86) and `HAMMOCK_C_EXP` (line 89) are consumed via **static pattern rules** keyed directly off the variable (lines 194, 210) — only the variable needs the edit, no separate rule-declaration line exists to duplicate it. The internal `sed`-based rate/type/expA extraction in the Mode-C recipes (lines 197-198, 213-215) uses unanchored `.*` trailing patterns — verified these still work correctly with `_full` appended (traced manually, `.*` absorbs it). Add `--metrics` to all four hammock invocations (lines 172-176 Mode A, 186-190 Mode B, 202-206 Mode C-default, 219-223 Mode C-expA). The grouped analysis rule (line 226, `$(FIGURE_TARGETS) $(SUMMARY_TARGETS) &: $(HAMMOCK_CSVS) ...`) needs no separate edit — it depends on `$(HAMMOCK_CSVS)` (line 90), which inherits the fix from `HAMMOCK_AB`/`HAMMOCK_C_DEF`/`HAMMOCK_C_EXP`.
   - Bonus effect: renaming the target patterns (not just adding the flag) also resolves the mtime-staleness hazard round-1 flagged — Make will no longer find the 84 old untagged files satisfying the new `_full`-tagged target names, so it can't silently skip regenerating. Combined with 3a's rename, existing files simply become the up-to-date `_full`-tagged targets directly — no regeneration needed, no staleness risk either.
6. **`experiments/synthetic_evolution/code/analyze.R`** — `parse_meta`'s regex (line 102): append `_full` before the closing `\\.csv$` anchor. Exact verified diff:
   ```r
   # before
   "^synth_r([0-9.]+)_([ABC])_hll_p(\\d+)_jacc([ABC])(?:_expA([0-9.]+))?(?:_B([0-9.]+))?\\.csv$"
   # after
   "^synth_r([0-9.]+)_([ABC])_hll_p(\\d+)_jacc([ABC])(?:_expA([0-9.]+))?(?:_B([0-9.]+))?_full\\.csv$"
   ```
   Same capture-group indices, no other code changes. Verified against concrete example filenames (both old-shape-turned-`_full` and Mode-C-with-expA-and-subB cases) — matches correctly, rejects untagged names (which no longer exist post-3a rename anyway).
7. **`experiments/maurano_dhs_validation/analyze.R`** — three edits, all verified:
   - `modeB_path` literal (line 143): `hammock_hll_p21_jaccB.csv` → `hammock_hll_p21_jaccB_full.csv`.
   - `parse_abc_name`'s regex (line 229): append `_full` before the closing `$` (operates on extension-stripped stem, so no `.csv` in the pattern) — `...(?:_B(\\d+\\.\\d+))?$` → `...(?:_B(\\d+\\.\\d+))?_full$`.
   - `parse_d_name`'s regex (line 240): `"^hammock_mnmzr_p(\\d+)_jaccD_k(\\d+)_w(\\d+)$"` → `"^hammock_mnmzr_p(\\d+)_jaccD_k(\\d+)_w(\\d+)_full$"`.
   - Confirmed no fourth hardcoded hammock-filename literal/regex exists elsewhere in this file (full-file grep for `.csv` literals and every `grepl`/`regmatches`/`sub`/`gsub`/`str_match*` call performed).
8. **`experiments/bedtools_benchmark/estimator_compare.py`** — add `"--metrics"` to `cmd` (lines 113-114). Path construction (line 127): `f"{prefix}_hll_p{precision}_jaccB.csv"` → `f"{prefix}_hll_p{precision}_jaccB_full.csv"`. (Considered a glob-based alternative matching `run_ie_subb.py`'s pattern; hardcoding is simpler and this function only ever runs one metrics mode, so no future-flexibility need — going with the hardcode.)
9. **`experiments/subB_mixed_stride/run_ie_subb.py`** — add `"--metrics"` to `cmd` (lines 45-47). No filename-tag fix needed — already glob-based (`glob.glob(prefix + "*")`), confirmed tag-agnostic.
10. **`experiments/primate-phylogeny/scripts/precision_probe.sh`** — add `"--metrics"` to the `sbatch --wrap=` invocation (line 47). Filename-selection fix (line 68): `ls "${cell_dir}"/probe*.csv | head -1` → `ls "${cell_dir}"/probe*_full.csv | head -1` (also fixes the ASCII-sort stale-pick hazard, since narrowing the glob leaves at most one match regardless of sort order — the awk block only reads `jaccard_similarity`, present in `_full`). The separate archived-baseline path read at lines 97-119 (a pre-existing, untouched `/vast/...` literal, predates this restructure entirely) is out of scope — different issue, not a live-run bug.

## Category 4 — (removed; folded into Category 5 per spec-fidelity finding + user confirmation)

`measure_cli_overhead.py`/`measure_threads_isolation.py` no longer get an
automatic `--metrics` fix. See Category 5.

## Category 5 — `# PART9:` markers (mark precisely; do not retarget flags or re-run, per user's repeated confirmation this session)

Every marker states the exact target concretely — **the arm that measures
inclusion-exclusion should become the bare/no-flag default**, not a vague
placeholder, per the user's explicit spec this session:

1. `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py` — `--metrics-arm`/`--metrics-all` (argparse ~1131-1144, `arms_for()`/`ie_tool_name_for_subb()` ~599-636). Marker + the Category 2 docstring rewrite happen together (help text currently says "instead of --no-metrics," which needs both a flag-name fix and a "this whole mechanism should switch to bare default" note).
2. `experiments/bedtools_benchmark/sweep.py` — the untimed second `--metrics` pass on the precision axis (~409-440, `if args.ie:` block). Verified cleanest case: reads only `jaccard_similarity_ie` (line 429-430, `parse_hammock_csv(..., column="jaccard_similarity_ie")`), nothing else.
3. `experiments/bedtools_benchmark/sweep.py` — `--metrics-arm` for threads/intervals axes (argparse ~859-866; actual `run_hammock(..., metrics=use_metrics)` call sites are ~535/~581, the `IE_ARM_TOOL` symbol itself is defined at ~512/~559 — corrected citation from v1, which conflated the two).
4. `experiments/bedtools_benchmark/measure_cli_overhead.py` and `measure_threads_isolation.py` — **moved here from v1's Category 4**. Marker at both `cpp_cmd`'s existing `--metrics` and `cli_cmd`'s bare invocation, explaining the now-mismatched arms and that Part 9 should decide the correct pairing (likely both at bare/register-equality, not both at `--metrics`) after a fresh measurement — per the seed doc's own item-4 carve-out for timing scripts, not a Part-6 default-policy `--metrics` fix.
5. `experiments/subB_mixed_stride/run_sweep.py` — at line 173 / the `--metrics` argparse definition (~264-269). **Caveat resolved, not just noted**: this arm's CSV also captures `containment_AB`/`containment_BA` (lines ~104-130) in addition to `jaccard_ie` — the marker states this explicitly and instructs Part 9 to verify those are actually unused by any downstream consumer of `run_sweep.py`'s own output CSV before retargeting; if used, this site should NOT drop to bare default and the marker should be removed instead.

No SLURM re-runs, no flag value changes in this category — comment-only, as
confirmed three times by the user this session.

## Category 6 — Confirmed no action needed (unchanged from v1)

Same as v1: `paper/pairwise_scaling/*.R` (all references are historical prose
or the unrelated `--metrics-arm`/`--metrics-all` wrapper-flag name, no `_j3`-style
live filename dependency); `RESULTS.md` files (dated archival); `estimator_rank_by_precision.py`
(never invokes hammock).

## New: pre-existing bugs fixed (unrelated to the rename itself, per user's "fix them in Part 6" confirmation)

1. `experiments/bedtools_benchmark/estimator_compare.py` — covered above (3b.8), already broken today independent of this restructure (bare call already writes a tagged filename that doesn't match the hardcoded path).
2. `experiments/ref-comparison/estimator_ie_crossref.py` (lines 116-120) — `sorted(glob.glob(...))[0]` → narrow the glob itself: `sorted(glob.glob(os.path.join(cell_dir, "*_full.csv")))`. Reads `jaccard_similarity`/`containment_AB`/`containment_BA` (verified, lines 84-85) — needs `_full` specifically, so narrowing (not mtime-sorting) is correct: a newest-by-mtime file of the wrong shape would still crash.
3. `experiments/primate-phylogeny/estimator_ie_topology.py` (lines 143-146) — same fix, same reasoning: `glob.glob(...)` → `sorted(glob.glob(os.path.join(RESULTS, mark, cell, "*_full.csv")))`.
4. `experiments/mus-homo/estimator_ie_tissue.py` (lines 80, 84, two call sites) — same fix: `glob.glob(...)` → `sorted(glob.glob(os.path.join(args.results, cell, "*_full.csv")))`.

All three (2-4) read `containment_AB`/`containment_BA` alongside `jaccard_similarity` — confirmed via source read, `_full` is the only shape that has them, so narrowing the glob (not switching to mtime-descending sort) is the correct fix in every case: it makes the "wrong shape" case impossible rather than just less likely.

## Done-criteria checklist

- [ ] `grep -rn -- '--no-metrics'` under `experiments/`/`paper/` returns nothing outside dated-archival `RESULTS.md` prose.
- [ ] Every live script/pipeline reading named similarity/containment columns (Category 3) passes `--metrics` AND its downstream filename-matching logic expects `_full`.
- [ ] Every genuine timing/benchmark script (Category 5, incl. the two moved from v1's Category 4) is PART9-marked with a concrete target, not fixed with an automatic flag choice.
- [ ] `grep -rn 'PART9:'` shows one marker per Category 5 item (5 files, ~6 markers).
- [ ] The four pre-existing glob/path bugs are fixed regardless of the rename (Category 3b.8 + the 3 new items above).
- [ ] 348 stale result files renamed (`raw_abc` 29, `raw_d` 235, `synthetic_evolution` 84), zero left with the old untagged naming.
- [ ] `fusion_ab.py` still works against archival pre-0.8.0 binaries after the Category 1 fix (version-gated, not a blind rename) — spot check by construction, not just by reading, once implemented.
- [ ] Full `pytest tests/` still green (none of these edits touch `tests/`, `python/`, or `cpp/`, so this is a regression check).

## Open items for round-2 reviewers

1. Is the `_full`-required (not optional) choice for the two R regex fixes (3b.6, 3b.7) correctly justified given 3a's rename makes it safe? Any file this plan doesn't know about that could still produce an untagged CSV in these directories between now and when these regexes run?
2. Double-check the Category 1 `fusion_ab.py` version-gating design for any call site or edge case the design agent's investigation might have missed (e.g. does `benchmark_cpp_vs_bedtools.py`'s own single `run_hammock` call at line ~780 correctly get the NEW binary's vintage, given it's a different code path than `fusion_ab.py`'s loop?).
3. Verify the file-rename step (3a) doesn't collide with anything mid-flight — check no SLURM job or process external to this worktree currently reads `raw_abc/`, `raw_d/`, or `hammock_g100000000/` by their old untagged names (the seed doc's two previously-noted jobs, 29761450/29763124, have already finished per an earlier `squeue` check this session, but re-verify since time has passed).
4. Any remaining bare-call-site or glob[0]-style bug not caught across the three rounds of investigation this session? A fresh sweep specifically for `sorted(glob` / unsorted `glob(...)[0]` patterns repo-wide (not just the 4 already found) would close this out with more confidence.
