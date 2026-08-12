# Part 6 implementation plan v3 (final — consolidated from two adversarial review rounds; ready for implementation)

Supersedes `docs/part6-plan-v2.md`. Same scope: `docs/seed-metrics-column-restructure.md`
Part 6. Contract unchanged: no flag → `jaccard_similarity_ie` only (`_ie`);
`--register-equality`/`--re` → `jaccard_similarity` + `register_equality_similarity`
(`_re`); `--metrics` → full 8-column block (`_full`).

## What changed from v2 (round-2 review, three reviewers)

**Restored (v2 silently dropped these during consolidation from v1):**
1. `run_sweep_abc.py`/`run_sweep_d.py` — the seed doc's own explicitly-named
   minimum ("At minimum, audit and fix..."). Missing from v2 entirely.
2. `experiments/bedtools_benchmark/README.md` and `experiments/subB_mixed_stride/README.md`
   updates.
3. `harness_floor_ms()`'s second hardcoded `--no-metrics` literal (line 383,
   `benchmark_cpp_vs_bedtools.py`) — separate from the `run_hammock()` fix.

**Corrected:**
4. **`fusion_ab.py`'s version-gating design was broken against the actual
   current binary** — empirically confirmed by building and running it: this
   worktree's `hammock-cpp` reports `--version` `0.7.1` but already had
   `--no-metrics` removed in `89febf5` (Part 4), and nothing has bumped
   `pyproject.toml` past `0.7.1` yet (that's Part 7's job). `0.7.1` is
   ambiguous — it covers both pre- and post-restructure binaries — so no
   semver threshold can distinguish them. **New design: capability probe, not
   version comparison** — check the binary's own `--help` output (cached per
   path, same rationale as before) for the literal string `--register-equality`.
   Sidesteps the Part-7 sequencing dependency entirely and is robust to this
   exact class of "code changed, version number didn't" drift, which this
   repo has already hit once.
5. **Category 3a (the stale-file rename) is deferred to Part 8, not done in
   Part 6.** `experiments/maurano_dhs_validation/results` and
   `experiments/synthetic_evolution/results` are symlinks into shared `/vast`
   storage that `main`'s own checkout resolves to identically — renaming now
   would break `main`'s current (unfixed) scripts until Part 8 merges the
   code fixes too. User-confirmed. The Category 3b **code** fixes (regexes,
   Snakefile/Makefile patterns requiring `_full`) still land in Part 6 as
   planned — they just won't find any matching files in this worktree until
   Part 8 does the rename, which is expected and fine (nothing in Part 6 runs
   these pipelines against live data).
6. **Stale-file count corrected: ≥432, not 348.** A parallel, actively-read
   sibling directory, `experiments/synthetic_evolution/results/hammock_g1000000/`
   (a different `GENOME_RANGE` tag, same live format), has 84 more untagged
   files the count missed. `raw_d_buggy_pre_fix/` (209 files) is correctly
   excluded — documented elsewhere as dead/superseded, unread by any script.
7. Category 5 item 4 (`measure_cli_overhead.py`/`measure_threads_isolation.py`)
   de-hedged: neither script reads any CSV column (verified), so there's no
   real ambiguity for Part 9 to resolve — the marker states the concrete
   target (drop `--metrics` from `cpp_cmd`, matching `cli_cmd`'s already-bare
   state) instead of "likely... not necessarily."
8. Category 5 item 5 (`subB_mixed_stride/run_sweep.py`) resolved, not just
   caveated: traced its one live consumer, `summarize_synthetic_ie.py` — it
   reads neither `jaccard_ie` nor the containment columns (confirmed by full
   read). So this site genuinely qualifies for "drop to bare default" per the
   seed's own item-2 criterion; the marker states this as settled, not as an
   open question for Part 9 to re-derive.
9. Category 5 item 3 (`sweep.py`'s threads/intervals `--metrics-arm`) needs
   **two** markers, not one — `sweep_threads` (~line 535) and
   `sweep_intervals` (~line 581) are structurally identical, ~46 lines apart,
   and the docstrings describe them as deliberately parallel. A single marker
   near the argparse definition risks a future Part 9 pass fixing one twin
   and silently leaving the other on `--metrics`.
10. Two more live consumers of the (Part-8-deferred) renamed files, found by
    a fresh trace neither prior round caught: `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R`
    (generates **manuscript Figure 6**) and `experiments/maurano_dhs_validation/mode_c_interpolation.R`
    (actively documented and run per that experiment's README). Both hardcode
    old untagged filenames/glob patterns against `raw_d`/`raw_abc`. Added to
    Category 3b.
11. A 5th glob-then-`[0]` bug, same defect shape as the four already being
    fixed: `benchmark_cpp_vs_bedtools.py:543-544`'s own `output_csv` capture
    (`csvs = [f for f in glob.glob(out_prefix + "*") if f.endswith(".csv")]; r["output_csv"] = csvs[0] if csvs else None`).
    Added to the pre-existing-bugs list.
12. `fusion_ab.py`'s new unconditional `--version`/`--help` probe (where none
    existed before) needs a try/except around the subprocess call — a bad or
    non-executable `--pre-binary`/`--post-binary` path should raise a clean
    error, not a raw `OSError`/`PermissionError` from inside the new probe.

**Noted, not fixed in Part 6 (flagged for whoever executes Part 8):**
- `sbatch_maurano.sh`/`sbatch_synthetic_ie.sh` pin a specific prebuilt
  `hammock-cpp` binary path on `main`. Part 8's merge plan should explicitly
  rebuild that pinned binary as part of the merge — if the Python-side rename
  lands on `main` before that binary is rebuilt, the next live SLURM run
  sends it a flag it doesn't understand and fails mid-allocation.
- `lru_cache`'s per-path (not per-content) keying means an in-place binary
  rebuild mid-run would serve a stale capability reading for the rest of that
  run. Accepted risk, consistent with this codebase's existing tolerance for
  the same class of hazard elsewhere (`check_binary_version`'s own docstring
  is about a related stale-binary case) — not a new problem introduced here,
  documented in the code comment rather than solved with new invalidation
  machinery.

## Category 1 — Mechanical rename `--no-metrics` → `--register-equality`, capability-gated

1. `experiments/subB_mixed_stride/run_sweep.py` line 173: straight rename to
   `"--register-equality"`. Confirmed no archival-binary use case in this
   directory (`run_sweep.py` takes one `--binary` arg; `sbatch_maurano.sh`/
   `sbatch_synthetic_ie.sh` each pin one binary, never two) — safe as a blind
   rename, unlike `fusion_ab.py`.
2. `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py` — **two
   separate `--no-metrics` literals**, both need fixing:
   - `harness_floor_ms()` line 383: straight rename to `"--register-equality"`
     (this probe always runs against the current/resolved binary, never an
     archival one).
   - `run_hammock()` line ~506: capability-gated, corrected design:
     - Add `_probe_binary_help(binary)`, `lru_cache`-wrapped, running
       `<binary> --help` once per distinct path (wrap the subprocess call in
       try/except, raising a clean `RuntimeError` on a bad/non-executable
       path rather than letting `OSError`/`PermissionError` propagate raw).
     - Add `_metrics_off_flag(binary) -> str`: `"--register-equality"` if
       `"--register-equality"` appears in the probed help text, else
       `"--no-metrics"`.
     - `run_hammock()`: `cmd += ["--metrics"] if metrics else [_metrics_off_flag(binary)]`.
     - No signature change (`binary` is already `run_hammock`'s first
       positional param) — no call-site changes needed in `fusion_ab.py`,
       `sweep.py`, or `pairwise_cost_by_precision.py`.
     - `check_binary_version()` keeps its existing `--version`-based
       minimum-version gate (unrelated concern — refuses binaries `< 0.7.0`)
       — this is an *additional*, separate probe for the flag-vocabulary
       question specifically, not a replacement.
     - Verified: this correctly handles `fusion_ab.py`'s dual-binary loop
       (each of `--pre-binary`/`--post-binary` gets its own cached
       capability read) and `benchmark_cpp_vs_bedtools.py`'s/`sweep.py`'s/
       `pairwise_cost_by_precision.py`'s single-binary call sites (all
       resolve to the current build, correctly land on `--register-equality`
       once that build actually has the flag — true today for
       `--register-equality`'s presence check regardless of the `--version`
       string ambiguity, since we're now checking `--help` text directly,
       not a version number).

## Category 2 — Prose/docstring corrections (unchanged from v2)

Same list as v2: `benchmark_cpp_vs_bedtools.py` (version-gate error string,
`run_hammock()` docstring — now also documents the capability-probe
mechanism instead of a version check, `--metrics-arm`/`--metrics-all` help
text paired with its Category 5 marker), `sweep.py` (`parse_hammock_csv()`
docstring + error message with corrected column counts: `--register-equality`
is 4 total columns not 3, bare default is 3 total not 9, `--metrics` is 10
total not 9; `--metrics-arm` docstring/help), `fusion_ab.py` (prose naming
`--no-metrics` literally — not its internal `"no_metrics"` arm-label
strings), `pairwise_cost_by_precision.py`, `run_sweep.py` (subB_mixed_stride),
`sbatch_fig3_largeN.sh`/`sbatch_fig3_panelA.sh`, `sbatch_maurano.sh`
(rename + incidental pre-existing stale-comment fix), `measure_cli_overhead.py`
(honest mismatch description, paired with its Category 5 marker).

## Category 3 — Bare Python-CLI call sites: filename-tag-aware code fixes now; data rename deferred to Part 8

### 3a. Stale result files — DEFERRED TO PART 8 (not a Part 6 action item)

Scope for Part 8, corrected count: `raw_abc/` (29), `raw_d/` (235),
`synthetic_evolution/results/hammock_g100000000/` (84),
`synthetic_evolution/results/hammock_g1000000/` (84) — **432 total**.
`raw_d_buggy_pre_fix/` (209 files) excluded — dead/superseded, unread by any
current script (documented elsewhere in the repo). Verified column-safety of
the rename holds for all four directories, same reasoning as v2: `analyze.R`'s
ABC scan uses only the default `jcols = "jaccard_similarity"` (never
overridden — confirmed by re-grepping every `scan_dir(` call in the file, and
confirming no second scan of `abc_dir` exists); Mode D's `jaccard_similarity_ie`
fallback reconstructs from `containment_AB`/`containment_BA` when absent
(lines 166-171, corrected from v2's "~161-169" citation); `synthetic_evolution/code/analyze.R`
only reads `jaccard_similarity`/legacy `jaccard`.

Left for whoever executes Part 8: a concrete rename command, e.g.
```bash
for f in *.csv; do
  case "$f" in *_ie.csv|*_re.csv|*_full.csv) continue ;; esac
  mv -- "$f" "${f%.csv}_full.csv"
done
```
run once per directory, executed only after the Category 3b code fixes below
have also landed on `main` (same merge, not a separate step — avoids the
window where renamed files exist but `main`'s consumers don't expect them).

### 3b. Code fixes — land in Part 6 now (filename-tag-aware, flag added)

1. **`experiments/maurano_dhs_validation/run_sweep_abc.py`** and
   **`run_sweep_d.py`** — restored from v1 (seed doc's named minimum, dropped
   in v2 by mistake). Add `"--metrics"` to each script's `cmd = [...]`
   (`run_sweep_abc.py` ~154-161, `run_sweep_d.py` ~140-148). Also fix their
   own filename-construction helpers so the summary CSV they write records
   the correct path: `run_sweep_abc.py`'s `output_csv_name()` (~L72-78) and
   `run_sweep_d.py`'s equivalent inline construction (~L138) both need `_full`
   appended, matching what `outprefix.py` will actually produce once
   `--metrics` is passed. (Verify exact current text against source at
   implementation time — cited from an earlier investigation pass, not
   independently re-confirmed in this final consolidation.)
2. **`experiments/maurano_dhs_validation/workflow/Snakefile`** — five
   pattern-string locations need `_full` before `.csv`, kept in sync
   (Snakemake requires wildcarded `output:` to literally match the pattern
   used elsewhere): `AB_CSVS` (line 54), `C_EXPA_CSVS` (57), `C_SUBB_CSVS`
   (60), `D_CSVS_BASE` (64), `D_CSVS_EXT` (66) — AND the four rules' own
   `output:` declarations (`run_mode_ab` 123, `run_mode_c_expA` 144,
   `run_mode_c_subB` 166, `run_mode_d` 188). `rule analyze`'s `input:`
   references the variables, no separate edit needed. Add `--metrics` to all
   four `shell:` blocks.
3. **`experiments/primate-phylogeny/workflow/Snakefile`** — `rule phylo_hammock`'s
   `output.csv` (161) and `rule build_phylogeny`'s `input.matrix` (190), both
   need `_full`, kept in sync. Add `--metrics` to the `shell:` block
   (172-179).
4. **`experiments/mus-homo/workflow/Snakefile`** — `rule mus_homo_hammock`'s
   `output.csv` (118) and `rule cluster_and_plot`'s `input.matrix` (144),
   both need `_full`. Add `--metrics` to the `shell:` block (129-136).
5. **`experiments/ref-comparison/workflow/Snakefile`** — `rule exp_a_hammock`'s
   `output.csv` (141) and `rule exp_a_validate_and_plot`'s `input.matrix`
   (167), both need `_full`. Add `--metrics` to the `shell:` block
   (153-161).
6. **`experiments/synthetic_evolution/Makefile`** — `HAMMOCK_AB` (82-83, feeds
   the `%`-pattern rules at 165/179 — variable AND both rule-target
   declarations need `_full`, kept in sync); `HAMMOCK_C_DEF` (86)/`HAMMOCK_C_EXP`
   (89) feed static pattern rules keyed directly off the variable (194, 210)
   — only the variable needs the edit. Verified the `sed`-based rate/type/expA
   extraction (197-198, 213-215) still works correctly with `_full` appended
   (unanchored `.*` absorbs it — traced against concrete example filenames).
   Add `--metrics` to all four hammock invocations (172-176, 186-190,
   202-206, 219-223). Grouped analysis rule (226) needs no separate edit —
   inherits via `$(HAMMOCK_CSVS)`.
7. **`experiments/synthetic_evolution/code/analyze.R`** — `parse_meta`'s regex
   (line 102): append `_full` before `\\.csv$`. Re-verified via live
   `Rscript` execution against 3+ concrete filenames incl. a Mode-C
   expA+subB edge case — correct, old untagged names now correctly rejected.
8. **`experiments/maurano_dhs_validation/analyze.R`** — three edits, all
   re-verified via live `Rscript`: `modeB_path` literal (line 143) →
   `hammock_hll_p21_jaccB_full.csv`; `parse_abc_name`'s regex (line 229) →
   append `_full$`; `parse_d_name`'s regex (line 240) → append `_full$`.
   Confirmed no fourth hardcoded hammock-filename literal/regex exists
   elsewhere in this file.
9. **`experiments/bedtools_benchmark/estimator_compare.py`** — add
   `"--metrics"` to `cmd` (113-114); path construction (line 127) →
   `f"{prefix}_hll_p{precision}_jaccB_full.csv"`.
10. **`experiments/subB_mixed_stride/run_ie_subb.py`** — add `"--metrics"` to
    `cmd` (45-47). No filename-tag fix needed (already glob-based, tag-agnostic).
11. **`experiments/primate-phylogeny/scripts/precision_probe.sh`** — add
    `"--metrics"` to the `sbatch --wrap=` invocation (line 47); narrow the
    glob (line 68) to `probe*_full.csv` (fixes both the missing-flag issue
    and the ASCII-sort stale-pick hazard in one change).
12. **`paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R`**
    (new, round-2 finding) — `default_csv` hardcoded path (~lines 62-65,
    reads `raw_d/hammock_mnmzr_p24_jaccD_k10_w30.csv`) → append `_full`
    before `.csv`. This script generates manuscript Figure 6 — verify the
    fix at implementation time by actually regenerating the figure once the
    Part-8 rename has landed, not just by inspection.
13. **`experiments/maurano_dhs_validation/mode_c_interpolation.R`** (new,
    round-2 finding) — hardcoded `raw_abc/hammock_hll_p21_jaccA.csv`/`jaccB.csv`
    literals (~lines 47-49) → append `_full`; unanchored `list.files(pattern = "^hammock_hll_p21_jaccC_B[0-9.]+\\.csv$")`
    and its `expA` sibling (~lines 101-102, 114-115) → require `_full` before
    `.csv` in the pattern, matching `analyze.R`'s regex treatment (item 8) —
    this file was documented in `analyze.R`'s own comment as reading "the
    same raw CSVs" but wasn't caught when `analyze.R` was fixed; fix both
    together this time.

## Category 5 — `# PART9:` markers (mark precisely; no flag changes, no re-runs)

1. `benchmark_cpp_vs_bedtools.py` — `--metrics-arm`/`--metrics-all` (argparse
   ~1131-1144, `arms_for()`/`ie_tool_name_for_subb()` ~599-636). Target:
   drop to bare default.
2. `sweep.py` — untimed second `--metrics` pass, precision axis (~409-440,
   `if args.ie:` block). Reads only `jaccard_similarity_ie` (verified,
   line ~429-430). Target: drop to bare default.
3. `sweep.py` — **two markers**, not one: `sweep_threads` (~line 535) AND
   `sweep_intervals` (~line 581), both `run_hammock(..., metrics=use_metrics)`
   calls under the `--metrics-arm`/`IE_ARM_TOOL` mechanism (symbol defined
   ~512/~559). Both must be marked — they're deliberately parallel per their
   own docstrings, and marking only one risks a future Part 9 pass breaking
   that parity.
4. `measure_cli_overhead.py`/`measure_threads_isolation.py` — de-hedged
   target: drop `--metrics` from `cpp_cmd`, matching `cli_cmd`'s already-bare
   state. Neither script reads any CSV column (verified) — no real
   ambiguity, state the concrete target.
5. `subB_mixed_stride/run_sweep.py` (line 173 / argparse ~264-269) — resolved,
   not caveated: traced the one live consumer (`summarize_synthetic_ie.py`)
   and confirmed it reads neither `jaccard_ie` nor `containment_AB`/`containment_BA`
   — this site cleanly qualifies for "drop to bare default." State as
   settled in the marker text.

## Category 2 (docs) restoration — README updates

`experiments/bedtools_benchmark/README.md` and `experiments/subB_mixed_stride/README.md`
"Conventions"-style sections need the new three-shape contract (currently
describe the old two-shape one).

## New: pre-existing bugs fixed (unrelated to the rename, per user's "fix them in Part 6" confirmation) — now 5, not 4

1. `experiments/bedtools_benchmark/estimator_compare.py` — covered above
   (3b.9), already broken today independent of this restructure.
2. `experiments/ref-comparison/estimator_ie_crossref.py` (116-120) —
   `sorted(glob.glob(...))[0]` → narrow to `"*_full.csv"`.
3. `experiments/primate-phylogeny/estimator_ie_topology.py` (143-146) — same
   fix, same reasoning.
4. `experiments/mus-homo/estimator_ie_tissue.py` (80, 84) — same fix, two
   call sites.
5. **New**: `experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py:543-544`
   — `run_hammock()`'s own `output_csv` capture has the identical unsorted-
   glob-then-`[0]` defect shape as items 2-4. Narrow similarly: since this
   function's `metrics` parameter determines the actual tag written, glob for
   the tag matching that call's own `metrics` value (`_full` when
   `metrics=True`, else whatever `_metrics_off_flag`/bare-default produces)
   rather than a single hardcoded suffix — this one needs to stay
   flag-aware, unlike items 2-4 which always want `_full` specifically.

(`precision_probe.sh`'s glob-narrowing fix, item 3b.11 above, is the same bug
class but was already counted under Category 3b rather than this list —
noting for inventory accuracy, not a separate action.)

## Done-criteria checklist

- [x] `grep -rn -- '--no-metrics'` under `experiments/`/`paper/` returns nothing outside dated-archival `RESULTS.md`/`draft.md` prose and explanatory "was `--no-metrics`"/"renamed from `--no-metrics`" context (verified 2026-08-11; every remaining hit is one of those two, plus the intentional `_metrics_off_flag` fallback string itself in `benchmark_cpp_vs_bedtools.py`).
- [x] Every live script/pipeline reading named similarity/containment columns (Category 3b) passes `--metrics` AND its filename-matching logic expects `_full` (all 13 items landed: `run_sweep_abc.py`/`run_sweep_d.py`, 4 Snakefiles, the Makefile, `analyze.R` x2, `estimator_compare.py`, `run_ie_subb.py`, `precision_probe.sh`, `plot_sequence_tissue_clustering.R`, `mode_c_interpolation.R`).
- [x] Every genuine timing/benchmark script (Category 5) is PART9-marked with a concrete, non-hedged target — 5 files (`benchmark_cpp_vs_bedtools.py`, `sweep.py`, `measure_cli_overhead.py`, `measure_threads_isolation.py`, `run_sweep.py`), item 3 has its required two markers (`sweep_threads`/`sweep_intervals`).
- [x] Five pre-existing glob/path bugs fixed, not four (`estimator_compare.py`, `estimator_ie_crossref.py`, `estimator_ie_topology.py`, `estimator_ie_tissue.py`, `benchmark_cpp_vs_bedtools.py`'s own `output_csv` capture — the last one resolved as "already safe, `sorted()` added defensively" rather than needing a flag-aware glob, since `out_prefix` is unique per call via `tempfile.NamedTemporaryFile`).
- [x] `fusion_ab.py`'s capability-probe fix verified against the actual current binary (not just by reading the design) — confirmed empirically during design: this worktree's `hammock-cpp` reports `--version 0.7.1` but already lacks `--no-metrics`, which is exactly the ambiguity the capability-probe design (vs. the rejected semver design) exists to handle.
- [x] Category 3a (432-file rename) explicitly NOT done in Part 6 — confirmed deferred to Part 8, noted in the Part 8 handoff below.
- [x] Full `pytest tests/` still green: 264 passed, 8 skipped (worktree baseline was 244p/8s at v0.7.1; growth is from Parts 2-5 landing more tests, not from Part 6 — none of Part 6's edits touch `tests/`, `python/`, or `cpp/`).

## Handoff note for Part 8

Two items beyond the deferred rename: (1) rebuild the pinned `hammock-cpp`
binary that `sbatch_maurano.sh`/`sbatch_synthetic_ie.sh` reference on `main`,
as part of the same merge that lands the Python-side flag rename — not a
separate step, to avoid a live SLURM job hitting a stale binary. (2) The
`lru_cache` staleness note above — no action needed unless Part 8's own
workflow rebuilds a binary at a fixed path mid-process.
