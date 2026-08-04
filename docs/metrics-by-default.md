# Making the two front-ends agree (v0.7.0)

Evidence and consequences for the change that flipped `hammock-cpp`'s defaults,
removed BagMinHash, and fixed two filename collisions. Companion to divergence
#9 in `CLAUDE.md`.

## What was wrong

`hammock-cpp` and the Python CLI disagreed about what a run produces, in two
ways that both failed silently.

**Columns.** The Python CLI has emitted `jaccard_similarity_ie` unconditionally
since v0.4.0. The binary emitted three columns unless `--metrics` was passed,
justified in `CLAUDE.md` as keeping benchmark timings comparable. That is a
developer's reason living in a user's default, and it points away from the
column the documentation recommends: `CLAUDE.md` §"Pick by what you need" says
to read `jaccard_similarity_ie`, because `jaccard_similarity` is a
register-equality statistic carrying a chance-agreement floor and is not
rank-faithful across pairs of differing size.

**Mode.** `hammock_cli.cpp` defaulted to `--mode A` while `cli.py`'s
`_autodetect_mode` returns **B** for BED input. The same two lists through the
two front-ends therefore produced a column *of the same name* holding numbers
from different algorithms.

Both are now B / 9 columns.

### Honest weighting

This was not an emergency. The binary is not on `$PATH` (it installs to
`<site-packages>/bin/`), `CMakeLists.txt` calls it the "benchmark binary", and
`README.md` already told you to pass `--metrics` "when the output is going to be
analyzed rather than timed". What survives is narrower: defaults should be
correct rather than fast, and the timing rationale belongs in the harnesses —
which is where it now lives, as an explicit `--no-metrics` on every timed path.

## Why the filename tag was inverted

The suffix now marks the **reduced** shape (`_j3`), not the metrics shape.

The invariant worth preserving is *filename ↔ column count*: a 9-column file and
a 3-column file must never collide on one path. Both directions satisfy it. The
tiebreaker is that tagging the metrics shape would have moved the default output
path for every existing invocation, while tagging the reduced shape leaves the
default exactly where every pre-0.7.0 default run already wrote.

That is safe here for a reason specific to this repo: **no raw `hammock-cpp` TSV
is archived anywhere.** Every one is written to a `NamedTemporaryFile` prefix and
deleted by its caller. A search across `experiments/` (with symlinks followed
into `/vast`), `docs/`, `paper/` and `tests/` for a file whose first line begins
`query` returns nothing. And all three parsers key on header names, not column
positions, so neither naming could cause a mis-parse.

One thing to know before "fixing" it: the C++ suffix grammar now has a `_j3`
component with **no counterpart in `outprefix.py`**, because the Python CLI has
only one output shape. That asymmetry is intentional.

Note also that `docs/data/hammock_hll_p{18,21,23}_jaccB.csv` *look* like
counterexamples — they match the C++ naming pattern — but they are
comma-separated Python-CLI output keyed `file1`/`file2`, and are unaffected.

## Measured cost of the block

`--mode B`, subB=1.0, N=64 files per side (4,096 pairs), 10k intervals/file,
`--threads 16`, 5 runs per cell with the two arms rotated by run index; medians
below. Generator `experiments/bedtools_benchmark/pairwise_cost_by_precision.py`,
data `docs/data/pairwise_cost_by_precision_20260804_164807.csv`.

µs/pair is a phase time divided by n·m, so it is a throughput figure **at that
thread count**, not a serial per-pair cost, and it does not transfer across
thread counts.

The two phases behave completely differently, which is the whole reason
`--verbose` splits them:

| p | `pair_time` µs/pair | | ×  | `write_time` (ms, whole run) | | Δ |
|---|---|---|---|---|---|---|
| | base | +block | | base | +block | |
| 12 | 0.28 | 0.72 | 2.54 | 1.7 | 9.9 | +8.2 |
| 14 | 1.05 | 2.66 | 2.53 | 1.8 | 10.0 | +8.2 |
| 16 | 4.85 | 10.77 | 2.22 | 1.8 | 10.1 | +8.3 |
| 18 | 16.74 | 40.88 | 2.44 | 1.9 | 10.1 | +8.3 |
| 20 | 64.20 | 164.50 | 2.56 | 1.9 | 10.1 | +8.2 |
| 22 | 264.16 | 676.31 | 2.56 | 1.9 | 16.2 | +14.3 |
| 24 | 1073.15 | 2806.33 | 2.62 | 4.4 | 10.0 | +5.7 |

- **The estimator multiplier is flat at ≈2.5×** (2.22–2.62 across twelve
  doublings of the register array), and `pair_time` scales as Θ(2^p) almost
  exactly: p=14→24 is **1021×** against the 1024× the register count predicts.
  A union plus two cardinality estimates is a constant factor on top of the
  register-equality pass; nothing about it is precision-sensitive.
- **The write cost is a constant ≈ +8 ms**, independent of p — six extra `%.17g`
  fields per row over 4,096 rows. (The p=22 and p=24 deltas are the same
  quantity plus filesystem noise on a ~10 ms median.)

Those two facts together give the wall-visible behaviour:

| p | `comparison_time` = pair + write, µs/pair | | multiplier |
|---|---|---|---|
| | base | +block | |
| 12 | 0.71 | 3.14 | **4.43×** |
| 14 | 1.49 | 5.10 | **3.43×** |
| 16 | 5.30 | 13.23 | 2.50× |
| 18 | 17.20 | 43.35 | 2.52× |
| 20 | 64.65 | 166.96 | 2.58× |
| 22 | 264.63 | 680.28 | 2.57× |
| 24 | 1073.62 | 2810.97 | 2.62× |

- **The multiplier is largest at *low* precision** — and it is entirely an
  output-formatting effect, not an estimator effect. At p=12 the +8 ms of
  `fprintf` is 2.9× the whole base comparison phase; by p=16 it is under 8% of
  it and the ratio has settled onto the flat 2.5×. Without the `Pairwise:` /
  `Write:` split a figure would charge that serial `fprintf` to the estimator
  arithmetic and report a spurious precision dependence.
- **It is small in wall-time terms everywhere the paper runs.** The pairwise
  phase is 0.61% of wall at N=512, p=14 (0.330 s of 53.721 s), so the block costs
  about **+1.5%** there. At N=64, p=24 it is 4–5% of wall base (4.4 s against
  98.5 s of sketching) and ~11% with the block. Pairwise does not overtake
  sketching until roughly N≈550 with the block on.

## Two silent-overwrite bugs fixed alongside

Both in `outprefix_with_suffix`, both found while adding the tag, both latent
(no in-repo caller triggers them today):

- **`--subA` never reached the filename.** `outprefix.py:26` emits
  `_A{subA:.2f}`; the C++ had no such branch, so two Mode C runs differing only
  in `--subA` wrote to the same path and the second silently won.
- **`subB` used an unconditional `%.2f`**, so `--subB 0.001` and `--subB 0.005`
  both rendered `_B0.00`. `outprefix.py:33-34` carries a `.4f` rule below 0.01
  for exactly this reason. The boundary is strict, so `subB == 0.01` keeps its
  historical `_B0.01`.

## Rollback

`git revert` **is not sufficient on its own.** A revert restores sources with
older mtimes than the built binary, so `test_binary_is_not_stale` passes and the
flipped 0.7.0 binary keeps living in both `build/<tag>/` and
`<site-packages>/bin/`. The full procedure:

```bash
git revert <commit>
CR=$CONDA_PREFIX
CC=$CR/bin/x86_64-conda-linux-gnu-gcc CXX=$CR/bin/x86_64-conda-linux-gnu-g++ \
  pip install -e . --no-build-isolation
readelf -d build/*/…/_core*.so | grep RPATH   # $ORIGIN-relative first
```

The reinstall matters independently of the revert: `tests/` globs
`build/*/hammock-cpp` while six sbatch scripts point `HAMMOCK_CPP_BIN` at the
site-packages copy, so rebuilding without reinstalling leaves the cluster runs on
the old binary. That is why `find_hammock_cpp()`'s callers now probe `--version`
on the **resolved** path rather than only on the build-tree glob — the env-var
branch is the one where the mismatch can actually happen.

## Verification performed

- Mode B `--no-metrics` output **byte-identical** to the pre-flip binary, and the
  first three columns of the new 9-column default byte-identical to it as well —
  the frozen-`jaccard_similarity` gate, run with `--mode B` pinned on both sides
  so the mode-default change could not confound it.
- Modes A, B and C, plus `--metrics` and a `--subA`/`--subB` run, all
  byte-identical across the BagMinHash removal.
- `Pairwise + Write` reconstructs `Pairwise+write` to 1 µs.
- `HAMMOCK_REQUIRE_CPP=1 pytest tests/` — 124 passed, 8 skipped (132 collected),
  with `test_metrics_block_matches_python_bit_for_bit` still asserting exact
  equality.
