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

`--mode B`, subB=1.0, 10k intervals/file, `--threads 16`. µs/pair is the
pairwise phase divided by n·m, so it is a throughput figure at that thread
count, not a serial per-pair cost, and it does not transfer across thread counts.

| p | base µs/pair | with the block | multiplier |
|---|---|---|---|
| 14 | 1.40–1.46 | 4.88–5.00 | **3.3–3.6×** |
| 18 | 15.1–16.3 | 41.6–43.5 | 2.6–2.7× |
| 24 | 1050–1295 | 2829 | 2.2–2.5× |

Two things worth stating plainly, because both are the opposite of what you would
guess:

- **The multiplier is largest at *low* precision.** Six extra `%.17g` fields per
  row are precision-independent, so they dominate when the register work is
  cheap. This is why `--verbose` now reports `Pairwise:` and `Write:`
  separately — without the split, a figure would attribute serial `fprintf` to
  the estimator arithmetic.
- **It is small in wall-time terms everywhere the paper runs.** The pairwise
  phase is 0.61% of wall at N=512, p=14 (0.330 s of 53.721 s), so the block costs
  about **+1.5%** there. At N=64, p=24 the phase is 4–5% of wall base and ~10%
  with the block, against 98.5 s of sketching. Pairwise does not overtake
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
