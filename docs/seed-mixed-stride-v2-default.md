# Seed: promote mixed-stride v2 to the default implementation

## Objective

Make the current `MixedStrideV2` implementation the production/default mixed-stride behavior in `hammock`, while preserving the legacy v1 implementation under an explicit legacy name for reproducibility.

This is a code-repository task in `jessicabonnie/hammock`. Do **not** edit the `hammock_paper` repository in this task. The user wants to audit and update manuscript text separately after the code behavior is settled.

Branch for this work:

```text
mixed-stride-v2-default
```

Before editing locally:

```bash
git fetch origin
git switch mixed-stride-v2-default
git pull --ff-only
```

## Why this change is appropriate

The current code contains two implementations:

- `SubBMethod::MixedStride`: legacy v1. For non-integral reciprocal rates, one of `S0=floor(1/r)` or `S1=ceil(1/r)` is chosen for the entire chromosome.
- `SubBMethod::MixedStrideV2`: chromosome-anchored fractional-interval grid. For non-integral reciprocal rates, it mixes adjacent integer gap lengths *within* each chromosome using fixed-point Beatty/Bresenham-style arithmetic.

The manuscript-facing benchmark rates are `subB=0.1` and `subB=0.01`. At both rates `1/r` is integral, so v2 deliberately delegates to the legacy grid. The implementation in `cpp/src/stride.cpp` contains this guard explicitly, and `cpp/tests/test_stride.cpp` already asserts that v1 and v2 produce identical sample counts and **identical HLL register arrays** at both rates.

Therefore the numerical manuscript findings at `0.1` and `0.01` should be invariant to the default switch by construction. However, I did **not** find a clearly documented full Maurano v1-v2 end-to-end rerun after v2 was added. Do not merely assume that happened: perform a compact end-to-end confirmation in this task.

## Current state to inspect first

Relevant files include at least:

```text
cpp/src/stride.cpp
cpp/include/hammock/stride.hpp
cpp/src/processing_modes.cpp
cpp/app/hammock_cli.cpp
bindings/_core.cpp
python/hammock/cli.py
python/hammock/runner.py
cpp/tests/test_stride.cpp
README.md
```

Search the repository for every reference to:

```text
MixedStride
MixedStrideV2
mixed-stride
mixed-stride-v2
subB_method
--subB-method
```

Do not assume the list above is exhaustive.

## Desired user-facing semantics

After this change:

- `mixed-stride` is the recommended/default public method and uses the v2 algorithm.
- The old v1 behavior remains available explicitly for reproducibility, preferably as `mixed-stride-v1` or another clear legacy spelling.
- `mixed-stride-v2` may remain accepted as an alias if useful for backward compatibility, but it should no longer be described as experimental if the verification gates below pass.
- `hash-threshold` and `single-hash` semantics remain unchanged.
- `subB=1.0` continues to keep every base regardless of method.
- `subB=0` continues to retain no points.

Prefer a design in which the public/default spelling `mixed-stride` maps to v2 while the enum/internal naming makes the distinction unambiguous. Do **not** simply delete v1; retaining it is valuable for historical reproduction and direct regression tests.

## Important compatibility issue already identified

At the current branch point, `MixedStrideV2` exists in the C++ implementation and standalone C++ CLI, but the Python binding/CLI is not fully wired for it:

- `bindings/_core.cpp::parse_subB_method` accepts `hash-threshold`, `mixed-stride`, and `single-hash`, but not `mixed-stride-v2`.
- the pybind default for `subB_method` is currently `mixed-stride`.
- `python/hammock/cli.py` exposes only `hash-threshold`, `mixed-stride`, and `single-hash`, with `mixed-stride` as default.

The promotion must therefore be consistent across **both** front ends, not only `hammock-cpp`.

## Implementation requirements

1. Preserve both algorithms in the C++ core.
2. Make the normal/default `mixed-stride` path invoke v2.
3. Add an explicit public spelling for the legacy v1 behavior.
4. Ensure both Python and standalone C++ CLIs accept the same relevant method names.
5. Ensure the Python binding parser accepts every method the Python CLI can emit.
6. Update help strings and README text so they describe v2 as the default behavior accurately.
7. Do not silently alter unrelated estimator/output semantics.
8. Do not regenerate manuscript figures or edit manuscript prose in this branch.
9. Audit experiment scripts that explicitly pass `--subB-method mixed-stride` before changing them. Do not mechanically rewrite them all:
   - if a script exists to reproduce an archived v1 result at a non-integral rate, pin it to the explicit legacy name;
   - if it is intended to exercise the current/default method, allow it to follow the new v2 default;
   - at manuscript rates `0.1` and `0.01`, v1 and v2 must remain byte-identical anyway.

## Verification gates

### A. Existing focused tests

Run the mixed-stride tests first. The current test suite includes checks that:

- v1 is deterministic per chromosome;
- v1 is interval-partition invariant on the normal path;
- v2 mixes gaps within a chromosome at `r=0.3`;
- v2 preserves legacy integral grids at `0.1` and `0.01`;
- v2 is partition invariant for long chromosome names;
- `subB=0` retains no points.

After the rename/default changes, adapt test names/spellings as needed but preserve these behavioral assertions.

### B. Add explicit default-behavior tests

Add tests that prove the public/default method now means v2, not merely that v2 exists. At minimum:

- at a non-integral rate such as `0.3`, default/`mixed-stride` matches explicit v2 behavior and differs from explicit legacy-v1 behavior on a chromosome where v1 would select only one gap length;
- at `0.1` and `0.01`, default/`mixed-stride`, explicit v2 alias (if retained), and explicit legacy v1 produce identical HLL register arrays;
- Python and standalone C++ front ends resolve the same method spelling to the same algorithm.

### C. End-to-end manuscript-rate parity check

Run a compact real-data confirmation on the Maurano 20-file corpus used for Figure 3B, using the manuscript configuration:

```text
mode B / interval
p = 18
threads = 8
subB = 1.0, 0.1, 0.01
jaccard_similarity_ie
```

At `0.1` and `0.01`, compare explicit legacy v1 against the new default/v2. Require exact equality of the pairwise numerical outputs, not just similar MAE.

Because the implementation delegates integral reciprocal rates to the legacy grid, any difference in pairwise values at these rates is a blocker and indicates a wiring/regression problem.

For runtime, exact equality is not required. Confirm only that the default switch does not introduce an obvious material regression at the manuscript rates.

Do not overwrite tracked canonical result CSVs or figures for this verification. Use `/tmp` or another scratch location.

### D. Non-integral behavior check

Use at least one non-integral rate (`0.3` is already represented in unit tests; optionally also `0.15` or `0.07`) and confirm:

- v2 uses both adjacent gap lengths within a chromosome;
- realized density approaches the requested rate on a sufficiently long span;
- interval partitioning does not change the sampled-coordinate sketch;
- changing thread count does not change the result;
- long chromosome names take the same v2 logic rather than the legacy fallback behavior.

### E. Build/test suite

Run the cleanest available full verification for this checkout. Prefer:

```bash
cmake --build build -j
ctest --test-dir build
pytest tests/
```

If the repository's normal build directory/configuration differs, follow the repo's current documented build procedure instead. Record any pre-existing failures separately from failures introduced by this branch.

## Documentation updates in the code repository

Update current user-facing documentation that describes the default method, including at least `README.md` and CLI help strings.

The accurate high-level description of v2 is approximately:

> deterministic chromosome-anchored fractional-interval sampling that mixes the two adjacent integer gap lengths needed to realize the requested sampling rate, using a chromosome-specific phase and no per-position sampling hash.

Avoid saying that v2 randomly chooses one stride for an entire chromosome; that describes v1.

Avoid claiming `mixed-stride` is a standardized literature term. It is hammock-specific terminology.

## Manuscript consequence to report, but do not edit here

Once v2 is the default, the currently visible manuscript Methods subsection should be re-audited against v2. The old formula

```text
r = q(1/S1) + (1-q)(1/S0)
```

and chromosome-level probabilistic stride choice describe v1, not v2.

For v2, if `S0 != S1`, the fraction `f` of gaps of length `S1` satisfies the mean-gap equation

```text
(1-f) S0 + f S1 = 1/r
f = (1/r - S0) / (S1 - S0)
```

The terminal task should finish by reporting the exact manuscript statements that need correction, but **must not edit `hammock_paper`** unless the user separately authorizes those edits.

## Scope safeguards

- Work only on the `mixed-stride-v2-default` branch.
- Do not push directly to `main`.
- Do not edit the manuscript repository.
- Preserve unrelated user changes.
- Do not delete archived experiment outputs.
- Do not reinterpret historical benchmark results at non-integral rates without checking which implementation generated them.

## Completion report

At the end, report:

1. exact code files changed;
2. exact public method names and which implementation each maps to;
3. focused mixed-stride test results;
4. full test/build results;
5. Maurano v1-v2 equality result at `subB=0.1` and `0.01`;
6. any runtime difference observed at those rates;
7. a non-integral-rate demonstration showing the new default really is v2;
8. any experiment scripts intentionally pinned to legacy v1 and why;
9. the specific visible manuscript Methods claims that now need updating;
10. commit SHA(s) on this branch.

Do not merge the branch automatically.