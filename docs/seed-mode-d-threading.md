# Seed: Mode D threading is a GIL convoy — the default makes it slower

Handoff note for a fresh chat. Written 2026-08-03. Nothing here is a decision;
it is the evidence gathered so far plus what still needs establishing.

## The finding

`--threads` makes Mode D **slower**, and the CLI defaults to
`min(8, cpu_count())` (`python/hammock/cli.py:345-346`). Every Mode D run in
this project's history — every sweep in `experiments/*/results/` — was run at
`--threads 8` and paid this.

Interval mode (A/B/C) is **not** affected: its work is inside the C++ extension,
which does release the GIL. This is a Mode D problem only.

## What was confirmed on 2026-08-03

All on a 48-core Rockfish node, `claude-ref-comparison` env, hammock v0.6.0.

**End-to-end CLI, 4 real Maurano DHS FASTAs, k=8 w=40 p=20:**

| threads | wall | maxRSS |
|---|---|---|
| 1 | **31.9 s** | 74.2 MB |
| 8 | **68.3 s** (2.1× slower) | 75.2 MB |

**Independently, 8 full-size Maurano FASTAs (864 MB/side, 1.73 Gbp sketched):**

| cell | v0.5.0 T1 | v0.5.0 T8 | v0.6.0 T1 | v0.6.0 T8 |
|---|---|---|---|---|
| k=8, w=40, p=20 | 95.3 s | 358.7 s | 61.8 s | 172.9 s |
| k=20, w=100, p=24 | 117.3 s | 586.1 s | 47.9 s | 129.2 s |

So the convoy costs 2.8–4.5× on v0.6.0 and 3.8–5.0× on v0.5.0 — i.e. it is
*larger* than the entire `_with_ends` removal win (divergence #8), which bought
1.5–2.5× single-threaded.

**In-process scaling, private per-thread sequence lists** (so it is not refcount
contention on shared input), k=8 w=40 p=20, 60k records/task:

```
threads=1  per-task 0.87s  1.00x
threads=2  per-task 1.43s  0.61x
threads=4  per-task 1.01s  0.86x
threads=8  per-task 2.20s  0.40x
```

**It is the threading model, not the machine.** Eight concurrent *processes* on
the identical workload scale ~7.1–7.6× (near-linear).

**Corroborating archived evidence:** every 8-thread Mode D cell in
`experiments/maurano_dhs_validation/results/sweep_d_*.csv` records
`cpu_s / wall_s ≈ 0.8`. A genuinely parallel 8-thread run would be near 8.
That ratio was sitting in the repo the whole time.

**Why.** `digest.window_minimizer` does not usefully release the GIL at this
call pattern — ~11 µs per call on a ~466 bp record — and alone tops out at 1.4×,
degrading past T4. The surrounding per-record Python in
`MinimizerSketch.add_string` holds the GIL outright. At k=8 w=40, `digest` is
only ~63% of v0.6.0's sketch time; the rest is Python.

## Where the code is

- `python/hammock/runner.py:168-184` — `_parallel_map`, the
  `ThreadPoolExecutor(max_workers=threads)`. Correctly falls through to a serial
  path when `threads <= 1 or n <= 1`.
- `python/hammock/runner.py:187-211` — `_sketch_many`, one task per input file.
  Note the unit of parallelism is the **file**, so a run over few large files
  has little to parallelize even if threading worked.
- `python/hammock/cli.py:137-138, 345-346` — the `--threads` flag and the
  `min(8, cpu_count())` default.
- `python/hammock/modes/sequence.py` — `add_string` is the hot Python.

Interval mode reaches the same `_parallel_map` but its per-file work is a C++
call that releases the GIL, so leave that path alone.

## What still needs establishing

1. **Whether interval mode's threading actually benefits**, measured rather than
   assumed. `_parallel_map` is shared. Before changing a global default, confirm
   A/B/C really do scale — nobody has published a threads-axis Mode B number
   against this claim, though `experiments/bedtools_benchmark/sweep.py` has a
   `threads` axis that could answer it directly.
2. **Where the Mode D crossover is.** At what (file count, file size) does the
   convoy start costing? With one input file there is one task and no contention;
   the 4-file run above already shows 2.1×.
3. ~~Whether a process pool is viable.~~ **Answered 2026-08-03: it is blocked,
   and unblocking it is a bindings change.** `_core.HLLSketch` is **not
   picklable** (`pickle.dumps` → `TypeError`) and the binding surface exposes
   **no register access at all** — the full public API is `add_hash64`,
   `add_string`, `clear`, `estimate_cardinality`, `estimate_intersection`,
   `estimate_jaccard`, `precision`. So a worker process cannot return a sketch
   today. First step for this route is adding `__getstate__`/`__setstate__` (or
   a register-array getter/setter) to `bindings/_core.cpp`. That is a small,
   self-contained change and is independently useful — it would also let sketches
   be cached to disk between runs.
4. **Whether the win survives process startup + IPC.** Each worker imports
   BioPython and `digest`; sketches are 2^p bytes to ship back (16.7 MB at p=24,
   × N files). Measure against the ~7.1× process-level ceiling above.
5. **Whether `digest` can be made to release the GIL.** It is a third-party
   pybind11 extension (VeryAmazed). If `window_minimizer` released the GIL and
   the per-record Python were cut, threads might work as intended. This is the
   principled fix but depends on an upstream package.
6. **Whether to just batch differently.** Parallelism is per-file; a coarser
   option is to leave threading alone and document `--threads 1` for Mode D,
   letting users get parallelism from SLURM array jobs across files — which is
   how the sweeps are already structured.

## Option space (none chosen)

- **Document only.** Tell users to pass `--threads 1` for Mode D. Zero risk,
  zero code change; already done in `CLAUDE.md`'s Architecture section. Leaves
  a footgun armed by default.
- **Change the Mode D default to 1** while leaving A/B/C at
  `min(8, cpu_count())`. Smallest real fix, and mechanically easy: `cli.py:343`
  already sets `args.mode = _autodetect_mode(args)` **two lines before** the
  `args.threads` default at `:345-346`, so the mode is known — verified. Makes
  every future Mode D run 2–4.5× faster with no accuracy change. **Changes
  wall-time numbers, not values** — output is byte-identical across thread
  counts (`tests/test_mode_d.py:81` runs `--threads 4` against serial and
  asserts byte equality; verified present).
- **Process pool for Mode D.** The real fix; ~7.1× available. Blocked on item 3
  — needs pickling support on `HLLSketch` first.
- **Upstream GIL release in `digest`.** Principled, slowest to land, and would
  still leave the per-record Python.

## Constraints carried over

- **Output must stay byte-identical across thread counts.**
  `tests/test_mode_d.py` asserts `--threads 4` output equals single-threaded
  byte-for-byte. Any change here must keep that test green.
- **`jaccard_similarity` is frozen.** Threading changes must not touch values —
  they are a wall-time-only concern.
- Mode D parity vs orig is structural, not byte-equal (CLAUDE.md divergence #6).
- Archived Mode D wall/RSS figures are inflated by this convoy *and* by the
  pre-v0.6.0 `_with_ends` work. Do not compare new timings against them; re-baseline.
- `experiments/bedtools_benchmark/` is Mode B only and is unaffected either way.

## Reproducing the measurements above

```bash
cd /home/jbonnie1/interval_sketch/hammock_claude
H=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock
ls experiments/maurano_dhs_validation/data/fastas/*.fa | head -4 > /tmp/fl4.txt
for T in 1 8; do
  /usr/bin/time -f "threads=$T wall=%e maxRSS=%MkB" \
    $H /tmp/fl4.txt /tmp/fl4.txt --mode D -k 8 -w 40 -p 20 --threads $T -o /tmp/th$T
done
```

**Two environment hazards that cost a reviewer real time here:**

- This shell's `PYTHONPATH` is the literal string
  `PYTHONPATH:/home/jbonnie1/interval_sketch/` (a missing `$`), which puts the
  *original* hammock's parent directory on `sys.path`.
- The conda env's editable install uses a meta-path finder
  (`_hammock_editable.pth`) that overrides `sys.path` entirely. A
  `--system-site-packages` venv will silently resolve `import hammock` to the
  repo working tree despite a correct venv install. To compare two versions,
  build real wheels into fully isolated venvs and assert `hammock.__version__`
  before trusting a number.
