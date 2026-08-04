# Seed: Mode D threading is a GIL convoy — PARTIALLY RESOLVED

Handoff note for a fresh chat. Written 2026-08-03; updated 2026-08-04.

**Status: the footgun is disarmed, the opportunity is still open.**

- **Decided and shipped (v0.6.1, 2026-08-04):** Mode D now defaults to
  `--threads 1`; interval modes keep `min(8, cpu_count())`. See "What was
  decided" below.
- **Still open:** Mode D runs on **one core**. Processes scale ~7.1× on this
  workload, so a process pool is the real prize — and its stated blocker turned
  out to be wrong (see item 3, corrected). Also still open: whether interval-mode
  threading benefits *on the Python path* (item 1, sharpened).

Everything below the decision section is evidence, not decisions.

## The finding

`--threads` makes Mode D **slower**, and the CLI used to default to
`min(8, cpu_count())` for it. Every Mode D run in this project's history —
every sweep in `experiments/*/results/` — was run at `--threads 8` and paid
this.

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

## What was decided (2026-08-04, v0.6.1)

The **default** flipped, nothing else. Chosen over the process pool because it
captures the measured 2–4.5× immediately at near-zero risk, and it does not
foreclose the pool.

- `cli._default_threads(mode)` returns `1` for `"D"`, `min(8, cpu_count())`
  otherwise. `args.mode` is already resolved two lines earlier in `main()`.
- An explicit `--threads > 1` in Mode D is **still honored**, with a one-line
  stderr note. That keeps the byte-equality tests exercising the thread path
  (`tests/test_mode_d.py:81`, `tests/test_parity_against_original.py:162`) and
  keeps the door open for anyone benchmarking the convoy.
- **`bed2fasta` extraction is exempt.** `--ref` runs are Mode D but spend their
  first phase in `b2f.convert_list` → `_parallel_map` over `bedtools getfasta`
  **subprocesses**, which parallelize for real. Clamping that to 1 would have
  been a regression, so `main()` resolves a second budget, `args.io_threads`
  (user value if given, else `min(8, cpu_count())`), read by
  `runner._run_bed2fasta`. **Keep this split if you rework threading** — it is
  the one place where "Mode D" and "GIL-bound" come apart.
- `_sketch_many` / `_parallel_map` were deliberately **not** touched: the clamp
  lives entirely in default resolution, so it stays overridable.
- Verified: 4 Maurano FASTAs, k=8 w=40 p=20 → 31.9 s default vs 63.1 s at
  `--threads 8` (1.98×), CSVs byte-identical, full suite green.

## Where the code is

- `python/hammock/runner.py:168-184` — `_parallel_map`, the
  `ThreadPoolExecutor(max_workers=threads)`. Correctly falls through to a serial
  path when `threads <= 1 or n <= 1`.
- `python/hammock/runner.py:187-211` — `_sketch_many`, one task per input file.
  Note the unit of parallelism is the **file**, so a run over few large files
  has little to parallelize even if threading worked.
- `python/hammock/cli.py` — the `--threads` flag, `_default_threads(mode)`, and
  the `args.threads` / `args.io_threads` resolution in `main()`.
- `python/hammock/runner.py` (`_run_bed2fasta`) — the `args.io_threads` consumer;
  `python/hammock/bed2fasta.py:148-164` (`convert_list`) is the subprocess fan-out.
- `python/hammock/modes/sequence.py` — `add_string` is the hot Python.

Interval mode reaches the same `_parallel_map` but its per-file work is a C++
call that releases the GIL, so leave that path alone.

## What still needs establishing

1. **Whether interval mode's threading actually benefits *on the Python path*.**
   Still open, and sharper than first written. `experiments/bedtools_benchmark/
   sweep.py --axis threads` does exist, and its results are on disk
   (`results/sweep_threads_20260512_150458.csv`, tabulated in
   `bedtools_benchmark/RESULTS.md`): 64 files × 10k intervals, p=14, Mode B wall
   **72.5 s @ t=1 → 6.6 s @ t=16 (≈10.9×)**. But that harness drives the
   standalone `hammock-cpp` binary's **OpenMP** path — *not*
   `runner._parallel_map`. So it establishes that the C++ Mode B kernel scales;
   it says nothing about the Python CLI's per-file thread pool, which is the code
   Mode D shares. That has never been measured.

   Cheapest route: add a threads loop to
   `experiments/maurano_dhs_validation/run_sweep_abc.py`, which already drives the
   **Python** CLI over the 20 Maurano BEDs and records `wall_s`/`cpu_s`/
   `max_rss_mb` (`--threads` is there but fixed, not swept — ~10 lines). SLURM
   template: `experiments/bedtools_benchmark/sbatch_threads.sh`
   (`--partition=shared --account=blangme2 --cpus-per-task=16`).

   Caveat when computing convoy ratios from archived sweeps: `run_sweep_d.py`
   parses `cpu_s` from `/usr/bin/time -v`'s **"User time" only**, dropping system
   time, so the `cpu_s/wall_s ≈ 0.8` figures quoted above are if anything
   *understated*.
2. **Where the Mode D crossover is.** At what (file count, file size) does the
   convoy start costing? With one input file there is one task and no contention;
   the 4-file run above already shows 2.1×.
3. ~~Whether a process pool is viable.~~ **Answered — and the 2026-08-03 answer
   overstated the blocker. Corrected 2026-08-04.**

   True: `_core.HLLSketch` is **not picklable** (`pickle.dumps` → `TypeError`),
   and the *binding* exposes no register access — the full bound API is
   `add_hash64`, `add_string`, `clear`, `estimate_cardinality`,
   `estimate_intersection`, `estimate_jaccard`, `precision`
   (`bindings/_core.cpp:221-258`). So a worker process cannot return a sketch
   today.

   **False, and this was the load-bearing error:** that the C++ class lacks
   register access. `cpp/include/hammock/hll_sketch.hpp:55-59` already declares
   `precision()`, `hash_size_bits()`, `num_registers()`, and **both const and
   non-const `registers()`, all public**. So pickling is a **binding-only,
   strictly additive** change — one `py::pickle(...)` in `bindings/_core.cpp`,
   no header edit, no friend declaration. `bindings/_core.cpp` is its own TU
   (~8 s to compile; the link does LTO over `libhammock_core.a`, so budget tens
   of seconds for the incremental rebuild, but `hammock_core` itself does not
   rebuild).

   State is `(precision_, hash_size_, registers_)` and that is complete —
   `num_registers_` is derived (`1 << precision_`, `cpp/src/hll_sketch.cpp:33`).
   Three gotchas for whoever writes it:

   - The ctor binding is `py::init<size_t>()` only, so `hash_size` is unreachable
     from Python; `__setstate__` must use the **two-arg C++ ctor**
     `HLLSketch(precision, hash_size)`.
   - `validated_precision` (`cpp/src/hll_sketch.cpp:18-26`) enforces
     `hash_size ∈ {32,64}` and `4 ≤ precision ≤ 24`, but nothing checks the
     register buffer length — `__setstate__` must assert
     `len(registers) == 1 << precision` itself or silently accept a short buffer.
   - `pairwise_metrics_hll` takes `const std::vector<HLLSketch>&` and, via
     `pybind11/stl.h`, already **deep-copies every sketch by value** on each call.
     So the per-sketch copy cost pickling implies is already being paid today —
     the marginal cost of a process pool is IPC, not copying.

   Independently useful: this also enables caching sketches to disk between runs.
   There is no sketch serialization anywhere in the repo today (repo-wide grep for
   `pickle|__reduce__|__getstate__` hits only this file).
4. **Whether the win survives process startup + IPC.** Each worker imports
   BioPython and `digest`; sketches are 2^p bytes to ship back (16.7 MB at p=24,
   × N files). Measure against the ~7.1× process-level ceiling above.

   Two things that make the port easier than it looks, both verified 2026-08-04:
   `MinimizerSketch` needs **no** `__reduce__` — it holds only ints plus the HLL
   (`python/hammock/modes/sequence.py:49-102`), and nothing in `runner.py` reads
   anything off it but `.minimizer_hll` (`runner.py:255-256`), so workers can
   return the bare `HLLSketch` (or a `(precision, bytes)` tuple) and never pickle
   the wrapper. `sketch_fasta(path, args)` takes an `argparse.Namespace` of
   scalars, which is picklable as-is — though a module-level
   `(path, k, w, seed, precision)` function is cleaner. What *does* need porting:
   `_sketch_many`'s verbose progress counter is a `threading.Lock` over shared
   state (`runner.py:198-209`); with processes that becomes a parent-side
   callback.
5. **Whether `digest` can be made to release the GIL.** It is a third-party
   pybind11 extension (VeryAmazed). If `window_minimizer` released the GIL and
   the per-record Python were cut, threads might work as intended. This is the
   principled fix but depends on an upstream package.
6. **Whether to just batch differently.** Parallelism is per-file; a coarser
   option is to leave threading alone and document `--threads 1` for Mode D,
   letting users get parallelism from SLURM array jobs across files — which is
   how the sweeps are already structured.

## Option space

- ~~**Document only.**~~ Superseded — the footgun is disarmed by default now.
- ~~**Change the Mode D default to 1.**~~ **Taken, v0.6.1.** See "What was
  decided" above.
- **Process pool for Mode D.** ← *the remaining prize.* ~7.1× available, vs the
  ~1× Mode D gets today. Item 3's blocker is smaller than originally recorded:
  `py::pickle` in `bindings/_core.cpp`, no C++ class change. Then flip
  `_sketch_many` to a `ProcessPoolExecutor` for Mode D only (leave A/B/C on
  threads — the C++ path already releases the GIL and would gain nothing but IPC
  cost), and port the progress counter. Must clear item 4 (startup + IPC) to be
  worth it.
- **Upstream GIL release in `digest`.** Principled, slowest to land, and would
  still leave the per-record Python. Only pays off if the surrounding Python is
  cut too — at k=8 w=40, `digest` is just ~63% of v0.6.0's sketch time.
- **Leave it to SLURM.** Item 6: parallelism is per-file, and the sweeps already
  fan out across files as array jobs. This is what users get today by default,
  and it is a legitimate answer for batch workloads — it just does nothing for a
  single interactive multi-file run.

## Constraints carried over

- **Output must stay byte-identical across thread counts.**
  `tests/test_mode_d.py` asserts `--threads 4` output equals single-threaded
  byte-for-byte. Any change here must keep that test green. A process pool
  preserves this for free as long as results stay index-aligned:
  `_parallel_map` indexes by `futures[fut]`, register merge is order-independent
  elementwise max, and a pickle round-trip of the register array is exact.
- **Don't clamp `args.io_threads`.** It is the bed2fasta/`bedtools getfasta`
  budget and is deliberately not subject to the Mode D thread clamp.
- **`jaccard_similarity` is frozen.** Threading changes must not touch values —
  they are a wall-time-only concern.
- Mode D parity vs orig is structural, not byte-equal (CLAUDE.md divergence #6).
- Archived Mode D wall/RSS figures are inflated by this convoy *and* by the
  pre-v0.6.0 `_with_ends` work. Do not compare new timings against them; re-baseline.
- `experiments/bedtools_benchmark/` is Mode B only **and drives the standalone
  `hammock-cpp` binary, not the Python CLI**, so it is unaffected either way —
  and, per item 1, cannot answer whether `_parallel_map` scales.

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

Since v0.6.1 the `--threads 1` arm is also what you get by omitting `--threads`
entirely, and the `--threads 8` arm prints the convoy note. Re-run 2026-08-04 on
this exact command: 31.9 s vs 63.1 s, CSVs byte-identical.

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
