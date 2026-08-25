# Mode D parallelism: current branch and preferred direction

Mode D currently sketches each FASTA through BioPython and
`MinimizerSketch.add_string`. The per-record Python loop holds the GIL, and the
short `digest.window_minimizer` calls do not release it usefully at this call
size. A Python thread pool is therefore slower than one thread on the measured
Maurano workloads.

`codex/mode-d-process-repair` is a bounded spawn-process implementation. It is
useful as a working fallback and benchmark reference, but it is not the desired
long-term architecture. Processes add startup, IPC, HLL transport, aggregate
RSS, and shutdown/error-handling complexity. The branch must always propagate
`sequence_hll_hash` to workers; otherwise serial and multi-worker runs can use
different Mode D feature-hash contracts.

## Preferred C++ boundary

Move one complete FASTA-to-HLL operation below the Python boundary, for example:

```text
sketch_fasta_hll(path, k, w, precision, seed, sequence_hll_hash)
```

The binding should release the GIL for the operation. The existing Python
file-dispatch thread pool can then parallelize Mode D in the same way it
parallelizes interval modes, without transporting completed HLLs between
processes. File-level dispatch is the first target; parallelism within one FASTA
can be considered separately if single-file workloads justify it.

Adding C++ threads around the current Python/BioPython loop is insufficient.
FASTA parsing, minimizer iteration, selector rehashing, and HLL updates must all
execute below the GIL for threads to scale.

## Required parity gates

- Default `rehash-selector64` and explicit `legacy-selector32` produce exactly
  the same HLL register arrays as the Python implementation.
- Results are invariant across thread counts and repeated runs, and final CSVs
  remain byte-identical apart from timing or progress output.
- Plain, gzipped, multiline, mixed-case, empty, ambiguous-base, and short-record
  FASTA behavior matches the current parser and fallback rules.
- Input ordering, duplicate paths, verbose progress, and error paths remain
  stable.
- `--ref` extraction keeps its independent `io_threads` budget, and pairwise
  OpenMP begins only after file sketching has completed.
- Worker/thread counts respect CPU affinity and memory limits without nested
  oversubscription.

Until this boundary exists, Mode D should retain its one-worker default. The
process branch should not be used as evidence that C++ file-level threading has
been validated.
