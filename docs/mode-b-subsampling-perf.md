# Mode B subsampling: future perf ideas

## Context

After three rounds of Mode B optimization (committed):

1. Per-position string-buffer reuse (eliminate 2 heap allocations per genomic
   position).
2. Race-free OMP via thread-local HLLs that merge at end (also fixes the
   non-atomic max-update on `registers_`).
3. Inline HLL `add` with `__builtin_ctzll` + `final`-class devirtualization
   of the inner stride-loop call site.

Cumulative speedup vs original on a 16-file (8q × 8r), ~50M-base-each
benchmark: **~32% faster single-thread, ~25% less CPU time at all configs.**

The two remaining ideas below are deferred deliberately. Both **only matter
when subsampling is active** (`--subB < 1.0`); on the full-coverage path the
hash work is the floor and there's nothing more to skip. So neither belongs
on the critical path until we have a profiled subsampled workload that says
they're worth the parity-risk.

## Idea 1 — hand-rolled BED parser

### What

`parse_bed_line` (`cpp/src/bed_parser.cpp:32`) currently uses
`std::istringstream` per BED line. For every line it allocates a stream
object, runs locale-aware `operator>>` for chromosome + start + end, and
optionally walks tokens for the peak-height column. Profiling on the
50K-interval synthetic shows parsing is a small fraction of total time
(~50K parses vs ~50M point hashes); replacing it with a hand-rolled parser
using `strtoll` + `memchr` would only shift the needle when **the parse
loop runs many times relative to the inner point loop** — i.e. when the
inner point loop is short.

### When it matters

The two regimes where `intervals.size()` grows large relative to total
points sampled:

- **Heavy subsampling** (`--subB 0.01`–`--subB 0.1`): same intervals, but
  ~1/100 to 1/10 the inner-loop work per interval. Parser cost rises from
  noise to a meaningful fraction.
- **Many small intervals**: e.g. SNP-like BED files with millions of 1bp
  intervals. The interval count itself dominates total work.
- **Mode A**: parser cost is *all* of the per-interval work (no point loop
  at all). A faster parser is a direct win for Mode A regardless of
  subsampling.

### Sketch of the change

```cpp
bool parse_bed_line_fast(const char* line, size_t len, ...) {
    const char* p = line;
    const char* end = line + len;
    // chr: scan to first \t
    const char* chr_end = static_cast<const char*>(memchr(p, '\t', end - p));
    if (!chr_end) return false;
    chr.assign(p, chr_end - p);
    p = chr_end + 1;
    // start, end: strtoll
    char* tmp;
    start = std::strtoll(p, &tmp, 10);
    if (tmp == p) return false;
    p = tmp + 1;
    end_pos = std::strtoll(p, &tmp, 10);
    ...
}
```

Read the file with `mmap` (or a single `fread` into a buffer) and walk
line-by-line via `memchr('\n', ...)` rather than `std::getline`. That also
removes the per-line `std::string` allocation in the read loop.

### Risks

- **None for parity** — same BED-format contract, same edge cases (header
  lines, `track`/`browser` prefixes, blank lines). The orig hammock just
  parses BED.
- Adds ~80 lines of careful pointer arithmetic; tests (`cpp/tests/
  test_bed_parser.cpp`) need new cases for short-line / ragged-EOF / no
  trailing newline.

### Estimated win

Maybe 5–15% on heavy subsampling workloads, ~20%+ on Mode A. Not measured
yet — verify by profiling a real subsampled workload before committing.

## Idea 2 — combine the xxh32 + xxh64 paths under subsampling

### What

When `--subB < 1.0`, the inner stride loop hashes each point string
**twice**:

1. `xxhash::hash32(point_buf, seed=31337)` — the gate hash. If above
   threshold, skip.
2. `xxhash::hash64(point_buf, seed=hll_seed)` — only on accepted points,
   for HLL ingestion.

Both are full passes over a ~14-byte buffer (~10ns each). For a typical
`--subB 0.1`, of every 10 points: 10× xxh32 + 1× xxh64. At Mode B's
~50M-points-per-file scale that's a lot of redundant byte-mixing.

The shape of the win: **fold the gate decision into the same hash pass as
the ingestion hash, so accepted points pay one pass and rejected points
pay (almost) nothing.**

### Why it's blocked on parity

The original Python contract is *specifically* xxh32 with seed 31337 for
the gate. Two separate algorithms (xxh32 and xxh64 are different — not
just truncations of one another), two separate seeds. We can't change
either without breaking byte-equal CSV parity for any subsampled run on
the orig.

So this idea has three flavors, in order of feasibility:

#### Flavor A — short-circuit subsample-rejected points

**Don't compute xxh64 on rejected points.** That's already what we do —
`stride.cpp:73-77` does `continue` on rejection before the xxh64 line. So
this is already optimal *given* the xxh32 contract. No-op.

#### Flavor B — opt-in `--fast-subsample` flag

A new flag that lets callers say "I don't need orig parity on subsampled
runs." Under it: compute one xxh64, use the high 32 bits as the gate,
the low 64 as the ingestion hash. Single pass, ~halves hash time on
accepted points and saves the xxh64 on rejected ones.

This is a real divergence, like the existing `--mixed-stride` and the
`--subB`-honoring change documented in `CLAUDE.md`. It would need:

- A flag (likely `--fast-subsample` or `--fast-gate`).
- The flag prints a one-line stderr note when used (parity asterisk).
- A test that confirms it's a near-equivalent estimator (Jaccard within
  noise) but warns it's *not* byte-equal.
- A line in `CLAUDE.md`'s "Intentional divergences" section.

Estimated win: at `--subB 0.1`, roughly 40–50% of inner-loop time saved
(was 10× xxh32 + 1× xxh64, becomes 1× xxh64). At `--subB 0.5` it's smaller
(~15–25%). At `--subB 1.0` it's zero (the gate path isn't entered).

#### Flavor C — replace xxh32 with a 32-bit slice of xxh64 in the gate

A more drastic version of B that just changes the hashing convention to
`hash32(s) := xxh64(s, seed=31337) >> 32`. Same single-pass benefit,
applies even when the user *doesn't* opt in, but breaks parity for every
subsampled run. Probably not worth the parity damage — orig users do
care about reproducibility.

### Decision

Defer until someone has a profiled workload where `--subB < 1.0` is the
hot path **and** they're willing to take a parity divergence. Then ship
flavor B (the opt-in flag) — flavor C is too costly.

## When to revisit

A trigger for either idea: a benchmark sweep on heavy-subsampled real
workloads (e.g. `experiments/bedtools_benchmark/sweep.py` with low `subB`)
that shows Mode B inner-loop time has plateaued and parsing/hashing are
the visible bars on the flame graph. Until then, the current
full-coverage path is well below memory bandwidth and the next bottleneck
is unclear without measurement.
