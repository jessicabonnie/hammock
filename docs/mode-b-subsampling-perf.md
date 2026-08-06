# Mode B perf — what we shipped and what's left

**Question:** after seven rounds of optimization on Mode B and its subB gate,
where does the time actually go, and which of the remaining ideas are worth
doing?

**Measured:** callgrind + gprof on 16 files × ~50 Mbp (May 2026). Plain Mode B
is 4.8× faster than the original baseline single-threaded; the inner loop is now
~60% `hash64_short`, ~17% stride bookkeeping, ~15% `HLL::add`. BED parsing is
0% — under 1% of total.

**Answer:** the hash work is the floor for the parity-required xxh64 algorithm,
so nothing parity-safe on the menu is worth more than a few percent. Mixed-stride
is already at its floor at low subB.

**Open / unbuilt:** Ideas B, C and D below (est. ~4–8%, ~1–5%, and 20–30%
respectively, each with a cost). Idea A **shipped** — see its dated note.

## Context

Mode B and the Mode B/C subsampling gate got several rounds of optimization.
Committed changes, in order:

1. **Per-position string-buffer reuse** — eliminated 2 heap allocations per
   genomic position.
2. **Race-free OMP** via thread-local HLLs that merge at end; also fixes the
   non-atomic max-update race on `registers_`.
3. **Inline HLL `add`** with `__builtin_ctzll` + `final`-class devirtualization
   at the inner stride-loop call site.
4. **`hash64_short`** — specialized xxh64 for <32-byte inputs, inlined into
   the stride loop. Drops the wide-lane branch.
5. **`point_buf` elimination** — replaced the thread-local `std::string` with
   a 64-byte stack buffer; `std::to_chars` writes the integer suffix directly
   at `prefix_len`.
6. **Incremental ASCII counter** — initial `std::to_chars` once per interval,
   then in-place ASCII increment on each iteration. Only the
   non-mixed-stride path; mixed-stride still per-iteration formats.
   (`std::to_chars` was 55% of inner-loop instructions per callgrind.)
7. **`hash32_short`** — same specialization as `hash64_short` for the subB
   gate hash.

### Cumulative speedup vs original baseline (16 files, ~50M bases each)

| Config | Original | Now | Speedup |
|---|---|---|---|
| `--threads 1` + `OMP=1` (single-thread) | 42.81s | **8.94s** | **4.8×** |
| default `--threads` + default OMP | 1.35s / 53.1 CPU-s | **0.83s / 20.5 CPU-s** | **1.6× wall, 2.6× CPU** |

### subB hash-threshold baseline (post-`hash32_short`)

| subB | hash-threshold | mixed-stride |
|---|---|---|
| 1.0 | 9.1s | — |
| 0.5 | 16.2s ← worst case | 12.4s |
| 0.1 | 8.3s | 3.0s |
| 0.01 | 6.1s | **0.8s** |

Mixed-stride is fundamentally cheaper at low subB (strides directly, no
gate) and already at the floor. Hash-threshold is where remaining perf
opportunity lives.

## Where the remaining time goes (callgrind, May 2026 after #1–6)

Total inner-loop instructions are now ~60% in `hash64_short`, ~17% in
`stride.cpp` (ASCII counter + loop control), ~15% in `HLL::add`, ~8% other.
The hash work is the floor for the parity-required xxh64 algorithm.

## Deferred ideas (with rough payoff estimates)

These are deliberately deferred. Most either give small gains or need a
parity-divergence flag. They're recorded so we know the menu when a future
profiled workload says "Mode B / subB is the hot path again."

### Idea A — `--fast-subsample` flag (combine gate + ingestion hash)

**Biggest remaining lever for subB hash-threshold.** Currently every
genomic position computes `xxh32(point, seed=31337)` to gate, then accepted
points additionally compute `xxh64(point, seed=hll_seed)`. The original
Python contract is parity-locked to those two distinct algorithms+seeds.

Behind an opt-in flag, we could compute one `xxh64`, use the high 32 bits
as the gate decision, low 64 as the ingestion hash. One hash pass instead
of two, eliminating the worst-case-at-subB=0.5 double-cost.

**Estimated impact** (extrapolating from how much the gate dominates):
- `subB=0.5` hash-threshold: 16.2s → ~10s (≈−35%)
- `subB=0.1`: 8.3s → ~5s (−35–40%)
- `subB=0.01`: 6.1s → ~3.5s (−40%+)
- `subB=1.0`: 0 (gate isn't entered)

**Cost:** parity divergence — would need a flag (`--fast-subsample` or
similar), a stderr asterisk when used, and an entry in `CLAUDE.md`'s
"Intentional divergences" section. No way to keep byte-equal CSV parity
with the orig Python while doing this.

> **Superseded 2026-08-06 — Idea A shipped, and the estimate above did not
> hold.** It landed as `--subB-method=single-hash` (CLAUDE.md divergence #4,
> `python/hammock/cli.py:237`), with exactly the design sketched here: one
> `xxh64(point, seed=hll_seed)` whose high 32 bits drive the gate and whose full
> 64 bits drive HLL ingestion, plus the stderr note and the divergence entry.
>
> The predicted −35–40% did **not** materialize. Measured on the Maurano DHS
> corpus, 5 reps (`docs/data/maurano_subB_summary.csv`, `wall_median`),
> single-hash is *slower* than hash-threshold at the low subB values the idea
> was aimed at, and wins only in the middle of the range:
>
> | subB | hash-threshold | single-hash | |
> |---|---|---|---|
> | 0.01 | 7.26 s | 7.97 s | 9.9% slower |
> | 0.05 | 8.10 s | 8.53 s | 5.3% slower |
> | 0.10 | 8.95 s | 9.21 s | 3.0% slower |
> | 0.25 | 11.18 s | 11.06 s | 1.1% faster |
> | 0.50 | 14.67 s | 13.62 s | 7.2% faster |
> | 1.00 | 9.55 s | 9.55 s | identical (gate not entered) |
>
> The best case is subB=0.5, and then by 7%, not 35% — and note that at both
> 0.25 and 0.50 *neither* method beats not subsampling at all (`wall_nosub`
> 9.54 s), so the win sits inside a regime that is already a net loss. The
> reasoning error was treating "one hash instead of
> two" as a saving: it trades a cheap `xxh32` for a more expensive `xxh64` on
> *every* position, so at low subB (where almost everything is rejected and the
> second hash was never computed anyway) it is a pure cost. Keep it as a
> research/comparison flag, not a speed lever. The numbers in the parent section
> above are the pre-ship estimates and are left as written.

### Idea B — `switch` over trailing-bytes count in `hash64_short`

In `hash64_short`, after the 8-byte and 4-byte chunks, the remaining
0–3 bytes go through a `while (p < end)` loop. Callgrind shows the loop
overhead (cond + body) is ~15% of the program for the full-coverage
benchmark. Replacing the loop with a `switch (end - p)` over the four
cases removes the loop bookkeeping.

**Estimated impact:** ~4–8% on full-coverage Mode B; similar fraction
on subB hash-threshold's gate (`hash32_short` has the same shape).
**Parity-safe.** Small change. Hasn't been done because the absolute
gain isn't large after everything else we shipped.

### Idea C — hand-rolled BED parser

Documented in the previous version of this file. After today's work,
parsing is well under 1% of total time on full-coverage Mode B (gprof
confirmed: `parse_bed_line` shows 0%). The only regime where it could
matter is *very* heavy subsampling on *very* sparse BEDs, or Mode A
(where parsing is most of the per-interval work).

**Estimated impact:**
- Full-coverage Mode B: ~1% (noise).
- Heavy subB (e.g. 0.01) on dense BEDs: ~5%.
- Mode A: 20–30% (parser is most of the work).

**Parity-safe.** Only worth doing if Mode A becomes a perf priority.

### Idea D — batched / SIMD hashing across multiple points

Hash N points at a time using SIMD lanes. xxh3 (the modern xxhash family)
has SIMD variants; we could either pick up a SIMD-friendly hash or write
a custom lane-parallel xxh64. Major refactor.

**Estimated impact:** 20–30% on plain Mode B; harder to estimate on subB
where the gate rejects most points.

**Cost:** significant refactor; risk of subtle correctness bugs;
architecture-specific intrinsics needed if we want predictable behavior
across CPUs. Defer until a profiled workload says it's the right move.

## When to revisit

Trigger any of these when:

- A profiled real workload on `experiments/subB_mixed_stride/` shows that
  Mode B / subB is the hot path of the user's analysis.
- The user is willing to accept a parity divergence (Idea A).
- Mode A becomes a perf priority (Idea C).
- We hit a wall and the only remaining lever is algorithmic redesign
  (Idea D).

Until then, plain Mode B is ~5× faster than the original baseline single-
threaded, and the subB hash-threshold path is within a small constant of
the mixed-stride path. That's probably good enough until the next user
need.
