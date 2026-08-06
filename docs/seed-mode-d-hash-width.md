# Seed: Mode D minimizer hash width / distribution vs the HLL estimator

Handoff note for a fresh chat. Written 2026-07-31; status block added
2026-08-06. **Nothing here is decided.** Everything below is evidence gathered
plus what still needs establishing.

- **Question.** `digest` returns ≤32-bit minimizer hashes while `HLLSketch`
  assumes 64. Does that bias Mode D's cardinality-derived columns, and if so by
  how much and in which direction?
- **Measured, and settled:** the ≤32-bit premise is true (§"What I confirmed"),
  and the spike into the top rho bucket is real at ~3–6% of minimizers.
- **Measured, and it overturns the original report's explanation:** the
  mechanism is **minimizer small-value bias** (a minimizer is a *minimum* over a
  w-window, so its hash distribution is heavily skewed small), **not** "only 8
  bits left for rho". The correction flips the precision dependence: the spike
  is a **p=24 phenomenon** (≤0.04% at p ≤ 18), not "worst at low precision".
- **Open.** Sign and size of the resulting cardinality bias — the reported
  −0.5% to −8.3% has never been reproduced and its sign disagrees with the
  mechanism. Also open: whether to fix it, and how (§"What still needs
  establishing").

**Dated note, 2026-08-06.** This file was written before v0.6.0 removed the
`_with_ends` family (CLAUDE.md divergence #8). Passages below that speak of "the
two HLLs Mode D merges" describe the pre-v0.6.0 code and are kept as written;
Mode D now builds a single minimizer HLL. The practical consequence is recorded
in the update block under item 4: `hash_size=32` is now a single-sketch change
with nothing to merge against.

## The original report (from another chat, unverified when made)

> `digest` returns ≤32-bit minimizer hashes while the ends/HLL path uses 64-bit
> xxh64. At p=24 that leaves 8 bits for rho, spiking 4.2% of elements into the
> rho=41 bucket and biasing Mode D cardinality by −0.5% to −8.3% (worst at low
> precision). This affects every cardinality-derived Mode D column, not just
> with_ends.

A verification agent was tasked with this and **died on a session limit without
reporting**, so the numbers above were never independently checked.

## What I confirmed on 2026-07-31

**The premise is correct.** `digest.window_minimizer(..., include_hash=True)`
returns hashes strictly below 2³² (max observed bit length 32, zero values
≥ 2³² across ~110k minimizers at three (k,w) settings). The start/end path in
the same sketch used full-width `xxhash.xxh64`. So the two HLLs that Mode D
merged for `*_with_ends` were fed hashes on **different scales**. (That merge is
gone as of v0.6.0; the ≤32-bit minimizer hashes are not — they are still what
the surviving minimizer HLL ingests, which is why this seed is still live.)

**The ~4% figure reproduces**, but the reported *mechanism* is wrong, and the
correction matters because it changes which precisions are affected.

The relevant code is `HLLSketch::add`, `cpp/include/hammock/hll_sketch.hpp:31-45`:

```cpp
const size_t idx  = hash_val & (num_registers_ - 1);   // low p bits
const uint64_t rest = hash_val >> precision_;
const size_t max_rho = hash_size_ - precision_;        // 64 - p
if (rest == 0) pos = max_rho + 1;                      // <-- the spike
else           pos = min(ctz(rest) + 1, max_rho + 1);
```

For a 32-bit hash, `rest` is nonzero and below 2^(32−p), so `ctz(rest)` is at
most 31−p — nowhere near the `max_rho` cap. **The 8-bits-left-for-rho framing
does not describe what happens.** The only way a 32-bit hash reaches the
`rho = 65 − p` bucket is `rest == 0`, i.e. the whole hash below 2^p.

For a *uniform* 32-bit hash that is P = 2^(p−32) = **0.39%** at p=24 — not 4%.
The extra factor comes from somewhere else:

**Minimizer selector hashes are order statistics, not uniform draws.** A
minimizer is the *minimum* hash over a w-window, so its distribution is heavily
biased small. Measured mean, as a fraction of 2³² (uniform would be 0.5):

| k, w | mean/2³² | P(hash < 2²⁴) | uniform-32bit would be |
| ---- | -------- | ------------- | ---------------------- |
| 8, 20  | 0.068 | **4.49%** | 0.39% |
| 10, 30 | 0.046 | **6.24%** | 0.39% |
| 15, 15 | 0.088 | **3.12%** | 0.39% |

So the spike is real and 8–16× larger than truncation alone explains. At
p ≤ 18 the same quantity is ≤0.04%, so **the spike is a p=24 phenomenon**, not
"worst at low precision" as originally reported. That inversion is worth
resolving early — it points the rerun scope in opposite directions.

## What still needs establishing

1. **Sign and size of the cardinality bias.** A spike into a *large* rho bucket
   should make the sketch look sparser and push the estimate **up**; the report
   says down (−0.5% to −8.3%). One of these is wrong. Note Mode D at p=24 runs
   at λ = n/m ≪ 1 (16.7M registers, far fewer minimizers), which is the regime
   where Ertl's estimator leans on the zero-register count, so intuition from
   the dense regime may not transfer. Measure it directly: build an HLL from
   known-distinct minimizer hashes and compare `estimate_cardinality()` to the
   true distinct count, sweeping p.
2. **Whether non-uniformity breaks anything beyond the `rest == 0` spike.** HLL
   assumes uniformly distributed hashes. Minimizer hashes are not uniform in
   magnitude — but the register index is the *low* p bits and rho is a
   trailing-zero count, both of which stay near-uniform under a magnitude bias.
   Plausibly the spike is the whole effect. Worth confirming rather than
   assuming.
3. **Which columns move, and by how much.** Everything cardinality-derived:
   `containment_*`, `cosketch_*`, and `jaccard_similarity_ie` (v0.5.0+). NOT
   `jaccard_similarity`, which is a ratio of two register counts and never
   calls the estimator.
4. **Whether the fix is a fix or a divergence.** Options sketched, none chosen:
   re-hash each minimizer through xxh64 before ingestion (changes every Mode D
   number ever produced); pass `hash_size=32` for the minimizer HLL (the
   parameter already exists and is already guarded — but then the minimizer and
   start/end HLLs cannot be merged, and `merge_max`/`union_with`/
   `jaccard_similarity` all correctly refuse mismatched `hash_size`, which is
   what `*_with_ends` needs); or document it and leave the numbers alone.

   > **Update 2026-07-31 — the second option is now unblocked.** v0.6.0 removed
   > the `_with_ends` family (CLAUDE.md divergence #8), so there is no longer a
   > start/end HLL to merge with and nothing requires the two sketches to share
   > a `hash_size`. `hash_size=32` on the minimizer HLL is now a single-sketch
   > change. Note the whole-sequence fallback in `sequence.py` still ingests an
   > `xxh64` value, so a mixed-width sketch is still possible on corpora with
   > sub-threshold records — that path would need re-hashing too.

## Constraints carried over from the current work

- **`jaccard_similarity` is frozen.** Changing it forces a rerun of every
  experiment. It is not affected by this issue anyway (point 3).
- Mode D parity vs orig is **structural, not byte-equal** — see CLAUDE.md
  divergence #6. Orig has the same 32-bit minimizer hashes, so whatever this
  turns out to be, orig has it too.
- Mode D needs a working `digest`; if it silently returns J=0, see
  `memory/project_modeD_zero_rpath_digest.md` (RPATH shadowing libstdc++)
  before suspecting anything else.
- Build with the conda env's own compiler, not a spack gcc module (CLAUDE.md
  "Cluster compiler caveat").

## Reproducing the measurements above

```python
import random
from hammock.modes.sequence import window_minimizer
random.seed(11)
seq = "".join(random.choice("ACGT") for _ in range(400000))
hs = [h for _, h in window_minimizer(seq, k=10, w=30, include_hash=True)]
print(max(hs) < 2**32, sum(hs)/len(hs)/2**32,
      sum(1 for h in hs if (h >> 24) == 0)/len(hs))
```
