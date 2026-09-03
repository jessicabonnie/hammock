#ifndef HAMMOCK_STRIDE_HPP
#define HAMMOCK_STRIDE_HPP

#include "hammock/abstract_sketch.hpp"
#include "hammock/hll_sketch.hpp"
#include <cstdint>
#include <string>

// Selects how subB (< 1.0) sampling decides which positions to keep.
//
//   HashThreshold — orig-parity: hash each position with xxh32 seed=31337
//                   and keep those below subsample * UINT32_MAX. One gate
//                   hash + one ingestion hash per accepted point.
//   MixedStrideV1 — legacy deterministic stride S ≈ 1/p chosen per chromosome
//                   via xxh32 seed=31337 hash of the chr. No per-position
//                   gate; only accepted positions are hashed at all.
//   MixedStrideV2 — default chromosome-anchored fractional-interval grid.
//                   Non-integral reciprocal rates mix the adjacent integer
//                   gaps within each chromosome. Integral rates use the
//                   legacy grid exactly, preserving published results.
//   SingleHash    — one xxh64 per position with seed=hll_seed; high 32 bits
//                   decide the gate, full 64 are the HLL ingestion hash.
//                   Opt-in parity divergence (different accepted-position set
//                   than HashThreshold).
//
// When subsample == 1.0, all methods behave identically (every position
// is sampled and the gate path isn't entered).
enum class SubBMethod {
    HashThreshold,
    MixedStrideV1,
    MixedStrideV2,
    SingleHash,
};

// Reject invalid rates and mixed-stride gaps that cannot be represented by
// the int64_t stride implementation. A zero rate is valid and retains no
// points; non-mixed methods do not have the reciprocal-size restriction.
void validate_subB_rate(double subsample, SubBMethod method);

size_t add_interval_points_to_sketch(const std::string& chr,
                                     int64_t start, int64_t end,
                                     AbstractSketch& sketch,
                                     const std::string& separator,
                                     double subsample,
                                     SubBMethod method,
                                     uint64_t hll_seed,
                                     uint32_t gate_seed = 31337);

// HLLSketch overload — same behavior as the AbstractSketch version but takes
// a concrete final-class reference so the compiler can inline HLLSketch::add
// inside the inner ~50M-iteration loop. Mode B/C call this directly.
size_t add_interval_points_to_sketch(const std::string& chr,
                                     int64_t start, int64_t end,
                                     HLLSketch& sketch,
                                     const std::string& separator,
                                     double subsample,
                                     SubBMethod method,
                                     uint64_t hll_seed,
                                     uint32_t gate_seed = 31337);

#endif
