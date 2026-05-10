#ifndef HAMMOCK_STRIDE_HPP
#define HAMMOCK_STRIDE_HPP

#include "hammock/abstract_sketch.hpp"
#include "hammock/hll_sketch.hpp"
#include <cstdint>
#include <string>

// Sample point positions in [start, end) on the given chromosome, hash each
// one with xxh64(seed=hll_seed), and add it to the sketch. Returns the count
// of points actually added.
//
// When subsample == 1.0, every position is sampled. When subsample < 1.0:
//   mixed_stride=false → hash-threshold: hash each position with xxh32
//                       seed=31337 and keep those whose hash falls below
//                       subsample * UINT32_MAX.
//   mixed_stride=true  → mixed-stride deterministic: stride S ≈ 1/p chosen
//                       per chromosome via xxh32 seed=31337 hash of the chr.
//
// The subsampling seed (31337) is hardcoded to match the reference Python
// contract; only the HLL ingestion seed is user-controllable.
size_t add_interval_points_to_sketch(const std::string& chr,
                                     int64_t start, int64_t end,
                                     AbstractSketch& sketch,
                                     const std::string& separator,
                                     double subsample,
                                     bool mixed_stride,
                                     uint64_t hll_seed);

// HLLSketch overload — same behavior as the AbstractSketch version but takes
// a concrete final-class reference so the compiler can inline HLLSketch::add
// inside the inner ~50M-iteration loop. Mode B/C call this directly.
size_t add_interval_points_to_sketch(const std::string& chr,
                                     int64_t start, int64_t end,
                                     HLLSketch& sketch,
                                     const std::string& separator,
                                     double subsample,
                                     bool mixed_stride,
                                     uint64_t hll_seed);

#endif
