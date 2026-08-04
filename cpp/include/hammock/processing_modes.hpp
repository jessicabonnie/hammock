#ifndef HAMMOCK_PROCESSING_MODES_HPP
#define HAMMOCK_PROCESSING_MODES_HPP

#include "hammock/abstract_sketch.hpp"
#include "hammock/stride.hpp"
#include <cstdint>
#include <string>

// Build the per-interval hashed string used by Mode A (and the interval
// component of Mode C). Default separator '\t' matches the reference Python.
std::string create_interval_string(const std::string& chr, int64_t start, int64_t end,
                                   const std::string& separator = "\t");

// Build the per-position hashed string used by Mode B (and the point
// component of Mode C). Note: the Python reference appends NO suffix here,
// just chr + sep + position.
std::string create_point_string(const std::string& chr, int64_t pos,
                                const std::string& separator = "\t");

// Subsampling-decision hash: xxh32 with the gate seed (default 31337 to
// preserve orig parity; user-settable via --gate-seed).
uint32_t hash_for_subsample(const std::string& s, uint32_t gate_seed = 31337);

size_t process_bed_file_mode_a(const std::string& filepath, AbstractSketch& sketch,
                               uint64_t hll_seed = 42,
                               const std::string& separator = "\t",
                               bool verbose = false);

size_t process_bed_file_mode_b(const std::string& filepath, AbstractSketch& sketch,
                               uint64_t hll_seed = 42,
                               const std::string& separator = "\t",
                               double subB = 1.0,
                               SubBMethod method = SubBMethod::MixedStride,
                               uint32_t gate_seed = 31337,
                               bool verbose = false);

size_t process_bed_file_mode_c(const std::string& filepath, AbstractSketch& sketch,
                               uint64_t hll_seed = 42,
                               const std::string& separator = "\t",
                               double subA = 1.0,
                               double subB = 1.0,
                               double expA = 0.0,
                               SubBMethod method = SubBMethod::MixedStride,
                               uint32_t gate_seed = 31337,
                               bool verbose = false);

#endif
