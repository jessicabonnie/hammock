#include "hammock/stride.hpp"
#include "hammock/processing_modes.hpp"
#include "hammock/xxhash.hpp"
#include <charconv>
#include <cmath>
#include <string>

namespace {

// Templated implementation so the HLLSketch instantiation can inline
// HLLSketch::add (defined in the header) into the inner ~50M-iteration loop.
// The AbstractSketch instantiation goes through the vtable — same as before.
template <typename Sketch>
size_t add_points_impl(const std::string& chr,
                       int64_t start, int64_t end,
                       Sketch& sketch,
                       const std::string& separator,
                       double subsample,
                       bool mixed_stride,
                       uint64_t hll_seed) {
    const bool do_subsample = (subsample < 1.0);
    size_t sampled = 0;

    thread_local std::string point_buf;
    point_buf.assign(chr);
    point_buf.append(separator);
    const size_t prefix_len = point_buf.size();
    char int_buf[24];

    if (mixed_stride && do_subsample) {
        const double p = subsample;
        const double inv = 1.0 / p;
        int64_t s0 = static_cast<int64_t>(std::floor(inv));
        int64_t s1 = static_cast<int64_t>(std::ceil(inv));
        if (s0 < 1) s0 = 1;
        if (s1 < 1) s1 = 1;

        bool choose_s1 = false;
        if (s0 != s1) {
            const double denom = (1.0 / s1) - (1.0 / s0);
            const double q = (denom != 0.0) ? (p - 1.0 / s0) / denom : 0.0;
            const uint32_t h = hash_for_subsample(chr);
            const double rfloat = h / 4294967295.0;
            choose_s1 = (rfloat < q);
        }
        const int64_t S = (s0 != s1 && choose_s1) ? s1 : s0;

        const uint32_t h2 = hash_for_subsample(chr + "|res");
        const int64_t residue = h2 % S;

        int64_t offset = (residue - (start % S)) % S;
        if (offset < 0) offset += S;
        const int64_t x0 = start + offset;

        for (int64_t pos = x0; pos < end; pos += S) {
            auto r = std::to_chars(int_buf, int_buf + sizeof(int_buf), pos);
            point_buf.resize(prefix_len);
            point_buf.append(int_buf, static_cast<size_t>(r.ptr - int_buf));
            const uint64_t hash_val = xxhash::hash64_short(point_buf.data(),
                                                          point_buf.size(),
                                                          hll_seed);
            sketch.add(hash_val);
            sampled++;
        }
    } else {
        const uint32_t threshold = static_cast<uint32_t>(subsample * 4294967295.0);

        for (int64_t pos = start; pos < end; pos++) {
            auto r = std::to_chars(int_buf, int_buf + sizeof(int_buf), pos);
            point_buf.resize(prefix_len);
            point_buf.append(int_buf, static_cast<size_t>(r.ptr - int_buf));

            if (do_subsample) {
                const uint32_t point_hash = xxhash::hash32(point_buf.data(),
                                                          point_buf.size(),
                                                          31337);
                if (point_hash > threshold) {
                    continue;
                }
            }

            const uint64_t hash_val = xxhash::hash64_short(point_buf.data(),
                                                          point_buf.size(),
                                                          hll_seed);
            sketch.add(hash_val);
            sampled++;
        }
    }

    return sampled;
}

}  // namespace

size_t add_interval_points_to_sketch(const std::string& chr,
                                     int64_t start, int64_t end,
                                     AbstractSketch& sketch,
                                     const std::string& separator,
                                     double subsample,
                                     bool mixed_stride,
                                     uint64_t hll_seed) {
    return add_points_impl(chr, start, end, sketch, separator,
                           subsample, mixed_stride, hll_seed);
}

size_t add_interval_points_to_sketch(const std::string& chr,
                                     int64_t start, int64_t end,
                                     HLLSketch& sketch,
                                     const std::string& separator,
                                     double subsample,
                                     bool mixed_stride,
                                     uint64_t hll_seed) {
    return add_points_impl(chr, start, end, sketch, separator,
                           subsample, mixed_stride, hll_seed);
}
