#include "hammock/stride.hpp"
#include "hammock/processing_modes.hpp"
#include "hammock/xxhash.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

void validate_subB_rate(double subsample, SubBMethod method) {
    if (!std::isfinite(subsample) || subsample < 0.0 || subsample > 1.0) {
        throw std::invalid_argument("subB must be finite and in [0, 1]");
    }
    const bool mixed_stride =
        method == SubBMethod::MixedStrideV1 || method == SubBMethod::MixedStrideV2;
    if (mixed_stride && subsample > 0.0 &&
        1.0L / static_cast<long double>(subsample) >
            static_cast<long double>(std::numeric_limits<int64_t>::max())) {
        throw std::invalid_argument(
            "positive mixed-stride subB is too small: reciprocal exceeds INT64_MAX; use 0 to retain no points");
    }
}

namespace {

// Stack buffer used to assemble chr+sep+pos for hashing. Sized large enough
// to fit any realistic genomics chr name (BED spec doesn't cap chr length,
// but in practice ≤30 chars even on patched assemblies) + sep + int64
// (≤20 digits). If a caller exceeds this we'd write past the end during
// std::to_chars; we guard with a fallback path.
constexpr size_t POINT_BUF_SIZE = 64;

// Fixed-point denominator for the fractional part of the requested mean gap.
// This is deliberately an integer accumulator: floating point is used once to
// convert the user-provided double, never to accumulate grid positions.
constexpr uint64_t FRACTION_SCALE = uint64_t{1} << 48;

__int128 floor_div(__int128 numerator, __int128 denominator) {
    const __int128 quotient = numerator / denominator;
    const __int128 remainder = numerator % denominator;
    return (remainder != 0 && numerator < 0) ? quotient - 1 : quotient;
}

// Coordinate k on a Beatty/Bresenham grid.  Consecutive coordinates differ
// by s0 or s1, and the long-gap fraction is fractional_numerator/scale.
__int128 mixed_stride_v2_coordinate(__int128 k, int64_t offset, int64_t s0,
                                    uint64_t phase,
                                    uint64_t fractional_numerator) {
    const __int128 scaled = static_cast<__int128>(phase) +
                            static_cast<__int128>(k) * fractional_numerator;
    const __int128 coordinate = static_cast<__int128>(offset) +
                                static_cast<__int128>(k) * s0 +
                                floor_div(scaled, FRACTION_SCALE);
    return coordinate;
}

// Find the first retained coordinate at or after start without walking the
// preceding chromosome. The bounds use the two possible gap lengths, then a
// monotone binary search on the exact fixed-point grid.
__int128 mixed_stride_v2_first_index(int64_t start, int64_t offset,
                                     int64_t s0, int64_t s1, uint64_t phase,
                                     uint64_t fractional_numerator) {
    const __int128 relative = static_cast<__int128>(start) - offset;
    const __int128 a = floor_div(relative, s0);
    const __int128 b = floor_div(relative, s1);
    __int128 low = std::min(a, b) - 4;
    __int128 high = std::max(a, b) + 4;

    // The four-coordinate padding makes these bounds strict for either sign.
    // Keep the checks as assertions in behavior rather than relying on that
    // derivation if a future fixed-point scale changes.
    while (mixed_stride_v2_coordinate(low, offset, s0, phase,
                                      fractional_numerator) >= start) {
        --low;
    }
    while (mixed_stride_v2_coordinate(high, offset, s0, phase,
                                      fractional_numerator) < start) {
        ++high;
    }
    while (high - low > 1) {
        const __int128 mid = low + (high - low) / 2;
        if (mixed_stride_v2_coordinate(mid, offset, s0, phase,
                                      fractional_numerator) < start) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return high;
}

template <typename Sketch>
size_t add_mixed_stride_v2_points(const std::string& chr, int64_t start,
                                  int64_t end, Sketch& sketch,
                                  const std::string& separator,
                                  double subsample, uint64_t hll_seed,
                                  uint32_t gate_seed) {
    if (start >= end || subsample <= 0.0) return 0;

    const double inverse = 1.0 / subsample;
    const int64_t s0 = std::max<int64_t>(1, static_cast<int64_t>(std::floor(inverse)));
    const int64_t s1 = std::max<int64_t>(1, static_cast<int64_t>(std::ceil(inverse)));

    // The integral branch delegates to the legacy implementation below. This
    // guard makes the published 0.1 and 0.01 grids byte-identical.
    if (s0 == s1) return std::numeric_limits<size_t>::max();

    const long double fractional = static_cast<long double>(inverse) - s0;
    uint64_t fractional_numerator = static_cast<uint64_t>(
        std::llround(fractional * static_cast<long double>(FRACTION_SCALE)));
    if (fractional_numerator == 0 || fractional_numerator >= FRACTION_SCALE) {
        // A user double extremely close to an integer reciprocal quantizes to
        // an integral grid. The legacy branch remains the stable definition.
        return std::numeric_limits<size_t>::max();
    }

    const uint32_t hash = hash_for_subsample(chr + "|v2", gate_seed);
    const int64_t offset = hash % s1;
    const uint64_t phase = (static_cast<uint64_t>(hash_for_subsample(
        chr + "|v2-phase", gate_seed)) << 16) ^ hash;
    const __int128 first_index = mixed_stride_v2_first_index(
        start, offset, s0, s1, phase % FRACTION_SCALE, fractional_numerator);

    std::string point = chr + separator;
    const size_t prefix_len = point.size();
    point.reserve(prefix_len + 24);
    char coordinate[24];
    size_t sampled = 0;
    for (__int128 index = first_index;; ++index) {
        const __int128 coordinate_value = mixed_stride_v2_coordinate(
            index, offset, s0, phase % FRACTION_SCALE, fractional_numerator);
        if (coordinate_value >= end) break;
        const int64_t pos = static_cast<int64_t>(coordinate_value);
        const auto result = std::to_chars(coordinate, coordinate + sizeof(coordinate), pos);
        point.resize(prefix_len);
        point.append(coordinate, static_cast<size_t>(result.ptr - coordinate));
        sketch.add(xxhash::hash64(point.data(), point.size(), hll_seed));
        ++sampled;
    }
    return sampled;
}

// Templated implementation so the HLLSketch instantiation can inline
// HLLSketch::add (defined in the header) into the inner ~50M-iteration loop.
// The AbstractSketch instantiation goes through the vtable — same as before.
template <typename Sketch>
size_t add_points_impl(const std::string& chr,
                       int64_t start, int64_t end,
                       Sketch& sketch,
                       const std::string& separator,
                       double subsample,
                       SubBMethod method,
                       uint64_t hll_seed,
                       uint32_t gate_seed) {
    if (subsample <= 0.0) return 0;
    const bool do_subsample = (subsample < 1.0);
    size_t sampled = 0;

    // When subsample is 1.0 all methods are equivalent; collapse to the
    // hash-threshold ASCII-counter fast path. (SingleHash doesn't enter its
    // gate either since we'd just be wasting the high-bits comparison.)
    const bool use_mixed_stride =
        (method == SubBMethod::MixedStrideV1 || method == SubBMethod::MixedStrideV2) &&
        do_subsample;
    const bool use_mixed_stride_v2 = (method == SubBMethod::MixedStrideV2) && do_subsample;
    const bool use_single_hash = (method == SubBMethod::SingleHash) && do_subsample;

    if (use_mixed_stride_v2) {
        const size_t v2_sampled = add_mixed_stride_v2_points(
            chr, start, end, sketch, separator, subsample, hll_seed, gate_seed);
        if (v2_sampled != std::numeric_limits<size_t>::max()) return v2_sampled;
        // Exact integral reciprocal rates intentionally use v1 below.
    }

    // Assemble chr+sep once into a stack buffer; per-point work then only
    // formats the integer at prefix_len and hashes the buffer.
    const size_t chr_len = chr.size();
    const size_t sep_len = separator.size();
    const size_t prefix_len = chr_len + sep_len;

    // ~20 chars reserved for the integer suffix. If a chr name is too long
    // to fit in the stack buffer we fall back to a thread-local std::string.
    if (prefix_len + 20 > POINT_BUF_SIZE) {
        thread_local std::string fallback;
        fallback.assign(chr);
        fallback.append(separator);
        char int_buf[24];
        const uint32_t threshold32 = static_cast<uint32_t>(subsample * 4294967295.0);
        const int64_t stride = use_mixed_stride
                                   ? std::max<int64_t>(1, static_cast<int64_t>(std::round(1.0 / subsample)))
                                   : 1;
        for (int64_t pos = start; pos < end; pos += stride) {
            auto r = std::to_chars(int_buf, int_buf + sizeof(int_buf), pos);
            fallback.resize(prefix_len);
            fallback.append(int_buf, static_cast<size_t>(r.ptr - int_buf));
            if (use_single_hash) {
                const uint64_t h = xxhash::hash64(fallback.data(), fallback.size(), hll_seed);
                if (static_cast<uint32_t>(h >> 32) > threshold32) continue;
                sketch.add(h);
                sampled++;
            } else {
                if (!use_mixed_stride && do_subsample) {
                    const uint32_t h32 = xxhash::hash32(fallback.data(),
                                                       fallback.size(), gate_seed);
                    if (h32 > threshold32) continue;
                }
                const uint64_t hash_val = xxhash::hash64(fallback.data(),
                                                        fallback.size(),
                                                        hll_seed);
                sketch.add(hash_val);
                sampled++;
            }
        }
        return sampled;
    }

    char buf[POINT_BUF_SIZE];
    std::memcpy(buf, chr.data(), chr_len);
    std::memcpy(buf + chr_len, separator.data(), sep_len);
    char* const int_start = buf + prefix_len;
    char* const buf_end = buf + POINT_BUF_SIZE;

    if (use_mixed_stride) {
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
            const uint32_t h = hash_for_subsample(chr, gate_seed);
            const double rfloat = h / 4294967295.0;
            choose_s1 = (rfloat < q);
        }
        const int64_t S = (s0 != s1 && choose_s1) ? s1 : s0;

        const uint32_t h2 = hash_for_subsample(chr + "|res", gate_seed);
        const int64_t residue = h2 % S;

        int64_t offset = (residue - (start % S)) % S;
        if (offset < 0) offset += S;
        const int64_t x0 = start + offset;

        for (int64_t pos = x0; pos < end; pos += S) {
            auto r = std::to_chars(int_start, buf_end, pos);
            const size_t total_len = static_cast<size_t>(r.ptr - buf);
            const uint64_t hash_val = xxhash::hash64_short(buf, total_len, hll_seed);
            sketch.add(hash_val);
            sampled++;
        }
    } else {
        if (start >= end) return sampled;
        const uint32_t threshold = static_cast<uint32_t>(subsample * 4294967295.0);

        // Format the starting position once; subsequent iterations just bump
        // the ASCII digits in place. Most increments only touch the last
        // byte; carries to longer numbers happen only at powers of 10 (rare
        // within a single BED interval). This avoids ~55% of the inner-loop
        // instructions previously spent in std::to_chars (callgrind, May
        // 2026).
        auto r = std::to_chars(int_start, buf_end, start);
        size_t int_len = static_cast<size_t>(r.ptr - int_start);

        for (int64_t pos = start; pos < end; pos++) {
            const size_t total_len = prefix_len + int_len;

            if (use_single_hash) {
                // One xxh64 does both jobs — gate (high 32 bits) and HLL
                // ingestion (full 64). Diverges from orig parity: different
                // accepted-position set than HashThreshold.
                const uint64_t h = xxhash::hash64_short(buf, total_len, hll_seed);
                if (static_cast<uint32_t>(h >> 32) <= threshold) {
                    sketch.add(h);
                    sampled++;
                }
            } else if (do_subsample) {
                const uint32_t point_hash = xxhash::hash32_short(buf, total_len, gate_seed);
                if (point_hash <= threshold) {
                    const uint64_t hash_val = xxhash::hash64_short(buf, total_len, hll_seed);
                    sketch.add(hash_val);
                    sampled++;
                }
            } else {
                const uint64_t hash_val = xxhash::hash64_short(buf, total_len, hll_seed);
                sketch.add(hash_val);
                sampled++;
            }

            // In-place ASCII increment. The common case is a single byte
            // bump on the last digit; carries propagate left only when the
            // digit was '9'.
            char* p = int_start + int_len - 1;
            while (*p == '9') {
                *p = '0';
                if (p == int_start) {
                    // Roll-over (e.g. 999 -> 1000): grow length by 1.
                    std::memmove(int_start + 1, int_start, int_len);
                    *int_start = '1';
                    int_len++;
                    goto inc_done;
                }
                p--;
            }
            (*p)++;
        inc_done:;
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
                                     SubBMethod method,
                                     uint64_t hll_seed,
                                     uint32_t gate_seed) {
    validate_subB_rate(subsample, method);
    return add_points_impl(chr, start, end, sketch, separator,
                           subsample, method, hll_seed, gate_seed);
}

size_t add_interval_points_to_sketch(const std::string& chr,
                                     int64_t start, int64_t end,
                                     HLLSketch& sketch,
                                     const std::string& separator,
                                     double subsample,
                                     SubBMethod method,
                                     uint64_t hll_seed,
                                     uint32_t gate_seed) {
    validate_subB_rate(subsample, method);
    return add_points_impl(chr, start, end, sketch, separator,
                           subsample, method, hll_seed, gate_seed);
}
