#include <doctest/doctest.h>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "hammock/abstract_sketch.hpp"
#include "hammock/hll_sketch.hpp"
#include "hammock/stride.hpp"
#include "hammock/xxhash.hpp"

namespace {

// Records the exact values offered to a sketch. This tests sampling semantics
// directly rather than inferring them through HLL cardinality noise.
class RecordingSketch final : public AbstractSketch {
public:
    std::vector<uint64_t> values;

    void add(uint64_t hash_val) override { values.push_back(hash_val); }
    double reg_eq_similarity(const AbstractSketch&) const override { return 0.0; }
    double cardinality() const override { return static_cast<double>(values.size()); }
    double intersection_size(const AbstractSketch&) const override { return 0.0; }
    std::unique_ptr<AbstractSketch> union_with(const AbstractSketch&) const override {
        return std::make_unique<RecordingSketch>(*this);
    }
    void clear() override { values.clear(); }
};

std::vector<uint64_t> expected_threshold_hashes(
    const std::string& chr, int64_t start, int64_t end,
    const std::string& separator, double rate, uint64_t hll_seed,
    uint32_t gate_seed) {
    std::vector<uint64_t> expected;
    const uint32_t threshold = static_cast<uint32_t>(rate * 4294967295.0);
    for (int64_t pos = start; pos < end; ++pos) {
        const std::string point = chr + separator + std::to_string(pos);
        if (xxhash::hash32(point, gate_seed) <= threshold) {
            expected.push_back(xxhash::hash64(point, hll_seed));
        }
    }
    return expected;
}

std::vector<uint64_t> expected_single_hashes(
    const std::string& chr, int64_t start, int64_t end,
    const std::string& separator, double rate, uint64_t hll_seed) {
    std::vector<uint64_t> expected;
    const uint32_t threshold = static_cast<uint32_t>(rate * 4294967295.0);
    for (int64_t pos = start; pos < end; ++pos) {
        const std::string point = chr + separator + std::to_string(pos);
        const uint64_t hash = xxhash::hash64(point, hll_seed);
        if (static_cast<uint32_t>(hash >> 32) <= threshold) expected.push_back(hash);
    }
    return expected;
}

std::vector<uint64_t> record_points(
    const std::string& chr, int64_t start, int64_t end,
    const std::string& separator, double rate, SubBMethod method,
    uint64_t hll_seed, uint32_t gate_seed) {
    RecordingSketch sketch;
    const size_t count = add_interval_points_to_sketch(
        chr, start, end, sketch, separator, rate, method, hll_seed, gate_seed);
    if (count != sketch.values.size()) {
        throw std::runtime_error("sampler count disagrees with recorded values");
    }
    return sketch.values;
}

}  // namespace

TEST_CASE("add_interval_points_to_sketch: subsample=1.0 samples every position") {
    HLLSketch sk(12);
    const size_t sampled = add_interval_points_to_sketch("1", 100, 200, sk, "-",
                                                         /*subsample=*/1.0,
                                                         SubBMethod::HashThreshold,
                                                         /*seed=*/42);
    CHECK(sampled == 100);
}

TEST_CASE("add_interval_points_to_sketch: hash-threshold matches exact two-hash oracle") {
    // 995..1005 crosses a digit-width rollover in the optimized ASCII counter.
    // The long chromosome forces the std::string fallback path.
    for (const std::string chr : {std::string("chr7"), std::string(80, 'x')}) {
        for (const uint32_t gate_seed : {uint32_t{0}, uint32_t{31337}}) {
            for (const uint64_t hll_seed : {uint64_t{0}, uint64_t{42}}) {
                const auto actual = record_points(
                    chr, 995, 1006, "-", 0.3, SubBMethod::HashThreshold,
                    hll_seed, gate_seed);
                const auto expected = expected_threshold_hashes(
                    chr, 995, 1006, "-", 0.3, hll_seed, gate_seed);
                CHECK(actual == expected);
            }
        }
    }
}

TEST_CASE("add_interval_points_to_sketch: hash-threshold matches python-xxhash vector") {
    // Generated independently with python-xxhash 3.5.0.
    const std::vector<uint64_t> expected = {
        0x358d00f9fd29cc30ULL, 0x13a48a90dccb85f8ULL,
        0x2211aae0ec685360ULL, 0x7416f7610e2bb0fcULL,
    };
    CHECK(record_points("chr7", 995, 1006, "-", 0.3,
                        SubBMethod::HashThreshold, 42, 31337) == expected);
}

TEST_CASE("add_interval_points_to_sketch: hash-threshold separates gate and HLL seeds") {
    const auto gate_a = record_points("chr2", 0, 256, "\t", 0.3,
                                      SubBMethod::HashThreshold, 42, 7);
    const auto gate_b = record_points("chr2", 0, 256, "\t", 0.3,
                                      SubBMethod::HashThreshold, 42, 8);
    const auto hash_b = record_points("chr2", 0, 256, "\t", 0.3,
                                      SubBMethod::HashThreshold, 43, 7);
    CHECK(gate_a != gate_b);
    CHECK(gate_a.size() == hash_b.size());
    CHECK(gate_a != hash_b);
}

TEST_CASE("add_interval_points_to_sketch: mixed-stride-v1 is deterministic per chromosome") {
    HLLSketch a(12), b(12);
    const size_t sa = add_interval_points_to_sketch("1", 100, 1100, a, "-", 0.25,
                                                    SubBMethod::MixedStrideV1, 42);
    const size_t sb = add_interval_points_to_sketch("1", 100, 1100, b, "-", 0.25,
                                                    SubBMethod::MixedStrideV1, 42);
    CHECK(sa == sb);
    // 1000 positions at p=0.25 → ~250 sampled.
    CHECK(sa > 200);
    CHECK(sa < 300);
}

TEST_CASE("add_interval_points_to_sketch: mixed-stride-v1 chooses one chromosome-wide gap") {
    // At p=0.3, S0=floor(1/p)=3 and S1=ceil(1/p)=4. With the default
    // gate seed, chromosome 1 selects S0 and chrX selects S1. An interval
    // whose length is divisible by both strides makes the expected counts
    // exact regardless of the chromosome-specific residue.
    HLLSketch s0_sketch(12), s1_sketch(12);
    const size_t sampled_s0 = add_interval_points_to_sketch(
        "1", 0, 1200, s0_sketch, "-", 0.3,
        SubBMethod::MixedStrideV1, 42);
    const size_t sampled_s1 = add_interval_points_to_sketch(
        "chrX", 0, 1200, s1_sketch, "-", 0.3,
        SubBMethod::MixedStrideV1, 42);

    CHECK(sampled_s0 == 400);
    CHECK(sampled_s1 == 300);
}

TEST_CASE("add_interval_points_to_sketch: mixed-stride-v1 phase is independent of interval boundaries") {
    // Splitting a chromosome span into adjacent BED intervals must retain the
    // same genomic coordinates as sketching that span as one interval.
    HLLSketch whole(12), partitioned(12);
    const size_t whole_count = add_interval_points_to_sketch(
        "chrX", 0, 1200, whole, "-", 0.3,
        SubBMethod::MixedStrideV1, 42);

    size_t partitioned_count = 0;
    partitioned_count += add_interval_points_to_sketch(
        "chrX", 0, 317, partitioned, "-", 0.3,
        SubBMethod::MixedStrideV1, 42);
    partitioned_count += add_interval_points_to_sketch(
        "chrX", 317, 811, partitioned, "-", 0.3,
        SubBMethod::MixedStrideV1, 42);
    partitioned_count += add_interval_points_to_sketch(
        "chrX", 811, 1200, partitioned, "-", 0.3,
        SubBMethod::MixedStrideV1, 42);

    CHECK(partitioned_count == whole_count);
    CHECK(partitioned.registers() == whole.registers());
}

TEST_CASE("add_interval_points_to_sketch: mixed-stride-v2 mixes gaps within a chromosome") {
    HLLSketch legacy(12), v2(12);
    const size_t legacy_count = add_interval_points_to_sketch(
        "chrX", 0, 3600, legacy, "-", 0.3, SubBMethod::MixedStrideV1, 42);
    const size_t v2_count = add_interval_points_to_sketch(
        "chrX", 0, 3600, v2, "-", 0.3, SubBMethod::MixedStrideV2, 42);

    // v1 uses either 3 or 4 everywhere (1200 or 900 points). v2's mean
    // gap is 10/3, so it retains 1080 points up to phase-boundary error.
    CHECK(v2_count >= 1078);
    CHECK(v2_count <= 1082);
    CHECK(v2_count != legacy_count);
}

TEST_CASE("add_interval_points_to_sketch: mixed-stride-v2 preserves legacy integral grids") {
    for (const double rate : {0.1, 0.01}) {
        HLLSketch legacy(12), v2(12);
        const size_t legacy_count = add_interval_points_to_sketch(
            "chr7", -135, 3017, legacy, "-", rate, SubBMethod::MixedStrideV1, 42);
        const size_t v2_count = add_interval_points_to_sketch(
            "chr7", -135, 3017, v2, "-", rate, SubBMethod::MixedStrideV2, 42);
        CHECK(v2_count == legacy_count);
        CHECK(v2.registers() == legacy.registers());
    }
}

TEST_CASE("add_interval_points_to_sketch: mixed-stride-v2 is partition invariant for long names") {
    const std::string long_chr(80, 'x');
    HLLSketch whole(12), partitioned(12);
    const size_t whole_count = add_interval_points_to_sketch(
        long_chr, -51, 8000, whole, "-", 0.15, SubBMethod::MixedStrideV2, 42);
    const size_t split_count =
        add_interval_points_to_sketch(long_chr, -51, 1777, partitioned, "-", 0.15,
                                      SubBMethod::MixedStrideV2, 42) +
        add_interval_points_to_sketch(long_chr, 1777, 5291, partitioned, "-", 0.15,
                                      SubBMethod::MixedStrideV2, 42) +
        add_interval_points_to_sketch(long_chr, 5291, 8000, partitioned, "-", 0.15,
                                      SubBMethod::MixedStrideV2, 42);
    CHECK(split_count == whole_count);
    CHECK(partitioned.registers() == whole.registers());
}

TEST_CASE("add_interval_points_to_sketch: zero subB retains no points") {
    HLLSketch sketch(12);
    CHECK(add_interval_points_to_sketch("chr1", 0, 100, sketch, "-", 0.0,
                                        SubBMethod::MixedStrideV2, 42) == 0);
    CHECK(sketch.registers() == std::vector<uint8_t>(sketch.registers().size(), 0));
}

TEST_CASE("add_interval_points_to_sketch: single-hash matches exact high-bits oracle") {
    for (const std::string chr : {std::string("chr7"), std::string(80, 'x')}) {
        for (const uint64_t hll_seed : {uint64_t{0}, uint64_t{42}}) {
            const auto actual = record_points(
                chr, 995, 1006, "-", 0.3, SubBMethod::SingleHash,
                hll_seed, 31337);
            const auto expected = expected_single_hashes(
                chr, 995, 1006, "-", 0.3, hll_seed);
            CHECK(actual == expected);
        }
    }
}

TEST_CASE("add_interval_points_to_sketch: single-hash matches python-xxhash vector") {
    // Generated independently with python-xxhash 3.5.0.
    const std::vector<uint64_t> expected = {
        0x358d00f9fd29cc30ULL, 0x13a48a90dccb85f8ULL,
        0x2211aae0ec685360ULL, 0x264fdf2951e05adcULL,
        0x156f2b7f080deb24ULL,
    };
    CHECK(record_points("chr7", 995, 1006, "-", 0.3,
                        SubBMethod::SingleHash, 42, 31337) == expected);
}

TEST_CASE("add_interval_points_to_sketch: single-hash uses HLL seed and ignores gate seed") {
    const auto baseline = record_points("chr2", 0, 256, "\t", 0.3,
                                        SubBMethod::SingleHash, 42, 7);
    const auto other_gate = record_points("chr2", 0, 256, "\t", 0.3,
                                          SubBMethod::SingleHash, 42, 999);
    const auto other_hash = record_points("chr2", 0, 256, "\t", 0.3,
                                          SubBMethod::SingleHash, 43, 7);
    CHECK(baseline == other_gate);
    CHECK(baseline != other_hash);
}
