#include <doctest/doctest.h>
#include <string>
#include <vector>
#include "hammock/hll_sketch.hpp"
#include "hammock/stride.hpp"

TEST_CASE("add_interval_points_to_sketch: subsample=1.0 samples every position") {
    HLLSketch sk(12);
    const size_t sampled = add_interval_points_to_sketch("1", 100, 200, sk, "-",
                                                         /*subsample=*/1.0,
                                                         SubBMethod::HashThreshold,
                                                         /*seed=*/42);
    CHECK(sampled == 100);
}

TEST_CASE("add_interval_points_to_sketch: hash-threshold subsample~=0.5 keeps roughly half") {
    HLLSketch sk(14);
    const size_t sampled = add_interval_points_to_sketch("1", 0, 10000, sk, "-",
                                                         /*subsample=*/0.5,
                                                         SubBMethod::HashThreshold,
                                                         /*seed=*/42);
    // For 10k positions with p=0.5, expect ~5000. Allow ±5% slack.
    CHECK(sampled > 4500);
    CHECK(sampled < 5500);
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

TEST_CASE("add_interval_points_to_sketch: single-hash gates with high 32 bits") {
    HLLSketch sk(14);
    const size_t sampled = add_interval_points_to_sketch("1", 0, 10000, sk, "-",
                                                         /*subsample=*/0.25,
                                                         SubBMethod::SingleHash,
                                                         /*seed=*/42);
    // Random uniform gate at p=0.25, 10k positions → ~2500. Allow ±10% slack.
    CHECK(sampled > 2200);
    CHECK(sampled < 2800);
}
