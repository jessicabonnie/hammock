#include <doctest/doctest.h>
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

TEST_CASE("add_interval_points_to_sketch: mixed-stride is deterministic per chromosome") {
    HLLSketch a(12), b(12);
    const size_t sa = add_interval_points_to_sketch("1", 100, 1100, a, "-", 0.25,
                                                    SubBMethod::MixedStride, 42);
    const size_t sb = add_interval_points_to_sketch("1", 100, 1100, b, "-", 0.25,
                                                    SubBMethod::MixedStride, 42);
    CHECK(sa == sb);
    // 1000 positions at p=0.25 → ~250 sampled.
    CHECK(sa > 200);
    CHECK(sa < 300);
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
