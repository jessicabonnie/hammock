#include <doctest/doctest.h>
#include "hammock/hll_sketch.hpp"
#include "hammock/xxhash.hpp"

TEST_CASE("HLLSketch reports zero cardinality when empty") {
    HLLSketch s(14);
    CHECK(s.cardinality() == doctest::Approx(0.0).epsilon(0.01));
}

TEST_CASE("HLLSketch identical inputs => Jaccard 1.0") {
    HLLSketch a(14), b(14);
    for (int i = 0; i < 1000; i++) {
        const uint64_t h = xxhash::hash64(std::to_string(i));
        a.add(h);
        b.add(h);
    }
    CHECK(a.jaccard_similarity(b) == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("HLLSketch disjoint inputs => low Jaccard") {
    HLLSketch a(14), b(14);
    for (int i = 0; i < 1000; i++) {
        a.add(xxhash::hash64("A_" + std::to_string(i)));
        b.add(xxhash::hash64("B_" + std::to_string(i)));
    }
    // Jaccard for two disjoint sets should be near zero. With p=14 and
    // 1000 elements each, give plenty of slack.
    CHECK(a.jaccard_similarity(b) < 0.1);
}

TEST_CASE("HLLSketch determinism: same inputs in different order produce identical state") {
    HLLSketch a(12), b(12);
    for (int i = 0; i < 500; i++) {
        a.add(xxhash::hash64("k_" + std::to_string(i)));
    }
    for (int i = 499; i >= 0; i--) {
        b.add(xxhash::hash64("k_" + std::to_string(i)));
    }
    CHECK(a.cardinality() == doctest::Approx(b.cardinality()));
    CHECK(a.jaccard_similarity(b) == doctest::Approx(1.0).epsilon(1e-9));
}
