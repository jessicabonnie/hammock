#include <doctest/doctest.h>
#include "hammock/hll_sketch.hpp"
#include "hammock/xxhash.hpp"

#include <stdexcept>
#include <string>

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

// ---------------------------------------------------------------------------
// intersection_size / union_with. Prior to these there was no test anywhere
// that checked a containment value against a known ground-truth intersection —
// the Python-side self-pair assertions pass by construction (the union
// registers are bit-identical to the operand's, so the two Ertl estimates
// cancel exactly) and so exercise column wiring, not the estimator.
// ---------------------------------------------------------------------------

namespace {
// Build a sketch over the half-open key range [lo, hi).
HLLSketch range_sketch(size_t precision, int lo, int hi) {
    HLLSketch s(precision);
    for (int i = lo; i < hi; i++) {
        s.add(xxhash::hash64("key_" + std::to_string(i)));
    }
    return s;
}
}  // namespace

TEST_CASE("HLLSketch intersection_size recovers a known overlap") {
    // |A| = 20000 over [0,20000), |B| = 20000 over [10000,30000):
    // ground-truth |A n B| = 10000, so containment_AB = containment_BA = 0.5.
    const HLLSketch a = range_sketch(16, 0, 20000);
    const HLLSketch b = range_sketch(16, 10000, 30000);

    const double inter = a.intersection_size(b);
    CHECK(inter == doctest::Approx(10000.0).epsilon(0.10));

    const double c_ab = inter / a.cardinality();
    const double c_ba = inter / b.cardinality();
    CHECK(c_ab == doctest::Approx(0.5).epsilon(0.10));
    CHECK(c_ba == doctest::Approx(0.5).epsilon(0.10));

    // The reconstruction documented in docs/jaccard-definitional-gap.md:
    // J = 1/(1/C_AB + 1/C_BA - 1). True set Jaccard here is 10000/30000.
    const double j_ie = 1.0 / (1.0 / c_ab + 1.0 / c_ba - 1.0);
    CHECK(j_ie == doctest::Approx(1.0 / 3.0).epsilon(0.15));
}

TEST_CASE("HLLSketch containment is asymmetric when the sets are") {
    // A strictly contains B: |A| = 40000, |B| = 10000, |A n B| = |B|.
    // containment_BA (= |A n B| / |B|) is exactly 1.0 by register dominance;
    // containment_AB (= |A n B| / |A|) should land near 0.25.
    const HLLSketch a = range_sketch(16, 0, 40000);
    const HLLSketch b = range_sketch(16, 0, 10000);

    const double inter = a.intersection_size(b);
    const double c_ab = inter / a.cardinality();
    const double c_ba = inter / b.cardinality();

    CHECK(c_ba == doctest::Approx(1.0).epsilon(0.02));
    CHECK(c_ab == doctest::Approx(0.25).epsilon(0.15));
    CHECK(c_ab < c_ba);
}

TEST_CASE("HLLSketch disjoint sets give containment near zero, unlike register-equality Jaccard") {
    // The two estimators are on different scales: inclusion-exclusion has no
    // chance-agreement floor, register-equality does. This is the measured
    // divergence documented in docs/jaccard-definitional-gap.md.
    const HLLSketch a = range_sketch(14, 0, 200000);
    const HLLSketch b = range_sketch(14, 200000, 400000);

    const double c_ab = a.intersection_size(b) / a.cardinality();
    CHECK(c_ab < 0.02);
    // ... while the register-equality Jaccard sits on its floor, well above 0.
    CHECK(a.jaccard_similarity(b) > 0.10);
}

TEST_CASE("HLLSketch intersection_size never returns a negative estimate") {
    // Inclusion-exclusion is a difference of large estimates and goes negative
    // for disjoint inputs about half the time; the >= 0 clamp must hold.
    for (int trial = 0; trial < 8; trial++) {
        const int base = trial * 100000;
        const HLLSketch a = range_sketch(12, base, base + 30000);
        const HLLSketch b = range_sketch(12, base + 50000, base + 80000);
        CHECK(a.intersection_size(b) >= 0.0);
    }
}

TEST_CASE("HLLSketch mismatched precision throws instead of overreading") {
    // union_with indexes the operand's registers over *this* sketch's count,
    // so an unguarded p=18 vs p=8 pair is a heap overread.
    HLLSketch big(18), small(8);
    for (int i = 0; i < 100; i++) {
        big.add(xxhash::hash64("x_" + std::to_string(i)));
        small.add(xxhash::hash64("x_" + std::to_string(i)));
    }
    CHECK_THROWS_AS(big.union_with(small), std::runtime_error);
    CHECK_THROWS_AS(big.intersection_size(small), std::runtime_error);
    CHECK_THROWS_AS(big.jaccard_similarity(small), std::runtime_error);
    // ...and in the other direction, where the overread would be a write.
    CHECK_THROWS_AS(small.union_with(big), std::runtime_error);
    CHECK_THROWS_AS(small.intersection_size(big), std::runtime_error);
}

TEST_CASE("HLLSketch union_with cardinality matches a directly-built union") {
    const HLLSketch a = range_sketch(16, 0, 20000);
    const HLLSketch b = range_sketch(16, 10000, 30000);
    const HLLSketch direct = range_sketch(16, 0, 30000);

    auto u = a.union_with(b);
    CHECK(u->cardinality() == doctest::Approx(direct.cardinality()).epsilon(1e-9));
}
