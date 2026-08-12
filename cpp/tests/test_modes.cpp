#include <doctest/doctest.h>
#include "hammock/hll_sketch.hpp"
#include "hammock/processing_modes.hpp"

#include <cstdio>
#include <fstream>
#include <string>

namespace {
std::string write_temp_bed(const char* contents) {
    char path_buf[L_tmpnam];
    std::tmpnam(path_buf);
    std::ofstream f(path_buf);
    f << contents;
    return std::string(path_buf);
}
}

TEST_CASE("Mode A: same BED file ⇒ Jaccard 1.0") {
    const std::string p = write_temp_bed(
        "chr1\t100\t200\n"
        "chr1\t300\t400\n"
        "chr2\t50\t75\n"
    );
    HLLSketch a(14), b(14);
    process_bed_file_mode_a(p, a);
    process_bed_file_mode_a(p, b);
    CHECK(a.reg_eq_similarity(b) == doctest::Approx(1.0).epsilon(1e-9));
    std::remove(p.c_str());
}

TEST_CASE("Mode A: chr/non-chr prefixes normalize identically") {
    const std::string p1 = write_temp_bed(
        "chr1\t100\t200\n"
        "chr2\t300\t400\n"
    );
    const std::string p2 = write_temp_bed(
        "1\t100\t200\n"
        "2\t300\t400\n"
    );
    HLLSketch a(14), b(14);
    process_bed_file_mode_a(p1, a);
    process_bed_file_mode_a(p2, b);
    CHECK(a.reg_eq_similarity(b) == doctest::Approx(1.0).epsilon(1e-9));
    std::remove(p1.c_str());
    std::remove(p2.c_str());
}

TEST_CASE("Mode B: short interval, every position sampled, low cardinality") {
    const std::string p = write_temp_bed(
        "chr1\t0\t100\n"
    );
    HLLSketch s(14);
    process_bed_file_mode_b(p, s, /*hll_seed=*/42, "-", /*subB=*/1.0,
                            SubBMethod::HashThreshold);
    // 100 unique points; HLL approximates cardinality.
    CHECK(s.cardinality() == doctest::Approx(100.0).epsilon(0.1));
    std::remove(p.c_str());
}

TEST_CASE("Mode C with expA=0 ≈ Mode A ∪ Mode B") {
    const std::string p = write_temp_bed(
        "chr1\t0\t100\n"
        "chr1\t500\t600\n"
    );
    HLLSketch a(14), b(14), c(14);
    process_bed_file_mode_a(p, a);
    process_bed_file_mode_b(p, b, /*hll_seed=*/42, "-", /*subB=*/1.0,
                            SubBMethod::HashThreshold);
    process_bed_file_mode_c(p, c, /*hll_seed=*/42, "-",
                            /*subA=*/1.0, /*subB=*/1.0, /*expA=*/0.0,
                            SubBMethod::HashThreshold);
    // Mode C with expA=0 adds intervals (×1) plus all points: should approx
    // equal merging A and B sketches.
    auto union_ab = a.union_with(b);
    HLLSketch* hll_union = dynamic_cast<HLLSketch*>(union_ab.get());
    REQUIRE(hll_union != nullptr);
    CHECK(c.cardinality() == doctest::Approx(hll_union->cardinality()).epsilon(0.05));
    std::remove(p.c_str());
}
