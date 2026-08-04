#include <doctest/doctest.h>
#include "hammock/bed_parser.hpp"

TEST_CASE("normalize_chromosome strips chr/CHR/Chr prefix") {
    CHECK(normalize_chromosome("chr1") == "1");
    CHECK(normalize_chromosome("CHR22") == "22");
    CHECK(normalize_chromosome("Chr X") == " X");
    CHECK(normalize_chromosome("1") == "1");
    CHECK(normalize_chromosome("chrX") == "X");
}

TEST_CASE("is_header_or_blank flags empties, comments, track/browser lines") {
    CHECK(is_header_or_blank(""));
    CHECK(is_header_or_blank("# comment"));
    CHECK(is_header_or_blank("@SQ"));
    CHECK(is_header_or_blank("track name=foo"));
    CHECK(is_header_or_blank("browser position chr1"));
    CHECK_FALSE(is_header_or_blank("chr1\t100\t200"));
}

TEST_CASE("parse_bed_line: minimal 3-column line") {
    std::string chr;
    int64_t start = 0, end = 0;
    REQUIRE(parse_bed_line("chr1\t100\t200", chr, start, end));
    // The raw name is kept, prefix and all -- bed_parser.cpp documents why:
    // hashing is byte-sensitive, so normalizing here would change the HLL
    // register state relative to the reference Python. (This assertion read
    // `== "1"` until 0.7.0, which is what normalize_chromosome does, not what
    // parse_bed_line does; the C++ suite is off by default so it never ran.)
    CHECK(chr == "chr1");
    CHECK(start == 100);
    CHECK(end == 200);
}

TEST_CASE("parse_bed_line: rejects malformed coordinates") {
    std::string chr;
    int64_t start = 0, end = 0;
    CHECK_FALSE(parse_bed_line("chr1\t-1\t100", chr, start, end));   // negative start
    CHECK_FALSE(parse_bed_line("chr1\t100\t100", chr, start, end));  // empty interval
    CHECK_FALSE(parse_bed_line("chr1\t200\t100", chr, start, end));  // inverted
    CHECK_FALSE(parse_bed_line("not_a_bed_line", chr, start, end));
}
