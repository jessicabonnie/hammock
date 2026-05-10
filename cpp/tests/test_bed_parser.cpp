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
    int64_t start = 0, end = 0, count = 0;
    REQUIRE(parse_bed_line("chr1\t100\t200", chr, start, end, count));
    CHECK(chr == "1");
    CHECK(start == 100);
    CHECK(end == 200);
    CHECK(count == 1);
}

TEST_CASE("parse_bed_line: rejects malformed coordinates") {
    std::string chr;
    int64_t start = 0, end = 0, count = 0;
    CHECK_FALSE(parse_bed_line("chr1\t-1\t100", chr, start, end, count));   // negative start
    CHECK_FALSE(parse_bed_line("chr1\t100\t100", chr, start, end, count));  // empty interval
    CHECK_FALSE(parse_bed_line("chr1\t200\t100", chr, start, end, count));  // inverted
    CHECK_FALSE(parse_bed_line("not_a_bed_line", chr, start, end, count));
}

TEST_CASE("parse_bed_line: reads count column") {
    std::string chr;
    int64_t start = 0, end = 0, count = 0;
    REQUIRE(parse_bed_line("chr2\t10\t20\tname\t42", chr, start, end, count, 5));
    CHECK(count == 42);
    REQUIRE(parse_bed_line("chr2\t10\t20\tname\t-3", chr, start, end, count, 5));
    CHECK(count == 0);   // negative clamped to 0
}
