#ifndef HAMMOCK_BED_PARSER_HPP
#define HAMMOCK_BED_PARSER_HPP

#include <cstdint>
#include <string>
#include <vector>

std::string normalize_chromosome(const std::string& chr);

bool is_header_or_blank(const std::string& line);

bool parse_bed_line(const std::string& line, std::string& chr, int64_t& start, int64_t& end,
                    int64_t& count, int peak_height_column = -1);

std::vector<std::string> read_filepath_list(const std::string& list_file);

#endif
