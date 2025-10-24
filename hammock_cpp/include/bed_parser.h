#ifndef BED_PARSER_H
#define BED_PARSER_H

#include <string>
#include <vector>
#include <cstdint>

// Helper function to remove "chr" prefix from chromosome name
std::string normalize_chromosome(const std::string& chr);

// Check if a line is a header or blank
bool is_header_or_blank(const std::string& line);

// Parse a BED line and extract chr, start, end, and optionally count
bool parse_bed_line(const std::string& line, std::string& chr, int64_t& start, int64_t& end, 
                    int64_t& count, int count_column = -1);

// Read list of filepaths from a text file
std::vector<std::string> read_filepath_list(const std::string& list_file);

#endif // BED_PARSER_H
