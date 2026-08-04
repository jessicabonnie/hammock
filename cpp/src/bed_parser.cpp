#include "hammock/bed_parser.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::string normalize_chromosome(const std::string& chr) {
    if (chr.length() >= 3 &&
        (chr[0] == 'c' || chr[0] == 'C') &&
        (chr[1] == 'h' || chr[1] == 'H') &&
        (chr[2] == 'r' || chr[2] == 'R')) {
        return chr.substr(3);
    }
    return chr;
}

bool is_header_or_blank(const std::string& line) {
    if (line.empty()) return true;
    if (line[0] == '#') return true;
    if (line[0] == '@') return true;
    if (line.find("track") == 0) return true;
    if (line.find("browser") == 0) return true;
    return false;
}

bool parse_bed_line(const std::string& line, std::string& chr, int64_t& start, int64_t& end) {
    if (is_header_or_blank(line)) {
        return false;
    }

    std::istringstream iss(line);
    std::string chr_raw;

    if (!(iss >> chr_raw >> start >> end)) {
        return false;
    }

    if (start < 0 || end <= start) {
        return false;
    }

    // Match the reference Python contract: keep the raw chromosome name
    // (including a "chr"/"CHR"/"Chr" prefix if the BED file has one).
    // Hashing is sensitive to byte content so any normalization changes the
    // resulting HLL register state.
    chr = chr_raw;

    return true;
}

std::vector<std::string> read_filepath_list(const std::string& list_file) {
    std::vector<std::string> filepaths;
    std::ifstream file(list_file);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open filepath list: " + list_file);
    }

    std::string line;
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

        if (!line.empty() && line[0] != '#') {
            filepaths.push_back(line);
        }
    }

    return filepaths;
}
