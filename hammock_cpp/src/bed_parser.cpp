#include "../include/bed_parser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
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

bool parse_bed_line(const std::string& line, std::string& chr, int64_t& start, int64_t& end, 
                    int64_t& count, int count_column) {
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
    
    chr = normalize_chromosome(chr_raw);
    
    // Extract count from specified column if requested
    if (count_column > 0) {
        count = 1;  // Default count if column doesn't exist or is invalid
        std::string token;
        int current_column = 3;  // We've already read chr, start, end (columns 1-3)
        
        while (iss >> token && current_column < count_column) {
            current_column++;
        }
        
        if (current_column == count_column) {
            try {
                count = std::stoll(token);
                if (count < 0) {
                    count = 0;  // Negative counts don't make sense
                }
            } catch (const std::exception&) {
                count = 1;  // Default to 1 if parsing fails
            }
        }
    } else {
        count = 1;  // Default count when no count column specified
    }
    
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
    
    file.close();
    return filepaths;
}
