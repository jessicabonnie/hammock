#include "../hll/hll.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// XXHash64 - Fast, high-quality 64-bit hash function
// ============================================================================

namespace xxhash {
    constexpr uint64_t PRIME64_1 = 0x9E3779B185EBCA87ULL;
    constexpr uint64_t PRIME64_2 = 0xC2B2AE3D27D4EB4FULL;
    constexpr uint64_t PRIME64_3 = 0x165667B19E3779F9ULL;
    constexpr uint64_t PRIME64_4 = 0x85EBCA77C2B2AE63ULL;
    constexpr uint64_t PRIME64_5 = 0x27D4EB2F165667C5ULL;

    inline uint64_t rotl64(uint64_t x, int r) {
        return (x << r) | (x >> (64 - r));
    }

    inline uint64_t avalanche(uint64_t h64) {
        h64 ^= h64 >> 33;
        h64 *= PRIME64_2;
        h64 ^= h64 >> 29;
        h64 *= PRIME64_3;
        h64 ^= h64 >> 32;
        return h64;
    }

    uint64_t hash64(const void* input, size_t len, uint64_t seed = 0) {
        const uint8_t* p = (const uint8_t*)input;
        const uint8_t* const end = p + len;
        uint64_t h64;

        if (len >= 32) {
            const uint8_t* const limit = end - 32;
            uint64_t v1 = seed + PRIME64_1 + PRIME64_2;
            uint64_t v2 = seed + PRIME64_2;
            uint64_t v3 = seed + 0;
            uint64_t v4 = seed - PRIME64_1;

            do {
                uint64_t k1, k2, k3, k4;
                memcpy(&k1, p, 8); p += 8;
                memcpy(&k2, p, 8); p += 8;
                memcpy(&k3, p, 8); p += 8;
                memcpy(&k4, p, 8); p += 8;

                v1 = rotl64(v1 + k1 * PRIME64_2, 31) * PRIME64_1;
                v2 = rotl64(v2 + k2 * PRIME64_2, 31) * PRIME64_1;
                v3 = rotl64(v3 + k3 * PRIME64_2, 31) * PRIME64_1;
                v4 = rotl64(v4 + k4 * PRIME64_2, 31) * PRIME64_1;
            } while (p <= limit);

            h64 = rotl64(v1, 1) + rotl64(v2, 7) + rotl64(v3, 12) + rotl64(v4, 18);

            v1 *= PRIME64_2; v1 = rotl64(v1, 31); v1 *= PRIME64_1;
            h64 ^= v1; h64 = h64 * PRIME64_1 + PRIME64_4;

            v2 *= PRIME64_2; v2 = rotl64(v2, 31); v2 *= PRIME64_1;
            h64 ^= v2; h64 = h64 * PRIME64_1 + PRIME64_4;

            v3 *= PRIME64_2; v3 = rotl64(v3, 31); v3 *= PRIME64_1;
            h64 ^= v3; h64 = h64 * PRIME64_1 + PRIME64_4;

            v4 *= PRIME64_2; v4 = rotl64(v4, 31); v4 *= PRIME64_1;
            h64 ^= v4; h64 = h64 * PRIME64_1 + PRIME64_4;
        } else {
            h64 = seed + PRIME64_5;
        }

        h64 += (uint64_t)len;

        while (p + 8 <= end) {
            uint64_t k1;
            memcpy(&k1, p, 8);
            k1 *= PRIME64_2;
            k1 = rotl64(k1, 31);
            k1 *= PRIME64_1;
            h64 ^= k1;
            h64 = rotl64(h64, 27) * PRIME64_1 + PRIME64_4;
            p += 8;
        }

        if (p + 4 <= end) {
            uint32_t k1;
            memcpy(&k1, p, 4);
            h64 ^= (uint64_t)k1 * PRIME64_1;
            h64 = rotl64(h64, 23) * PRIME64_2 + PRIME64_3;
            p += 4;
        }

        while (p < end) {
            h64 ^= (*p++) * PRIME64_5;
            h64 = rotl64(h64, 11) * PRIME64_1;
        }

        return avalanche(h64);
    }

    inline uint64_t hash64(const std::string& str, uint64_t seed = 0) {
        return hash64(str.data(), str.size(), seed);
    }
}  // namespace xxhash

// ============================================================================
// THREAD MANAGEMENT
// ============================================================================

// Intelligently choose number of threads based on workload
int choose_optimal_threads(size_t num_files, size_t num_primary_files, int mode) {
#ifdef _OPENMP
    // Get system defaults
    int max_threads = omp_get_max_threads();
    
    // If user has set OMP_NUM_THREADS, respect it
    char* omp_env = std::getenv("OMP_NUM_THREADS");
    if (omp_env != nullptr) {
        return max_threads;  // User has explicitly set thread count
    }
    
    // Calculate total number of comparisons
    size_t total_comparisons = num_files * num_primary_files;
    size_t total_files = num_files + num_primary_files;
    
    // Heuristics for choosing thread count:
    // Mode A: Scale primarily with file count (fast per file)
    // Mode B/C: Use more threads even for few files (slow per file, parallelizes within file)
    
    int optimal_threads;
    
    if (mode == 0) {  // MODE_A - interval-based, fast
        if (total_comparisons <= 4) {
            optimal_threads = std::min(4, max_threads);
        } else if (total_comparisons <= 16) {
            optimal_threads = std::min(8, max_threads);
        } else if (total_comparisons <= 64) {
            optimal_threads = std::min(16, max_threads);
        } else {
            optimal_threads = std::min(24, max_threads);
        }
    } else {  // MODE_B or MODE_C - point-based, slow, benefits from more parallelism
        // Use more threads for Mode B/C since they have:
        // 1. File-level parallelism (processing multiple files)
        // 2. Interval-level parallelism (processing intervals within each file)
        if (total_files <= 2) {
            // Even with few files, use more threads for interval parallelism
            optimal_threads = std::min(16, max_threads);
        } else if (total_files <= 8) {
            optimal_threads = std::min(24, max_threads);
        } else if (total_files <= 16) {
            // Medium file counts: still 24 (will process sequentially for good internal parallelism)
            optimal_threads = std::min(24, max_threads);
        } else {
            // Many files: increase threads for file-level parallelism
            optimal_threads = std::min(32, max_threads);
        }
    }
    
    return optimal_threads;
#else
    return 1;  // No OpenMP support
#endif
}

// ============================================================================
// SHARED UTILITY FUNCTIONS (used by both Mode A and Mode B)
// ============================================================================

// Helper function to remove "chr" prefix from chromosome name
std::string normalize_chromosome(const std::string& chr) {
    if (chr.length() >= 3 && 
        (chr[0] == 'c' || chr[0] == 'C') &&
        (chr[1] == 'h' || chr[1] == 'H') &&
        (chr[2] == 'r' || chr[2] == 'R')) {
        return chr.substr(3);
    }
    return chr;
}

// Check if a line is a header or blank
bool is_header_or_blank(const std::string& line) {
    if (line.empty()) return true;
    if (line[0] == '#') return true;
    if (line[0] == '@') return true;
    if (line.find("track") == 0) return true;
    if (line.find("browser") == 0) return true;
    return false;
}

// Parse a BED line and extract chr, start, end
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
    
    chr = normalize_chromosome(chr_raw);
    return true;
}

// Read list of filepaths from a text file
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

// Calculate Jaccard similarity between two sketches
// Use register-based method which is empirically most accurate for HLL
// This avoids the fundamental errors in HLL cardinality-based intersection/union estimation
double calculate_jaccard(hll::hll_t& sketch1, hll::hll_t& sketch2) {
    // Register-based Jaccard compares register patterns directly
    // This method has been empirically validated and avoids cardinality estimation errors
    return sketch1.jaccard_similarity_registers(sketch2);
}

// ============================================================================
// MODE A: INTERVAL-BASED COMPARISON
// ============================================================================

std::string create_interval_string(const std::string& chr, int64_t start, int64_t end, 
                                   const std::string& separator = "-") {
    return chr + separator + std::to_string(start) + separator + std::to_string(end) + separator + "A";
}

size_t process_bed_file_mode_a(const std::string& filepath, hll::hll_t& sketch, 
                               const std::string& separator = "-", bool verbose = false) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    size_t interval_count = 0;
    std::string line;
    std::string chr;
    int64_t start, end;
    
    while (std::getline(file, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        
        if (parse_bed_line(line, chr, start, end)) {
            std::string interval_str = create_interval_string(chr, start, end, separator);
            
            // Hash the interval string using XXHash for proper distribution
            uint64_t hash_val = xxhash::hash64(interval_str);
            
            // Use add() directly to avoid double-hashing
            sketch.add(hash_val);
            interval_count++;
            
            if (verbose && interval_count % 10000 == 0) {
                std::cerr << "Processed " << interval_count << " intervals...\r" << std::flush;
            }
        }
    }
    
    if (verbose && interval_count > 0) {
        std::cerr << "Processed " << interval_count << " intervals total.       " << std::endl;
    }
    
    file.close();
    return interval_count;
}

// ============================================================================
// MODE B: POINT-BASED COMPARISON
// ============================================================================

std::string create_point_string(const std::string& chr, int64_t pos, 
                                const std::string& separator = "-") {
    return chr + separator + std::to_string(pos) + separator + "B";
}

// Hash a string to a 32-bit value for subsampling decisions
uint32_t hash_string_32(const std::string& s) {
    // Simple but effective 32-bit hash
    uint32_t hash = 2166136261u;  // FNV offset basis
    for (char c : s) {
        hash ^= static_cast<uint32_t>(static_cast<unsigned char>(c));
        hash *= 16777619u;  // FNV prime
    }
    return hash;
}

size_t process_bed_file_mode_b(const std::string& filepath, hll::hll_t& sketch,
                               const std::string& separator = "-", 
                               double subsample = 1.0,
                               bool verbose = false) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    // Calculate threshold for subsampling
    uint32_t threshold = static_cast<uint32_t>(subsample * 4294967295.0);  // 2^32 - 1
    bool do_subsample = (subsample < 1.0);
    
    // Read all intervals first
    struct Interval {
        std::string chr;
        int64_t start, end;
    };
    std::vector<Interval> intervals;
    
    std::string line;
    std::string chr;
    int64_t start, end;
    
    while (std::getline(file, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if (parse_bed_line(line, chr, start, end)) {
            intervals.push_back({chr, start, end});
        }
    }
    file.close();
    
    size_t total_points = 0;
    size_t sampled_points = 0;
    
    // Parallel processing of intervals
    // Use static scheduling for deterministic results
    // NOTE: Using critical section around add() because HLL operator+= has bugs
#pragma omp parallel for schedule(static) reduction(+:total_points,sampled_points)
    for (size_t i = 0; i < intervals.size(); i++) {
        const auto& interval = intervals[i];
        
        // Add each point in the interval [start, end)
        for (int64_t pos = interval.start; pos < interval.end; pos++) {
            std::string point_str = create_point_string(interval.chr, pos, separator);
            total_points++;
            
            // Deterministic subsampling: hash the point to decide if we keep it
            if (do_subsample) {
                uint32_t point_hash = hash_string_32(point_str);
                if (point_hash > threshold) {
                    continue;  // Skip this point
                }
            }
            
            // Hash the point string for HLL using XXHash for proper distribution
            uint64_t hash_val = xxhash::hash64(point_str);
            
            // CRITICAL SECTION: HLL add() is thread-safe with -DTHREADSAFE
            sketch.add(hash_val);
            sampled_points++;
        }
    }
    
    if (verbose && intervals.size() > 0) {
        if (do_subsample) {
            std::cerr << "Processed " << intervals.size() << " intervals, " 
                     << sampled_points << "/" << total_points 
                     << " points (subB=" << subsample << ").       " << std::endl;
        } else {
            std::cerr << "Processed " << intervals.size() << " intervals, " 
                     << total_points << " points total.       " << std::endl;
        }
    }
    
    return sampled_points;
}

// ============================================================================
// MODE C: COMBINED INTERVAL + POINT COMPARISON
// ============================================================================

size_t process_bed_file_mode_c(const std::string& filepath, hll::hll_t& sketch,
                               const std::string& separator = "-",
                               double subsample = 1.0,
                               double expA = 0.0,
                               bool verbose = false) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    // Calculate threshold for point subsampling
    uint32_t threshold = static_cast<uint32_t>(subsample * 4294967295.0);
    bool do_subsample = (subsample < 1.0);
    
    // Calculate interval expansion factor: 10^expA
    size_t interval_expansions = 1;
    if (expA > 0) {
        interval_expansions = static_cast<size_t>(std::pow(10.0, expA));
        if (interval_expansions < 1) {
            interval_expansions = 1;
        }
    }
    
    // Read all intervals first
    struct Interval {
        std::string chr;
        int64_t start, end;
    };
    std::vector<Interval> intervals;
    
    std::string line;
    std::string chr;
    int64_t start, end;
    
    while (std::getline(file, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if (parse_bed_line(line, chr, start, end)) {
            intervals.push_back({chr, start, end});
        }
    }
    file.close();
    
    size_t total_interval_elements = 0;
    size_t total_points = 0;
    size_t sampled_points = 0;
    
    // Parallel processing of intervals
    // Use static scheduling for deterministic results
    // NOTE: Using thread-safe add() because HLL operator+= has bugs
#pragma omp parallel for schedule(static) reduction(+:total_interval_elements,total_points,sampled_points)
    for (size_t i = 0; i < intervals.size(); i++) {
        const auto& interval = intervals[i];
        
        // Add the interval with expansions (Mode A component)
        std::string base_interval_str = create_interval_string(interval.chr, interval.start, interval.end, separator);
        
        // Hash intervals using XXHash for proper distribution
        for (size_t exp = 0; exp < interval_expansions; exp++) {
            std::string interval_str = base_interval_str + separator + std::to_string(exp);
            uint64_t interval_hash = xxhash::hash64(interval_str);
            // CRITICAL SECTION: HLL add() is thread-safe with -DTHREADSAFE
            sketch.add(interval_hash);
            total_interval_elements++;
        }
        
        // Add points (Mode B component)
        for (int64_t pos = interval.start; pos < interval.end; pos++) {
            std::string point_str = create_point_string(interval.chr, pos, separator);
            total_points++;
            
            // Deterministic subsampling for points
            if (do_subsample) {
                uint32_t point_hash = hash_string_32(point_str);
                if (point_hash > threshold) {
                    continue;
                }
            }
            
            // Hash the point string for HLL using XXHash for proper distribution
            uint64_t hash_val = xxhash::hash64(point_str);
            
            // CRITICAL SECTION: HLL add() is thread-safe with -DTHREADSAFE
            sketch.add(hash_val);
            sampled_points++;
        }
    }
    
    if (verbose && intervals.size() > 0) {
        std::cerr << "Processed " << intervals.size() << " intervals";
        if (expA > 0) {
            std::cerr << " (" << total_interval_elements << " expanded elements, expA=" << expA << ")";
        }
        if (do_subsample) {
            std::cerr << ", " << sampled_points << "/" << total_points 
                     << " points (subB=" << subsample << ")";
        } else {
            std::cerr << ", " << total_points << " points";
        }
        std::cerr << ".       " << std::endl;
    }
    
    return total_interval_elements + sampled_points;  // Total elements added
}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

#ifndef HAMMOCK_TEST_MODE

enum Mode {
    MODE_A,  // Interval-based
    MODE_B,  // Point-based
    MODE_C   // Combined (intervals + points)
};

void print_usage(const char* program_name) {
    std::cerr << "Hammock - BED File Jaccard Similarity Calculator\n\n";
    std::cerr << "Usage:\n";
    std::cerr << "  " << program_name << " <files_list> <primary_list> [options]\n\n";
    std::cerr << "Arguments:\n";
    std::cerr << "  files_list    : Text file containing paths to BED files (one per line)\n";
    std::cerr << "  primary_list  : Text file containing paths to primary BED files\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  --mode <A|B|C>         : Comparison mode (default: A)\n";
    std::cerr << "                           A = interval-based (chr-start-end)\n";
    std::cerr << "                           B = point-based (all bases in intervals)\n";
    std::cerr << "                           C = combined (intervals + points)\n";
    std::cerr << "  --subB <float>         : Subsampling rate for points (0.0-1.0, default: 1.0)\n";
    std::cerr << "                           Used in Mode B and Mode C\n";
    std::cerr << "  --expA <float>         : Interval expansion exponent (default: 0)\n";
    std::cerr << "                           Adds 10^expA versions of each interval (Mode C only)\n";
    std::cerr << "                           Increases interval weight relative to points\n";
    std::cerr << "                           Examples: 0.5 = 3x, 1.0 = 10x, 2.0 = 100x\n";
    std::cerr << "  --precision, -p <int>  : HyperLogLog precision (default: 18)\n";
    std::cerr << "  --separator, -s <str>  : Separator for strings (default: \"-\")\n";
    std::cerr << "  --output, -o <prefix>  : Output file prefix (default: hammock)\n";
    std::cerr << "                           Output: {prefix}_hll_p{precision}_jacc{mode}[_expA{expA}][_B{subB}].csv\n";
    std::cerr << "  --verbose, -v          : Print verbose progress information\n";
    std::cerr << "  --help, -h             : Show this help message\n\n";
    std::cerr << "Examples:\n";
    std::cerr << "  # Mode A: Compare intervals only\n";
    std::cerr << "  " << program_name << " files.txt primary.txt --mode A -o results\n";
    std::cerr << "    → results_hll_p18_jaccA.csv\n\n";
    std::cerr << "  # Mode B: Compare points only\n";
    std::cerr << "  " << program_name << " files.txt primary.txt --mode B -o results\n";
    std::cerr << "    → results_hll_p18_jaccB.csv\n\n";
    std::cerr << "  # Mode C: Compare both intervals and points\n";
    std::cerr << "  " << program_name << " files.txt primary.txt --mode C -o results\n";
    std::cerr << "    → results_hll_p18_jaccC.csv\n\n";
    std::cerr << "  # Mode C with interval expansion (expA=2.5 → 316x)\n";
    std::cerr << "  " << program_name << " files.txt primary.txt --mode C --expA 2.5 -o results\n";
    std::cerr << "    → results_hll_p18_jaccC_expA2.50.csv\n\n";
    std::cerr << "  # Mode C with both expansion and subsampling\n";
    std::cerr << "  " << program_name << " files.txt primary.txt --mode C --expA 2 --subB 0.1 -o results\n";
    std::cerr << "    → results_hll_p18_jaccC_expA2.00_B0.10.csv\n";
}

int main(int argc, char* argv[]) {
    // Default parameters
    std::string files_list;
    std::string primary_list;
    Mode mode = MODE_A;
    size_t precision = 18;
    std::string separator = "-";
    std::string output_prefix = "hammock";
    double subB = 1.0;
    double expA = 0.0; // Default interval expansion exponent
    bool verbose = false;
    
    // Parse command line arguments
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    files_list = argv[1];
    primary_list = argv[2];
    
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = true;
        } else if (arg == "--mode" && i + 1 < argc) {
            std::string mode_str = argv[++i];
            if (mode_str == "A" || mode_str == "a") {
                mode = MODE_A;
            } else if (mode_str == "B" || mode_str == "b") {
                mode = MODE_B;
            } else if (mode_str == "C" || mode_str == "c") {
                mode = MODE_C;
            } else {
                std::cerr << "Error: Invalid mode '" << mode_str << "'. Must be A, B, or C." << std::endl;
                return 1;
            }
        } else if (arg == "--subB" && i + 1 < argc) {
            subB = std::stod(argv[++i]);
            if (subB < 0.0 || subB > 1.0) {
                std::cerr << "Error: --subB must be between 0.0 and 1.0" << std::endl;
                return 1;
            }
        } else if (arg == "--expA" && i + 1 < argc) {
            expA = std::stod(argv[++i]);
            if (expA < 0.0) {
                std::cerr << "Error: --expA must be >= 0.0" << std::endl;
                return 1;
            }
        } else if ((arg == "--precision" || arg == "-p") && i + 1 < argc) {
            precision = std::stoul(argv[++i]);
        } else if ((arg == "--separator" || arg == "-s") && i + 1 < argc) {
            separator = argv[++i];
        } else if ((arg == "--output" || arg == "-o" || arg == "--outprefix") && i + 1 < argc) {
            output_prefix = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    try {
        // Read file lists
        if (verbose) {
            std::cerr << "Reading file lists..." << std::endl;
            if (mode == MODE_A) {
                std::cerr << "Mode: A (interval-based)" << std::endl;
            } else if (mode == MODE_B) {
                std::cerr << "Mode: B (point-based)" << std::endl;
            } else {
                std::cerr << "Mode: C (combined intervals + points)" << std::endl;
            }
            if ((mode == MODE_B || mode == MODE_C) && subB < 1.0) {
                std::cerr << "Subsampling: " << (subB * 100.0) << "% of points" << std::endl;
            }
            if (mode == MODE_C && expA > 0.0) {
                size_t expansions = static_cast<size_t>(std::pow(10.0, expA));
                std::cerr << "Interval expansion: 10^" << expA << " = " << expansions 
                         << " versions per interval" << std::endl;
            }
        }
        
        std::vector<std::string> files = read_filepath_list(files_list);
        std::vector<std::string> primary_files = read_filepath_list(primary_list);
        
        if (files.empty() || primary_files.empty()) {
            std::cerr << "Error: One or both file lists are empty" << std::endl;
            return 1;
        }
        
        // Intelligently set number of threads based on workload and mode
#ifdef _OPENMP
        int max_available_threads = omp_get_max_threads();
        int optimal_threads = choose_optimal_threads(files.size(), primary_files.size(), mode);
        omp_set_num_threads(optimal_threads);
        
        // Enable nested parallelism for Mode B/C (file-level + interval-level)
        omp_set_max_active_levels(2);
#else
        int optimal_threads = 1;
#endif
        
        if (verbose) {
            std::cerr << "Found " << files.size() << " files and " 
                     << primary_files.size() << " primary files" << std::endl;
            std::cerr << "Using precision: " << precision << std::endl;
        }
        
        // Create sketches for all files
        std::vector<hll::hll_t> sketches(files.size(), hll::hll_t(precision));
        std::vector<hll::hll_t> primary_sketches;
        
        if (verbose) {
            std::cerr << "\nProcessing files..." << std::endl;
#ifdef _OPENMP
            std::cerr << "Using " << optimal_threads << " threads";
            if (max_available_threads != optimal_threads) {
                std::cerr << " (auto-selected from " << max_available_threads << " available)";
            }
            std::cerr << std::endl;
#endif
        }
        
        // Start timing for file processing
        auto file_proc_start = std::chrono::high_resolution_clock::now();
        
        // Parallel file processing strategy:
        // - Mode A: Always parallelize at file level (fast per file)
        // - Mode B/C with < 8 files: Sequential (maximize internal parallelism)
        // - Mode B/C with 8-24 files: Parallel only if 3+ threads per file available
        // - Mode B/C with > 24 files: Always parallel (sequential too slow)
        bool parallelize_files = (mode == MODE_A) || 
                                 (files.size() > 24) ||
                                 (files.size() >= 8 && optimal_threads >= (int)files.size() * 3);
        
        if (verbose) {
            std::cerr << "File processing strategy: " 
                     << (parallelize_files ? "PARALLEL" : "SEQUENTIAL")
                     << " (files=" << files.size() 
                     << ", threads=" << optimal_threads
                     << ", ratio=" << (double)optimal_threads / files.size() << ")" << std::endl;
        }
        
        if (parallelize_files) {
#pragma omp parallel for schedule(dynamic)
            for (size_t i = 0; i < files.size(); i++) {
                if (verbose) {
#pragma omp critical
                    {
                        std::cerr << "Processing file " << (i + 1) << "/" << files.size() 
                                 << ": " << files[i] << std::endl;
                    }
                }
                
                size_t count;
                
                if (mode == MODE_A) {
                    count = process_bed_file_mode_a(files[i], sketches[i], separator, false);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(files[i], sketches[i], separator, subB, false);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(files[i], sketches[i], separator, subB, expA, false);
                }
                
                if (verbose) {
#pragma omp critical
                    {
                        std::string count_desc;
                        if (mode == MODE_A) {
                            count_desc = "Intervals";
                        } else if (mode == MODE_B) {
                            count_desc = "Points";
                        } else {
                            count_desc = "Elements (intervals + points)";
                        }
                        std::cerr << "  " << count_desc << ": " << count 
                                 << ", Estimated cardinality: " << sketches[i].report_ertl_improved() << std::endl;
                    }
                }
            }
        } else {
            // Process files sequentially, allowing full parallelism within each file
            for (size_t i = 0; i < files.size(); i++) {
                if (verbose) {
                    std::cerr << "Processing file " << (i + 1) << "/" << files.size() 
                             << ": " << files[i] << std::endl;
                }
                
                size_t count;
                
                if (mode == MODE_A) {
                    count = process_bed_file_mode_a(files[i], sketches[i], separator, verbose);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(files[i], sketches[i], separator, subB, verbose);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(files[i], sketches[i], separator, subB, expA, verbose);
                }
                
                if (verbose) {
                    std::string count_desc;
                    if (mode == MODE_A) {
                        count_desc = "Intervals";
                    } else if (mode == MODE_B) {
                        count_desc = "Points";
                    } else {
                        count_desc = "Elements (intervals + points)";
                    }
                    std::cerr << "  " << count_desc << ": " << count 
                             << ", Estimated cardinality: " << sketches[i].report_ertl_improved() << std::endl;
                }
            }
        }
        
        // End timing for file processing
        auto file_proc_end = std::chrono::high_resolution_clock::now();
        
        if (verbose) {
            std::cerr << "\nProcessing primary files..." << std::endl;
        }
        
        // Pre-allocate primary sketches vector
        primary_sketches.resize(primary_files.size(), hll::hll_t(precision));
        
        // Start timing for primary file processing
        auto primary_proc_start = std::chrono::high_resolution_clock::now();
        
        // Parallel primary file processing (same strategy as main files)
        bool parallelize_primary = (mode == MODE_A) || 
                                   (primary_files.size() > 24) ||
                                   (primary_files.size() >= 8 && optimal_threads >= (int)primary_files.size() * 3);
        
        if (verbose) {
            std::cerr << "Primary file processing strategy: " 
                     << (parallelize_primary ? "PARALLEL" : "SEQUENTIAL")
                     << " (primary_files=" << primary_files.size() 
                     << ", threads=" << optimal_threads
                     << ", ratio=" << (double)optimal_threads / primary_files.size() << ")" << std::endl;
        }
        
        if (parallelize_primary) {
#pragma omp parallel for schedule(dynamic)
            for (size_t i = 0; i < primary_files.size(); i++) {
                if (verbose) {
#pragma omp critical
                    {
                        std::cerr << "Processing primary file " << (i + 1) << "/" 
                                 << primary_files.size() << ": " << primary_files[i] << std::endl;
                    }
                }
                
                size_t count;
                
                if (mode == MODE_A) {
                    count = process_bed_file_mode_a(primary_files[i], primary_sketches[i], separator, false);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(primary_files[i], primary_sketches[i], separator, subB, false);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(primary_files[i], primary_sketches[i], separator, subB, expA, false);
                }
                
                if (verbose) {
#pragma omp critical
                    {
                        std::string count_desc;
                        if (mode == MODE_A) {
                            count_desc = "Intervals";
                        } else if (mode == MODE_B) {
                            count_desc = "Points";
                        } else {
                            count_desc = "Elements (intervals + points)";
                        }
                        std::cerr << "  " << count_desc << ": " << count 
                                 << ", Estimated cardinality: " << primary_sketches[i].report_ertl_improved() << std::endl;
                    }
                }
            }
        } else {
            // Process primary files sequentially for better interval parallelism
            for (size_t i = 0; i < primary_files.size(); i++) {
                if (verbose) {
                    std::cerr << "Processing primary file " << (i + 1) << "/" 
                             << primary_files.size() << ": " << primary_files[i] << std::endl;
                }
                
                size_t count;
                
                if (mode == MODE_A) {
                    count = process_bed_file_mode_a(primary_files[i], primary_sketches[i], separator, verbose);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(primary_files[i], primary_sketches[i], separator, subB, verbose);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(primary_files[i], primary_sketches[i], separator, subB, expA, verbose);
                }
                
                if (verbose) {
                    std::string count_desc;
                    if (mode == MODE_A) {
                        count_desc = "Intervals";
                    } else if (mode == MODE_B) {
                        count_desc = "Points";
                    } else {
                        count_desc = "Elements (intervals + points)";
                    }
                    std::cerr << "  " << count_desc << ": " << count 
                             << ", Estimated cardinality: " << primary_sketches[i].report_ertl_improved() << std::endl;
                }
            }
        }
        
        // End timing for primary file processing
        auto primary_proc_end = std::chrono::high_resolution_clock::now();
        
        // Calculate all pairwise Jaccard similarities
        if (verbose) {
            std::cerr << "\nCalculating Jaccard similarities..." << std::endl;
        }
        
        // Start timing for comparisons
        auto comparison_start = std::chrono::high_resolution_clock::now();
        
        // Generate mode string
        std::string mode_str;
        if (mode == MODE_A) mode_str = "A";
        else if (mode == MODE_B) mode_str = "B";
        else mode_str = "C";
        
        // Generate output filename using Python hammock convention
        std::string output_file = output_prefix + "_hll_p" + std::to_string(precision) + "_jacc" + mode_str;
        
        // Add expA if present (and not default)
        if (mode == MODE_C && expA > 0.0) {
            std::ostringstream expA_stream;
            expA_stream << std::fixed << std::setprecision(2) << expA;
            output_file += "_expA" + expA_stream.str();
        }
        
        // Add subB if not default
        if ((mode == MODE_B || mode == MODE_C) && subB < 1.0) {
            std::ostringstream subB_stream;
            subB_stream << std::fixed << std::setprecision(2) << subB;
            output_file += "_B" + subB_stream.str();
        }
        
        output_file += ".csv";
        
        std::ofstream outfile(output_file);
        if (!outfile.is_open()) {
            throw std::runtime_error("Cannot open output file: " + output_file);
        }
        
        // Write CSV header (matching Python format)
        outfile << "bed1,bed2,sketch_type,mode,precision,jaccard,intersection,union,cardinality1,cardinality2";
        if (mode == MODE_C) {
            outfile << ",expA,subB";
        }
        outfile << "\n";
        
        size_t total_comparisons = files.size() * primary_files.size();
        
        // Pre-compute all comparisons in parallel
        struct ComparisonResult {
            size_t i, j;
            double jaccard, card1, card2, intersection, union_size;
        };
        std::vector<ComparisonResult> results(total_comparisons);
        
#pragma omp parallel for schedule(dynamic) collapse(2)
        for (size_t i = 0; i < files.size(); i++) {
            for (size_t j = 0; j < primary_files.size(); j++) {
                size_t idx = i * primary_files.size() + j;
                
                results[idx].i = i;
                results[idx].j = j;
                
                // Use inclusion-exclusion Jaccard with MLE cardinalities (should be most accurate)
                results[idx].jaccard = calculate_jaccard(sketches[i], primary_sketches[j]);
                
                // Individual cardinalities using MLE for better accuracy
                results[idx].card1 = sketches[i].report_ertl_improved();
                results[idx].card2 = primary_sketches[j].report_ertl_improved();
                
                // Intersection via MIN of registers (higher error than Jaccard)
                results[idx].intersection = hll::intersection_size(sketches[i], primary_sketches[j]);
                
                // Union via MAX of registers using MLE for better accuracy
                // Using operator+ which creates a copy (operator+= now fixed but keeping for consistency)
                hll::hll_t union_sketch = hll::operator+(sketches[i], primary_sketches[j]);
                results[idx].union_size = union_sketch.report_ertl_improved();
                
                // NOTE: Due to different estimation methods, jaccard ≠ intersection/union
                // The register-based Jaccard is most accurate; use it as ground truth
                
                if (verbose && (idx + 1) % 10 == 0) {
#pragma omp critical
                    {
                        std::cerr << "Progress: " << (idx + 1) << "/" << total_comparisons 
                                 << " comparisons...\r" << std::flush;
                    }
                }
            }
        }
        
        // Write results sequentially to maintain order
        for (const auto& result : results) {
            outfile << files[result.i] << ","
                   << primary_files[result.j] << ","
                   << "hyperloglog," << mode_str << ","
                   << precision << ","
                   << result.jaccard << ","
                   << result.intersection << ","
                   << result.union_size << ","
                   << result.card1 << ","
                   << result.card2;
            
            // Add mode C specific parameters
            if (mode == MODE_C) {
                outfile << "," << expA << "," << subB;
            }
            
            outfile << "\n";
        }
        
        outfile.close();
        
        // End timing for comparisons
        auto comparison_end = std::chrono::high_resolution_clock::now();
        
        // Calculate durations in seconds
        double file_proc_time = std::chrono::duration<double>(file_proc_end - file_proc_start).count();
        double primary_proc_time = std::chrono::duration<double>(primary_proc_end - primary_proc_start).count();
        double sketch_creation_time = file_proc_time + primary_proc_time;
        double comparison_time = std::chrono::duration<double>(comparison_end - comparison_start).count();
        double total_time = sketch_creation_time + comparison_time;
        
        if (verbose) {
            std::cerr << "\nTiming breakdown:" << std::endl;
            std::cerr << "  File processing:         " << std::fixed << std::setprecision(4) << file_proc_time << " s" << std::endl;
            std::cerr << "  Primary file processing: " << std::fixed << std::setprecision(4) << primary_proc_time << " s" << std::endl;
            std::cerr << "  Sketch creation total:   " << std::fixed << std::setprecision(4) << sketch_creation_time << " s" << std::endl;
            std::cerr << "  Comparisons:             " << std::fixed << std::setprecision(4) << comparison_time << " s" << std::endl;
            std::cerr << "  Total time:              " << std::fixed << std::setprecision(4) << total_time << " s" << std::endl;
            std::cerr << "\nCompleted! Results written to: " << output_file 
                     << "                    " << std::endl;
        } else {
            std::cout << "Results written to: " << output_file << std::endl;
        }
        
        // Always output timing in parseable format for benchmarking (to stderr)
        std::cerr << "TIMING: sketch_creation=" << std::fixed << std::setprecision(6) << sketch_creation_time
                  << " comparison=" << comparison_time
                  << " total=" << total_time << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
#endif // HAMMOCK_TEST_MODE

