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
#include <memory>
#include <unordered_map>
#include <queue>

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
// ABSTRACT SKETCH INTERFACE
// ============================================================================

// Abstract base class for all sketch types
class AbstractSketch {
public:
    virtual ~AbstractSketch() = default;
    virtual void add(uint64_t hash_val) = 0;
    virtual double jaccard_similarity(const AbstractSketch& other) const = 0;
    virtual double cardinality() const = 0;
    virtual double intersection_size(const AbstractSketch& other) const = 0;
    virtual std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const = 0;
    virtual std::string get_sketch_type() const = 0;
    virtual void clear() = 0;
};

// ============================================================================
// HLL SKETCH WRAPPER
// ============================================================================

class HLLSketch : public AbstractSketch {
private:
    hll::hll_t sketch_;
    
public:
    explicit HLLSketch(size_t precision) : sketch_(precision) {}
    
    void add(uint64_t hash_val) override {
        sketch_.add(hash_val);
    }
    
    double jaccard_similarity(const AbstractSketch& other) const override {
        const HLLSketch* hll_other = dynamic_cast<const HLLSketch*>(&other);
        if (!hll_other) {
            throw std::runtime_error("Cannot compute Jaccard between different sketch types");
        }
        return sketch_.jaccard_similarity_registers(hll_other->sketch_);
    }
    
    double cardinality() const override {
        return const_cast<hll::hll_t&>(sketch_).report_ertl_improved();
    }
    
    double intersection_size(const AbstractSketch& other) const override {
        const HLLSketch* hll_other = dynamic_cast<const HLLSketch*>(&other);
        if (!hll_other) {
            throw std::runtime_error("Cannot compute intersection between different sketch types");
        }
        return hll::intersection_size(sketch_, hll_other->sketch_);
    }
    
    std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const override {
        const HLLSketch* hll_other = dynamic_cast<const HLLSketch*>(&other);
        if (!hll_other) {
            throw std::runtime_error("Cannot compute union between different sketch types");
        }
        auto result = std::make_unique<HLLSketch>(sketch_.get_np());
        result->sketch_ = hll::operator+(sketch_, hll_other->sketch_);
        return result;
    }
    
    std::string get_sketch_type() const override {
        return "hyperloglog";
    }
    
    void clear() override {
        sketch_.clear();
    }
};

// ============================================================================
// BAGMINHASH SKETCH IMPLEMENTATION
// ============================================================================

class BagMinHashSketch : public AbstractSketch {
private:
    size_t num_hashes_;
    std::vector<std::pair<uint64_t, uint64_t>> min_hashes_;  // (hash, count) pairs
    uint64_t seed_;
    
    // Hash function for generating multiple hash values
    uint64_t hash_with_seed(uint64_t value, uint64_t seed) const {
        // Convert uint64_t to string for xxhash
        std::string value_str = std::to_string(value);
        return xxhash::hash64(value_str, seed);
    }
    
public:
    explicit BagMinHashSketch(size_t num_hashes, uint64_t seed = 0) 
        : num_hashes_(num_hashes), min_hashes_(num_hashes, {UINT64_MAX, 0}), seed_(seed) {}
    
    void add(uint64_t hash_val) override {
        // For BagMinHash, we need to track counts, but since we're getting individual hash values,
        // we'll treat each add as a count of 1. In a real implementation, you'd want to batch
        // adds by the same hash value to properly handle counts.
        for (size_t i = 0; i < num_hashes_; ++i) {
            uint64_t hash_seed = seed_ + i;
            uint64_t hashed_value = hash_with_seed(hash_val, hash_seed);
            
            if (hashed_value < min_hashes_[i].first) {
                min_hashes_[i].first = hashed_value;
                min_hashes_[i].second = 1;  // Reset count when new minimum found
            } else if (hashed_value == min_hashes_[i].first) {
                min_hashes_[i].second++;  // Increment count for same minimum
            }
        }
    }
    
    // Add with explicit count (for use with --count argument)
    void add_with_count(uint64_t hash_val, uint64_t count) {
        for (size_t i = 0; i < num_hashes_; ++i) {
            uint64_t hash_seed = seed_ + i;
            uint64_t hashed_value = hash_with_seed(hash_val, hash_seed);
            
            if (hashed_value < min_hashes_[i].first) {
                min_hashes_[i].first = hashed_value;
                min_hashes_[i].second = count;
            } else if (hashed_value == min_hashes_[i].first) {
                min_hashes_[i].second += count;
            }
        }
    }
    
    double jaccard_similarity(const AbstractSketch& other) const override {
        const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
        if (!bmh_other) {
            throw std::runtime_error("Cannot compute Jaccard between different sketch types");
        }
        
        if (min_hashes_.size() != bmh_other->min_hashes_.size()) {
            throw std::runtime_error("BagMinHash sketches must have same number of hashes");
        }
        
        double intersection_sum = 0.0;
        double union_sum = 0.0;
        
        for (size_t i = 0; i < min_hashes_.size(); ++i) {
            if (min_hashes_[i].first == bmh_other->min_hashes_[i].first) {
                // Same minimum hash - intersection is min of counts, union is max
                intersection_sum += std::min(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
                union_sum += std::max(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
            } else {
                // Different minimum hashes - no intersection, union is sum of counts
                union_sum += min_hashes_[i].second + bmh_other->min_hashes_[i].second;
            }
        }
        
        return (union_sum > 0) ? (intersection_sum / union_sum) : 0.0;
    }
    
    double cardinality() const override {
        // For BagMinHash, cardinality estimation is more complex
        // This is a simplified version - in practice you'd want a more sophisticated estimator
        double total_count = 0.0;
        for (const auto& pair : min_hashes_) {
            total_count += pair.second;
        }
        return total_count / num_hashes_;
    }
    
    double intersection_size(const AbstractSketch& other) const override {
        const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
        if (!bmh_other) {
            throw std::runtime_error("Cannot compute intersection between different sketch types");
        }
        
        double intersection_sum = 0.0;
        for (size_t i = 0; i < min_hashes_.size(); ++i) {
            if (min_hashes_[i].first == bmh_other->min_hashes_[i].first) {
                intersection_sum += std::min(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
            }
        }
        return intersection_sum;
    }
    
    std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const override {
        const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
        if (!bmh_other) {
            throw std::runtime_error("Cannot compute union between different sketch types");
        }
        
        auto result = std::make_unique<BagMinHashSketch>(num_hashes_, seed_);
        
        for (size_t i = 0; i < min_hashes_.size(); ++i) {
            if (min_hashes_[i].first < bmh_other->min_hashes_[i].first) {
                result->min_hashes_[i] = min_hashes_[i];
            } else if (bmh_other->min_hashes_[i].first < min_hashes_[i].first) {
                result->min_hashes_[i] = bmh_other->min_hashes_[i];
            } else {
                // Same minimum hash - take the maximum count
                result->min_hashes_[i] = {min_hashes_[i].first, 
                                         std::max(min_hashes_[i].second, bmh_other->min_hashes_[i].second)};
            }
        }
        
        return result;
    }
    
    std::string get_sketch_type() const override {
        return "bagminhash";
    }
    
    void clear() override {
        std::fill(min_hashes_.begin(), min_hashes_.end(), std::make_pair(UINT64_MAX, 0));
    }
};

// ============================================================================
// THREAD MANAGEMENT
// ============================================================================

// Intelligently choose number of threads based on workload
int choose_optimal_threads(size_t num_files, size_t num_primary_files, int mode, int manual_threads = -1) {
#ifdef _OPENMP
    // Get system defaults
    int max_threads = omp_get_max_threads();
    
    // Priority 1: Command-line argument (highest priority)
    if (manual_threads > 0) {
        return std::min(manual_threads, max_threads);  // User-specified via --threads
    }
    
    // Priority 2: OMP_NUM_THREADS environment variable
    char* omp_env = std::getenv("OMP_NUM_THREADS");
    if (omp_env != nullptr) {
        return max_threads;  // User has explicitly set thread count via environment
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

// Parse a BED line and extract chr, start, end, and optionally count
bool parse_bed_line(const std::string& line, std::string& chr, int64_t& start, int64_t& end, 
                    int64_t& count, int count_column = -1) {
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
// Use polymorphic method that works with any sketch type
double calculate_jaccard(const AbstractSketch& sketch1, const AbstractSketch& sketch2) {
    return sketch1.jaccard_similarity(sketch2);
}

// ============================================================================
// MODE A: INTERVAL-BASED COMPARISON
// ============================================================================

std::string create_interval_string(const std::string& chr, int64_t start, int64_t end, 
                                   const std::string& separator = "-") {
    return chr + separator + std::to_string(start) + separator + std::to_string(end) + separator + "A";
}

size_t process_bed_file_mode_a(const std::string& filepath, AbstractSketch& sketch, 
                               const std::string& separator = "-", int count_column = -1, bool verbose = false) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    size_t interval_count = 0;
    std::string line;
    std::string chr;
    int64_t start, end, count;
    
    while (std::getline(file, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        
        if (parse_bed_line(line, chr, start, end, count, count_column)) {
            std::string interval_str = create_interval_string(chr, start, end, separator);
            
            // Hash the interval string using XXHash for proper distribution
            uint64_t hash_val = xxhash::hash64(interval_str);
            
            // Use count-aware adding if count column is specified
            if (count_column > 0) {
                // Try to cast to BagMinHashSketch to use count-aware adding
                BagMinHashSketch* bmh_sketch = dynamic_cast<BagMinHashSketch*>(&sketch);
                if (bmh_sketch) {
                    bmh_sketch->add_with_count(hash_val, count);
                } else {
                    // Fall back to regular add for HLL
                    for (int64_t i = 0; i < count; i++) {
                        sketch.add(hash_val);
                    }
                }
            } else {
                // Use add() directly to avoid double-hashing
                sketch.add(hash_val);
            }
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
uint32_t hash_string_32(const std::string& s, uint64_t seed = 0) {
    // Use XXHash32 for consistency with Python version
    uint64_t hash64 = xxhash::hash64(s, seed);
    return static_cast<uint32_t>(hash64 & 0xFFFFFFFF);
}

size_t process_bed_file_mode_b(const std::string& filepath, AbstractSketch& sketch,
                               const std::string& separator = "-", 
                               double subsample = 1.0,
                               bool mixed_stride = false,
                               uint64_t seed = 0,
                               int count_column = -1,
                               bool verbose = false) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    // Calculate threshold for hash-threshold subsampling
    uint32_t threshold = static_cast<uint32_t>(subsample * 4294967295.0);  // 2^32 - 1
    bool do_subsample = (subsample < 1.0);
    
    // Read all intervals first
    struct Interval {
        std::string chr;
        int64_t start, end, count;
    };
    std::vector<Interval> intervals;
    
    std::string line;
    std::string chr;
    int64_t start, end, count;
    
    while (std::getline(file, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if (parse_bed_line(line, chr, start, end, count, count_column)) {
            intervals.push_back({chr, start, end, count});
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
        
        if (mixed_stride && do_subsample) {
            // Mixed-stride deterministic subsampling
            double p = subsample;
            double inv = 1.0 / p;
            int64_t s0 = static_cast<int64_t>(std::floor(inv));
            int64_t s1 = static_cast<int64_t>(std::ceil(inv));
            if (s0 < 1) s0 = 1;
            if (s1 < 1) s1 = 1;
            
            bool choose_s1 = false;
            if (s0 != s1) {
                double denom = (1.0 / s1) - (1.0 / s0);
                double q = (denom != 0.0) ? (p - 1.0 / s0) / denom : 0.0;
                // Deterministic choice per chromosome
                uint32_t h = hash_string_32(interval.chr, seed);
                double rfloat = h / 4294967295.0;
                choose_s1 = (rfloat < q);
            }
            int64_t S = (s0 != s1 && choose_s1) ? s1 : s0;
            
            // Residue per chromosome
            uint32_t h2 = hash_string_32(interval.chr + "|res", seed);
            int64_t residue = h2 % S;
            
            // First hit in [start, end)
            int64_t offset = (residue - (interval.start % S)) % S;
            if (offset < 0) offset += S;  // Handle negative modulo
            int64_t x0 = interval.start + offset;
            
            for (int64_t pos = x0; pos < interval.end; pos += S) {
                std::string point_str = create_point_string(interval.chr, pos, separator);
                total_points++;  // Count all potential points for reporting
                
                // Hash the point string for HLL using XXHash for proper distribution
                uint64_t hash_val = xxhash::hash64(point_str);
                
                // CRITICAL SECTION: HLL add() is thread-safe with -DTHREADSAFE
                sketch.add(hash_val);
                sampled_points++;
            }
            
            // Add unsampled points to total count for accurate reporting
            total_points += (interval.end - interval.start) - ((interval.end - x0 + S - 1) / S);
            
        } else {
            // Hash-threshold subsampling (original method)
            for (int64_t pos = interval.start; pos < interval.end; pos++) {
                std::string point_str = create_point_string(interval.chr, pos, separator);
                total_points++;
                
                // Deterministic subsampling: hash the point to decide if we keep it
                if (do_subsample) {
                    uint32_t point_hash = hash_string_32(point_str, seed);
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
    }
    
    if (verbose && intervals.size() > 0) {
        if (do_subsample) {
            std::string method = mixed_stride ? "mixed-stride" : "hash-threshold";
            std::cerr << "Processed " << intervals.size() << " intervals, " 
                     << sampled_points << "/" << total_points 
                     << " points (subB=" << subsample << ", method=" << method << ").       " << std::endl;
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

size_t process_bed_file_mode_c(const std::string& filepath, AbstractSketch& sketch,
                               const std::string& separator = "-",
                               double subsample = 1.0,
                               double expA = 0.0,
                               bool mixed_stride = false,
                               uint64_t seed = 0,
                               int count_column = -1,
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
        int64_t start, end, count;
    };
    std::vector<Interval> intervals;
    
    std::string line;
    std::string chr;
    int64_t start, end, count;
    
    while (std::getline(file, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if (parse_bed_line(line, chr, start, end, count, count_column)) {
            intervals.push_back({chr, start, end, count});
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
        if (mixed_stride && do_subsample) {
            // Mixed-stride deterministic subsampling
            double p = subsample;
            double inv = 1.0 / p;
            int64_t s0 = static_cast<int64_t>(std::floor(inv));
            int64_t s1 = static_cast<int64_t>(std::ceil(inv));
            if (s0 < 1) s0 = 1;
            if (s1 < 1) s1 = 1;
            
            bool choose_s1 = false;
            if (s0 != s1) {
                double denom = (1.0 / s1) - (1.0 / s0);
                double q = (denom != 0.0) ? (p - 1.0 / s0) / denom : 0.0;
                // Deterministic choice per chromosome
                uint32_t h = hash_string_32(interval.chr, seed);
                double rfloat = h / 4294967295.0;
                choose_s1 = (rfloat < q);
            }
            int64_t S = (s0 != s1 && choose_s1) ? s1 : s0;
            
            // Residue per chromosome
            uint32_t h2 = hash_string_32(interval.chr + "|res", seed);
            int64_t residue = h2 % S;
            
            // First hit in [start, end)
            int64_t offset = (residue - (interval.start % S)) % S;
            if (offset < 0) offset += S;  // Handle negative modulo
            int64_t x0 = interval.start + offset;
            
            for (int64_t pos = x0; pos < interval.end; pos += S) {
                std::string point_str = create_point_string(interval.chr, pos, separator);
                total_points++;  // Count all potential points for reporting
                
                // Hash the point string for HLL using XXHash for proper distribution
                uint64_t hash_val = xxhash::hash64(point_str);
                
                // CRITICAL SECTION: HLL add() is thread-safe with -DTHREADSAFE
                sketch.add(hash_val);
                sampled_points++;
            }
            
            // Add unsampled points to total count for accurate reporting
            total_points += (interval.end - interval.start) - ((interval.end - x0 + S - 1) / S);
            
        } else {
            // Hash-threshold subsampling (original method)
            for (int64_t pos = interval.start; pos < interval.end; pos++) {
                std::string point_str = create_point_string(interval.chr, pos, separator);
                total_points++;
                
                // Deterministic subsampling for points
                if (do_subsample) {
                    uint32_t point_hash = hash_string_32(point_str, seed);
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
    }
    
    if (verbose && intervals.size() > 0) {
        std::cerr << "Processed " << intervals.size() << " intervals";
        if (expA > 0) {
            std::cerr << " (" << total_interval_elements << " expanded elements, expA=" << expA << ")";
        }
        if (do_subsample) {
            std::string method = mixed_stride ? "mixed-stride" : "hash-threshold";
            std::cerr << ", " << sampled_points << "/" << total_points 
                     << " points (subB=" << subsample << ", method=" << method << ")";
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
    std::cerr << "  --mixed-stride         : Use mixed-stride subsampling strategy (default: hash-threshold)\n";
    std::cerr << "                           Mixed-stride avoids per-point hashing for better performance\n";
    std::cerr << "                           at low subB rates. Both methods are deterministic.\n";
    std::cerr << "  --seed <int>           : Random seed for hashing (default: 0)\n";
    std::cerr << "  --expA <float>         : Interval expansion exponent (default: 0)\n";
    std::cerr << "                           Adds 10^expA versions of each interval (Mode C only)\n";
    std::cerr << "                           Increases interval weight relative to points\n";
    std::cerr << "                           Examples: 0.5 = 3x, 1.0 = 10x, 2.0 = 100x\n";
    std::cerr << "  --precision, -p <int>  : HyperLogLog precision (default: 18)\n";
    std::cerr << "  --separator, -s <str>  : Separator for strings (default: \"-\")\n";
    std::cerr << "  --output, -o <prefix>  : Output file prefix (default: hammock)\n";
    std::cerr << "                           Output: {prefix}_hll_p{precision}_jacc{mode}[_expA{expA}][_B{subB}].csv\n";
    std::cerr << "  --threads, -t <int>    : Number of threads to use (default: auto-select based on workload)\n";
    std::cerr << "                           Overrides OMP_NUM_THREADS environment variable\n";
    std::cerr << "  --count <int>          : Column index (1-based) containing precomputed read counts\n";
    std::cerr << "                           When specified, uses BagMinHash sketching with count values\n";
    std::cerr << "                           When not specified, uses HyperLogLog sketching (default)\n";
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
    std::cerr << "    → results_hll_p18_jaccC_expA2.00_B0.10.csv\n\n";
    std::cerr << "  # Using BagMinHash with count column (column 4)\n";
    std::cerr << "  " << program_name << " files.txt primary.txt --mode A --count 4 -o results\n";
    std::cerr << "    → results_bagminhash_p18_jaccA.csv\n";
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
    bool mixed_stride = false;  // Default to hash-threshold method
    uint64_t seed = 0;  // Default seed
    int manual_threads = -1;  // -1 means auto-select
    int count_column = -1;  // -1 means no count column specified
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
        } else if (arg == "--mixed-stride") {
            mixed_stride = true;
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
        } else if (arg == "--seed" && i + 1 < argc) {
            seed = std::stoull(argv[++i]);
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
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            manual_threads = std::stoi(argv[++i]);
            if (manual_threads < 1) {
                std::cerr << "Error: --threads must be >= 1" << std::endl;
                return 1;
            }
        } else if (arg == "--count" && i + 1 < argc) {
            count_column = std::stoi(argv[++i]);
            if (count_column < 1) {
                std::cerr << "Error: --count must be >= 1 (1-based column index)" << std::endl;
                return 1;
            }
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
                std::string method = mixed_stride ? "mixed-stride" : "hash-threshold";
                std::cerr << "Subsampling: " << (subB * 100.0) << "% of points (method: " << method << ")" << std::endl;
            }
            if (mode == MODE_C && expA > 0.0) {
                size_t expansions = static_cast<size_t>(std::pow(10.0, expA));
                std::cerr << "Interval expansion: 10^" << expA << " = " << expansions 
                         << " versions per interval" << std::endl;
            }
            if (seed != 0) {
                std::cerr << "Random seed: " << seed << std::endl;
            }
        }
        
        std::vector<std::string> files = read_filepath_list(files_list);
        std::vector<std::string> primary_files = read_filepath_list(primary_list);
        
        if (files.empty() || primary_files.empty()) {
            std::cerr << "Error: One or both file lists are empty" << std::endl;
            return 1;
        }
        
        // Optimization: Check if same file list provided for both inputs
        bool same_file_list = (files_list == primary_list);
        if (same_file_list && verbose) {
            std::cerr << "Note: Same file list detected - will reuse sketches for efficiency" << std::endl;
        }
        
        // Intelligently set number of threads based on workload and mode
#ifdef _OPENMP
        int max_available_threads = omp_get_max_threads();
        int optimal_threads = choose_optimal_threads(files.size(), primary_files.size(), mode, manual_threads);
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
        std::vector<std::unique_ptr<AbstractSketch>> sketches;
        std::vector<std::unique_ptr<AbstractSketch>> primary_sketches;
        
        // Create appropriate sketch type based on count column
        if (count_column > 0) {
            // Use BagMinHash when count column is specified
            for (size_t i = 0; i < files.size(); i++) {
                sketches.push_back(std::make_unique<BagMinHashSketch>(precision, seed));
            }
            for (size_t i = 0; i < primary_files.size(); i++) {
                primary_sketches.push_back(std::make_unique<BagMinHashSketch>(precision, seed));
            }
        } else {
            // Use HLL when no count column is specified
            for (size_t i = 0; i < files.size(); i++) {
                sketches.push_back(std::make_unique<HLLSketch>(precision));
            }
            for (size_t i = 0; i < primary_files.size(); i++) {
                primary_sketches.push_back(std::make_unique<HLLSketch>(precision));
            }
        }
        
        if (verbose) {
            std::cerr << "\nProcessing files..." << std::endl;
#ifdef _OPENMP
            std::cerr << "Using " << optimal_threads << " threads";
            if (manual_threads > 0) {
                std::cerr << " (user-specified via --threads)";
            } else {
                char* omp_env = std::getenv("OMP_NUM_THREADS");
                if (omp_env != nullptr) {
                    std::cerr << " (set via OMP_NUM_THREADS)";
                } else if (max_available_threads != optimal_threads) {
                    std::cerr << " (auto-selected from " << max_available_threads << " available)";
                }
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
                    count = process_bed_file_mode_a(files[i], *sketches[i], separator, count_column, false);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(files[i], *sketches[i], separator, subB, mixed_stride, seed, count_column, false);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(files[i], *sketches[i], separator, subB, expA, mixed_stride, seed, count_column, false);
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
                                 << ", Estimated cardinality: " << sketches[i]->cardinality() << std::endl;
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
                    count = process_bed_file_mode_a(files[i], *sketches[i], separator, count_column, verbose);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(files[i], *sketches[i], separator, subB, mixed_stride, seed, count_column, verbose);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(files[i], *sketches[i], separator, subB, expA, mixed_stride, seed, count_column, verbose);
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
                             << ", Estimated cardinality: " << sketches[i]->cardinality() << std::endl;
                }
            }
        }
        
        // End timing for file processing
        auto file_proc_end = std::chrono::high_resolution_clock::now();
        
        // Start timing for primary file processing
        auto primary_proc_start = std::chrono::high_resolution_clock::now();
        
        // Optimization: Skip if same file list (reuse sketches)
        if (same_file_list) {
            if (verbose) {
                std::cerr << "\nSkipping primary file processing (reusing sketches from main files)" << std::endl;
            }
            // Point primary_sketches reference to sketches - but we can't do that with vector
            // Instead, we'll handle this in the comparison loop
        } else {
            if (verbose) {
                std::cerr << "\nProcessing primary files..." << std::endl;
            }
            
            // Pre-allocate primary sketches vector
            primary_sketches.clear();
            if (count_column > 0) {
                for (size_t i = 0; i < primary_files.size(); i++) {
                    primary_sketches.push_back(std::make_unique<BagMinHashSketch>(precision, seed));
                }
            } else {
                for (size_t i = 0; i < primary_files.size(); i++) {
                    primary_sketches.push_back(std::make_unique<HLLSketch>(precision));
                }
            }
        
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
                    count = process_bed_file_mode_a(primary_files[i], *primary_sketches[i], separator, count_column, false);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(primary_files[i], *primary_sketches[i], separator, subB, mixed_stride, seed, count_column, false);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(primary_files[i], *primary_sketches[i], separator, subB, expA, mixed_stride, seed, count_column, false);
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
                                 << ", Estimated cardinality: " << primary_sketches[i]->cardinality() << std::endl;
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
                    count = process_bed_file_mode_a(primary_files[i], *primary_sketches[i], separator, count_column, verbose);
                } else if (mode == MODE_B) {
                    count = process_bed_file_mode_b(primary_files[i], *primary_sketches[i], separator, subB, mixed_stride, seed, count_column, verbose);
                } else {  // MODE_C
                    count = process_bed_file_mode_c(primary_files[i], *primary_sketches[i], separator, subB, expA, mixed_stride, seed, count_column, verbose);
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
                             << ", Estimated cardinality: " << primary_sketches[i]->cardinality() << std::endl;
                }
            }
        }
        }  // End of else block for primary file processing
        
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
                
                // Get the appropriate sketch for comparison
                // If same file list, use sketches[j]; otherwise use primary_sketches[j]
                AbstractSketch& sketch_j = same_file_list ? *sketches[j] : *primary_sketches[j];
                
                // Use polymorphic Jaccard calculation
                results[idx].jaccard = calculate_jaccard(*sketches[i], sketch_j);
                
                // Individual cardinalities using polymorphic method
                results[idx].card1 = sketches[i]->cardinality();
                results[idx].card2 = sketch_j.cardinality();
                
                // Intersection using polymorphic method
                results[idx].intersection = sketches[i]->intersection_size(sketch_j);
                
                // Union using polymorphic method
                auto union_sketch = sketches[i]->union_with(sketch_j);
                results[idx].union_size = union_sketch->cardinality();
                
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
                   << sketches[result.i]->get_sketch_type() << "," << mode_str << ","
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

