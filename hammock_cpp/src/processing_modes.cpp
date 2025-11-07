#include "../include/processing_modes.h"
#include "../include/xxhash.h"
#include "../include/bagminhash_sketch.h"
#include "../include/bed_parser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

std::string create_interval_string(const std::string& chr, int64_t start, int64_t end, 
                                   const std::string& separator) {
    return chr + separator + std::to_string(start) + separator + std::to_string(end) + separator + "A";
}

std::string create_point_string(const std::string& chr, int64_t pos, 
                                const std::string& separator) {
    return chr + separator + std::to_string(pos) + separator + "B";
}

uint32_t hash_string_32(const std::string& s, uint64_t seed) {
    // Use XXHash32 for consistency with Python version
    uint64_t hash64 = xxhash::hash64(s, seed);
    return static_cast<uint32_t>(hash64 & 0xFFFFFFFF);
}

double calculate_jaccard(const AbstractSketch& sketch1, const AbstractSketch& sketch2) {
    return sketch1.jaccard_similarity(sketch2);
}

size_t process_bed_file_mode_a(const std::string& filepath, AbstractSketch& sketch, 
                               const std::string& separator, int peak_height_column, bool verbose) {
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
        
        if (parse_bed_line(line, chr, start, end, count, peak_height_column)) {
            std::string interval_str = create_interval_string(chr, start, end, separator);
            
            // Hash the interval string using XXHash for proper distribution
            uint64_t hash_val = xxhash::hash64(interval_str);
            
            // Use count-aware adding if count column is specified
            if (peak_height_column > 0) {
                // Try to cast to BagMinHashSketch to use count-aware adding
                BagMinHashSketch* bmh_sketch = dynamic_cast<BagMinHashSketch*>(&sketch);
                if (bmh_sketch) {
                    // Use normalized weights for BagMinHash to ensure scale-comparable weights
                    bmh_sketch->add_with_normalized_count(hash_val, count);
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

size_t process_bed_file_mode_b(const std::string& filepath, AbstractSketch& sketch,
                               const std::string& separator, 
                               double subsample,
                               bool mixed_stride,
                               uint64_t seed,
                               int peak_height_column,
                               bool verbose) {
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
        if (parse_bed_line(line, chr, start, end, count, peak_height_column)) {
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

size_t process_bed_file_mode_c(const std::string& filepath, AbstractSketch& sketch,
                               const std::string& separator,
                               double subsample,
                               double expA,
                               bool mixed_stride,
                               uint64_t seed,
                               int peak_height_column,
                               bool verbose) {
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
        if (parse_bed_line(line, chr, start, end, count, peak_height_column)) {
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
