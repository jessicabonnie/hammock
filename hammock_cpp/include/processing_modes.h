#ifndef PROCESSING_MODES_H
#define PROCESSING_MODES_H

#include "abstract_sketch.h"
#include <string>
#include <cstdint>

// String creation functions
std::string create_interval_string(const std::string& chr, int64_t start, int64_t end, 
                                   const std::string& separator = "-");
std::string create_point_string(const std::string& chr, int64_t pos, 
                                const std::string& separator = "-");

// Hash function for subsampling decisions
uint32_t hash_string_32(const std::string& s, uint64_t seed = 0);

// Processing functions for different modes
size_t process_bed_file_mode_a(const std::string& filepath, AbstractSketch& sketch, 
                               const std::string& separator = "-", int peak_height_column = -1, bool verbose = false);

size_t process_bed_file_mode_b(const std::string& filepath, AbstractSketch& sketch,
                               const std::string& separator = "-", 
                               double subsample = 1.0,
                               bool mixed_stride = false,
                               uint64_t seed = 0,
                               int peak_height_column = -1,
                               bool verbose = false);

size_t process_bed_file_mode_c(const std::string& filepath, AbstractSketch& sketch,
                               const std::string& separator = "-",
                               double subsample = 1.0,
                               double expA = 0.0,
                               bool mixed_stride = false,
                               uint64_t seed = 0,
                               int peak_height_column = -1,
                               bool verbose = false);

// Jaccard calculation
double calculate_jaccard(const AbstractSketch& sketch1, const AbstractSketch& sketch2);

#endif // PROCESSING_MODES_H
