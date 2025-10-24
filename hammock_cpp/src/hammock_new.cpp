#include "../include/abstract_sketch.h"
#include "../include/hll_sketch.h"
#include "../include/bagminhash_sketch.h"
#include "../include/bed_parser.h"
#include "../include/processing_modes.h"
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

#ifdef _OPENMP
#include <omp.h>
#endif

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
    std::cerr << "                           Output: {prefix}_{hll|bmh}_p{precision}_jacc{mode}[_expA{expA}][_B{subB}].csv\n";
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
    std::cerr << "  # BagMinHash with count column (column 4)\n";
    std::cerr << "  " << program_name << " files.txt primary.txt --mode A --count 4 -o results\n";
    std::cerr << "    → results_bmh_p18_jaccA.csv\n";
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
        std::string sketch_type = (count_column >= 0) ? "bmh" : "hll";
        std::string output_file = output_prefix + "_" + sketch_type + "_p" + std::to_string(precision) + "_jacc" + mode_str;
        
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
