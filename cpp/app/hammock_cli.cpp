// hammock-cpp: standalone benchmark binary. Same algorithm as the Python
// `hammock`, no Python in the loop. Modes A/B/C only — Mode D requires the
// digest library and lives on the Python side.
#include "hammock/bagminhash_sketch.hpp"
#include "hammock/bed_parser.hpp"
#include "hammock/hll_sketch.hpp"
#include "hammock/processing_modes.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Args {
    std::string filepaths_file;
    std::string primary_file;
    std::string mode = "A";
    std::string outprefix = "hammock";
    std::string separator = "\t";   // Python parity default
    size_t precision = 18;
    double subA = 1.0;
    double subB = 1.0;
    double expA = 0.0;
    SubBMethod subB_method = SubBMethod::HashThreshold;
    uint64_t seed = 42;
    uint32_t gate_seed = 31337;
    int peak_height_column = -1;
    int threads = 0;
    bool verbose = false;
};

void print_help(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " <filepaths_file> <primary_file> [options]\n\n"
        "Options:\n"
        "  --mode <A|B|C>          Comparison mode (default: A)\n"
        "  --subA <float>          Subsampling rate for intervals (Mode C, 0..1, default: 1.0)\n"
        "  --subB <float>          Subsampling rate for points (0..1, default: 1.0)\n"
        "  --subB-method <name>    Point-sampling method (default: hash-threshold)\n"
        "                          Values: hash-threshold | mixed-stride | single-hash\n"
        "  --mixed-stride          [deprecated] alias for --subB-method=mixed-stride\n"
        "  --seed <int>            HLL hashing seed for xxh64 (default: 42)\n"
        "  --gate-seed <int>       Seed for the subB gate hash xxh32 and the\n"
        "                          mixed-stride chr->stride hash (default: 31337,\n"
        "                          matches orig hammock). Ignored for single-hash.\n"
        "  --expA <float>          Interval expansion exponent (default: 0)\n"
        "  --precision, -p <int>   HyperLogLog precision 4..24 (default: 18)\n"
        "  --separator, -s <str>   Separator for hashed strings (default: \"-\")\n"
        "  --output, -o <prefix>   Output filename prefix (default: hammock)\n"
        "  --threads, -t <int>     Thread count (default: OpenMP auto)\n"
        "  --peak-height <int>     1-based column index for count weights (Mode A → BagMinHash)\n"
        "  --verbose, -v           Print progress to stderr\n"
        "  --help, -h              Show this help\n";
}

bool parse_args(int argc, char** argv, Args& out) {
    std::vector<std::string> positional;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            print_help(argv[0]);
            std::exit(0);
        } else if (a == "--verbose" || a == "-v") {
            out.verbose = true;
        } else if (a == "--mixed-stride") {
            out.subB_method = SubBMethod::MixedStride;
        } else if (a == "--subB-method" && i + 1 < argc) {
            std::string m = argv[++i];
            if (m == "hash-threshold")      out.subB_method = SubBMethod::HashThreshold;
            else if (m == "mixed-stride")   out.subB_method = SubBMethod::MixedStride;
            else if (m == "single-hash") out.subB_method = SubBMethod::SingleHash;
            else {
                std::cerr << "Error: --subB-method must be one of hash-threshold, "
                             "mixed-stride, single-hash (got '" << m << "')\n";
                return false;
            }
        } else if (a == "--mode" && i + 1 < argc) {
            out.mode = argv[++i];
        } else if (a == "--subA" && i + 1 < argc) {
            out.subA = std::stod(argv[++i]);
        } else if (a == "--subB" && i + 1 < argc) {
            out.subB = std::stod(argv[++i]);
        } else if (a == "--seed" && i + 1 < argc) {
            out.seed = std::stoull(argv[++i]);
        } else if (a == "--gate-seed" && i + 1 < argc) {
            out.gate_seed = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (a == "--expA" && i + 1 < argc) {
            out.expA = std::stod(argv[++i]);
        } else if ((a == "--precision" || a == "-p") && i + 1 < argc) {
            out.precision = std::stoul(argv[++i]);
        } else if ((a == "--separator" || a == "-s") && i + 1 < argc) {
            out.separator = argv[++i];
        } else if ((a == "--output" || a == "-o" || a == "--outprefix") && i + 1 < argc) {
            out.outprefix = argv[++i];
        } else if ((a == "--threads" || a == "-t") && i + 1 < argc) {
            out.threads = std::stoi(argv[++i]);
        } else if (a == "--peak-height" && i + 1 < argc) {
            out.peak_height_column = std::stoi(argv[++i]);
        } else if (!a.empty() && a[0] != '-') {
            positional.push_back(a);
        } else {
            std::cerr << "Error: unknown argument '" << a << "'\n";
            return false;
        }
    }
    if (positional.size() != 2) {
        std::cerr << "Error: expected 2 positional arguments (filepaths_file, primary_file)\n";
        return false;
    }
    out.filepaths_file = positional[0];
    out.primary_file = positional[1];

    if (out.mode != "A" && out.mode != "B" && out.mode != "C") {
        std::cerr << "Error: --mode must be A, B, or C (got '" << out.mode << "')\n";
        return false;
    }
    if (out.subA < 0.0 || out.subA > 1.0) {
        std::cerr << "Error: --subA must be in [0, 1]\n";
        return false;
    }
    if (out.subB < 0.0 || out.subB > 1.0) {
        std::cerr << "Error: --subB must be in [0, 1]\n";
        return false;
    }
    if (out.precision < 4 || out.precision > 24) {
        std::cerr << "Error: --precision must be in [4, 24]\n";
        return false;
    }
    return true;
}

std::unique_ptr<AbstractSketch> make_sketch(const Args& a) {
    if (a.peak_height_column > 0) {
        // BagMinHash for count-weighted similarity. Default 128 hashes
        // matches the Python CLI's --num_hashes default.
        return std::make_unique<BagMinHashSketch>(128, a.seed);
    }
    return std::make_unique<HLLSketch>(a.precision);
}

void process_one(const std::string& path, AbstractSketch& sketch, const Args& a) {
    if (a.mode == "A") {
        process_bed_file_mode_a(path, sketch, a.seed, a.separator, a.peak_height_column, a.verbose);
    } else if (a.mode == "B") {
        process_bed_file_mode_b(path, sketch, a.seed, a.separator, a.subB, a.subB_method,
                                a.gate_seed, a.peak_height_column, a.verbose);
    } else {  // "C"
        process_bed_file_mode_c(path, sketch, a.seed, a.separator, a.subA, a.subB, a.expA,
                                a.subB_method, a.gate_seed, a.peak_height_column, a.verbose);
    }
}

std::string outprefix_with_suffix(const Args& a) {
    std::string out = a.outprefix;
    out += (a.peak_height_column > 0) ? "_bmh_jacc" : "_hll_p" + std::to_string(a.precision) + "_jacc";
    out += a.mode;
    if (a.expA > 0) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", a.expA);
        out += "_expA";
        out += buf;
    }
    if (a.subB != 1.0) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", a.subB);
        out += "_B";
        out += buf;
    }
    return out + ".csv";
}

std::string basename_of(const std::string& path) {
    const size_t s = path.find_last_of('/');
    return (s == std::string::npos) ? path : path.substr(s + 1);
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args)) {
        return 2;
    }

#ifdef _OPENMP
    if (args.threads > 0) {
        omp_set_num_threads(args.threads);
    }
    if (args.verbose) {
        std::cerr << "OpenMP threads: " << omp_get_max_threads() << "\n";
    }
#endif

    const std::vector<std::string> queries = read_filepath_list(args.filepaths_file);
    const std::vector<std::string> refs = read_filepath_list(args.primary_file);
    if (queries.empty() || refs.empty()) {
        std::cerr << "Error: empty filepath list(s)\n";
        return 2;
    }

    if (args.verbose) {
        std::cerr << "Sketching " << queries.size() << " query and " << refs.size()
                  << " primary files (mode " << args.mode << ").\n";
    }

    auto t0 = std::chrono::steady_clock::now();

    std::vector<std::unique_ptr<AbstractSketch>> qsk, rsk;
    qsk.reserve(queries.size());
    rsk.reserve(refs.size());
    for (const auto& p : queries) {
        auto s = make_sketch(args);
        process_one(p, *s, args);
        qsk.push_back(std::move(s));
    }
    for (const auto& p : refs) {
        auto s = make_sketch(args);
        process_one(p, *s, args);
        rsk.push_back(std::move(s));
    }

    auto t1 = std::chrono::steady_clock::now();

    const std::string out_path = outprefix_with_suffix(args);
    FILE* fp = std::fopen(out_path.c_str(), "w");
    if (!fp) {
        std::cerr << "Error: could not open '" << out_path << "' for writing\n";
        return 1;
    }
    std::fprintf(fp, "query\treference\tjaccard_similarity\n");

    const size_t n = queries.size();
    const size_t m = refs.size();
    std::vector<double> matrix(n * m);
#pragma omp parallel for collapse(2) schedule(static)
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            matrix[i * m + j] = qsk[i]->jaccard_similarity(*rsk[j]);
        }
    }

    for (size_t i = 0; i < n; i++) {
        const std::string ql = basename_of(queries[i]);
        for (size_t j = 0; j < m; j++) {
            std::fprintf(fp, "%s\t%s\t%.17g\n",
                         ql.c_str(), basename_of(refs[j]).c_str(), matrix[i * m + j]);
        }
    }
    std::fclose(fp);

    auto t2 = std::chrono::steady_clock::now();

    if (args.verbose) {
        const auto sketch_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        const auto pair_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cerr << "Sketching: " << sketch_ms << " ms\n"
                  << "Pairwise+write: " << pair_ms << " ms\n"
                  << "Wrote " << out_path << "\n";
    }
    return 0;
}
