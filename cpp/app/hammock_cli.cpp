// hammock-cpp: standalone benchmark binary. Same algorithm as the Python
// `hammock`, no Python in the loop. Modes A/B/C only — Mode D requires the
// digest library and lives on the Python side.
#include "hammock/bed_parser.hpp"
#include "hammock/hll_sketch.hpp"
#include "hammock/processing_modes.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <type_traits>
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
    std::string mode = "B";
    std::string outprefix = "hammock";
    std::string separator = "\t";   // Python parity default
    size_t precision = 18;
    double subA = 1.0;
    double subB = 1.0;
    double expA = 0.0;
    SubBMethod subB_method = SubBMethod::MixedStride;
    uint64_t seed = 42;
    uint32_t gate_seed = 31337;
    int threads = 0;
    bool verbose = false;
    bool metrics = true;
};

// Mirror of runner.py's _jaccard_ie_from_containments. Kept expression-for-
// expression identical (same clamp, same operand order, same divisions) so the
// two programs agree bit-for-bit rather than to a couple of ulp -- the
// cross-tool test asserts exact equality on this.
double jaccard_ie_from_containments(double c_ab, double c_ba) {
    c_ab = std::min(c_ab, 1.0);
    c_ba = std::min(c_ba, 1.0);
    if (!(c_ab > 0.0 && c_ba > 0.0)) return 0.0;
    return 1.0 / (1.0 / c_ab + 1.0 / c_ba - 1.0);
}

// Guarded numeric parsing. Bare std::sto* throws out of main on a
// non-numeric argument -- `hammock-cpp -p abc` used to SIGABRT -- and silently
// truncates trailing garbage (`-p 4.9` parsed as 4). Reject both.
template <typename T>
bool parse_num(const char* name, const char* text, T& out) {
    try {
        size_t pos = 0;
        if constexpr (std::is_floating_point_v<T>) {
            const double v = std::stod(text, &pos);
            if (text[pos] != '\0') { std::cerr << "Error: " << name << " expects a number (got '" << text << "')\n"; return false; }
            out = static_cast<T>(v);
        } else {
            const long long v = std::stoll(text, &pos);
            if (text[pos] != '\0') { std::cerr << "Error: " << name << " expects an integer (got '" << text << "')\n"; return false; }
            if (v < 0 && std::is_unsigned_v<T>) { std::cerr << "Error: " << name << " must not be negative (got " << v << ")\n"; return false; }
            out = static_cast<T>(v);
        }
        return true;
    } catch (const std::exception&) {
        std::cerr << "Error: " << name << " expects a number (got '" << text << "')\n";
        return false;
    }
}

void print_help(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " <filepaths_file> <primary_file> [options]\n\n"
        "Standalone benchmark front-end for hammock: same sketching algorithm as the\n"
        "Python CLI, no Python in the loop. Interval modes A/B/C only -- sequence mode\n"
        "(D) needs the digest library and lives on the Python side.\n\n"
        "Positional arguments (both required):\n"
        "  <filepaths_file>        Text file listing one plain-text BED path per line --\n"
        "                          NOT a BED file itself. This is the query side (the\n"
        "                          'query' and containment_AB columns).\n"
        "  <primary_file>          Same, for the reference side (the 'reference' and\n"
        "                          containment_BA columns). May be the same file, for\n"
        "                          an all-vs-all run.\n"
        "                          Neither .gz nor BigBed is decoded; such a file\n"
        "                          sketches to nothing and scores 0 against everything.\n\n"
        "Options:\n"
        "  --mode <A|B|C>          Comparison mode (default: B, matching the Python CLI's\n"
        "                          autodetect for BED input). A hashes whole interval\n"
        "                          records; B hashes every covered base position;\n"
        "                          C does both and is the only mode reading --subA/--expA.\n"
        "  --subA <float>          Mode C only: fraction of interval records kept as\n"
        "                          whole-interval elements, 0..1 (default: 1.0). A record\n"
        "                          that fails the gate still contributes its base points.\n"
        "  --subB <float>          Modes B and C: fraction of base positions sketched,\n"
        "                          0..1 (default: 1.0 = every base).\n"
        "  --subB-method <name>    How --subB picks positions (default: mixed-stride; no\n"
        "                          effect at --subB 1.0, where all three keep every base).\n"
        "                          mixed-stride   deterministic ~1/subB stride, stride and\n"
        "                                         phase keyed by chromosome name; fastest.\n"
        "                          hash-threshold keep when xxh32(point, --gate-seed) is\n"
        "                                         under the rate; orig-hammock parity.\n"
        "                          single-hash    one xxh64 is both gate and HLL input;\n"
        "                                         opt-in divergence, not orig-comparable.\n"
        "  --seed <int>            Seed of the xxh64 hash ingested by the HLL (default: 42).\n"
        "                          Under single-hash it also decides which points survive.\n"
        "  --gate-seed <int>       Seed of the xxh32 subsampling gates: the --subB point\n"
        "                          gate under hash-threshold, the mixed-stride\n"
        "                          chr->stride/offset hash, and the Mode C --subA interval\n"
        "                          gate (default: 31337, matches orig hammock). Under\n"
        "                          single-hash the point gate uses --seed instead, but the\n"
        "                          --subA gate still uses this seed.\n"
        "  --expA <float>          Mode C only: give each kept interval record 10**expA\n"
        "                          distinct elements instead of one (default: 0, i.e. no\n"
        "                          expansion), up-weighting intervals against points.\n"
        "  --precision, -p <int>   HyperLogLog precision: 2**p one-byte registers per\n"
        "                          sketch, p in 4..24 (default: 18).\n"
        "  --separator, -s <str>   Field separator used *inside* the strings that get\n"
        "                          hashed (default: tab). This is not the output\n"
        "                          delimiter; changing it changes every sketch and breaks\n"
        "                          comparability with the Python CLI.\n"
        "  --metrics               Accepted for compatibility; this is now the default.\n"
        "  --no-metrics            Emit only query/reference/jaccard_similarity, and\n"
        "                          tag the output filename '_j3'. Skips a union plus\n"
        "                          two cardinality estimates per pair -- use it for\n"
        "                          timing runs that must stay comparable with results\n"
        "                          published before 0.7.0.\n"
        "  --output, -o <prefix>   Output filename prefix (default: hammock); --outprefix\n"
        "                          is accepted as an alias. Suffixes are appended, see\n"
        "                          Output below.\n"
        "  --threads, -t <int>     OpenMP thread count for both per-file sketching and the\n"
        "                          pairwise loop (default: whatever OpenMP picks, i.e.\n"
        "                          $OMP_NUM_THREADS or all cores).\n"
        "  --verbose, -v           Progress and per-phase timings (sketching, pairwise,\n"
        "                          write, in microseconds) on stderr.\n"
        "  --help, -h              Show this help\n"
        "  --version               Print the version and exit\n\n"
        "Output:\n"
        "  A TAB-delimited file (despite the .csv extension) named\n"
        "    <prefix>_hll_p<precision>_jacc<MODE>[_expA<e>|_A<subA>][_B<subB>][_j3].csv\n"
        "  with columns query, reference, jaccard_similarity, jaccard_similarity_ie,\n"
        "  containment_AB, containment_BA, cosketch_geom, cosketch_arith, cosketch_max\n"
        "  (only the first three under --no-metrics). Inputs are always named by\n"
        "  basename; there is no --full-paths here. Note the metric columns match the\n"
        "  Python CLI bit-for-bit, but the row layout does not: the Python CSV is\n"
        "  comma-delimited and carries extra metadata columns.\n";
}

bool parse_args(int argc, char** argv, Args& out) {
    std::vector<std::string> positional;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            print_help(argv[0]);
            std::exit(0);
        } else if (a == "--version") {
            // stdout, unlike --help: a harness probes this to refuse a stale
            // binary, and parsing stderr would also catch warnings.
            std::cout << "hammock-cpp " << HAMMOCK_VERSION << "\n";
            std::exit(0);
        } else if (a == "--verbose" || a == "-v") {
            out.verbose = true;
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
            if (!parse_num("--subA", argv[++i], out.subA)) return false;
        } else if (a == "--subB" && i + 1 < argc) {
            if (!parse_num("--subB", argv[++i], out.subB)) return false;
        } else if (a == "--seed" && i + 1 < argc) {
            if (!parse_num("--seed", argv[++i], out.seed)) return false;
        } else if (a == "--gate-seed" && i + 1 < argc) {
            if (!parse_num("--gate-seed", argv[++i], out.gate_seed)) return false;
        } else if (a == "--expA" && i + 1 < argc) {
            if (!parse_num("--expA", argv[++i], out.expA)) return false;
        } else if ((a == "--precision" || a == "-p") && i + 1 < argc) {
            if (!parse_num("--precision", argv[++i], out.precision)) return false;
        } else if ((a == "--separator" || a == "-s") && i + 1 < argc) {
            out.separator = argv[++i];
        } else if ((a == "--output" || a == "-o" || a == "--outprefix") && i + 1 < argc) {
            out.outprefix = argv[++i];
        } else if ((a == "--threads" || a == "-t") && i + 1 < argc) {
            if (!parse_num("--threads", argv[++i], out.threads)) return false;
        } else if (a == "--metrics") {
            out.metrics = true;
        } else if (a == "--no-metrics") {
            out.metrics = false;
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
    return std::make_unique<HLLSketch>(a.precision);
}

void process_one(const std::string& path, AbstractSketch& sketch, const Args& a) {
    if (a.mode == "A") {
        process_bed_file_mode_a(path, sketch, a.seed, a.separator, a.verbose);
    } else if (a.mode == "B") {
        process_bed_file_mode_b(path, sketch, a.seed, a.separator, a.subB, a.subB_method,
                                a.gate_seed, a.verbose);
    } else {  // "C"
        process_bed_file_mode_c(path, sketch, a.seed, a.separator, a.subA, a.subB, a.expA,
                                a.subB_method, a.gate_seed, a.verbose);
    }
}

std::string outprefix_with_suffix(const Args& a) {
    std::string out = a.outprefix;
    out += "_hll_p" + std::to_string(a.precision) + "_jacc";
    out += a.mode;
    // Mirrors outprefix.py::get_new_prefix, including its elif: expA wins over
    // subA. Until 0.7.0 the subA branch was missing here entirely, so two Mode C
    // runs differing only in --subA wrote to the same path and the second
    // silently overwrote the first.
    char buf[32];
    if (a.expA > 0) {
        std::snprintf(buf, sizeof(buf), "%.2f", a.expA);
        out += "_expA";
        out += buf;
    } else if (a.subA != 1.0) {
        std::snprintf(buf, sizeof(buf), "%.2f", a.subA);
        out += "_A";
        out += buf;
    }
    if (a.subB != 1.0) {
        // Four decimals strictly below 0.01, matching outprefix.py: under %.2f
        // every subB in (0, 0.005) renders as "0.00" and collides. The boundary
        // is strict so subB == 0.01 keeps its historical "_B0.01".
        std::snprintf(buf, sizeof(buf), a.subB < 0.01 ? "%.4f" : "%.2f", a.subB);
        out += "_B";
        out += buf;
    }
    // Tag the *reduced* shape, so filename <-> column count stays one-to-one
    // while the default output keeps the path it has always had. Every reader
    // of this file parses by header name, so the tag is provenance for humans,
    // not a parsing hook.
    if (!a.metrics) out += "_j3";
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
    // process_one throws (e.g. "Cannot open file") and used to propagate out of
    // main -> std::terminate -> SIGABRT with a core dump. That is not just ugly
    // for a typo'd path: an open() can fail transiently under system-wide fd
    // exhaustion, so a long batch could abort mid-run on a file that is
    // perfectly readable a second later. Report and exit non-zero instead.
    try {
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
    } catch (const std::exception& e) {
        std::cerr << "Error while sketching: " << e.what() << "\n";
        return 3;
    } catch (...) {
        std::cerr << "Error while sketching: unknown exception\n";
        return 3;
    }

    auto t1 = std::chrono::steady_clock::now();

    const std::string out_path = outprefix_with_suffix(args);
    FILE* fp = std::fopen(out_path.c_str(), "w");
    if (!fp) {
        std::cerr << "Error: could not open '" << out_path << "' for writing\n";
        return 1;
    }
    if (args.metrics) {
        std::fprintf(fp, "query\treference\tjaccard_similarity\tjaccard_similarity_ie"
                         "\tcontainment_AB\tcontainment_BA"
                         "\tcosketch_geom\tcosketch_arith\tcosketch_max\n");
    } else {
        std::fprintf(fp, "query\treference\tjaccard_similarity\n");
    }

    const size_t n = queries.size();
    const size_t m = refs.size();
    const size_t stride = args.metrics ? 7 : 1;
    std::vector<double> matrix(n * m * stride);

    // Cardinalities are pair-invariant, so hoist them: calling
    // intersection_size() per pair would recompute both operands' Ertl
    // estimates every time. Only done under --metrics -- an unconditional
    // hoist would add n+m register passes to the default path.
    std::vector<double> qcard, rcard;
    if (args.metrics) {
        qcard.resize(n);
        rcard.resize(m);
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n; i++) qcard[i] = qsk[i]->cardinality();
#pragma omp parallel for schedule(static)
        for (size_t j = 0; j < m; j++) rcard[j] = rsk[j]->cardinality();
    }

    // The fused jaccard+union path below is a non-virtual HLLSketch method, so
    // downcast once here rather than per pair. HLLSketch is the only
    // AbstractSketch implementation, so this always succeeds; if a second
    // backend ever lands, a null entry falls back to the virtual route.
    std::vector<const HLLSketch*> qh(n, nullptr), rh(m, nullptr);
    if (args.metrics) {
        for (size_t i = 0; i < n; i++) qh[i] = dynamic_cast<const HLLSketch*>(qsk[i].get());
        for (size_t j = 0; j < m; j++) rh[j] = dynamic_cast<const HLLSketch*>(rsk[j].get());
    }

    // The guards inside the sketch methods throw. An exception escaping an
    // OpenMP structured block calls std::terminate, so catch and latch instead.
    // Fixed buffer, not std::string: a std::string assignment would allocate
    // inside the critical section -- throwing out of it and calling
    // std::terminate, i.e. exactly what the latch exists to prevent. snprintf
    // into stack storage cannot allocate.
    std::atomic<bool> failed{false};
    char first_error[256] = {0};
#pragma omp parallel for collapse(2) schedule(static)
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            if (failed.load(std::memory_order_relaxed)) continue;
            try {
                double* cell = &matrix[(i * m + j) * stride];
                if (!args.metrics) {
                    cell[0] = qsk[i]->jaccard_similarity(*rsk[j]);
                    continue;
                }
                double u;
                if (qh[i] && rh[j]) {
                    qh[i]->jaccard_and_union_cardinality(*rh[j], cell[0], u);
                } else {
                    cell[0] = qsk[i]->jaccard_similarity(*rsk[j]);
                    u = qsk[i]->union_with(*rsk[j])->cardinality();
                }
                const double inter = std::max(0.0, qcard[i] + rcard[j] - u);
                const double c_ab = (qcard[i] > 0.0) ? inter / qcard[i] : 0.0;
                const double c_ba = (rcard[j] > 0.0) ? inter / rcard[j] : 0.0;
                cell[1] = jaccard_ie_from_containments(c_ab, c_ba);
                cell[2] = c_ab;
                cell[3] = c_ba;
                // Cosketches use the *unclamped* containments, matching
                // runner._cosketch_from_containments.
                cell[4] = std::sqrt(std::max(c_ab * c_ba, 0.0));
                cell[5] = 0.5 * (c_ab + c_ba);
                cell[6] = std::max(c_ab, c_ba);
            } catch (const std::exception& e) {
#pragma omp critical
                {
                    if (!failed.exchange(true))
                        std::snprintf(first_error, sizeof(first_error), "%s", e.what());
                }
            } catch (...) {
                // Anything not derived from std::exception would otherwise
                // escape the parallel region and terminate.
#pragma omp critical
                {
                    if (!failed.exchange(true))
                        std::snprintf(first_error, sizeof(first_error),
                                      "unknown exception");
                }
            }
        }
    }
    if (failed.load()) {
        std::fclose(fp);
        // Remove the header-only file: downstream harnesses glob for the
        // output and would otherwise collect it as a valid zero-row result.
        std::remove(out_path.c_str());
        std::cerr << "Error computing pairwise metrics: " << first_error << "\n";
        return 1;
    }

    // Split the pairwise phase from the serial write. The write loop below is
    // single-threaded fprintf whose cost scales with the *column count*, so
    // folding it into one number would charge --metrics' six extra %.17g fields
    // to the estimator arithmetic. At low precision that formatting dominates.
    auto t_w0 = std::chrono::steady_clock::now();

    for (size_t i = 0; i < n; i++) {
        const std::string ql = basename_of(queries[i]);
        for (size_t j = 0; j < m; j++) {
            const double* cell = &matrix[(i * m + j) * stride];
            std::fprintf(fp, "%s\t%s\t%.17g",
                         ql.c_str(), basename_of(refs[j]).c_str(), cell[0]);
            for (size_t k = 1; k < stride; k++) std::fprintf(fp, "\t%.17g", cell[k]);
            std::fprintf(fp, "\n");
        }
    }
    std::fclose(fp);

    auto t2 = std::chrono::steady_clock::now();

    if (args.verbose) {
        // Microseconds, not milliseconds: at p=14 the whole pairwise phase is
        // ~1 ms per 1024 pairs, so the old integer-ms timer reported a flat 0
        // for every small collection and the harnesses had to floor it.
        // "Pairwise+write" keeps its historical meaning (pair loop + write) so
        // archived comparison_time values stay comparable; the two new lines
        // decompose it.
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        const auto sketch_us = duration_cast<microseconds>(t1 - t0).count();
        const auto pair_only_us = duration_cast<microseconds>(t_w0 - t1).count();
        const auto write_us = duration_cast<microseconds>(t2 - t_w0).count();
        const auto pair_us = duration_cast<microseconds>(t2 - t1).count();
        std::cerr << "Sketching: " << sketch_us << " us\n"
                  << "Pairwise: " << pair_only_us << " us\n"
                  << "Write: " << write_us << " us\n"
                  << "Pairwise+write: " << pair_us << " us\n"
                  << "Wrote " << out_path << "\n";
    }
    return 0;
}
