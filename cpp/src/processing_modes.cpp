#include "hammock/processing_modes.hpp"
#include "hammock/bagminhash_sketch.hpp"
#include "hammock/bed_parser.hpp"
#include "hammock/hll_sketch.hpp"
#include "hammock/stride.hpp"
#include "hammock/xxhash.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

std::string create_interval_string(const std::string& chr, int64_t start, int64_t end,
                                   const std::string& separator) {
    return chr + separator + std::to_string(start) + separator + std::to_string(end) + separator + "A";
}

std::string create_point_string(const std::string& chr, int64_t pos,
                                const std::string& separator) {
    // Python parity: chr + sep + pos, no trailing "B".
    return chr + separator + std::to_string(pos);
}

uint32_t hash_for_subsample(const std::string& s, uint32_t gate_seed) {
    // Python parity uses xxh32 seed=31337; users can override via --gate-seed.
    return xxhash::hash32(s, gate_seed);
}

double calculate_jaccard(const AbstractSketch& sketch1, const AbstractSketch& sketch2) {
    return sketch1.jaccard_similarity(sketch2);
}

size_t process_bed_file_mode_a(const std::string& filepath, AbstractSketch& sketch,
                               uint64_t hll_seed,
                               const std::string& separator, int peak_height_column,
                               bool verbose) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    BagMinHashSketch* bmh = (peak_height_column > 0)
                                ? dynamic_cast<BagMinHashSketch*>(&sketch)
                                : nullptr;

    size_t interval_count = 0;
    std::string line;
    std::string chr;
    int64_t start, end, count;

    while (std::getline(file, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

        if (parse_bed_line(line, chr, start, end, count, peak_height_column)) {
            const std::string interval_str = create_interval_string(chr, start, end, separator);
            const uint64_t hash_val = xxhash::hash64(interval_str, hll_seed);

            if (bmh) {
                bmh->add_with_normalized_count(hash_val, count);
            } else if (peak_height_column > 0) {
                for (int64_t i = 0; i < count; i++) {
                    sketch.add(hash_val);
                }
            } else {
                sketch.add(hash_val);
            }
            interval_count++;

            if (verbose && interval_count % 10000 == 0) {
                std::cerr << "Processed " << interval_count << " intervals...\r" << std::flush;
            }
        }
    }

    if (verbose && interval_count > 0) {
        std::cerr << "Processed " << interval_count << " intervals total.       \n";
    }

    return interval_count;
}

namespace {
struct Interval {
    std::string chr;
    int64_t start;
    int64_t end;
    int64_t count;
};

std::vector<Interval> read_intervals(const std::string& filepath, int peak_height_column) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

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
    return intervals;
}
}  // namespace

namespace {
inline const char* method_label(SubBMethod m) {
    switch (m) {
        case SubBMethod::HashThreshold: return "hash-threshold";
        case SubBMethod::MixedStride:   return "mixed-stride";
        case SubBMethod::SingleHash: return "single-hash";
    }
    return "?";
}
}  // namespace

size_t process_bed_file_mode_b(const std::string& filepath, AbstractSketch& sketch,
                               uint64_t hll_seed,
                               const std::string& separator,
                               double subB,
                               SubBMethod method,
                               uint32_t gate_seed,
                               int peak_height_column,
                               bool verbose) {
    const std::vector<Interval> intervals = read_intervals(filepath, peak_height_column);

    size_t total_points = 0;
    size_t sampled_points = 0;

    // Parallelize via thread-local HLL accumulators that merge at the end.
    // The previous design parallelized straight into the shared sketch with
    // a non-atomic check-then-store on uint8_t registers, which was a real
    // data race that could lose max-updates. Thread-locals also avoid false
    // sharing on the registers vector. Non-HLL sketches (e.g. BagMinHash on
    // Mode A) don't go through this path; if one is ever passed here we fall
    // back to serial.
    HLLSketch* main_hll = dynamic_cast<HLLSketch*>(&sketch);
    if (main_hll) {
        const size_t precision = main_hll->precision();
        const size_t hash_size = main_hll->hash_size_bits();

#pragma omp parallel reduction(+:total_points,sampled_points)
        {
            HLLSketch local(precision, hash_size);
#pragma omp for schedule(static) nowait
            for (size_t i = 0; i < intervals.size(); i++) {
                const auto& iv = intervals[i];
                total_points += (iv.end - iv.start);
                sampled_points += add_interval_points_to_sketch(iv.chr, iv.start, iv.end,
                                                               local, separator, subB,
                                                               method, hll_seed, gate_seed);
            }
#pragma omp critical
            { main_hll->merge_max(local); }
        }
    } else {
        for (size_t i = 0; i < intervals.size(); i++) {
            const auto& iv = intervals[i];
            total_points += (iv.end - iv.start);
            sampled_points += add_interval_points_to_sketch(iv.chr, iv.start, iv.end,
                                                           sketch, separator, subB,
                                                           method, hll_seed, gate_seed);
        }
    }

    if (verbose && !intervals.empty()) {
        if (subB < 1.0) {
            std::cerr << "Processed " << intervals.size() << " intervals, "
                      << sampled_points << "/" << total_points
                      << " points (subB=" << subB << ", method=" << method_label(method) << ").       \n";
        } else {
            std::cerr << "Processed " << intervals.size() << " intervals, "
                      << total_points << " points total.       \n";
        }
    }

    return sampled_points;
}

size_t process_bed_file_mode_c(const std::string& filepath, AbstractSketch& sketch,
                               uint64_t hll_seed,
                               const std::string& separator,
                               double subA,
                               double subB,
                               double expA,
                               SubBMethod method,
                               uint32_t gate_seed,
                               int peak_height_column,
                               bool verbose) {
    const size_t mult = (expA > 0)
                            ? std::max<size_t>(1, static_cast<size_t>(std::pow(10.0, expA)))
                            : 1;
    const bool do_subA = (subA < 1.0);
    const uint32_t subA_threshold = static_cast<uint32_t>(subA * 4294967295.0);

    const std::vector<Interval> intervals = read_intervals(filepath, peak_height_column);

    size_t total_interval_elements = 0;
    size_t kept_intervals = 0;
    size_t total_points = 0;
    size_t sampled_points = 0;

    // See process_bed_file_mode_b for the rationale: thread-local HLL
    // accumulators avoid the non-atomic max-update race on the shared
    // registers array.
    HLLSketch* main_hll = dynamic_cast<HLLSketch*>(&sketch);

    // Generic lambda so the HLLSketch instantiation picks the HLLSketch&
    // overload of add_interval_points_to_sketch (and resolves s.add() to the
    // inline final body).
    auto process_one = [&](const Interval& iv, auto& s,
                           size_t& te, size_t& kept,
                           size_t& tp, size_t& sp) {
        const std::string base = create_interval_string(iv.chr, iv.start, iv.end, separator);
        if (do_subA && hash_for_subsample(base, gate_seed) > subA_threshold) {
            tp += (iv.end - iv.start);
            sp += add_interval_points_to_sketch(iv.chr, iv.start, iv.end,
                                                s, separator, subB,
                                                method, hll_seed, gate_seed);
            return;
        }
        kept++;
        if (mult > 1) {
            // Python (intervals.py:344-348): for i in range(mult): add(interval + str(i+1))
            for (size_t k = 1; k <= mult; k++) {
                const std::string str = base + std::to_string(k);
                s.add(xxhash::hash64(str, hll_seed));
                te++;
            }
        } else {
            s.add(xxhash::hash64(base, hll_seed));
            te++;
        }
        tp += (iv.end - iv.start);
        sp += add_interval_points_to_sketch(iv.chr, iv.start, iv.end,
                                            s, separator, subB,
                                            method, hll_seed, gate_seed);
    };

    if (main_hll) {
        const size_t precision = main_hll->precision();
        const size_t hash_size = main_hll->hash_size_bits();

#pragma omp parallel reduction(+:total_interval_elements,kept_intervals,total_points,sampled_points)
        {
            HLLSketch local(precision, hash_size);
            size_t te = 0, kept = 0, tp = 0, sp = 0;
#pragma omp for schedule(static) nowait
            for (size_t i = 0; i < intervals.size(); i++) {
                process_one(intervals[i], local, te, kept, tp, sp);
            }
            total_interval_elements += te;
            kept_intervals += kept;
            total_points += tp;
            sampled_points += sp;
#pragma omp critical
            { main_hll->merge_max(local); }
        }
    } else {
        for (size_t i = 0; i < intervals.size(); i++) {
            process_one(intervals[i], sketch,
                        total_interval_elements, kept_intervals,
                        total_points, sampled_points);
        }
    }

    if (verbose && !intervals.empty()) {
        std::cerr << "Processed " << intervals.size() << " intervals";
        if (do_subA) {
            std::cerr << " (kept " << kept_intervals << " after subA=" << subA << ")";
        }
        if (expA > 0) {
            std::cerr << " (" << total_interval_elements << " expanded elements, expA=" << expA << ")";
        }
        if (subB < 1.0) {
            std::cerr << ", " << sampled_points << "/" << total_points
                      << " points (subB=" << subB << ", method=" << method_label(method) << ")";
        } else {
            std::cerr << ", " << total_points << " points";
        }
        std::cerr << ".       \n";
    }

    return total_interval_elements + sampled_points;
}
