#ifndef HAMMOCK_HLL_SKETCH_HPP
#define HAMMOCK_HLL_SKETCH_HPP

#include "hammock/abstract_sketch.hpp"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Hammock-parity HyperLogLog. Implements *exactly* the algorithm in
// hammock/lib/hyperloglog.py (the pure-Python parent class), not the
// optimized hll/ library's variant — the two use different bit conventions
// (idx low-bits vs high-bits, ctz vs clz for rho) and would not produce
// byte-equal CSV output even on identical inputs.
//
// `final`: lets the compiler devirtualize add() at any HLLSketch& call site,
// which matters because the Mode B/C inner loop calls add() ~50M times per
// large file and the body is small enough to inline (see add() below).
class HLLSketch final : public AbstractSketch {
public:
    explicit HLLSketch(size_t precision, size_t hash_size = 64);
    HLLSketch(const HLLSketch&) = default;
    HLLSketch& operator=(const HLLSketch&) = default;

    // Defined inline so the stride loop can inline the body. Equivalent to
    // hyperloglog.py::_process_kmer + _rho:
    //   idx = hash_val low precision bits
    //   pos = trailing-zero count of (hash_val >> precision) + 1, capped at
    //         hash_size - precision + 1
    //   registers[idx] = max(registers[idx], pos)
    void add(uint64_t hash_val) final {
        const size_t idx = hash_val & (num_registers_ - 1);
        const uint64_t rest = hash_val >> precision_;
        const size_t max_rho = hash_size_ - precision_;
        size_t pos;
        if (rest == 0) {
            pos = max_rho + 1;
        } else {
            const size_t ctz = static_cast<size_t>(__builtin_ctzll(rest));
            pos = (ctz < max_rho) ? ctz + 1 : max_rho + 1;
        }
        if (registers_[idx] < pos) {
            registers_[idx] = static_cast<uint8_t>(pos);
        }
    }
    // In-place elementwise max — used to combine thread-local accumulators.
    void merge_max(const HLLSketch& other);
    double reg_eq_similarity(const AbstractSketch& other) const override;
    double cardinality() const override;

    // One pass over both register arrays, yielding the register-equality
    // Jaccard *and* the cardinality of the union — without materializing the
    // union sketch. This is what the pairwise metric loops want: the old route
    // (reg_eq_similarity + union_with()->cardinality(), or worse
    // intersection_size()) walks the registers three to five times per pair and
    // heap-allocates 16 MiB at p=24 for a union that is read once and thrown
    // away.
    //
    // Bit-identical to the separate calls by construction, not by measurement:
    // ertl_improved_estimate() consumes the registers *only* as the integer
    // histogram built at the top of it, and `counts[max(a[i], b[i])]++` is the
    // same multiset as `counts[union[i]]++`. Identical integers in, identical
    // doubles out of the same tau/sigma code. The one place that reads the
    // registers directly is the z < 1e-10 Flajolet fallback, which is
    // summation-order sensitive; it is reproduced here over max(a[i], b[i]) in
    // index order.
    void reg_eq_and_union_cardinality(const HLLSketch& other,
                                      double& reg_eq,
                                      double& union_cardinality) const;

    double intersection_size(const AbstractSketch& other) const override;
    std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const override;
    void clear() override;

    size_t precision() const { return precision_; }
    size_t hash_size_bits() const { return hash_size_; }
    size_t num_registers() const { return num_registers_; }
    const std::vector<uint8_t>& registers() const { return registers_; }
    std::vector<uint8_t>& registers() { return registers_; }

private:
    size_t precision_;
    size_t hash_size_;
    size_t num_registers_;
    std::vector<uint8_t> registers_;

    double get_alpha() const;
    double get_tau(double x) const;
    double get_sigma(double x) const;
    double ertl_improved_estimate() const;

    // Ertl 2017 Eq. 13 evaluated from a register-value histogram. `counts` must
    // be sized as ertl_improved_estimate() sizes it. Returns false without
    // touching `out` when z < 1e-10, i.e. when the caller must fall back to the
    // Flajolet estimate — that path needs the registers themselves in index
    // order, which a histogram cannot reconstruct, so its owner supplies it.
    bool ertl_from_counts(const std::vector<size_t>& counts, double& out) const;
};

#endif
