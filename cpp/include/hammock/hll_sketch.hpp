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
class HLLSketch : public AbstractSketch {
public:
    explicit HLLSketch(size_t precision, size_t hash_size = 64);
    HLLSketch(const HLLSketch&) = default;
    HLLSketch& operator=(const HLLSketch&) = default;

    void add(uint64_t hash_val) override;
    // In-place elementwise max — used to combine thread-local accumulators.
    void merge_max(const HLLSketch& other);
    double jaccard_similarity(const AbstractSketch& other) const override;
    double cardinality() const override;
    double intersection_size(const AbstractSketch& other) const override;
    std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const override;
    std::string get_sketch_type() const override;
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
};

#endif
