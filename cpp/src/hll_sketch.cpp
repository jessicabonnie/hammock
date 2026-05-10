#include "hammock/hll_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

HLLSketch::HLLSketch(size_t precision, size_t hash_size)
    : precision_(precision),
      hash_size_(hash_size),
      num_registers_(static_cast<size_t>(1) << precision),
      registers_(num_registers_, 0) {
    if (precision < 4) {
        throw std::invalid_argument("HLL precision must be >= 4");
    }
    if (hash_size != 32 && hash_size != 64) {
        throw std::invalid_argument("hash_size must be 32 or 64");
    }
}

// HLLSketch::add is defined inline in the header.

void HLLSketch::merge_max(const HLLSketch& other) {
    if (precision_ != other.precision_ || hash_size_ != other.hash_size_) {
        throw std::runtime_error("HLLs must match for merge_max");
    }
    for (size_t i = 0; i < num_registers_; ++i) {
        if (other.registers_[i] > registers_[i]) {
            registers_[i] = other.registers_[i];
        }
    }
}

double HLLSketch::jaccard_similarity(const AbstractSketch& other) const {
    const HLLSketch* o = dynamic_cast<const HLLSketch*>(&other);
    if (!o) {
        throw std::runtime_error("Cannot compute Jaccard between different sketch types");
    }
    if (precision_ != o->precision_) {
        throw std::runtime_error("HLLs must have same precision for Jaccard");
    }
    // Python parity: estimate_jaccard_registers
    //   active = (R_a != 0) | (R_b != 0)
    //   J = |{i in active: R_a[i] == R_b[i]}| / |active|
    size_t matching = 0, active = 0;
    for (size_t i = 0; i < num_registers_; ++i) {
        const uint8_t a = registers_[i];
        const uint8_t b = o->registers_[i];
        if (a != 0 || b != 0) {
            ++active;
            if (a == b) ++matching;
        }
    }
    return (active == 0) ? 0.0 : static_cast<double>(matching) / static_cast<double>(active);
}

double HLLSketch::get_alpha() const {
    const double m = static_cast<double>(num_registers_);
    if (hash_size_ == 32) {
        if (num_registers_ == 16384) return 0.673;
        if (num_registers_ == 32768) return 0.697;
        if (num_registers_ == 65536) return 0.709;
        return 0.7213 / (1.0 + 1.079 / m);
    }
    // hash_size_ == 64
    if (num_registers_ == 16384) return 0.351;
    if (num_registers_ == 32768) return 0.532;
    if (num_registers_ == 65536) return 0.625;
    return 0.7213 / (1.0 + 1.079 / m) * 0.85;  // Python's 64-bit scaling
}

double HLLSketch::get_tau(double x) const {
    if (x == 0.0 || x == 1.0) return 0.0;
    double z = 1.0 - x;
    double y = 1.0;
    double zprev = x;
    while (zprev != z) {
        x = std::sqrt(x);
        zprev = z;
        y *= 0.5;
        const double tmp = 1.0 - x;
        z -= tmp * tmp * y;
    }
    return z / 3.0;
}

double HLLSketch::get_sigma(double x) const {
    if (x == 1.0) return std::numeric_limits<double>::infinity();
    double z = x;
    double zprev = 0.0;
    double y = 1.0;
    while (z != zprev) {
        x = x * x;
        zprev = z;
        z += x * y;
        y += y;
    }
    return z;
}

double HLLSketch::ertl_improved_estimate() const {
    // Ertl 2017 "New cardinality estimation algorithms for HyperLogLog sketches", Eq. 13:
    //   N̂ = α_∞ · m² / Z
    //   Z = m · τ(1 − C[q+1]/m) · 2^(−q)
    //       + Σ_{k=1}^q C[k] · 2^(−k)
    //       + m · σ(C[0]/m)
    // where q = hash_size − precision, and α_∞ = 1 / (2·ln 2).
    // The geometric loop below produces the m·τ·2^(−q) term and the C[k]·2^(−k)
    // sum in one pass. τ and σ are the iterative series from Ertl Eq. 11–12.
    const double alpha_inf = 1.0 / (2.0 * std::log(2.0));
    const double m = static_cast<double>(num_registers_);
    const size_t max_register_value = hash_size_ - precision_;

    // register_counts[k] = number of registers with value k. Size must
    // accommodate at least max_register_value + 1.
    std::vector<size_t> counts(std::max<size_t>(max_register_value + 2, 64), 0);
    for (size_t v : registers_) counts[v]++;

    const double n_maxed = static_cast<double>(counts[max_register_value + 1]);
    const double non_maxreg_frac = (m - n_maxed) / m;
    const double zero_reg_frac = static_cast<double>(counts[0]) / m;

    const double tau = get_tau(non_maxreg_frac);
    double z = m * tau;
    for (size_t i = max_register_value; i >= 1; --i) {
        if (i < counts.size()) z += counts[i];
        z *= 0.5;
    }
    const double sigma = get_sigma(zero_reg_frac);
    z += m * sigma;
    if (z < 1e-10) {
        // Fall back to the original Flajolet estimate. Python does this too.
        const double alpha = get_alpha();
        double sum = 0.0;
        for (size_t v : registers_) sum += std::pow(2.0, -static_cast<int>(v));
        return alpha * m * m / sum;
    }
    return alpha_inf * m * m / z;
}

double HLLSketch::cardinality() const {
    bool all_zero = true;
    for (size_t v : registers_) if (v != 0) { all_zero = false; break; }
    if (all_zero) return 0.0;
    return ertl_improved_estimate();
}

double HLLSketch::intersection_size(const AbstractSketch& other) const {
    // Python's `estimate_intersection` uses inclusion-exclusion:
    //   |A| + |B| - |A ∪ B|, with cardinality via Ertl-improved on each.
    const HLLSketch* o = dynamic_cast<const HLLSketch*>(&other);
    if (!o) {
        throw std::runtime_error("Cannot compute intersection between different sketch types");
    }
    auto u = union_with(other);
    HLLSketch* hu = dynamic_cast<HLLSketch*>(u.get());
    if (!hu) throw std::runtime_error("union_with returned wrong type");
    const double inter = cardinality() + o->cardinality() - hu->cardinality();
    return std::max(0.0, inter);
}

std::unique_ptr<AbstractSketch> HLLSketch::union_with(const AbstractSketch& other) const {
    const HLLSketch* o = dynamic_cast<const HLLSketch*>(&other);
    if (!o) {
        throw std::runtime_error("Cannot compute union between different sketch types");
    }
    auto result = std::make_unique<HLLSketch>(precision_, hash_size_);
    for (size_t i = 0; i < num_registers_; ++i) {
        result->registers_[i] = std::max(registers_[i], o->registers_[i]);
    }
    return result;
}

std::string HLLSketch::get_sketch_type() const {
    return "hyperloglog";
}

void HLLSketch::clear() {
    std::fill(registers_.begin(), registers_.end(), static_cast<uint8_t>(0));
}
