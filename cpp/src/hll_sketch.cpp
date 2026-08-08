#include "hammock/hll_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

// Validate before any member initializer runs. Member initializers execute in
// declaration order and precede the constructor body, so a check written in
// the body is unreachable for exactly the inputs that break: `1 << precision`
// is UB at precision >= 64, and `registers_(1 << 32, 0)` would try to allocate
// 4 GiB before the body could reject it. `precision_` is declared first
// (hll_sketch.hpp:62), so routing the check through it runs first of all.
//
// The upper bound also keeps `hash_size_ - precision_` (max rho, computed as
// size_t at hll_sketch.hpp:34 and .cpp:110) from underflowing: at precision 65
// it wrapped to SIZE_MAX and turned the estimator loop into ~1.8e19 iterations.
size_t validated_precision(size_t precision, size_t hash_size) {
    if (hash_size != 32 && hash_size != 64) {
        throw std::invalid_argument("hash_size must be 32 or 64");
    }
    if (precision < 4 || precision > 24) {
        throw std::invalid_argument("HLL precision must be in 4..24");
    }
    return precision;
}

}  // namespace

HLLSketch::HLLSketch(size_t precision, size_t hash_size)
    : precision_(validated_precision(precision, hash_size)),
      hash_size_(hash_size),
      num_registers_(static_cast<size_t>(1) << precision_),
      registers_(num_registers_, 0) {}

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
    // hash_size must match as well as precision: equal precision alone gives
    // equal register *counts* (so no overread) but the rho values are on
    // different scales, which would silently return a meaningless J. merge_max
    // and union_with both guard this; this one used to check precision only.
    // Unreachable from Python -- the pybind ctor exposes only `precision`, so
    // every sketch there is hash_size 64 -- so adding it moves no value.
    if (precision_ != o->precision_ || hash_size_ != o->hash_size_) {
        throw std::runtime_error("HLLs must have same precision and hash size for Jaccard");
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
    // The geometric loop in ertl_from_counts produces the m·τ·2^(−q) term and
    // the C[k]·2^(−k) sum in one pass. τ and σ are the iterative series from
    // Ertl Eq. 11–12.
    const double m = static_cast<double>(num_registers_);
    const size_t max_register_value = hash_size_ - precision_;

    // register_counts[k] = number of registers with value k. Size must
    // accommodate at least max_register_value + 1.
    std::vector<size_t> counts(std::max<size_t>(max_register_value + 2, 64), 0);
    for (size_t v : registers_) counts[v]++;

    double out;
    if (ertl_from_counts(counts, out)) return out;

    // Fall back to the original Flajolet estimate. Python does this too.
    // Order-sensitive: the sum runs over the registers in index order.
    const double alpha = get_alpha();
    double sum = 0.0;
    for (size_t v : registers_) sum += std::pow(2.0, -static_cast<int>(v));
    return alpha * m * m / sum;
}

bool HLLSketch::ertl_from_counts(const std::vector<size_t>& counts, double& out) const {
    const double alpha_inf = 1.0 / (2.0 * std::log(2.0));
    const double m = static_cast<double>(num_registers_);
    const size_t max_register_value = hash_size_ - precision_;

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
    if (z < 1e-10) return false;
    out = alpha_inf * m * m / z;
    return true;
}

void HLLSketch::jaccard_and_union_cardinality(const HLLSketch& other,
                                              double& jaccard,
                                              double& union_cardinality) const {
    // Same guards as jaccard_similarity() and union_with(): equal precision
    // *and* hash size. Unequal precision would overread the operand's
    // registers; equal precision with unequal hash size gives rho values on
    // different scales.
    if (precision_ != other.precision_ || hash_size_ != other.hash_size_) {
        throw std::runtime_error(
            "HLLs must have same precision and hash size for fused jaccard/union");
    }

    const size_t max_register_value = hash_size_ - precision_;
    std::vector<size_t> counts(std::max<size_t>(max_register_value + 2, 64), 0);

    size_t matching = 0, active = 0;
    bool union_all_zero = true;
    for (size_t i = 0; i < num_registers_; ++i) {
        const uint8_t a = registers_[i];
        const uint8_t b = other.registers_[i];
        if (a != 0 || b != 0) {
            ++active;
            if (a == b) ++matching;
        }
        const uint8_t u = (a > b) ? a : b;
        counts[u]++;
        if (u != 0) union_all_zero = false;
    }

    jaccard = (active == 0) ? 0.0
                            : static_cast<double>(matching) / static_cast<double>(active);

    // Mirrors cardinality()'s all-zero short-circuit on the union sketch.
    if (union_all_zero) {
        union_cardinality = 0.0;
        return;
    }
    double out;
    if (ertl_from_counts(counts, out)) {
        union_cardinality = out;
        return;
    }
    // Flajolet fallback, over the union registers in index order.
    const double m = static_cast<double>(num_registers_);
    const double alpha = get_alpha();
    double sum = 0.0;
    for (size_t i = 0; i < num_registers_; ++i) {
        const uint8_t u = (registers_[i] > other.registers_[i]) ? registers_[i]
                                                                : other.registers_[i];
        sum += std::pow(2.0, -static_cast<int>(u));
    }
    union_cardinality = alpha * m * m / sum;
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
    if (precision_ != o->precision_) {
        throw std::runtime_error("HLLs must have same precision for intersection");
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
    // Without this the loop below reads o->registers_[i] for i < this->
    // num_registers_ — a heap overread (and, via the write, corruption) when
    // the operand has the smaller precision. jaccard_similarity and merge_max
    // both guard; this one did not.
    if (precision_ != o->precision_ || hash_size_ != o->hash_size_) {
        throw std::runtime_error("HLLs must have same precision for union");
    }
    auto result = std::make_unique<HLLSketch>(precision_, hash_size_);
    for (size_t i = 0; i < num_registers_; ++i) {
        result->registers_[i] = std::max(registers_[i], o->registers_[i]);
    }
    return result;
}


void HLLSketch::clear() {
    std::fill(registers_.begin(), registers_.end(), static_cast<uint8_t>(0));
}
