#include "hll.h"
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
namespace hll {

constexpr long double TWO_POW_32 = (1ull << 32) * 1.;

// Helper functions for Ertl's estimators
static double gen_tau(double x) {
    if(x == 0.0 || x == 1.0) return 0.0;
    
    double z = 1.0 - x;
    double tmp = 0.0;
    double y = 1.0;
    double zp = x;
    
    while(zp != z) {
        x = std::sqrt(x);
        zp = z;
        y *= 0.5;
        tmp = (1.0 - x);
        z -= tmp * tmp * y;
    }
    
    return z / 3.0;
}

static double gen_sigma(double x) {
    if(x == 1.0) return std::numeric_limits<double>::infinity();
    
    double z = x;
    for(double zp = 0.0, y = 1.0; z != zp;) {
        x *= x;
        zp = z;
        z += x * y;
        y += y;
        if(std::isnan(z)) {
            return zp;
        }
    }
    return z;
}

static double calculate_register_probability(double cardinality, int register_value) {
    // Probability that a register has value j given cardinality n
    // P(M[j] = j | n) = (1 - 2^(-j-1))^n - (1 - 2^(-j))^n
    
    if(register_value == 0) {
        // P(M[j] = 0 | n) = (1 - 1/2)^n = 2^(-n)
        return std::pow(0.5, cardinality);
    }
    
    const double p_j_minus_1 = 1.0 - std::pow(2.0, -(register_value - 1));
    const double p_j = 1.0 - std::pow(2.0, -register_value);
    
    return std::pow(p_j, cardinality) - std::pow(p_j_minus_1, cardinality);
}

void hll_t::sum() {
    sum_ = 0;
    for(unsigned i(0); i < m_; ++i) sum_ += 1. / (1ull << core_[i]);
    is_calculated_ = 1;
}

double hll_t::creport() const {
    if(!is_calculated_) throw std::runtime_error("Result must be calculated in order to report."
                                                 " Try the report() function.");
    const long double ret(alpha_ * m_ * m_ / sum_);
    // Small/large range corrections
    // See Flajolet, et al. HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm
    if(ret < small_range_correction_threshold()) {
        int t(0);
        for(const auto i: core_) t += i == 0;
        if(t) {
            LOG_DEBUG("Small value correction. Original estimate %lf. New estimate %lf.\n",
                      ret, m_ * std::log((double)m_ / t));
            return m_ * std::log((double)(m_) / t);
        }
    } else if(ret > LARGE_RANGE_CORRECTION_THRESHOLD) {
        const long double corr(-TWO_POW_32 * std::log(1. - ret / TWO_POW_32));
        if(!std::isnan(corr)) return corr;
        LOG_WARNING("Large range correction returned nan.\n");
    }
    return ret;
}

double hll_t::cest_err() const {
    if(!is_calculated_) throw std::runtime_error("Result must be calculated in order to report.");
    return relative_error_ * creport();
}

double hll_t::est_err() noexcept {
    if(!is_calculated_) sum();
    return cest_err();
}

double hll_t::report() noexcept {
    if(!is_calculated_) sum();
    return creport();
}

hll_t const &hll_t::operator+=(const hll_t &other) {
    if(other.np_ != np_) {
        char buf[256];
        sprintf(buf, "For operator +=: np_ (%zu) != other.np_ (%zu)\n", np_, other.np_);
        throw std::runtime_error(buf);
    }
    unsigned i;
#if HAS_AVX_512
    // Use unaligned loads/stores because AlignedAllocator only provides 16-byte alignment
    uint8_t *els_ptr = core_.data();
    const uint8_t *oels_ptr = other.core_.data();
    // Process full SIMD chunks (64 bytes per chunk)
    for(i = 0; i < m_ >> 6; ++i) {
        __m512i a = _mm512_loadu_si512(reinterpret_cast<const void*>(els_ptr + i * 64));
        __m512i b = _mm512_loadu_si512(reinterpret_cast<const void*>(oels_ptr + i * 64));
        __m512i result = _mm512_or_epi64(a, b);
        _mm512_storeu_si512(reinterpret_cast<void*>(els_ptr + i * 64), result);
    }
    // Process remaining bytes (start after last SIMD chunk)
    for(i = (m_ >> 6) << 6; i < m_; ++i) core_[i] |= other.core_[i];
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#elif __AVX2__
    // Use unaligned loads/stores because AlignedAllocator only provides 16-byte alignment
    uint8_t *els_ptr = core_.data();
    const uint8_t *oels_ptr = other.core_.data();
    // Process full SIMD chunks (32 bytes per chunk)
    for(i = 0; i < m_ >> 5; ++i) {
        __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(els_ptr + i * 32));
        __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(oels_ptr + i * 32));
        __m256i result = _mm256_or_si256(a, b);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(els_ptr + i * 32), result);
    }
    // Process remaining bytes (start after last SIMD chunk)
    for(i = (m_ >> 5) << 5; i < m_; ++i) core_[i] |= other.core_[i];
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#elif __SSE2__
    // Use unaligned loads/stores for safety
    uint8_t *els_ptr = core_.data();
    const uint8_t *oels_ptr = other.core_.data();
    // Process full SIMD chunks (16 bytes per chunk)
    for(i = 0; i < m_ >> 4; ++i) {
        __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i*>(els_ptr + i * 16));
        __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i*>(oels_ptr + i * 16));
        __m128i result = _mm_or_si128(a, b);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(els_ptr + i * 16), result);
    }
    // Process remaining bytes (start after last SIMD chunk)
    for(i = (m_ >> 4) << 4; i < m_; ++i) core_[i] |= other.core_[i];
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#else
    for(i = 0; i < m_; ++i) core_[i] |= other.core_[i];
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#endif
    return *this;
}

hll_t const &hll_t::operator&=(const hll_t &other) {
    if(other.np_ != np_) {
        char buf[256];
        sprintf(buf, "For operator &=: np_ (%zu) != other.np_ (%zu)\n", np_, other.np_);
        throw std::runtime_error(buf);
    }
    unsigned i;
#if HAS_AVX_512
    // Use unaligned loads/stores because AlignedAllocator only provides 16-byte alignment
    uint8_t *els_ptr = core_.data();
    const uint8_t *oels_ptr = other.core_.data();
    // Process full SIMD chunks (64 bytes per chunk)
    for(i = 0; i < m_ >> 6; ++i) {
        __m512i a = _mm512_loadu_si512(reinterpret_cast<const void*>(els_ptr + i * 64));
        __m512i b = _mm512_loadu_si512(reinterpret_cast<const void*>(oels_ptr + i * 64));
        __m512i result = _mm512_min_epu8(a, b);
        _mm512_storeu_si512(reinterpret_cast<void*>(els_ptr + i * 64), result);
    }
    // Process remaining bytes (start after last SIMD chunk)
    for(i = (m_ >> 6) << 6; i < m_; ++i) core_[i] = std::min(core_[i], other.core_[i]);
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#elif __AVX2__
    // Use unaligned loads/stores because AlignedAllocator only provides 16-byte alignment
    uint8_t *els_ptr = core_.data();
    const uint8_t *oels_ptr = other.core_.data();
    // Process full SIMD chunks (32 bytes per chunk)
    for(i = 0; i < m_ >> 5; ++i) {
        __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(els_ptr + i * 32));
        __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(oels_ptr + i * 32));
        __m256i result = _mm256_min_epu8(a, b);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(els_ptr + i * 32), result);
    }
    // Process remaining bytes (start after last SIMD chunk)
    for(i = (m_ >> 5) << 5; i < m_; ++i) core_[i] = std::min(core_[i], other.core_[i]);
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#elif __SSE2__
    // Use unaligned loads/stores for safety
    uint8_t *els_ptr = core_.data();
    const uint8_t *oels_ptr = other.core_.data();
    // Process full SIMD chunks (16 bytes per chunk)
    for(i = 0; i < m_ >> 4; ++i) {
        __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i*>(els_ptr + i * 16));
        __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i*>(oels_ptr + i * 16));
        __m128i result = _mm_min_epu8(a, b);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(els_ptr + i * 16), result);
    }
    // Process remaining bytes (start after last SIMD chunk)
    for(i = (m_ >> 4) << 4; i < m_; ++i) core_[i] = std::min(core_[i], other.core_[i]);
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#else
    for(i = 0; i < m_; ++i) core_[i] = std::min(core_[i], other.core_[i]);
    // Invalidate cached sum since we modified the registers
    is_calculated_ = 0;
#endif
    return *this;
}

// Returns the size of a symmetric set difference.
double operator^(hll_t &first, hll_t &other) {
    return 2*(hll_t(first) + other).report() - first.report() - other.report();
}

// Returns the set intersection
hll_t operator&(hll_t &first, hll_t &other) {
    hll_t tmp(first);
    return tmp &= other;
}

hll_t operator+(const hll_t &one, const hll_t &other) {
    if(other.get_np() != one.get_np())
        LOG_EXIT("np_ (%zu) != other.get_np() (%zu)\n", one.get_np(), other.get_np());
    hll_t ret(one);
    return ret += other;
}
// Returns the size of the set intersection
double intersection_size(const hll_t &first, const hll_t &other) {
    hll_t tmp(first);
    tmp &= other;
    return tmp.report_ertl_improved();
}

// Returns the size of the set intersection
double intersection_size(hll_t &first, hll_t &other) noexcept {
    hll_t tmp(first);
    tmp &= other;
    return tmp.report_ertl_improved();
}

// Clears, allows reuse with different np.
void hll_t::resize(std::size_t new_size) {
    new_size = roundup64(new_size);
    LOG_DEBUG("Resizing to %zu, with np = %zu\n", new_size, (std::size_t)std::log2(new_size));
    clear();
    core_.resize(new_size);
    np_ = (std::size_t)std::log2(new_size);
    m_ = new_size;
    alpha_ = make_alpha(m_);
    relative_error_ = 1.03896 / std::sqrt(m_);
}

void hll_t::clear() {
     std::fill(std::begin(core_), std::end(core_), 0u);
     sum_ = is_calculated_ = 0;
}

std::string hll_t::to_string() const {
    return is_calculated_ ? std::to_string(creport()) + ", +- " + std::to_string(cest_err())
                          : desc_string();
}

std::string hll_t::desc_string() const {
    char buf[1024];
    std::sprintf(buf, "Size: %zu. nb: %zu. error: %lf. Is calculated: %s. sum: %lf\n",
                 np_, m_, relative_error_, is_calculated_ ? "true": "false", sum_);
    return buf;
}

void hll_t::free() {
    core_.resize(0);
    core_.shrink_to_fit();
}

double hll_t::jaccard_hll_ratio(const hll_t& other) const {
    if (np_ != other.np_) {
        throw std::invalid_argument("Cannot compare HyperLogLog sketches with different precision");
    }

    hll_t intersection_sketch(*this);
    intersection_sketch &= other;
    double intersection_estimate = intersection_sketch.report_ertl_improved();

    hll_t union_sketch(*this);
    union_sketch += other;
    double union_estimate = union_sketch.report_ertl_improved();

    if (union_estimate <= 0.0) {
        return 0.0;
    }
    return intersection_estimate / union_estimate;
}

double hll_t::jaccard_register_match(const hll_t& other) const {
    if (np_ != other.np_) {
        throw std::invalid_argument("Cannot compare HyperLogLog sketches with different precision");
    }

    std::size_t matching_nonzero = 0;
    std::size_t total_active = 0;
    for (std::size_t i = 0; i < m_; ++i) {
        const uint8_t reg1 = core_[i];
        const uint8_t reg2 = other.core_[i];
        if (reg1 > 0 || reg2 > 0) {
            total_active++;
            if (reg1 == reg2) {
                matching_nonzero++;
            }
        }
    }
    return (total_active > 0)
        ? static_cast<double>(matching_nonzero) / static_cast<double>(total_active)
        : 0.0;
}

double hll_t::jaccard_harmonic_direct(const hll_t& other) const {
    if (np_ != other.np_) {
        throw std::invalid_argument("Cannot compare HyperLogLog sketches with different precision");
    }
    
    // Direct register-pair Jaccard estimator - no cardinality estimation!
    // Insight: For each register pair (A[i], B[i]), estimate local "contribution"
    // to intersection vs union based on the relationship between values.
    //
    // If A[i] = B[i] = k, both sets likely hashed similar rare elements → strong intersection signal
    // If A[i] << B[i], set A likely has many more elements in this bucket → weak intersection
    // 
    // Weight by inverse probability: elements with leading zeros k contribute ~2^k to cardinality
    
    double weighted_intersection = 0.0;
    double weighted_union = 0.0;
    
    for (std::size_t i = 0; i < m_; ++i) {
        const uint8_t a = core_[i];
        const uint8_t b = other.core_[i];
        
        if (a == 0 && b == 0) continue;  // No signal from empty registers
        
        const uint8_t min_val = (a < b) ? a : b;
        const uint8_t max_val = (a > b) ? a : b;
        
        // The "weight" of this bucket based on rarity (higher value = rarer = more weight)
        const double weight = std::exp2(static_cast<double>(max_val));
        
        //  Contribution to intersection: similarity between register values
        // If registers match → full contribution; if different → partial based on overlap
        const double similarity = (max_val == 0) ? 1.0 : static_cast<double>(min_val) / static_cast<double>(max_val);
        
        weighted_intersection += similarity * weight;
        weighted_union += weight;
    }
    
    return (weighted_union > 0.0) ? (weighted_intersection / weighted_union) : 0.0;
}

double hll_t::jaccard_super_registers(const hll_t& other, std::size_t super_size) const {
    if (np_ != other.np_) {
        throw std::invalid_argument("Cannot compare HyperLogLog sketches with different precision");
    }
    if (super_size == 0) {
        throw std::invalid_argument("super_size must be >= 1");
    }
    // Lightweight permutations: for m_ = 2^p, (i * odd + add) & (m_-1) is a bijection.
    static const std::uint32_t multipliers[] = {1u, 3u, 5u, 7u};
    static const std::size_t num_perms = sizeof(multipliers) / sizeof(multipliers[0]);
    const std::size_t mask = m_ - 1u;

    // Limit number of offsets to m_ to avoid degenerate cases when m_ < super_size
    const std::size_t num_offsets = std::min(super_size, m_);
    double avg_fraction = 0.0;
    std::size_t used_terms = 0; // permutations * offsets that had active groups

    for (std::size_t perm_idx = 0; perm_idx < num_perms; ++perm_idx) {
        const std::uint32_t mult = multipliers[perm_idx];
        const std::size_t add = (perm_idx * 2654435761u) & mask; // odd-ish rotation
        for (std::size_t offset = 0; offset < num_offsets; ++offset) {
            std::size_t active_groups = 0;
            std::size_t matching_groups = 0;
            // Walk groups starting at this offset on the permuted index space
            for (std::size_t i = offset; i < m_; i += super_size) {
                const std::size_t group_end = std::min(i + super_size, m_);
                uint8_t group_max_a = 0;
                uint8_t group_max_b = 0;
                for (std::size_t j = i; j < group_end; ++j) {
                    const std::size_t mapped = ((j * mult) + add) & mask;
                    const uint8_t a = core_[mapped];
                    const uint8_t b = other.core_[mapped];
                    if (a > group_max_a) group_max_a = a;
                    if (b > group_max_b) group_max_b = b;
                }
                const bool active = (group_max_a > 0) || (group_max_b > 0);
                if (active) {
                    ++active_groups;
                    if (group_max_a == group_max_b) ++matching_groups;
                }
            }
            if (active_groups > 0) {
                avg_fraction += static_cast<double>(matching_groups) / static_cast<double>(active_groups);
                ++used_terms;
            }
        }
    }
    if (used_terms == 0) return 0.0;
    return avg_fraction / static_cast<double>(used_terms);
}

double hll_t::jaccard_similarity_registers(const hll_t& other) const {
    double jaccard_hll = jaccard_hll_ratio(other);
    // Replace raw register matching with Dashing-style super-register matching (default size 4)
    double register_jaccard = jaccard_super_registers(other, 4);

    // Cardinality ratio guides how we weight the two estimators:
    //  - near-equal set sizes → rely on the HLL ratio
    //  - moderately imbalanced sizes → blend the two
    //  - highly imbalanced sizes → lean on the register match heuristic
    hll_t sketch_a(*this);
    hll_t sketch_b(other);
    double card_a = sketch_a.report_ertl_improved();
    double card_b = sketch_b.report_ertl_improved();

    double min_card = std::min(card_a, card_b);
    double max_card = std::max(card_a, card_b);
    double card_ratio = (max_card > 0.0) ? (min_card / max_card) : 0.0;

    double combined;
    if (card_ratio > 0.8) {
        combined = jaccard_hll;
    } else if (card_ratio > 0.5) {
        combined = 0.6 * jaccard_hll + 0.4 * register_jaccard;
    } else {
        combined = 0.3 * jaccard_hll + 0.7 * register_jaccard;
    }

    if (combined < 0.0) {
        combined = 0.0;
    } else if (combined > 1.0) {
        combined = 1.0;
    }

    return combined;
}

// MLE cardinality estimation methods
double hll_t::report_mle() noexcept {
    if(!is_calculated_) sum();
    return creport();  // For now, same as standard - can be enhanced
}

double hll_t::report_ertl_improved() noexcept {
    if(!is_calculated_) sum();
    
    // Ertl's improved estimator
    // Based on: "New cardinality estimation algorithms for HyperLogLog sketches" by Otmar Ertl
    
    const double m = static_cast<double>(m_);
    const double alpha = alpha_;
    
    // Count registers by value
    std::vector<int> counts(64 - np_ + 1, 0);
    for(const auto& reg : core_) {
        if(reg < counts.size()) {
            counts[reg]++;
        }
    }
    
    // Calculate tau and sigma corrections
    const double divinv = 1.0 / (2.0 * std::log(2.0));
    
    // Tau correction for maxed registers
    const double maxed_registers = counts[64 - np_];
    const double tau = gen_tau((m - maxed_registers) / m);
    
    // Calculate harmonic sum with tau correction
    double z = m * tau;
    for(int i = 64 - np_; i > 0; --i) {
        z += counts[i];
        z *= 0.5;
    }
    
    // Sigma correction for zero registers
    const double zeros = counts[0];
    const double sigma = gen_sigma(zeros / m);
    z += m * sigma;
    
    return m * divinv * m / z;
}

double hll_t::report_ertl_mle() noexcept {
    if(!is_calculated_) sum();
    
    // Ertl's MLE estimator
    // Based on maximum likelihood estimation of HyperLogLog sketches
    
    const double m = static_cast<double>(m_);
    const double alpha = alpha_;
    
    // Count registers by value
    std::vector<int> counts(64 - np_ + 1, 0);
    for(const auto& reg : core_) {
        if(reg < counts.size()) {
            counts[reg]++;
        }
    }
    
    // Initial estimate using standard HLL
    double estimate = creport();
    
    // Iterative MLE refinement
    const double relative_error = 1e-4;
    const int max_iterations = 10;
    
    for(int iter = 0; iter < max_iterations; ++iter) {
        double prev_estimate = estimate;
        
        // Calculate likelihood gradient
        double gradient = 0.0;
        double hessian = 0.0;
        
        for(int j = 0; j < static_cast<int>(counts.size()); ++j) {
            if(counts[j] == 0) continue;
            
            const double prob = calculate_register_probability(estimate, j);
            if(prob > 0) {
                gradient += counts[j] / prob;
                hessian += counts[j] / (prob * prob);
            }
        }
        
        // Newton-Raphson update
        if(hessian > 0) {
            estimate = estimate - gradient / hessian;
            estimate = std::max(estimate, 1.0);  // Ensure positive
        }
        
        // Check convergence
        if(std::abs(estimate - prev_estimate) / estimate < relative_error) {
            break;
        }
    }
    
    return estimate;
}

} // namespace hll
