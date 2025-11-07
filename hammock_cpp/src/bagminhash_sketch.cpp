#include "../include/bagminhash_sketch.h"
#include "../include/xxhash.h"
#include <stdexcept>
#include <algorithm>
#include <climits>
#include <cmath>

/**
 * BagMinHashSketch Constructor
 * 
 * BagMinHash is a sketching technique that extends MinHash to handle weighted sets (bags)
 * where elements can appear multiple times. Instead of just tracking the minimum hash
 * for each hash function, it tracks both the minimum hash value and the count of how
 * many times that minimum value appears.
 * 
 * @param num_hashes Number of hash functions to use (determines sketch size)
 * @param seed Random seed for hash function generation
 */
BagMinHashSketch::BagMinHashSketch(size_t num_hashes, uint64_t seed) 
    : num_hashes_(num_hashes), min_hashes_(num_hashes, {UINT64_MAX, 0}), seed_(seed) {}

/**
 * Hash a value with a specific seed using XXHash64
 * 
 * Each hash function in BagMinHash uses a different seed to ensure independence.
 * This creates multiple independent hash functions from a single hash algorithm.
 * 
 * @param value The input value to hash
 * @param seed The seed for this particular hash function
 * @return 64-bit hash value
 */
uint64_t BagMinHashSketch::hash_with_seed(uint64_t value, uint64_t seed) const {
    // Convert uint64_t to string for xxhash (xxhash expects string input)
    std::string value_str = std::to_string(value);
    return xxhash::hash64(value_str, seed);
}

/**
 * Add a single hash value to the sketch (treats as count of 1)
 * 
 * This method is used when we don't have explicit count information.
 * Each call to add() is treated as adding one occurrence of the hash value.
 * 
 * @param hash_val The hash value to add to the sketch
 */
void BagMinHashSketch::add(uint64_t hash_val) {
    // For BagMinHash, we need to track counts, but since we're getting individual hash values,
    // we'll treat each add as a count of 1. In a real implementation, you'd want to batch
    // adds by the same hash value to properly handle counts.
    
    // Apply each hash function to the input value
    for (size_t i = 0; i < num_hashes_; ++i) {
        uint64_t hash_seed = seed_ + i;  // Each hash function uses a different seed
        uint64_t hashed_value = hash_with_seed(hash_val, hash_seed);
        
        if (hashed_value < min_hashes_[i].first) {
            // Found a new minimum hash value - update both value and count
            min_hashes_[i].first = hashed_value;
            min_hashes_[i].second = 1;  // Reset count when new minimum found
        } else if (hashed_value == min_hashes_[i].first) {
            // Same minimum hash value - increment the count
            min_hashes_[i].second++;
        }
        // If hashed_value > min_hashes_[i].first, we ignore it (not a new minimum)
    }
}


/**
 * Compute Jaccard similarity between two BagMinHash sketches
 * 
 * For BagMinHash, Jaccard similarity is computed as:
 * J(A,B) = |A ∩ B| / |A ∪ B|
 * 
 * Where:
 * - |A ∩ B| = sum of min(counts) for hash functions where both sketches have the same minimum hash
 * - |A ∪ B| = sum of max(counts) for hash functions where both sketches have the same minimum hash
 *            + sum of counts for hash functions where sketches have different minimum hashes
 * 
 * @param other The other BagMinHash sketch to compare with
 * @return Jaccard similarity (0.0 to 1.0)
 */
double BagMinHashSketch::jaccard_similarity(const AbstractSketch& other) const {
    const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
    if (!bmh_other) {
        throw std::runtime_error("Cannot compute Jaccard between different sketch types");
    }
    
    if (min_hashes_.size() != bmh_other->min_hashes_.size()) {
        throw std::runtime_error("BagMinHash sketches must have same number of hashes");
    }
    
    double intersection_sum = 0.0;
    double union_sum = 0.0;
    
    // Compare each hash function's minimum hash and count
    for (size_t i = 0; i < min_hashes_.size(); ++i) {
        if (min_hashes_[i].first == bmh_other->min_hashes_[i].first) {
            // Same minimum hash - intersection is min of counts, union is max of counts
            intersection_sum += std::min(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
            union_sum += std::max(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
        } else {
            // Different minimum hashes - no intersection, union is sum of both counts
            union_sum += min_hashes_[i].second + bmh_other->min_hashes_[i].second;
        }
    }
    
    return (union_sum > 0) ? (intersection_sum / union_sum) : 0.0;
}

/**
 * Estimate the cardinality (total count) of the sketch
 * 
 * This is a simplified cardinality estimator for BagMinHash. In practice,
 * more sophisticated estimators could be used that take into account the
 * distribution of minimum hash values and their counts.
 * 
 * Current implementation: average of all counts across hash functions
 * 
 * @return Estimated cardinality of the sketched set
 */
double BagMinHashSketch::cardinality() const {
    // For BagMinHash, cardinality estimation is more complex
    // This is a simplified version - in practice you'd want a more sophisticated estimator
    double total_count = 0.0;
    for (const auto& pair : min_hashes_) {
        total_count += pair.second;
    }
    return total_count / num_hashes_;
}

/**
 * Compute the estimated intersection size between two BagMinHash sketches
 * 
 * The intersection size is computed as the sum of minimum counts for hash functions
 * where both sketches have the same minimum hash value.
 * 
 * @param other The other BagMinHash sketch to compare with
 * @return Estimated intersection size
 */
double BagMinHashSketch::intersection_size(const AbstractSketch& other) const {
    const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
    if (!bmh_other) {
        throw std::runtime_error("Cannot compute intersection between different sketch types");
    }
    
    double intersection_sum = 0.0;
    for (size_t i = 0; i < min_hashes_.size(); ++i) {
        if (min_hashes_[i].first == bmh_other->min_hashes_[i].first) {
            // Same minimum hash - intersection size is the minimum of the two counts
            intersection_sum += std::min(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
        }
        // If different minimum hashes, intersection is 0 for this hash function
    }
    return intersection_sum;
}

/**
 * Compute the union of two BagMinHash sketches
 * 
 * For each hash function, the union takes:
 * - The minimum hash value from whichever sketch has the smaller minimum
 * - If both sketches have the same minimum hash value, take the maximum count
 * 
 * @param other The other BagMinHash sketch to union with
 * @return A new BagMinHash sketch representing the union
 */
std::unique_ptr<AbstractSketch> BagMinHashSketch::union_with(const AbstractSketch& other) const {
    const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
    if (!bmh_other) {
        throw std::runtime_error("Cannot compute union between different sketch types");
    }
    
    auto result = std::make_unique<BagMinHashSketch>(num_hashes_, seed_);
    
    for (size_t i = 0; i < min_hashes_.size(); ++i) {
        if (min_hashes_[i].first < bmh_other->min_hashes_[i].first) {
            // This sketch has the smaller minimum hash - use its value and count
            result->min_hashes_[i] = min_hashes_[i];
        } else if (bmh_other->min_hashes_[i].first < min_hashes_[i].first) {
            // Other sketch has the smaller minimum hash - use its value and count
            result->min_hashes_[i] = bmh_other->min_hashes_[i];
        } else {
            // Same minimum hash - take the maximum count (union of bags)
            result->min_hashes_[i] = {min_hashes_[i].first, 
                                     std::max(min_hashes_[i].second, bmh_other->min_hashes_[i].second)};
        }
    }
    
    return result;
}

/**
 * Get the sketch type identifier
 * 
 * @return String identifier for this sketch type
 */
std::string BagMinHashSketch::get_sketch_type() const {
    return "bagminhash";
}

/**
 * Clear the sketch (reset all hash functions to initial state)
 * 
 * Resets all minimum hash values to UINT64_MAX and counts to 0,
 * effectively clearing the sketch back to its initial empty state.
 */
void BagMinHashSketch::clear() {
    std::fill(min_hashes_.begin(), min_hashes_.end(), std::make_pair(UINT64_MAX, 0));
}

/**
 * Add a hash value with normalized weight
 * 
 * This method normalizes the raw count to a scale-comparable weight.
 * The normalization uses a logarithmic scaling approach to compress the
 * range of counts while preserving relative differences.
 * 
 * Formula: weight = log2(count + 1) * scale_factor
 * Where scale_factor is chosen to keep weights in a reasonable range.
 * 
 * This is the primary method for adding weighted elements to BagMinHash sketches.
 * Raw counts are automatically normalized to ensure scale-comparable weights
 * across different datasets.
 * 
 * @param hash_val The hash value to add to the sketch
 * @param raw_count The raw count from the BED file
 * @param scale_factor Scaling factor for normalization (default: 100.0)
 */
void BagMinHashSketch::add_with_normalized_count(uint64_t hash_val, uint64_t raw_count, double scale_factor) {
    // Normalize the raw count using logarithmic scaling
    double normalized_weight = std::log2(raw_count + 1.0) * scale_factor;
    uint64_t weight = static_cast<uint64_t>(std::round(normalized_weight));
    
    // Ensure minimum weight of 1
    if (weight < 1) weight = 1;
    
    // Apply each hash function to the input value
    for (size_t i = 0; i < num_hashes_; ++i) {
        uint64_t hash_seed = seed_ + i;  // Each hash function uses a different seed
        uint64_t hashed_value = hash_with_seed(hash_val, hash_seed);
        
        if (hashed_value < min_hashes_[i].first) {
            // Found a new minimum hash value - update both value and count
            min_hashes_[i].first = hashed_value;
            min_hashes_[i].second = weight;  // Set count to the normalized weight
        } else if (hashed_value == min_hashes_[i].first) {
            // Same minimum hash value - add the weight to existing count
            min_hashes_[i].second += weight;
        }
        // If hashed_value > min_hashes_[i].first, we ignore it (not a new minimum)
    }
}
