#include "../include/bagminhash_sketch.h"
#include "../include/xxhash.h"
#include <stdexcept>
#include <algorithm>
#include <climits>

BagMinHashSketch::BagMinHashSketch(size_t num_hashes, uint64_t seed) 
    : num_hashes_(num_hashes), min_hashes_(num_hashes, {UINT64_MAX, 0}), seed_(seed) {}

uint64_t BagMinHashSketch::hash_with_seed(uint64_t value, uint64_t seed) const {
    // Convert uint64_t to string for xxhash
    std::string value_str = std::to_string(value);
    return xxhash::hash64(value_str, seed);
}

void BagMinHashSketch::add(uint64_t hash_val) {
    // For BagMinHash, we need to track counts, but since we're getting individual hash values,
    // we'll treat each add as a count of 1. In a real implementation, you'd want to batch
    // adds by the same hash value to properly handle counts.
    for (size_t i = 0; i < num_hashes_; ++i) {
        uint64_t hash_seed = seed_ + i;
        uint64_t hashed_value = hash_with_seed(hash_val, hash_seed);
        
        if (hashed_value < min_hashes_[i].first) {
            min_hashes_[i].first = hashed_value;
            min_hashes_[i].second = 1;  // Reset count when new minimum found
        } else if (hashed_value == min_hashes_[i].first) {
            min_hashes_[i].second++;  // Increment count for same minimum
        }
    }
}

void BagMinHashSketch::add_with_count(uint64_t hash_val, uint64_t count) {
    for (size_t i = 0; i < num_hashes_; ++i) {
        uint64_t hash_seed = seed_ + i;
        uint64_t hashed_value = hash_with_seed(hash_val, hash_seed);
        
        if (hashed_value < min_hashes_[i].first) {
            min_hashes_[i].first = hashed_value;
            min_hashes_[i].second = count;
        } else if (hashed_value == min_hashes_[i].first) {
            min_hashes_[i].second += count;
        }
    }
}

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
    
    for (size_t i = 0; i < min_hashes_.size(); ++i) {
        if (min_hashes_[i].first == bmh_other->min_hashes_[i].first) {
            // Same minimum hash - intersection is min of counts, union is max
            intersection_sum += std::min(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
            union_sum += std::max(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
        } else {
            // Different minimum hashes - no intersection, union is sum of counts
            union_sum += min_hashes_[i].second + bmh_other->min_hashes_[i].second;
        }
    }
    
    return (union_sum > 0) ? (intersection_sum / union_sum) : 0.0;
}

double BagMinHashSketch::cardinality() const {
    // For BagMinHash, cardinality estimation is more complex
    // This is a simplified version - in practice you'd want a more sophisticated estimator
    double total_count = 0.0;
    for (const auto& pair : min_hashes_) {
        total_count += pair.second;
    }
    return total_count / num_hashes_;
}

double BagMinHashSketch::intersection_size(const AbstractSketch& other) const {
    const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
    if (!bmh_other) {
        throw std::runtime_error("Cannot compute intersection between different sketch types");
    }
    
    double intersection_sum = 0.0;
    for (size_t i = 0; i < min_hashes_.size(); ++i) {
        if (min_hashes_[i].first == bmh_other->min_hashes_[i].first) {
            intersection_sum += std::min(min_hashes_[i].second, bmh_other->min_hashes_[i].second);
        }
    }
    return intersection_sum;
}

std::unique_ptr<AbstractSketch> BagMinHashSketch::union_with(const AbstractSketch& other) const {
    const BagMinHashSketch* bmh_other = dynamic_cast<const BagMinHashSketch*>(&other);
    if (!bmh_other) {
        throw std::runtime_error("Cannot compute union between different sketch types");
    }
    
    auto result = std::make_unique<BagMinHashSketch>(num_hashes_, seed_);
    
    for (size_t i = 0; i < min_hashes_.size(); ++i) {
        if (min_hashes_[i].first < bmh_other->min_hashes_[i].first) {
            result->min_hashes_[i] = min_hashes_[i];
        } else if (bmh_other->min_hashes_[i].first < min_hashes_[i].first) {
            result->min_hashes_[i] = bmh_other->min_hashes_[i];
        } else {
            // Same minimum hash - take the maximum count
            result->min_hashes_[i] = {min_hashes_[i].first, 
                                     std::max(min_hashes_[i].second, bmh_other->min_hashes_[i].second)};
        }
    }
    
    return result;
}

std::string BagMinHashSketch::get_sketch_type() const {
    return "bagminhash";
}

void BagMinHashSketch::clear() {
    std::fill(min_hashes_.begin(), min_hashes_.end(), std::make_pair(UINT64_MAX, 0));
}
