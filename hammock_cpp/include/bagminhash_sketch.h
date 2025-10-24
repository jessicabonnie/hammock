#ifndef BAGMINHASH_SKETCH_H
#define BAGMINHASH_SKETCH_H

#include "abstract_sketch.h"
#include <vector>
#include <string>
#include <cstdint>

class BagMinHashSketch : public AbstractSketch {
private:
    size_t num_hashes_;
    std::vector<std::pair<uint64_t, uint64_t>> min_hashes_;  // (hash, count) pairs
    uint64_t seed_;
    
    // Hash function for generating multiple hash values
    uint64_t hash_with_seed(uint64_t value, uint64_t seed) const;
    
public:
    explicit BagMinHashSketch(size_t num_hashes, uint64_t seed = 0);
    
    void add(uint64_t hash_val) override;
    void add_with_count(uint64_t hash_val, uint64_t count);
    double jaccard_similarity(const AbstractSketch& other) const override;
    double cardinality() const override;
    double intersection_size(const AbstractSketch& other) const override;
    std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const override;
    std::string get_sketch_type() const override;
    void clear() override;
};

#endif // BAGMINHASH_SKETCH_H
