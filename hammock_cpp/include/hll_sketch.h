#ifndef HLL_SKETCH_H
#define HLL_SKETCH_H

#include "abstract_sketch.h"
#include "../hll/hll.h"
#include <memory>
#include <string>

class HLLSketch : public AbstractSketch {
private:
    hll::hll_t sketch_;
    
public:
    explicit HLLSketch(size_t precision);
    
    void add(uint64_t hash_val) override;
    double jaccard_similarity(const AbstractSketch& other) const override;
    double cardinality() const override;
    double intersection_size(const AbstractSketch& other) const override;
    std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const override;
    std::string get_sketch_type() const override;
    void clear() override;
};

#endif // HLL_SKETCH_H
