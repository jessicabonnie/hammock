#include "../include/hll_sketch.h"
#include <stdexcept>

HLLSketch::HLLSketch(size_t precision) : sketch_(precision) {}

void HLLSketch::add(uint64_t hash_val) {
    sketch_.add(hash_val);
}

double HLLSketch::jaccard_similarity(const AbstractSketch& other) const {
    const HLLSketch* hll_other = dynamic_cast<const HLLSketch*>(&other);
    if (!hll_other) {
        throw std::runtime_error("Cannot compute Jaccard between different sketch types");
    }
    return sketch_.jaccard_similarity_registers(hll_other->sketch_);
}

double HLLSketch::cardinality() const {
    return const_cast<hll::hll_t&>(sketch_).report_ertl_improved();
}

double HLLSketch::intersection_size(const AbstractSketch& other) const {
    const HLLSketch* hll_other = dynamic_cast<const HLLSketch*>(&other);
    if (!hll_other) {
        throw std::runtime_error("Cannot compute intersection between different sketch types");
    }
    return hll::intersection_size(sketch_, hll_other->sketch_);
}

std::unique_ptr<AbstractSketch> HLLSketch::union_with(const AbstractSketch& other) const {
    const HLLSketch* hll_other = dynamic_cast<const HLLSketch*>(&other);
    if (!hll_other) {
        throw std::runtime_error("Cannot compute union between different sketch types");
    }
    auto result = std::make_unique<HLLSketch>(sketch_.get_np());
    result->sketch_ = hll::operator+(sketch_, hll_other->sketch_);
    return result;
}

std::string HLLSketch::get_sketch_type() const {
    return "hyperloglog";
}

void HLLSketch::clear() {
    sketch_.clear();
}
