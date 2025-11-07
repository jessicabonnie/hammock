#ifndef ABSTRACT_SKETCH_H
#define ABSTRACT_SKETCH_H

#include <memory>
#include <string>
#include <cstdint>

// Abstract base class for all sketch types
class AbstractSketch {
public:
    virtual ~AbstractSketch() = default;
    virtual void add(uint64_t hash_val) = 0;
    virtual double jaccard_similarity(const AbstractSketch& other) const = 0;
    virtual double cardinality() const = 0;
    virtual double intersection_size(const AbstractSketch& other) const = 0;
    virtual std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const = 0;
    virtual std::string get_sketch_type() const = 0;
    virtual void clear() = 0;
};

#endif // ABSTRACT_SKETCH_H
