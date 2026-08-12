#ifndef HAMMOCK_ABSTRACT_SKETCH_HPP
#define HAMMOCK_ABSTRACT_SKETCH_HPP

#include <cstdint>
#include <memory>
#include <string>

class AbstractSketch {
public:
    virtual ~AbstractSketch() = default;
    virtual void add(uint64_t hash_val) = 0;
    virtual double reg_eq_similarity(const AbstractSketch& other) const = 0;
    virtual double cardinality() const = 0;
    virtual double intersection_size(const AbstractSketch& other) const = 0;
    virtual std::unique_ptr<AbstractSketch> union_with(const AbstractSketch& other) const = 0;
    virtual void clear() = 0;
};

#endif
