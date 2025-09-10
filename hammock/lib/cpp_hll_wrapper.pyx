# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import numpy as np
cimport numpy as np
from libc.stdint cimport uint64_t, uint32_t, uint8_t, int8_t
from libc.string cimport memcpy
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

# Import the C++ hyperloglog implementation
cdef extern from "hll/hll.h" namespace "hll":
    cdef cppclass hll_t:
        hll_t(size_t np) except +
        void add(uint64_t hashval)
        void addh(uint64_t element)
        double report()
        double creport() const
        double est_err()
        double cest_err() const
        void clear()
        void sum()
        size_t get_np() const
        bool is_ready() const
        string to_string() const
        string desc_string() const
        void resize(size_t new_size)
        void free()
        bool within_bounds(uint64_t actual_size) const
        double jaccard_similarity_registers(const hll_t& other) const
    
    # Global functions for set operations
    double operator^(hll_t& first, hll_t& other)
    hll_t operator&(hll_t& first, hll_t& other)
    double intersection_size(hll_t& first, hll_t& other)
    double intersection_size(const hll_t& first, const hll_t& other)
    hll_t operator+(const hll_t& one, const hll_t& other)

# Import xxhash for hashing
import xxhash

cdef class CppHyperLogLog:
    """C++ HyperLogLog implementation wrapped for Python."""
    
    cdef:
        hll_t* _hll
        int precision
        int hash_size
        uint64_t seed
        bool _initialized
    
    def __init__(self, int precision=12, int hash_size=64, uint64_t seed=42):
        """Initialize C++ HyperLogLog sketch.
        
        Args:
            precision: Number of bits for register indexing (4-24)
            hash_size: Size of hash in bits (32 or 64)
            seed: Random seed for hashing
        """
        if precision < 4:
            raise ValueError("Precision must be at least 4")
        # Note: Upper precision limit removed for better flexibility
        # if precision > 24:
        #     raise ValueError("Precision must be at most 24")
        if hash_size not in [32, 64]:
            raise ValueError("hash_size must be 32 or 64")
            
        self.precision = precision
        self.hash_size = hash_size
        self.seed = seed
        self._initialized = False
        
        # Create C++ HyperLogLog instance
        self._hll = new hll_t(precision)
        self._initialized = True
    
    def __dealloc__(self):
        """Clean up C++ object."""
        if self._initialized and self._hll != NULL:
            self._hll.free()
            del self._hll
    
    def add(self, uint64_t element):
        """Add an element to the sketch."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        self._hll.addh(element)
    
    def add_hash(self, uint64_t hashval):
        """Add a pre-computed hash value to the sketch."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        self._hll.add(hashval)
    
    def add_string(self, str text):
        """Add a string to the sketch by hashing it."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        
        # Hash the string using xxhash
        text_bytes = text.encode('utf-8')
        if self.hash_size == 32:
            hash_val = xxhash.xxh32(text_bytes, seed=self.seed).intdigest()
        else:
            hash_val = xxhash.xxh64(text_bytes, seed=self.seed).intdigest()
        
        self._hll.add(hash_val)
    
    def add_bytes(self, bytes data):
        """Add bytes to the sketch by hashing them."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        
        # Hash the bytes using xxhash
        if self.hash_size == 32:
            hash_val = xxhash.xxh32(data, seed=self.seed).intdigest()
        else:
            hash_val = xxhash.xxh64(data, seed=self.seed).intdigest()
        
        self._hll.add(hash_val)
    
    def add_batch(self, list elements):
        """Add multiple elements to the sketch efficiently."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        
        for element in elements:
            if isinstance(element, str):
                self.add_string(element)
            elif isinstance(element, bytes):
                self.add_bytes(element)
            elif isinstance(element, int):
                self.add(element)
            else:
                # Convert to string and hash
                self.add_string(str(element))
    
    def cardinality(self):
        """Get the estimated cardinality."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.report()
    
    def cardinality_const(self):
        """Get the estimated cardinality (const version)."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.creport()
    
    def error_estimate(self):
        """Get the error estimate."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.est_err()
    
    def error_estimate_const(self):
        """Get the error estimate (const version)."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.cest_err()
    
    def clear(self):
        """Clear the sketch."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        self._hll.clear()
    
    def is_ready(self):
        """Check if the sketch is ready for reporting."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.is_ready()
    
    def to_string(self):
        """Get string representation of the sketch."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.to_string().decode('utf-8')
    
    def description(self):
        """Get descriptive string of the sketch."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.desc_string().decode('utf-8')
    
    def within_bounds(self, uint64_t actual_size):
        """Check if estimate is within error bounds of actual size."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll.within_bounds(actual_size)
    
    def resize(self, size_t new_size):
        """Resize the sketch to a new precision."""
        if not self._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        self._hll.resize(new_size)
        self.precision = new_size
    
    def union_(self, CppHyperLogLog other):
        """Create union of two sketches."""
        if not self._initialized or not other._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        result = CppHyperLogLog(self.precision, self.hash_size, self.seed)
        result._hll[0] = self._hll[0] + other._hll[0]
        return result
    
    def intersection(self, CppHyperLogLog other):
        """Create intersection of two sketches."""
        if not self._initialized or not other._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        result = CppHyperLogLog(self.precision, self.hash_size, self.seed)
        result._hll[0] = self._hll[0] & other._hll[0]
        return result
    
    def symmetric_difference(self, CppHyperLogLog other):
        """Calculate symmetric difference size."""
        if not self._initialized or not other._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return self._hll[0] ^ other._hll[0]
    
    def intersection_size(self, CppHyperLogLog other):
        """Calculate intersection size."""
        if not self._initialized or not other._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        return intersection_size(self._hll[0], other._hll[0])
    
    def jaccard_similarity(self, CppHyperLogLog other):
        """Calculate Jaccard similarity."""
        if not self._initialized or not other._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        
        # Calculate intersection size
        intersection = intersection_size(self._hll[0], other._hll[0])
        
        # Calculate union size using the union method
        union_sketch = self.union_(other)
        union_size = union_sketch.cardinality()
        
        if union_size == 0:
            return 0.0
        
        return intersection / union_size
    
    def jaccard_similarity_registers(self, CppHyperLogLog other):
        """Calculate Jaccard similarity using register-based comparison."""
        if not self._initialized or not other._initialized:
            raise RuntimeError("HyperLogLog not initialized")
        
        return self._hll[0].jaccard_similarity_registers(other._hll[0])
    
    @property
    def precision(self):
        """Get the precision of the sketch."""
        return self.precision
    
    @property
    def num_registers(self):
        """Get the number of registers."""
        return 1 << self.precision
