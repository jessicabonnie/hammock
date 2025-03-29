#!/usr/bin/env python3
"""
Compatibility wrapper for native RustHLL to match the API expected by benchmark scripts.
"""

try:
    import rust_hll
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False

class RustHLLWrapper:
    """Wrapper for native RustHLL to make it API-compatible with the benchmark scripts."""
    
    def __init__(self, precision=12, use_threading=True, min_thread_batch=50000):
        """Initialize a new RustHLL instance."""
        if not RUST_AVAILABLE:
            raise ImportError("RustHLL module is not available")
        
        self._sketch = rust_hll.RustHLL(precision, use_threading, min_thread_batch)
    
    def add(self, value):
        """Add a value to the sketch."""
        self._sketch.add_value(str(value))
    
    def add_string(self, value):
        """Add a string value to the sketch."""
        self._sketch.add_value(value)
    
    def add_batch(self, values):
        """Add a batch of values to the sketch."""
        # Convert all values to strings for consistent handling
        string_values = [str(value) for value in values]
        self._sketch.add_batch(string_values)
    
    def cardinality(self):
        """Estimate the cardinality."""
        return self._sketch.estimate_cardinality()
    
    def estimate_cardinality(self, method=None):
        """Estimate the cardinality (method parameter is ignored)."""
        return self._sketch.estimate_cardinality()
    
    def merge(self, other):
        """Merge with another sketch."""
        if isinstance(other, RustHLLWrapper):
            self._sketch.merge(other._sketch)
        else:
            # Try to access the underlying native sketch if possible
            try:
                self._sketch.merge(other._sketch)
            except (AttributeError, TypeError):
                raise TypeError("Cannot merge with incompatible sketch type")
    
    def jaccard(self, other):
        """Calculate Jaccard similarity with another sketch."""
        if isinstance(other, RustHLLWrapper):
            return self._sketch.jaccard(other._sketch)
        else:
            # Try to access the underlying native sketch if possible
            try:
                return self._sketch.jaccard(other._sketch)
            except (AttributeError, TypeError):
                raise TypeError("Cannot calculate Jaccard with incompatible sketch type")
    
    def __str__(self):
        """Return a string representation."""
        return f"RustHLLWrapper(precision={self.precision})"
    
    @property
    def precision(self):
        """Get the sketch precision."""
        # This attribute is needed for the benchmark
        # We don't have direct access to it, so use debug_info to extract it
        debug_info = self._sketch.debug_info()
        precision_str = debug_info.split("precision=")[1].split(",")[0]
        return int(precision_str) 