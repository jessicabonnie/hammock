"""
C++ HyperLogLog implementation for hammock.

This module provides a Python interface to the C++ HyperLogLog implementation,
integrating with the existing hammock sketch architecture.
"""

import numpy as np
from typing import Dict, Optional, List, Union, Any
from hammock.lib.abstractsketch import AbstractSketch

try:
    from .cpp_hll_wrapper import CppHyperLogLog
    CPP_HLL_AVAILABLE = True
except ImportError:
    CPP_HLL_AVAILABLE = False
    CppHyperLogLog = None

class CppHyperLogLogSketch(AbstractSketch):
    """
    C++ HyperLogLog sketch implementation for hammock.
    
    This class provides a Python interface to the C++ HyperLogLog implementation,
    maintaining compatibility with the existing hammock sketch architecture.
    """
    
    def __init__(self, 
                 precision: int = 12, 
                 kmer_size: int = 0, 
                 window_size: int = 0, 
                 seed: Optional[int] = None,
                 debug: bool = False,
                 expected_cardinality: Optional[int] = None,
                 hash_size: int = 64):
        """Initialize C++ HyperLogLog sketch.
        
        Args:
            precision: Number of bits for register indexing (4-24)
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            debug: Whether to print debug information
            expected_cardinality: Expected number of unique items
            hash_size: Size of hash in bits (32 or 64)
        """
        if not CPP_HLL_AVAILABLE:
            raise ImportError("C++ HyperLogLog implementation not available. "
                            "Please ensure the C++ code is compiled and Cython extension is built.")
        
        super().__init__()
        
        # Validate parameters
        if precision < 4:
            raise ValueError("Precision must be at least 4")
        # Note: Upper precision limit removed for better flexibility
        # if precision > 24:
        #     raise ValueError("Precision must be at most 24")
        if hash_size not in [32, 64]:
            raise ValueError("hash_size must be 32 or 64")
        if kmer_size < 0:
            raise ValueError("k-mer size must be non-negative")
        
        # Adjust precision based on expected cardinality if provided
        if expected_cardinality is not None:
            if expected_cardinality < 1000:
                precision = 8
            elif expected_cardinality < 10000:
                precision = 12
            elif expected_cardinality < 100000:
                precision = 18
            else:
                precision = 22
        
        self.precision = precision
        self.num_registers = 1 << precision
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.debug = debug
        self.seed = seed if seed is not None else 42
        self.hash_size = hash_size
        
        # Initialize C++ HyperLogLog instance
        self._cpp_hll = CppHyperLogLog(precision, hash_size, self.seed)
        
        # Track item count for compatibility
        self.item_count = 0
        
        if debug:
            print(f"C++ HyperLogLog initialized: precision={precision}, "
                  f"hash_size={hash_size}, registers={self.num_registers}")
    
    def add_string(self, s: str) -> None:
        """Add a string to the sketch.
        
        Args:
            s: String to add
        """
        self._cpp_hll.add_string(s)
        self.item_count += 1
    
    def add(self, item: Union[str, bytes, int]) -> None:
        """Add an item to the sketch.
        
        Args:
            item: Item to add (string, bytes, or integer)
        """
        if isinstance(item, str):
            self._cpp_hll.add_string(item)
        elif isinstance(item, bytes):
            self._cpp_hll.add_bytes(item)
        elif isinstance(item, int):
            self._cpp_hll.add(item)
        else:
            # Convert to string and add
            self._cpp_hll.add_string(str(item))
        
        self.item_count += 1
    
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        for s in strings:
            self._cpp_hll.add_string(s)
        self.item_count += len(strings)
    
    def add_batch_mixed(self, items: List[Union[str, bytes, int]]) -> None:
        """Add multiple items to the sketch efficiently.
        
        Args:
            items: List of items to add
        """
        self._cpp_hll.add_batch(items)
        self.item_count += len(items)
    
    def add_kmers(self, sequence: str) -> None:
        """Add k-mers from a sequence to the sketch.
        
        Args:
            sequence: DNA/protein sequence to process
        """
        if self.kmer_size == 0:
            # Whole sequence mode
            self.add(sequence)
        else:
            # K-mer mode
            if len(sequence) >= self.kmer_size:
                kmers = []
                for i in range(len(sequence) - self.kmer_size + 1):
                    kmer = sequence[i:i + self.kmer_size]
                    kmers.append(kmer)
                self.add_batch(kmers)
    
    def add_kmers_with_window(self, sequence: str) -> None:
        """Add k-mers from a sequence using sliding window approach.
        
        Args:
            sequence: DNA/protein sequence to process
        """
        if self.kmer_size == 0:
            # Whole sequence mode
            self.add(sequence)
        else:
            # Sliding window mode
            if len(sequence) >= self.window_size:
                for i in range(len(sequence) - self.window_size + 1):
                    window = sequence[i:i + self.window_size]
                    # Extract k-mers from window
                    for j in range(len(window) - self.kmer_size + 1):
                        kmer = window[j:j + self.kmer_size]
                        self.add(kmer)
    
    def cardinality(self) -> float:
        """Get the estimated cardinality.
        
        Returns:
            Estimated number of unique items
        """
        return self._cpp_hll.cardinality()
    
    def error_estimate(self) -> float:
        """Get the error estimate.
        
        Returns:
            Estimated error in cardinality estimation
        """
        return self._cpp_hll.error_estimate()
    
    def clear(self) -> None:
        """Clear the sketch."""
        self._cpp_hll.clear()
        self.item_count = 0
    
    def similarity_values(self, other: 'CppHyperLogLogSketch') -> Dict[str, float]:
        """Calculate similarity values between two sketches.
        
        Args:
            other: Another C++ HyperLogLog sketch
            
        Returns:
            Dictionary containing similarity measures
        """
        if not isinstance(other, CppHyperLogLogSketch):
            raise TypeError("Can only compare with other C++ HyperLogLog sketches")
        
        # Calculate Jaccard similarity
        jaccard = self._cpp_hll.jaccard_similarity(other._cpp_hll)
        
        # Calculate intersection size
        intersection = self._cpp_hll.intersection_size(other._cpp_hll)
        
        # Calculate union size (approximate)
        union_sketch = self._cpp_hll.union_(other._cpp_hll)
        union_size = union_sketch.cardinality()
        
        # Calculate symmetric difference
        sym_diff = self._cpp_hll.symmetric_difference(other._cpp_hll)
        
        return {
            'jaccard': jaccard,
            'intersection': intersection,
            'union': union_size,
            'symmetric_difference': sym_diff
        }
    
    def union_(self, other: 'CppHyperLogLogSketch') -> 'CppHyperLogLogSketch':
        """Create union of two sketches.
        
        Args:
            other: Another C++ HyperLogLog sketch
            
        Returns:
            New sketch representing the union
        """
        if not isinstance(other, CppHyperLogLogSketch):
            raise TypeError("Can only union with other C++ HyperLogLog sketches")
        
        result = CppHyperLogLogSketch(
            precision=self.precision,
            kmer_size=self.kmer_size,
            window_size=self.window_size,
            seed=self.seed,
            debug=self.debug,
            hash_size=self.hash_size
        )
        
        result._cpp_hll = self._cpp_hll.union_(other._cpp_hll)
        return result
    
    def intersection(self, other: 'CppHyperLogLogSketch') -> 'CppHyperLogLogSketch':
        """Create intersection of two sketches.
        
        Args:
            other: Another C++ HyperLogLog sketch
            
        Returns:
            New sketch representing the intersection
        """
        if not isinstance(other, CppHyperLogLogSketch):
            raise TypeError("Can only intersect with other C++ HyperLogLog sketches")
        
        result = CppHyperLogLogSketch(
            precision=self.precision,
            kmer_size=self.kmer_size,
            window_size=self.window_size,
            seed=self.seed,
            debug=self.debug,
            hash_size=self.hash_size
        )
        
        result._cpp_hll = self._cpp_hll.intersection(other._cpp_hll)
        return result
    
    def within_bounds(self, actual_size: int) -> bool:
        """Check if estimate is within error bounds of actual size.
        
        Args:
            actual_size: Actual number of unique items
            
        Returns:
            True if estimate is within error bounds
        """
        return self._cpp_hll.within_bounds(actual_size)
    
    def to_string(self) -> str:
        """Get string representation of the sketch.
        
        Returns:
            String representation
        """
        return self._cpp_hll.to_string()
    
    def description(self) -> str:
        """Get descriptive string of the sketch.
        
        Returns:
            Descriptive string
        """
        return self._cpp_hll.description()
    
    def get_registers(self) -> np.ndarray:
        """Get the register values as a numpy array.
        
        Note: This is not directly available from the C++ implementation,
        but we can reconstruct it by querying the C++ object.
        
        Returns:
            Numpy array of register values
        """
        # This would require additional C++ methods to expose register data
        # For now, return a placeholder
        raise NotImplementedError("Register access not yet implemented in C++ wrapper")
    
    def save(self, filename: str) -> None:
        """Save the sketch to a file.
        
        Args:
            filename: Output file path
        """
        # This would require serialization methods in the C++ implementation
        raise NotImplementedError("Save functionality not yet implemented")
    
    @classmethod
    def load(cls, filename: str) -> 'CppHyperLogLogSketch':
        """Load a sketch from a file.
        
        Args:
            filename: Input file path
            
        Returns:
            Loaded sketch
        """
        # This would require deserialization methods in the C++ implementation
        raise NotImplementedError("Load functionality not yet implemented")
    
    def __repr__(self) -> str:
        """String representation of the sketch."""
        return (f"CppHyperLogLogSketch(precision={self.precision}, "
                f"hash_size={self.hash_size}, cardinality={self.cardinality():.2f})")
    
    def __str__(self) -> str:
        """String representation of the sketch."""
        return self.__repr__()


# Factory function for easy creation
def create_cpp_hll_sketch(precision: int = 12, 
                         kmer_size: int = 0,
                         window_size: int = 0,
                         seed: Optional[int] = None,
                         debug: bool = False,
                         expected_cardinality: Optional[int] = None,
                         hash_size: int = 64) -> CppHyperLogLogSketch:
    """Create a C++ HyperLogLog sketch.
    
    Args:
        precision: Number of bits for register indexing (4-24)
        kmer_size: Size of k-mers (0 for whole string mode)
        window_size: Size of sliding window (0 or == kmer_size for no windowing)
        seed: Random seed for hashing
        debug: Whether to print debug information
        expected_cardinality: Expected number of unique items
        hash_size: Size of hash in bits (32 or 64)
        
    Returns:
        C++ HyperLogLog sketch instance
    """
    return CppHyperLogLogSketch(
        precision=precision,
        kmer_size=kmer_size,
        window_size=window_size,
        seed=seed,
        debug=debug,
        expected_cardinality=expected_cardinality,
        hash_size=hash_size
    )
