"""
Fast HyperLogLog implementation using Rust.
"""
from __future__ import annotations
from typing import Optional, Dict, Any, List, Set, Tuple, Union
import sys
import os

from hammock.lib.abstractsketch import AbstractSketch
import numpy as np  # type: ignore
import xxhash  # type: ignore

# Try to import the Rust extension
try:
    # Add the rust_hll directory to the Python path
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "rust_hll"))
    # Also add the target/release directory
    rust_target_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "rust_hll", "target", "release")
    sys.path.append(rust_target_dir)
    
    # Import the Rust module
    import rust_hll
    RUST_AVAILABLE = True
    print(f"Found Rust HyperLogLog module")
except ImportError as e:
    print(f"Error importing Rust HyperLogLog module: {e}")
    RUST_AVAILABLE = False
    rust_hll = None

# If the Rust extension is not available, use the Python implementation
if not RUST_AVAILABLE:
    from .hyperloglog import HyperLogLog as PyHyperLogLog
    print("Rust HyperLogLog extension not found. Using Python implementation.")


class RustHyperLogLog(AbstractSketch):
    """Rust implementation of HyperLogLog."""
    
    def __init__(self, 
                 precision: int = 12,
                 kmer_size: int = 0,
                 window_size: int = 0,
                 seed: Optional[int] = None,
                 expected_cardinality: Optional[int] = None,
                 hash_size: int = 32,
                 debug: bool = False):
        """Initialize the sketch.
        
        Args:
            precision: Number of bits to use for register addressing.
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            expected_cardinality: Expected number of unique items. If provided, precision will be adjusted.
            hash_size: Size of hash in bits (32 or 64)
            debug: Debug flag (not used, provided for backward compatibility)
        """
        super().__init__()
        
        if hash_size not in [32, 64]:
            raise ValueError("hash_size must be 32 or 64")
            
        if expected_cardinality is not None:
            # Adjust precision based on expected cardinality
            if expected_cardinality < 1000:
                precision = 8  # For small sets, use lower precision
            elif expected_cardinality < 10000:
                precision = 12  # For medium sets
            elif expected_cardinality < 100000:
                precision = 18  # For larger sets
            else:
                precision = 22  # For very large sets
        
        # Enforce only theoretical limits
        if precision < 4:
            raise ValueError(f"Precision must be at least 4")
        if precision >= hash_size:
            raise ValueError(f"Precision must be less than hash_size ({hash_size})")
            
        if kmer_size < 0:
            raise ValueError("k-mer size must be non-negative")
            
        if window_size and window_size < kmer_size:
            raise ValueError("Window size must be >= kmer size")
        
        self.precision = precision
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.seed = seed if seed is not None else 42
        self.hash_size = hash_size
        
        if RUST_AVAILABLE:
            try:
                # Try to create the Rust HLL - no arbitrary precision limit
                self._sketch = rust_hll.RustHLL(precision, None, None, hash_size)
                self._using_rust = True
            except ValueError as e:
                # Fallback to Python implementation if Rust fails
                print(f"Warning: Rust implementation rejected precision={precision}: {e}")
                print(f"Falling back to Python implementation")
                self._using_rust = False
                self.registers = np.zeros(1 << precision, dtype=np.float64)
                
                # Create a wrapper object that redirects calls to fallback methods
                class FallbackWrapper:
                    def __init__(self, parent):
                        self.parent = parent
                    
                    def add_value(self, value):
                        return self.parent.add_value(value)
                    
                    def add_batch(self, values):
                        return self.parent.add_batch_fallback(values)
                    
                    def merge(self, other):
                        if hasattr(other, '_sketch'):
                            other = other._sketch
                        return self.parent.merge_fallback(other)
                    
                    def estimate_cardinality(self):
                        return self.parent.estimate_cardinality_fallback()
                    
                    def is_empty(self):
                        return self.parent.is_empty_fallback()
                    
                    def debug_info(self):
                        return f"FallbackWrapper(precision={self.parent.precision}, using_rust=False)"
                
                self._sketch = FallbackWrapper(self)
                return
        else:
            print(f"Rust HLL not available, using Python implementation with precision={precision}")
            self.registers = np.zeros(1 << precision, dtype=np.float64)
            self._using_rust = False
            
            # Create a wrapper object that redirects calls to fallback methods
            class FallbackWrapper:
                def __init__(self, parent):
                    self.parent = parent
                
                def add_value(self, value):
                    return self.parent.add_value(value)
                
                def add_batch(self, values):
                    return self.parent.add_batch_fallback(values)
                
                def merge(self, other):
                    if hasattr(other, '_sketch'):
                        other = other._sketch
                    return self.parent.merge_fallback(other)
                
                def estimate_cardinality(self):
                    return self.parent.estimate_cardinality_fallback()
                
                def is_empty(self):
                    return self.parent.is_empty_fallback()
                
                def debug_info(self):
                    return f"FallbackWrapper(precision={self.parent.precision}, using_rust=False)"
            
            self._sketch = FallbackWrapper(self)
        
    def is_using_rust(self) -> bool:
        """Check if using Rust implementation."""
        return self._using_rust
        
    def add(self, value: Union[str, int, float]) -> None:
        """Add a value to the sketch."""
        # Convert to string and use add for Rust implementation
        if isinstance(value, (int, float)):
            value = str(value)
        self._sketch.add_value(value)
        
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self._sketch.add_value(s)
        
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch."""
        # Convert any non-string values to strings
        string_values = [str(s) if not isinstance(s, str) else s for s in strings]
        self._sketch.add_batch(string_values)
        
    def estimate_cardinality(self, method: str = "fast_mle") -> float:
        """Estimate the cardinality of the set.
        
        Args:
            method: The estimation method to use. One of:
                - "original": Original HyperLogLog estimation
                - "ertl_mle": ERTL Maximum Likelihood Estimation
                - "fast_mle": Fast MLE implementation (default)
        """
        if method not in ["original", "ertl_mle", "fast_mle"]:
            raise ValueError("Method must be one of: original, ertl_mle, fast_mle")
        return float(self._sketch.estimate_cardinality())
        
    def estimate_intersection(self, other: 'RustHyperLogLog') -> float:
        """Estimate intersection cardinality with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compute intersection of HLLs with different precision")
        
        # Basic inclusion-exclusion
        a = self.estimate_cardinality()
        b = other.estimate_cardinality()
        union = self.estimate_union(other)
        intersection = max(0.0, a + b - union)
        
        return float(intersection)
        
    def estimate_union(self, other: 'RustHyperLogLog') -> float:
        """Estimate union cardinality with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compute union of HLLs with different precision")
        
        # Create a temporary sketch for the union
        temp_sketch = rust_hll.RustHLL(self.precision, None, None, self.hash_size)
        temp_sketch.merge(self._sketch)
        temp_sketch.merge(other._sketch)
        return float(temp_sketch.estimate_cardinality())
        
    def estimate_jaccard(self, other: 'RustHyperLogLog') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compare HLLs with different precision")
        
        if self.is_empty() or other.is_empty():
            return 0.0
        
        intersection = self.estimate_intersection(other)
        union = self.estimate_union(other)
        
        if union == 0:
            return 0.0
            
        return float(intersection / union)
        
    def is_empty(self) -> bool:
        """Check if the sketch is empty."""
        return self._sketch.is_empty()
        
    def merge(self, other: 'RustHyperLogLog') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only merge with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        self._sketch.merge(other._sketch)
        
    def save(self, filepath: str) -> None:
        """Write sketch to file."""
        # Convert to absolute path and ensure directory exists
        filepath = os.path.abspath(filepath)
        dir_path = os.path.dirname(filepath)
        if dir_path:
            os.makedirs(dir_path, exist_ok=True)
        try:
            self._sketch.write(filepath)
        except Exception as e:
            print(f"Error saving sketch to {filepath}: {e}", file=sys.stderr)
            raise
        
    @classmethod
    def load(cls, filepath: str) -> 'RustHyperLogLog':
        """Load sketch from file."""
        # Convert to absolute path
        filepath = os.path.abspath(filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Sketch file not found: {filepath}")
        try:
            sketch = cls()
            sketch._sketch = rust_hll.RustHLL.load(filepath)
            return sketch
        except Exception as e:
            print(f"Error loading sketch from {filepath}: {e}", file=sys.stderr)
            raise
    
    def cardinality(self) -> float:
        """Get the estimated cardinality (number of unique values)."""
        return float(self._sketch.estimate_cardinality())
    
    def jaccard(self, other: 'RustHyperLogLog') -> float:
        """Calculate the Jaccard similarity with another sketch."""
        return float(self.estimate_jaccard(other))
    
    def similarity_values(self, other: 'RustHyperLogLog') -> Dict[str, float]:
        """Get similarity values between this sketch and another."""
        return {"jaccard": float(self.jaccard(other))}
    
    def __str__(self) -> str:
        """Get string representation."""
        return f"RustHyperLogLog(precision={self.precision}, kmer_size={self.kmer_size}, window_size={self.window_size}, seed={self.seed}, hash_size={self.hash_size})"
    
    def __repr__(self) -> str:
        """Get representation string."""
        return self.__str__()

    def add_value(self, value: str) -> None:
        """Add a value to the sketch (Python fallback)."""
        if self._using_rust:
            raise RuntimeError("This method should not be called when using Rust")
        # Convert value to bytes if needed
        if isinstance(value, str):
            value = value.encode() if hasattr(value, 'encode') else str(value).encode()
        hash_val = self.hash_str(value)
        idx = hash_val & ((1 << self.precision) - 1)
        # Get rank using _rho function
        rank = self._rho(hash_val)
        
        # Update register
        self.registers[idx] = max(self.registers[idx], rank)
    
    def hash_str(self, s: bytes) -> int:
        """Hash a string for Python fallback."""
        import xxhash  # type: ignore
        hasher = xxhash.xxh64(seed=self.seed) if self.hash_size == 64 else xxhash.xxh32(seed=self.seed)
        hasher.update(s)
        return hasher.intdigest()
    
    def _get_alpha(self, m: float) -> float:
        """Get alpha correction factor for HyperLogLog."""
        if m == 16:
            return 0.673
        elif m == 32:
            return 0.697
        elif m == 64:
            return 0.709
        else:
            return 0.7213 / (1.0 + 1.079 / m)
            
    def add_batch_fallback(self, values: List[str]) -> None:
        """Add multiple values to the sketch (Python fallback)."""
        if self._using_rust:
            raise RuntimeError("This method should not be called when using Rust")
        for value in values:
            self.add_value(value)
    
    def merge_fallback(self, other) -> None:
        """Merge another sketch (Python fallback)."""
        if self._using_rust:
            raise RuntimeError("This method should not be called when using Rust")
        if hasattr(other, 'registers'):
            # Take element-wise maximum
            np.maximum(self.registers, other.registers, out=self.registers)
    
    def estimate_cardinality_fallback(self) -> float:
        """Estimate cardinality (Python fallback)."""
        if self._using_rust:
            raise RuntimeError("This method should not be called when using Rust")
        # Implement a basic HyperLogLog cardinality estimation
        m = float(1 << self.precision)
        alpha = self._get_alpha(m)
        
        # Calculate the sum of 2^(-register[i])
        register_harmonics = np.power(2.0, -self.registers)
        estimate = alpha * m * m / np.sum(register_harmonics)
        
        # Apply corrections
        if estimate <= 2.5 * m:
            # Small range correction
            zero_count = np.sum(self.registers == 0)
            if zero_count > 0:
                estimate = m * np.log(m / zero_count)
        
        # Large range correction omitted for simplicity
        
        return float(estimate)
    
    def is_empty_fallback(self) -> bool:
        """Check if sketch is empty (Python fallback)."""
        if self._using_rust:
            raise RuntimeError("This method should not be called when using Rust")
        return np.all(self.registers == 0)

    def _rho(self, value: int) -> int:
        """Calculate the position of the leftmost 1-bit."""
        if value == 0:
            return self.hash_size - self.precision + 1
        # Shift right by precision bits to get the remaining hash bits
        value >>= self.precision
        # Count leading zeros in the remaining bits
        rank = 1
        while value & 1 == 0 and rank <= self.hash_size - self.precision:
            rank += 1
            value >>= 1
        return rank

# Alias for backward compatibility
RustHLL = RustHyperLogLog 