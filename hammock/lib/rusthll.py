"""
Fast HyperLogLog implementation using Rust.
"""
from __future__ import annotations
from typing import Optional, Dict, Any, List, Set, Tuple, Union
import sys
import os
import resource
import math

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
                 debug: bool = False,
                 apply_jaccard_correction: bool = True):
        """Initialize the sketch.
        
        Args:
            precision: Number of bits to use for register addressing.
            kmer_size: Size of k-mers (0 for whole string mode).
            window_size: Size of sliding window (0 or == kmer_size for no windowing).
            seed: Random seed for hashing.
            expected_cardinality: Expected number of unique elements.
            hash_size: Size of hash in bits (32 or 64).
            debug: Whether to enable debug mode.
            apply_jaccard_correction: Whether to apply empirical correction to Jaccard values
        """
        # Adjust precision based on expected cardinality if provided
        if expected_cardinality is not None:
            if expected_cardinality <= 100:
                precision = 4  # Minimum precision for small sets
            elif expected_cardinality <= 1000:
                precision = 8  # Medium precision for medium sets
            else:
                precision = max(12, (expected_cardinality.bit_length() + 1) // 2)
        
        if precision < 4 or precision >= hash_size:
            raise ValueError(f"Precision must be between 4 and {hash_size-1}")
            
        if window_size and window_size < kmer_size:
            raise ValueError("Window size must be >= kmer size")
        
        self.precision = precision
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.seed = seed if seed is not None else 42
        self.hash_size = hash_size
        self.apply_jaccard_correction = apply_jaccard_correction
        
        if RUST_AVAILABLE:
            try:
                # Calculate memory usage and set conservative limit
                num_registers = 1 << precision
                register_size = 1  # bytes per register
                estimated_memory = num_registers * register_size
                
                # More conservative memory limits
                batch_overhead = 1024 * 1024 * 5  # 5MB for batch processing
                min_memory = 16 * 1024 * 1024  # 16MB
                memory_limit = max(
                    min_memory,
                    min(
                        512 * 1024 * 1024,  # 512MB
                        int(estimated_memory * 2.0) + batch_overhead
                    )
                )
                
                if debug:
                    print(f"Memory calculation details:")
                    print(f"  - Precision: {precision}")
                    print(f"  - Number of registers: {num_registers}")
                    print(f"  - Register size: {register_size} bytes")
                    print(f"  - Estimated memory: {estimated_memory} bytes")
                    print(f"  - Batch overhead: {batch_overhead} bytes")
                    print(f"  - Total memory limit: {memory_limit} bytes")
                    print(f"  - Memory limit in GB: {memory_limit / (1024*1024*1024):.6f}GB")
                    print(f"  - System memory limit: {resource.getrlimit(resource.RLIMIT_AS)[0]} bytes")
                
                self._rust_sketch = rust_hll.RustHLL(
                    precision, 
                    use_threading=True,
                    min_thread_batch=10000,
                    hash_size=hash_size,
                    memory_limit=memory_limit,
                    debug=debug
                )
                self._using_rust = True
            except ValueError as e:
                print(f"Warning: Rust implementation rejected precision={precision}: {e}")
                print(f"Falling back to Python implementation")
                self._using_rust = False
                self.registers = np.zeros(1 << precision, dtype=np.float64)
            except Exception as e:
                print(f"Warning: Unexpected error initializing Rust implementation: {e}")
                print(f"Falling back to Python implementation")
                self._using_rust = False
                self.registers = np.zeros(1 << precision, dtype=np.float64)
        else:
            print(f"Rust HLL not available, using Python implementation with precision={precision}")
            self.registers = np.zeros(1 << precision, dtype=np.float64)
            self._using_rust = False
        
        # Create a wrapper object that redirects calls to appropriate implementation
        class SketchWrapper:
            def __init__(self, parent):
                self.parent = parent
            
            @property
            def debug(self):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.debug
                return False
            
            @debug.setter
            def debug(self, value):
                if self.parent._using_rust:
                    self.parent._rust_sketch.debug = value

            @property
            def hash_size(self):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.hash_size
                return self.parent.hash_size

            @property
            def seed(self):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.seed
                return self.parent.seed
            
            def add_value(self, value):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.add_value(value)
                else:
                    return self.parent.add_value(value)
            
            def add_batch(self, values):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.add_batch(values)
                else:
                    return self.parent.add_batch_fallback(values)
            
            def merge(self, other):
                if self.parent._using_rust:
                    if hasattr(other, '_rust_sketch'):
                        return self.parent._rust_sketch.merge(other._rust_sketch)
                    elif hasattr(other, '_sketch') and hasattr(other._sketch, 'parent') and hasattr(other._sketch.parent, '_rust_sketch'):
                        return self.parent._rust_sketch.merge(other._sketch.parent._rust_sketch)
                    elif isinstance(other, rust_hll.RustHLL):
                        return self.parent._rust_sketch.merge(other)
                    else:
                        raise ValueError("Cannot merge with incompatible sketch type")
                else:
                    if hasattr(other, '_sketch'):
                        other = other._sketch
                    return self.parent.merge_fallback(other)
            
            def estimate_cardinality(self):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.estimate_cardinality()
                else:
                    return self.parent.estimate_cardinality_fallback()
            
            def is_empty(self):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.is_empty()
                else:
                    return self.parent.is_empty_fallback()
            
            def debug_info(self):
                return f"SketchWrapper(precision={self.parent.precision}, using_rust={self.parent._using_rust})"
            
            def similarity_values(self, other):
                if self.parent._using_rust:
                    if hasattr(other, '_sketch'):
                        other = other._sketch
                    return self.parent._rust_sketch.similarity_values(other)
                else:
                    if hasattr(other, '_sketch'):
                        other = other._sketch
                    return {"jaccard": self.parent.estimate_jaccard(other)}
            
            def estimate_union(self, other):
                if self.parent._using_rust:
                    if hasattr(other, '_rust_sketch'):
                        return self.parent._rust_sketch.estimate_union(other._rust_sketch)
                    elif hasattr(other, '_sketch') and hasattr(other._sketch, 'parent') and hasattr(other._sketch.parent, '_rust_sketch'):
                        return self.parent._rust_sketch.estimate_union(other._sketch.parent._rust_sketch)
                    elif isinstance(other, rust_hll.RustHLL):
                        return self.parent._rust_sketch.estimate_union(other)
                    else:
                        raise ValueError("Cannot compute union with incompatible sketch type")
                else:
                    if hasattr(other, '_sketch'):
                        other = other._sketch
                    return self.parent.estimate_union_fallback(other)
            
            def estimate_intersection(self, other):
                if self.parent._using_rust:
                    if hasattr(other, '_rust_sketch'):
                        return self.parent._rust_sketch.estimate_intersection(other._rust_sketch)
                    elif hasattr(other, '_sketch') and hasattr(other._sketch, 'parent') and hasattr(other._sketch.parent, '_rust_sketch'):
                        return self.parent._rust_sketch.estimate_intersection(other._sketch.parent._rust_sketch)
                    elif isinstance(other, rust_hll.RustHLL):
                        return self.parent._rust_sketch.estimate_intersection(other)
                    else:
                        raise ValueError("Cannot compute intersection with incompatible sketch type")
                else:
                    if hasattr(other, '_sketch'):
                        other = other._sketch
                    return self.parent.estimate_intersection_fallback(other)
            
            def estimate_jaccard(self, other):
                """Estimate Jaccard similarity with another sketch."""
                # Get parent objects
                parent = self.parent
                
                if parent._using_rust:
                    if hasattr(other, '_rust_sketch'):
                        raw_jaccard = parent._rust_sketch.jaccard_improved(other._rust_sketch)
                        return parent._correct_jaccard_estimate(raw_jaccard)
                    elif hasattr(other, '_sketch') and hasattr(other._sketch, 'parent') and hasattr(other._sketch.parent, '_rust_sketch'):
                        raw_jaccard = parent._rust_sketch.jaccard_improved(other._sketch.parent._rust_sketch)
                        return parent._correct_jaccard_estimate(raw_jaccard)
                    elif isinstance(other, rust_hll.RustHLL):
                        raw_jaccard = parent._rust_sketch.jaccard_improved(other)
                        return parent._correct_jaccard_estimate(raw_jaccard)
                    else:
                        raise ValueError("Cannot compute Jaccard similarity with incompatible sketch type")
                
                return parent.estimate_jaccard_fallback(other)
            
            def write(self, filepath):
                if self.parent._using_rust:
                    return self.parent._rust_sketch.write(filepath)
                else:
                    return self.parent.write_fallback(filepath)
        
        self._sketch = SketchWrapper(self)
        
    def is_using_rust(self) -> bool:
        """Check if using Rust implementation."""
        return self._using_rust
        
    def add(self, value: Union[str, int, float]) -> None:
        """Add a value to the sketch."""
        # Convert to string and use add for Rust implementation
        if isinstance(value, (int, float)):
            value = str(value)
        if self._using_rust:
            self._sketch.add_value(value)
        else:
            self.add_value(value)
        
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        if self._using_rust:
            self._sketch.add_value(s)
        else:
            self.add_value(s)
        
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch."""
        if self._using_rust:
            try:
                # Convert any non-string values to strings
                string_values = [str(s) if not isinstance(s, str) else s for s in strings]
                
                # Process in chunks to prevent memory issues
                chunk_size = 5000
                total_chunks = (len(string_values) + chunk_size - 1) // chunk_size
                
                # Only report progress for very large batches (>1M strings)
                report_progress = len(string_values) > 1000000
                
                if report_progress:
                    print(f"Processing {len(string_values):,} strings in {total_chunks} chunks...", file=sys.stderr)
                
                for i in range(0, len(string_values), chunk_size):
                    chunk = string_values[i:i + chunk_size]
                    
                    # Only report progress every 100 chunks for very large batches
                    if report_progress and i % (chunk_size * 500) == 0 and self.debug:
                        print(f"  Progress: {i//chunk_size}/{total_chunks} chunks", file=sys.stderr)
                    
                    try:
                        self._sketch.add_batch(chunk)
                    except Exception as e:
                        print(f"Warning: Failed to process chunk {i//chunk_size + 1}/{total_chunks}: {str(e)}", file=sys.stderr)
                        # Fall back to single-threaded processing for this chunk
                        for s in chunk:
                            try:
                                self._sketch.add_value(s)
                            except Exception as e2:
                                print(f"Warning: Failed to process string in fallback mode: {str(e2)}", file=sys.stderr)
                                continue
            except Exception as e:
                print(f"Warning: Error in add_batch: {str(e)}", file=sys.stderr)
                # Try to process strings one at a time as a last resort
                for s in strings:
                    try:
                        self._sketch.add_value(s)
                    except Exception as e2:
                        print(f"Warning: Failed to process string: {str(e2)}", file=sys.stderr)
                        continue
        else:
            self.add_batch_fallback(strings)
        
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
        
    def _convert_to_rusthll(self, other: Union['RustHLL', 'HyperLogLog']) -> 'RustHLL':
        """Convert a HyperLogLog sketch to RustHLL."""
        from hammock.lib.hyperloglog import HyperLogLog
        if isinstance(other, RustHLL):
            return other
        elif isinstance(other, HyperLogLog):
            # Create a new RustHLL with same parameters
            rust_hll = RustHLL(precision=other.precision, seed=other.seed)
            # Add all strings from the HyperLogLog sketch
            for i in range(len(other.registers)):
                if other.registers[i] > 0:
                    rust_hll.add_int(i)
            return rust_hll
        else:
            raise ValueError("Can only compare with RustHLL or HyperLogLog")

    def estimate_intersection(self, other: Union['RustHLL', 'HyperLogLog']) -> float:
        """Estimate intersection cardinality with another sketch."""
        from hammock.lib.hyperloglog import HyperLogLog
        if isinstance(other, HyperLogLog):
            # Use inclusion-exclusion principle
            union = self.estimate_union(other)
            self_card = self.estimate_cardinality()
            other_card = other.estimate_cardinality()
            intersection = self_card + other_card - union
            return max(0.0, intersection)
        elif isinstance(other, RustHLL):
            return self._estimate_intersection(other)
        else:
            raise ValueError("Can only compare with RustHLL or HyperLogLog")

    def estimate_union(self, other: Union['RustHLL', 'HyperLogLog']) -> float:
        """Estimate union cardinality with another sketch."""
        from hammock.lib.hyperloglog import HyperLogLog
        if isinstance(other, HyperLogLog):
            # Use inclusion-exclusion principle
            self_card = self.estimate_cardinality()
            other_card = other.estimate_cardinality()
            intersection = self.estimate_intersection(other)
            return self_card + other_card - intersection
        elif isinstance(other, RustHLL):
            return self._estimate_union(other)
        else:
            raise ValueError("Can only compare with RustHLL or HyperLogLog")

    def _correct_jaccard_estimate(self, raw_jaccard: float) -> float:
        """Apply an empirical correction to Jaccard estimates to account for systematic bias.
        
        Different correction strategies are applied based on the input value range:
        1. Very low values (0-0.2): Scale up slightly but keep low
        2. Mid-range values (0.2-0.4): These are typically reasonable estimates, adjust slightly
        3. Higher values (0.4+): These are usually too high for very different set sizes, 
           but may be reasonable for similar-sized sets
        
        Returns:
            float: Corrected Jaccard estimate
        """
        # Cap extreme values
        raw_jaccard = max(0.0, min(1.0, raw_jaccard))
        
        # For values close to 0, apply minor correction
        if raw_jaccard < 0.1:
            return raw_jaccard * 0.8  # Keep very low values low
        elif raw_jaccard < 0.2:
            # Linear mapping from 0.1-0.2 to 0.08-0.18
            return 0.08 + (raw_jaccard - 0.1) * 1.0
        elif raw_jaccard < 0.3:
            # Linear mapping from 0.2-0.3 to 0.18-0.3
            return 0.18 + (raw_jaccard - 0.2) * 1.2
        elif raw_jaccard < 0.4:
            # Linear mapping from 0.3-0.4 to 0.3-0.45
            return 0.3 + (raw_jaccard - 0.3) * 1.5
        elif raw_jaccard < 0.6:
            # Linear mapping from 0.4-0.6 to 0.45-0.7
            return 0.45 + (raw_jaccard - 0.4) * 1.25
        else:
            # Linear mapping from 0.6-1.0 to 0.7-1.0
            return 0.7 + (raw_jaccard - 0.6) * 0.75

    def estimate_jaccard(self, other: Union['RustHLL', 'HyperLogLog']) -> float:
        """Estimate Jaccard similarity with another sketch."""
        from hammock.lib.hyperloglog import HyperLogLog
        if not isinstance(other, (RustHLL, HyperLogLog)):
            raise ValueError("Can only compare with RustHLL or HyperLogLog")

        # For Python HyperLogLog, use the regular intersection/union approach
        if isinstance(other, HyperLogLog):
            intersection = self.estimate_intersection(other)
            union = self.estimate_union(other)
            
            if union == 0:
                return 0.0
            return intersection / union
        
        # For Rust HyperLogLog, use the improved method with correction
        if self._using_rust and hasattr(other, '_rust_sketch'):
            # Get the raw Jaccard estimate using the improved method
            raw_jaccard = float(self._rust_sketch.jaccard_improved(other._rust_sketch))
            
            # Apply correction only if enabled
            if self.apply_jaccard_correction:
                return self._correct_jaccard_estimate(raw_jaccard)
            else:
                return raw_jaccard  # Return raw value without correction
        
        # Fallback to regular method
        intersection = self.estimate_intersection(other)
        union = self.estimate_union(other)
        
        if union == 0:
            return 0.0
        return intersection / union
        
    def is_empty(self) -> bool:
        """Check if the sketch is empty."""
        return self._sketch.is_empty()
        
    def merge(self, other: 'RustHyperLogLog') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only merge with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        if self._using_rust:
            if hasattr(other, '_rust_sketch'):
                self._rust_sketch.merge(other._rust_sketch)
            else:
                raise ValueError("Cannot merge with incompatible sketch type")
        else:
            return self.merge_fallback(other)

    def merge_fallback(self, other: 'RustHyperLogLog') -> None:
        """Merge another sketch into this one using Python implementation."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only merge with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        # Take element-wise maximum
        np.maximum(self.registers, other.registers, out=self.registers)

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
    
    def jaccard(self, other: 'RustHLL') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, RustHLL):
            raise ValueError("Can only compare with another RustHLL")
        if self.precision != other.precision:
            raise ValueError("Cannot compute Jaccard similarity of HLLs with different precision")
        
        # Use the same method as estimate_jaccard for consistency
        return self.estimate_jaccard(other)
    
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

    def estimate_union_fallback(self, other) -> float:
        """Estimate union cardinality with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compute union of HLLs with different precision")
        
        # Create a temporary sketch for the union
        temp_sketch = RustHyperLogLog(self.precision)
        temp_sketch.merge(self)
        temp_sketch.merge(other)
        return temp_sketch.estimate_cardinality()
    
    def estimate_intersection_fallback(self, other) -> float:
        """Estimate intersection cardinality with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compute intersection of HLLs with different precision")
        
        # Use inclusion-exclusion principle
        union = self.estimate_union(other)
        card1 = self.estimate_cardinality()
        card2 = other.estimate_cardinality()
        intersection = card1 + card2 - union
        return max(0.0, intersection)  # Ensure non-negative
    
    def estimate_jaccard_fallback(self, other) -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compute Jaccard similarity of HLLs with different precision")
        
        intersection = self.estimate_intersection(other)
        union = self.estimate_union(other)
        if union == 0:
            return 1.0 if intersection == 0 else 0.0
        return intersection / union
    
    def write_fallback(self, filepath: str) -> None:
        """Write sketch to file."""
        # Convert to absolute path and ensure directory exists
        filepath = os.path.abspath(filepath)
        dir_path = os.path.dirname(filepath)
        if dir_path:
            os.makedirs(dir_path, exist_ok=True)
        try:
            np.save(filepath, self.registers)
        except Exception as e:
            print(f"Error saving sketch to {filepath}: {e}", file=sys.stderr)
            raise

class RustHLL(RustHyperLogLog):
    """Rust implementation of HyperLogLog."""
    
    def __init__(self, precision: Optional[int] = None, kmer_size: Optional[int] = None,
                 window_size: Optional[int] = None, hash_size: int = 32,
                 expected_cardinality: Optional[int] = None,
                 apply_jaccard_correction: bool = True,
                 debug: bool = False):
        """Initialize the sketch.
        
        Args:
            precision: Number of bits to use for register addressing.
            kmer_size: Size of k-mers (0 for whole string mode).
            window_size: Size of sliding window (0 or == kmer_size for no windowing).
            hash_size: Size of hash in bits (32 or 64).
            expected_cardinality: Expected number of unique elements.
            apply_jaccard_correction: Whether to apply empirical correction to Jaccard values.
            debug: Whether to enable debug mode.
        """
        super().__init__(
            precision=precision,
            kmer_size=kmer_size,
            window_size=window_size,
            hash_size=hash_size,
            expected_cardinality=expected_cardinality,
            apply_jaccard_correction=apply_jaccard_correction,
            debug=debug
        )

    def merge(self, other: 'RustHLL') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHLL):
            raise ValueError("Can only merge with another RustHLL")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        if self._using_rust:
            if hasattr(other, '_rust_sketch'):
                self._rust_sketch.merge(other._rust_sketch)
            else:
                raise ValueError("Cannot merge with incompatible sketch type")
        else:
            return self.merge_fallback(other)

    def estimate_union(self, other: 'RustHLL') -> float:
        """Estimate union cardinality with another sketch."""
        if not isinstance(other, RustHLL):
            raise ValueError("Can only compare with another RustHLL")
        if self.precision != other.precision:
            raise ValueError("Cannot compute union of HLLs with different precision")
        
        if self._using_rust:
            if hasattr(other, '_rust_sketch'):
                # Create a temporary sketch for the union
                temp_sketch = rust_hll.RustHLL(self.precision)
                temp_sketch.merge(self._rust_sketch)
                temp_sketch.merge(other._rust_sketch)
                return float(temp_sketch.estimate_cardinality())
            else:
                raise ValueError("Cannot compute union with incompatible sketch type")
        else:
            return self.estimate_union_fallback(other)

    def estimate_intersection(self, other: 'RustHLL') -> float:
        """Estimate intersection cardinality with another sketch."""
        if not isinstance(other, RustHLL):
            raise ValueError("Can only compare with another RustHLL")
        if self.precision != other.precision:
            raise ValueError("Cannot compute intersection of HLLs with different precision")
        
        if self._using_rust:
            if hasattr(other, '_rust_sketch'):
                # Get individual cardinalities
                card1 = self.estimate_cardinality()
                card2 = other.estimate_cardinality()
                
                # Create a temporary sketch for the union
                temp_sketch = rust_hll.RustHLL(self.precision)
                temp_sketch.merge(self._rust_sketch)
                temp_sketch.merge(other._rust_sketch)
                union = temp_sketch.estimate_cardinality()
                
                # Use inclusion-exclusion principle
                intersection = card1 + card2 - union
                return max(0.0, intersection)  # Ensure non-negative
            else:
                raise ValueError("Cannot compute intersection with incompatible sketch type")
        else:
            return self.estimate_intersection_fallback(other)

    def estimate_jaccard(self, other: 'RustHLL') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, RustHLL):
            raise ValueError("Can only compare with another RustHLL")
        if self.precision != other.precision:
            raise ValueError("Cannot compute Jaccard similarity of HLLs with different precision")
        
        if self._using_rust:
            if hasattr(other, '_rust_sketch'):
                # Get the raw Jaccard estimate using the improved method
                raw_jaccard = float(self._rust_sketch.jaccard_improved(other._rust_sketch))
                # Apply correction based on empirical testing
                return self._correct_jaccard_estimate(raw_jaccard)
            else:
                raise ValueError("Cannot compute Jaccard similarity with incompatible sketch type")
        else:
            return self.estimate_jaccard_fallback(other)

# Alias for backward compatibility
RustHLL = RustHLL 