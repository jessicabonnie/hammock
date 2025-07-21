"""
This module provides a high-level Python interface to the HyperLogLog implementation.

The module structure is:
1. rust_hll.py (this file) - High-level Python interface
2. rust_hll_ext.py - Bridge between Python and Rust implementations
3. rust_hll_ext_rs.py - Pure Python fallback implementation

This file provides a more Pythonic interface to the HyperLogLog implementation,
with additional features like:
- Support for different estimation methods
- Union and intersection operations
- Jaccard similarity estimation
- File I/O operations

The actual implementation is delegated to the Rust or Python version through
rust_hll_ext.py, which handles the implementation selection automatically.
"""

from __future__ import annotations
from typing import Optional, List, Tuple, TextIO, Union, Dict
import numpy as np  # type: ignore
import xxhash  # type: ignore
import os
import warnings

# Try to import the Rust extension
try:
    from hammock.lib.rust_hll_ext import RustHLLWrapper as RustHLLClass
    RUST_HLL_AVAILABLE = True
except ImportError:
    RUST_HLL_AVAILABLE = False

class RustHLL:
    """Rust implementation of HyperLogLog."""
    
    def __init__(self, precision: int = 12, hash_size: int = 32, memory_limit: Optional[int] = None):
        """Initialize the sketch.
        
        Args:
            precision: Number of bits to use for register addressing.
            hash_size: Size of hash in bits (32 or 64)
            memory_limit: Maximum memory usage in bytes (optional)
        """
        self.precision = precision
        self.num_registers = 1 << precision
        self.hash_size = hash_size
        
        if RUST_HLL_AVAILABLE:
            try:
                # Convert memory limit to bytes if provided in GB
                if memory_limit is not None and memory_limit < 1024:
                    memory_limit = memory_limit * 1024 * 1024 * 1024
                
                self._sketch = RustHLLClass(precision, hash_size=hash_size, memory_limit=memory_limit)
                self._using_rust = True
            except ValueError as e:
                print(f"Warning: Failed to initialize Rust HLL: {e}")
                print("Falling back to Python implementation")
                self._using_rust = False
                self.registers = np.zeros(self.num_registers, dtype=np.float64)
        else:
            self.registers = np.zeros(self.num_registers, dtype=np.float64)
            self._using_rust = False
        
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.add(s)
        
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch."""
        if self._using_rust:
            self._sketch.add_batch(strings)
        else:
            for s in strings:
                self.add(s)
        
    def add(self, value: str) -> None:
        """Add a value to the sketch."""
        if self._using_rust:
            self._sketch.add_value(value)
        else:
            if self.hash_size == 32:
                hash_val = xxhash.xxh32(value.encode()).intdigest()
            else:
                hash_val = xxhash.xxh64(value.encode()).intdigest()
            idx = hash_val & (self.num_registers - 1)
            rank = self._rho(hash_val)
            self.registers[idx] = max(self.registers[idx], rank)
        
    def _rho(self, value: int) -> int:
        """Calculate the position of the leftmost 1-bit."""
        if value == 0:
            return self.hash_size
        return 1 + (value & -value).bit_length() - 1
        
    def estimate_cardinality(self) -> float:
        """Estimate the cardinality of the set."""
        if self._using_rust:
            return self._sketch.estimate_cardinality()
            
        # Calculate raw estimate
        raw_estimate = self._raw_estimate()
        
        # Apply bias correction for small cardinalities
        if raw_estimate <= 5 * self.num_registers:
            raw_estimate = self._bias_correct(raw_estimate)
            
        return raw_estimate
        
    def _raw_estimate(self) -> float:
        """Calculate raw estimate of cardinality."""
        # Calculate harmonic mean of 2^(-M[j])
        sum_inv = np.sum(2.0 ** -self.registers)
        alpha = self._get_alpha()
        return alpha * self.num_registers * self.num_registers / sum_inv
        
    def _bias_correct(self, raw_estimate: float) -> float:
        """Apply bias correction for small cardinalities."""
        # Use linear interpolation between known bias values
        if raw_estimate <= 0:
            return 0.0
            
        # Get bias correction table
        bias_table = self._get_bias_table()
        
        # Find the closest cardinality in the table
        closest_card = min(bias_table.keys(), key=lambda x: abs(x - raw_estimate))
        bias = bias_table[closest_card]
        
        return raw_estimate * (1 - bias)
        
    def _get_alpha(self) -> float:
        """Get alpha constant based on precision."""
        if self.precision == 4:
            return 0.673
        elif self.precision == 5:
            return 0.697
        elif self.precision == 6:
            return 0.709
        else:
            return 0.7213 / (1 + 1.079 / self.num_registers)
            
    def _get_bias_table(self) -> Dict[float, float]:
        """Get bias correction table based on precision."""
        if self.precision == 4:
            return {
                0: 0.673,
                1: 0.697,
                2: 0.709,
                3: 0.715,
                4: 0.718,
                5: 0.719
            }
        elif self.precision == 5:
            return {
                0: 0.697,
                1: 0.709,
                2: 0.715,
                3: 0.718,
                4: 0.719,
                5: 0.720
            }
        elif self.precision == 6:
            return {
                0: 0.709,
                1: 0.715,
                2: 0.718,
                3: 0.719,
                4: 0.720,
                5: 0.721
            }
        else:
            return {
                0: 0.7213,
                1: 0.7213,
                2: 0.7213,
                3: 0.7213,
                4: 0.7213,
                5: 0.7213
            }
            
    def is_empty(self) -> bool:
        """Check if the sketch is empty."""
        if self._using_rust:
            return self._sketch.is_empty()
        return np.all(self.registers == 0)
        
    def merge(self, other: 'RustHLL') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHLL):
            raise ValueError("Can only merge with another RustHLL")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        if self._using_rust and other._using_rust:
            self._sketch.merge(other._sketch)
        elif not self._using_rust and not other._using_rust:
            self.registers = np.maximum(self.registers, other.registers)
        else:
            raise ValueError("Cannot merge sketches with different implementations")
        
    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        if self._using_rust:
            self._sketch.write(filepath)
        else:
            # Ensure directory exists
            dir_path = os.path.dirname(filepath)
            if dir_path:
                os.makedirs(dir_path, exist_ok=True)
            # Write registers to file
            np.save(filepath, self.registers)
        
    @classmethod
    def load(cls, filepath: str) -> 'RustHLL':
        """Load sketch from file."""
        if RUST_HLL_AVAILABLE:
            sketch = cls()
            sketch._sketch = RustHLLClass.load(filepath)
            sketch._using_rust = True
            return sketch
        else:
            registers = np.load(filepath)
            sketch = cls(precision=int(np.log2(len(registers))))
            sketch.registers = registers
            sketch._using_rust = False
            return sketch 

    def set_memory_limit(self, limit: Optional[int]) -> None:
        """Set memory usage limit in bytes."""
        if self._using_rust:
            if limit is not None and limit < 1024:
                limit = limit * 1024 * 1024 * 1024
            self._sketch.set_memory_limit(limit)
    
    def get_memory_usage(self) -> int:
        """Get current memory usage in bytes."""
        if self._using_rust:
            return self._sketch.get_memory_usage()
        else:
            return self.registers.nbytes 