"""
This module contains the pure Python implementation of HyperLogLog that serves as a fallback
when the Rust implementation is not available.

This file is not meant to be imported directly by users. Instead, users should import
from rust_hll_ext.py, which will automatically use this implementation if the Rust
version is not available.

The implementation here is a direct translation of the Rust algorithm to Python,
using numpy for efficient array operations. It provides the same interface and
functionality as the Rust version, but may be slower for large datasets.
"""

from __future__ import annotations
from typing import Optional, List, Tuple, TextIO, Union, Dict
import numpy as np  # type: ignore
import xxhash  # type: ignore

class RustHLLExtRS:
    """Rust implementation of HyperLogLog with fallback to pure Python."""
    
    def __init__(self, precision: int = 12, hash_size: int = 32):
        """Initialize the sketch.
        
        Args:
            precision: Number of bits to use for register addressing.
            hash_size: Size of hash in bits (32 or 64)
        """
        self.precision = precision
        self.num_registers = 1 << precision
        self.hash_size = hash_size
        
        if RUST_AVAILABLE:
            self._sketch = RustHLLClass(precision, hash_size=hash_size)
            self._using_rust = True
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
        return np.all(self.registers == 0)
        
    def merge(self, other: 'RustHLLExtRS') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHLLExtRS):
            raise ValueError("Can only merge with another RustHLLExtRS")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        self.registers = np.maximum(self.registers, other.registers)
        
    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        np.save(filepath, self.registers)
        
    @classmethod
    def load(cls, filepath: str) -> 'RustHLLExtRS':
        """Load sketch from file."""
        registers = np.load(filepath)
        sketch = cls(precision=int(np.log2(len(registers))))
        sketch.registers = registers
        return sketch 