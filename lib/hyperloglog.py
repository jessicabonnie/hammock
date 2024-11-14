import math
from typing import Optional, List, Union
import numpy as np
import xxhash

class HyperLogLog:
    def __init__(self, 
                 precision: int = 14, 
                 kmer_size: int = 0, 
                 window_size: int = 0, 
                 seed: int = 0):
        """Initialize HyperLogLog sketch.
        
        Args:
            precision: Number of bits for register indexing (4-16)
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
        """
        if not 4 <= precision <= 16:
            raise ValueError("Precision must be between 4 and 16")
        
        if window_size and window_size < kmer_size:
            raise ValueError("Window size must be >= kmer size")
        
        self.precision = precision
        self.num_registers = 1 << precision
        self.seed = seed
        self.registers = np.zeros(self.num_registers, dtype=np.uint8)
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        
        # Calculate alpha_mm (bias correction factor)
        if self.num_registers == 16:
            self.alpha_mm = 0.673
        elif self.num_registers == 32:
            self.alpha_mm = 0.697
        elif self.num_registers == 64:
            self.alpha_mm = 0.709
        else:
            self.alpha_mm = 0.7213 / (1 + 1.079 / self.num_registers)

    def _hash64(self, x: int) -> int:
        """64-bit hash function."""
        hasher = xxhash.xxh64(seed=self.seed)
        hasher.update(x.to_bytes(8, byteorder='little'))
        return hasher.intdigest()
    # def _hash64(self, x: int) -> int:
        # """64-bit hash function."""
        # x = (x ^ self.seed) & 0xFFFFFFFFFFFFFFFF
        # x ^= (x >> 33) & 0xFFFFFFFFFFFFFFFF
        # x = (x * 0xff51afd7ed558ccd) & 0xFFFFFFFFFFFFFFFF
        # x ^= (x >> 33) & 0xFFFFFFFFFFFFFFFF
        # x = (x * 0xc4ceb9fe1a85ec53) & 0xFFFFFFFFFFFFFFFF
        # x ^= (x >> 33) & 0xFFFFFFFFFFFFFFFF
        # return x

    @staticmethod
    def _hash_str(s: bytes, seed: int = 0) -> int:
        """Hash a string using xxhash."""
        hasher = xxhash.xxh64(seed=seed)
        hasher.update(s)
        return hasher.intdigest()

    def _rho(self, hash_val: int) -> int:
        """Calculate position of leftmost 1-bit."""
        hash_val >>= self.precision
        pos = 1
        while (hash_val & 1) == 0 and pos <= 64:
            pos += 1
            hash_val >>= 1
        return pos

    def _process_kmer(self, kmer: str) -> None:
        """Process a single k-mer or string."""
        hash_val = self._hash_str(kmer.encode())
        idx = hash_val & (self.num_registers - 1)
        rank = self._rho(hash_val)
        self.registers[idx] = max(self.registers[idx], rank)

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        if len(s) < self.kmer_size:
            return

        # Whole string mode
        if self.kmer_size == 0:
            self._process_kmer(s)
            return

        # K-mer mode without windowing
        if self.window_size == self.kmer_size:
            for i in range(len(s) - self.kmer_size + 1):
                self._process_kmer(s[i:i + self.kmer_size])
            return

        # Windowing mode
        for i in range(len(s) - self.window_size + 1):
            window = s[i:i + self.window_size]
            min_hash = float('inf')
            min_pos = 0
            
            # Find minimum hash in this window
            for j in range(self.window_size - self.kmer_size + 1):
                kmer = window[j:j + self.kmer_size]
                h = self._hash_str(kmer.encode())
                if h < min_hash:
                    min_hash = h
                    min_pos = j
            
            # Process the minimizer
            self._process_kmer(window[min_pos:min_pos + self.kmer_size])

    def add_int(self, value: int) -> None:
        """Add an integer to the sketch."""
        hash_val = self._hash64(value)
        idx = hash_val & (self.num_registers - 1)
        rank = self._rho(hash_val)
        self.registers[idx] = max(self.registers[idx], rank)

    def estimate_cardinality(self) -> float:
        """Estimate cardinality of the set."""
        sum_inv = np.sum(2.0 ** -self.registers)
        num_zeros = np.sum(self.registers == 0)
        
        # Initial estimate in 64-bit space
        estimate = self.alpha_mm * (2**64) * self.num_registers ** 2 / sum_inv

        # Small range correction
        if estimate <= 2.5 * self.num_registers:
            if num_zeros > 0:
                estimate = self.num_registers * math.log(self.num_registers / num_zeros)
                # Scale small range correction to 64-bit space
                estimate = estimate * (2**64) / self.num_registers
        # Large range correction
        elif estimate > 2**64 / 2:
            estimate = -2**64 * math.log(1 - estimate / 2**64)

        return estimate

    def merge(self, other: 'HyperLogLog') -> None:
        """Merge another HLL sketch into this one."""
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        self.registers = np.maximum(self.registers, other.registers)

    def estimate_intersection(self, other: 'HyperLogLog') -> float:
        """Estimate intersection cardinality with another HLL."""
        if self.precision != other.precision:
            raise ValueError("Cannot compute intersection of HLLs with different precision")

        a = self.estimate_cardinality()
        b = other.estimate_cardinality()
        c = self.estimate_union(other)

        return max(0.0, a + b - c)

    def estimate_jaccard(self, other: 'HyperLogLog') -> float:
        """Estimate Jaccard similarity with another HyperLogLog sketch.
        
        Returns:
            Float between 0 and 1 representing Jaccard similarity
        """
        if self.precision != other.precision:
            raise ValueError("Cannot compare HLLs with different precision")
            
        intersection = self.estimate_intersection(other)
        union = self.estimate_cardinality() + other.estimate_cardinality() - intersection
        
        return intersection / union if union > 0 else 0.0

    def estimate_cardinality(self) -> float:
        """Estimate cardinality of the set."""
        sum_inv = np.sum(2.0 ** -self.registers)
        num_zeros = np.sum(self.registers == 0)
        
        # Initial estimate in 64-bit space
        estimate = self.alpha_mm * (2**64) * self.num_registers ** 2 / sum_inv

        # Small range correction
        if estimate <= 2.5 * self.num_registers:
            if num_zeros > 0:
                estimate = self.num_registers * math.log(self.num_registers / num_zeros)
                # Scale small range correction to 64-bit space
                estimate = estimate * (2**64) / self.num_registers
        # Large range correction
        elif estimate > 2**64 / 2:
            estimate = -2**64 * math.log(1 - estimate / 2**64)

        return estimate

    def estimate_union(self, other: 'HyperLogLog') -> float:
        """Estimate union cardinality with another HyperLogLog sketch."""
        if self.precision != other.precision:
            raise ValueError("Cannot compute union of HLLs with different precision")
            
        merged = HyperLogLog(self.precision, self.kmer_size, self.window_size, self.seed)
        merged.registers = np.maximum(self.registers, other.registers)
        
        # Use the same scaling as cardinality estimation
        sum_inv = np.sum(2.0 ** -merged.registers)
        num_zeros = np.sum(merged.registers == 0)
        
        estimate = merged.alpha_mm * (2**64) * merged.num_registers ** 2 / sum_inv

        if estimate <= 2.5 * merged.num_registers:
            if num_zeros > 0:
                estimate = merged.num_registers * math.log(merged.num_registers / num_zeros)
                estimate = estimate * (2**64) / merged.num_registers
        elif estimate > 2**64 / 2:
            estimate = -2**64 * math.log(1 - estimate / 2**64)

        return estimate

    def estimate_intersection(self, other: 'HyperLogLog') -> float:
        """Estimate intersection cardinality with another HLL."""
        if self.precision != other.precision:
            raise ValueError("Cannot compute intersection of HLLs with different precision")

        a = self.estimate_cardinality()
        b = other.estimate_cardinality()
        c = self.estimate_union(other)

        return max(0.0, a + b - c)
