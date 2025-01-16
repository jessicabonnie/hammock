from __future__ import annotations
import numpy as np # type: ignore
import xxhash # type: ignore
from typing import Optional, List, Union
from hammock.lib.abstractsketch import AbstractSketch

class MinHash(AbstractSketch):
    def __init__(self, 
                 num_hashes: int = 128,
                 kmer_size: int = 0,
                 window_size: int = 0,
                 seed: int = 0,
                 debug: bool = False):
        """Initialize MinHash sketch.
        
        Args:
            num_hashes: Number of hash functions
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            debug: Whether to print debug information
        """
        self.num_hashes = num_hashes
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.seed = seed
        self.signatures = np.full(num_hashes, np.inf)
        
        # Create xxhash objects for each hash function
        self.hashers = [xxhash.xxh64(seed=seed + i) for i in range(num_hashes)]
    
    def _hash_str(self, s: str) -> np.ndarray:
        """Generate all hash values for a string."""
        hashes = np.zeros(self.num_hashes)
        for i, hasher in enumerate(self.hashers):
            hasher.reset()  # Reset the hasher state
            hasher.update(s.encode())  # Hash the string
            hashes[i] = hasher.intdigest()  # Don't normalize to [0,1]
        return hashes
    
    def _process_kmer(self, kmer: str) -> None:
        """Process a single k-mer or string."""
        hashes = self._hash_str(kmer)
        self.signatures = np.minimum(self.signatures, hashes)

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        if len(s) < self.kmer_size:
            return

        if self.kmer_size == 0:
            self._process_kmer(s)
            return

        if self.window_size == self.kmer_size:
            for i in range(len(s) - self.kmer_size + 1):
                self._process_kmer(s[i:i + self.kmer_size])
            return

        # Windowing mode
        for i in range(len(s) - self.window_size + 1):
            window = s[i:i + self.window_size]
            min_hash = float('inf')
            min_kmer = ''
            
            for j in range(self.window_size - self.kmer_size + 1):
                kmer = window[j:j + self.kmer_size]
                h = hash(kmer)
                if h < min_hash:
                    min_hash = h
                    min_kmer = kmer
            
            self._process_kmer(min_kmer)

    def merge(self, other: 'MinHash') -> None:
        """Merge another MinHash sketch into this one."""
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot merge MinHash sketches with different sizes")
        self.signatures = np.minimum(self.signatures, other.signatures)

    def estimate_cardinality(self) -> float:
        """Estimate cardinality of the set using MinHash signatures."""
        mask = ~np.isinf(self.signatures)
        if not np.any(mask):
            return 0.0
            
        k = max(1, int(0.7 * self.num_hashes))
        kth_min = np.partition(self.signatures[mask], k)[k]
        
        return (k / kth_min) * (2**64) if kth_min > 0 else 0.0

    def estimate_intersection(self, other: 'MinHash') -> float:
        """Estimate intersection cardinality with another MinHash sketch."""
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot compute intersection of MinHash sketches with different sizes")
        
        jaccard = self.estimate_jaccard(other)
        union = self.estimate_union(other)
        
        return max(0.0, jaccard * union)

    def estimate_jaccard(self, other: 'MinHash') -> float:
        """Estimate Jaccard similarity with another MinHash sketch."""
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot compare MinHash sketches with different sizes")
        
        mask = (~np.isinf(self.signatures)) & (~np.isinf(other.signatures))
        if not np.any(mask):
            return 0.0
            
        matches = np.sum((self.signatures == other.signatures) & mask)
        total = np.sum(mask)
        
        return matches / total if total > 0 else 0.0

    def estimate_union(self, other: 'MinHash') -> float:
        """Estimate union cardinality with another MinHash sketch."""
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot compute union of MinHash sketches with different sizes")
        
        combined = np.minimum(self.signatures, other.signatures)
        mask = ~np.isinf(combined)
        if not np.any(mask):
            return 0.0
            
        k = max(1, int(0.7 * self.num_hashes))
        kth_min = np.partition(combined[mask], k)[k]
        
        return (k / kth_min) * (2**64) if kth_min > 0 else 0.0

    def write(self, filepath: str) -> None:
        """Write sketch to file in binary format."""
        np.savez_compressed(
            filepath,
            signatures=self.signatures,
            num_hashes=np.array([self.num_hashes]),
            kmer_size=np.array([self.kmer_size]),
            window_size=np.array([self.window_size]),
            seed=np.array([self.seed])
        )

    @classmethod
    def load(cls, filepath: str) -> 'MinHash':
        """Load sketch from file in binary format."""
        data = np.load(filepath)
        sketch = cls(
            num_hashes=int(data['num_hashes'][0]),
            kmer_size=int(data['kmer_size'][0]),
            window_size=int(data['window_size'][0]),
            seed=int(data['seed'][0])
        )
        sketch.signatures = data['signatures']
        return sketch

    def add_int(self, value: int) -> None:
        hashes = np.array([h.intdigest() for h in self.hashers])
        self.signatures = np.minimum(self.signatures, hashes)
