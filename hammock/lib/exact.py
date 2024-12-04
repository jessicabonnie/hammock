from typing import Optional, Set
import xxhash

class ExactCounter:
    def __init__(self,
                 precision: int = 14,  # Unused but kept for interface compatibility
                 kmer_size: int = 0,
                 window_size: int = 0,
                 seed: int = 0,
                 debug: bool = False):
        """Initialize exact counter.
        
        Args match HyperLogLog for compatibility, but only kmer_size,
        window_size, and seed are used.
        """
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.seed = seed
        self.elements: Set[str] = set()
        self.debug = debug
    def add_string(self, s: str) -> None:
        """Add a string to the counter."""
        if len(s) < self.kmer_size:
            return

        if self.kmer_size == 0:
            self.elements.add(s)
            return

        if self.window_size == self.kmer_size:
            for i in range(len(s) - self.kmer_size + 1):
                self.elements.add(s[i:i + self.kmer_size])
            return

        # Windowing mode
        for i in range(len(s) - self.window_size + 1):
            window = s[i:i + self.window_size]
            min_hash = float('inf')
            min_kmer = ''
            
            for j in range(self.window_size - self.kmer_size + 1):
                kmer = window[j:j + self.kmer_size]
                hasher = xxhash.xxh64(seed=self.seed)
                hasher.update(kmer.encode())
                h = hasher.intdigest()
                if h < min_hash:
                    min_hash = h
                    min_kmer = kmer
            
            self.elements.add(min_kmer)

    def merge(self, other: 'ExactCounter') -> None:
        """Merge another counter into this one."""
        self.elements.update(other.elements)

    def cardinality(self) -> float:
        """Return exact cardinality."""
        return float(len(self.elements))

    def intersection(self, other: 'ExactCounter') -> float:
        """Return exact intersection cardinality."""
        return float(len(self.elements & other.elements))

    def union(self, other: 'ExactCounter') -> float:
        """Return exact union cardinality."""
        return float(len(self.elements | other.elements))

    def jaccard(self, other: 'ExactCounter') -> float:
        """Return exact Jaccard similarity.
        
        Returns:
            Float between 0 and 1 representing Jaccard similarity
        """
        intersection = len(self.elements & other.elements)
        union = len(self.elements | other.elements)
        return float(intersection / union) if union > 0 else 0.0

    # Keep estimate_* methods for interface compatibility
    def estimate_cardinality(self) -> float:
        """Alias for cardinality()."""
        return self.cardinality()

    def estimate_intersection(self, other: 'ExactCounter') -> float:
        """Alias for intersection()."""
        return self.intersection(other)

    def estimate_union(self, other: 'ExactCounter') -> float:
        """Alias for union()."""
        return self.union(other)

    def estimate_jaccard(self, other: 'ExactCounter') -> float:
        """Alias for jaccard()."""
        return self.jaccard(other)