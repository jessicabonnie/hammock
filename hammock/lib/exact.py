from typing import Optional, Set, Union
import xxhash # type: ignore
from hammock.lib.abstractsketch import AbstractSketch

class ExactCounter(AbstractSketch):
    """Exact set counter that implements the sketch interface."""
    
    def __init__(self,
                 kmer_size: int = 0,
                 seed: int = 0):
        """Initialize exact counter.
        
        Args:
            kmer_size: Size of kmers (0 for whole string mode)
            seed: Random seed for hashing
        """
        self.kmer_size = kmer_size
        self.seed = seed
        self.elements: Set[str] = set()

    def add_string(self, s: str) -> None:
        """Add a string to the counter."""
        if len(s) < self.kmer_size:
            return

        if self.kmer_size == 0:
            self.elements.add(s)
            return

        # Process k-mers
        for i in range(len(s) - self.kmer_size + 1):
            self.elements.add(s[i:i + self.kmer_size])

    def estimate_cardinality(self) -> float:
        """Return exact cardinality."""
        return float(len(self.elements))

    def estimate_jaccard(self, other: 'AbstractSketch') -> float:
        """Return exact Jaccard similarity."""
        if not isinstance(other, ExactCounter):
            raise ValueError("Can only compare with another ExactCounter")
        intersection = len(self.elements & other.elements)
        union = len(self.elements | other.elements)
        return float(intersection) / union if union > 0 else 0.0

    def merge(self, other: 'AbstractSketch') -> None:
        """Merge another counter into this one."""
        if not isinstance(other, ExactCounter):
            raise ValueError("Can only merge with another ExactCounter")
        self.elements.update(other.elements)