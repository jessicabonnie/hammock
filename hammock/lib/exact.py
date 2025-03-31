from __future__ import annotations
# import json
from typing import Set, Optional, Union, List, Dict, Tuple
from hammock.lib.abstractsketch import AbstractSketch

class ExactCounter(AbstractSketch):
    """Exact counter for k-mers."""
    
    def __init__(self, kmer_size: int = 8, seed: int = 42):
        """Initialize exact counter."""
        super().__init__()
        self.kmer_size = kmer_size
        self.seed = seed
        self.elements: Set[str] = set()
        
    def add_string(self, s: str) -> None:
        """Add a string to the counter."""
        self.elements.add(s)
        
    def estimate_cardinality(self) -> float:
        """Return exact cardinality."""
        return float(len(self.elements))
        
    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Compute exact Jaccard similarity.
        
        Returns:
            Dictionary containing 'jaccard_similarity'
        """
        if not isinstance(other, ExactCounter):
            raise ValueError("Can only compare with another ExactCounter")
            
        if not self.elements and not other.elements:
            return {'jaccard_similarity': 0.0}
            
        intersection = len(self.elements & other.elements)
        union = len(self.elements | other.elements)
        
        jaccard = intersection / union if union > 0 else 0.0
        return {'jaccard_similarity': jaccard}
        
    def merge(self, other: 'ExactCounter') -> None:
        """Merge another counter into this one."""
        if not isinstance(other, ExactCounter):
            raise TypeError("Can only merge with another ExactCounter")
        self.elements.update(other.elements)

    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the counter.
        
        Args:
            strings: List of strings to add to the counter
        """
        self.elements.update(strings)

    def write(self, filepath: str) -> None:
        """Write operation not supported for ExactCounter."""
        raise NotImplementedError("ExactCounter does not support writing to file")
        # """Write sketch to file."""
        # with open(filepath, 'w') as f:
        #     for element in self.elements:
        #         f.write(f"{element}\n")

    @classmethod
    def load(cls, filepath: str) -> 'ExactCounter':
        """Load operation not supported for ExactCounter."""
        raise NotImplementedError("ExactCounter does not support loading from file")
        # """Load sketch from file."""
        # sketch = cls()
        # with open(filepath, 'r') as f:
        #     for line in f:
        #         sketch.add_string(line.strip())
        # return sketch