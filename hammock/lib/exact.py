import json
from typing import Set, Optional
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
        
    def estimate_jaccard(self, other: 'ExactCounter') -> float:
        """Calculate exact Jaccard similarity."""
        if not isinstance(other, ExactCounter):
            raise TypeError("Can only compare with another ExactCounter")
        intersection = len(self.elements & other.elements)
        union = len(self.elements | other.elements)
        return intersection / union if union > 0 else 0.0
        
    def merge(self, other: 'ExactCounter') -> None:
        """Merge another counter into this one."""
        if not isinstance(other, ExactCounter):
            raise TypeError("Can only merge with another ExactCounter")
        self.elements.update(other.elements)

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