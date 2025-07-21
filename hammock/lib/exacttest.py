from __future__ import annotations
from typing import Optional, Set, Tuple, Dict, Any
from .abstractsketch import AbstractSketch

class ExactTest(AbstractSketch):
    """Exact counter for testing purposes only."""
    
    def __init__(self, seed: int = 0):
        """Initialize exact counter."""
        self.elements: Set[str] = set()
        self.seed = seed
    
    def add_string(self, s: str) -> None:
        """Add string to counter."""
        self.elements.add(s)
    
    def add_batch(self, strings: list[str]) -> None:
        """Add multiple strings to counter."""
        for s in strings:
            self.add_string(s)
    
    def estimate_cardinality(self) -> float:
        """Return exact cardinality."""
        return len(self.elements)
    
    def similarity_values(self, other: 'AbstractSketch') -> dict[str, float]:
        """Calculate similarity values.
        
        Returns:
            Dictionary containing 'jaccard_similarity'
        """
        if not isinstance(other, ExactTest):
            raise ValueError("Can only compare with another ExactTest")
        return {'jaccard_similarity': self.estimate_jaccard(other)}
    
    def estimate_jaccard(self, other: 'ExactTest') -> float:
        """Return exact Jaccard similarity."""
        intersection = len(self.elements & other.elements)
        union = len(self.elements | other.elements)
        return intersection / union if union > 0 else 0.0
    
    # _hash_str method now inherited from AbstractSketch base class
    
    def merge(self, other: 'ExactTest') -> None:
        """Merge another ExactTest into this one."""
        self.elements.update(other.elements)
    
    def write(self, path: str) -> None:
        """Write sketch to file."""
        with open(path, 'w') as f:
            for element in sorted(self.elements):
                f.write(element + '\n')
    
    @classmethod
    def load(cls, path: str) -> 'ExactTest':
        """Load sketch from file."""
        sketch = cls()
        with open(path, 'r') as f:
            for line in f:
                sketch.elements.add(line.strip())
        return sketch 