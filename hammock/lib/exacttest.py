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
    
    def estimate_cardinality(self) -> float:
        """Return exact cardinality."""
        return len(self.elements)
    
    def estimate_jaccard(self, other: 'ExactTest') -> float:
        """Return exact Jaccard similarity."""
        intersection = len(self.elements & other.elements)
        union = len(self.elements | other.elements)
        return intersection / union if union > 0 else 0.0
    
    def _hash_str(self, s: bytes, seed: Optional[int] = None) -> int:
        """Hash bytes using Python's built-in hash function."""
        if seed is None:
            seed = self.seed
        return hash(s + str(seed).encode('utf-8')) % (2**32)
    
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