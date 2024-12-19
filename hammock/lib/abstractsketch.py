from abc import ABC, abstractmethod
from typing import Union, Optional

class AbstractSketch(ABC):
    """Abstract base class defining the sketch interface."""
    
    @abstractmethod
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        pass
    
    @abstractmethod
    def add_int(self, value: int) -> None:
        """Add an integer to the sketch."""
        pass
    
    @abstractmethod
    def estimate_cardinality(self) -> float:
        """Estimate cardinality of the sketch."""
        pass
    
    @abstractmethod
    def estimate_jaccard(self, other: 'AbstractSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        pass
    
    @abstractmethod
    def _hash_str(self, s: Union[str, bytes], seed: int = 0) -> int:
        """Hash a string or bytes object."""
        pass
        
    @abstractmethod
    def write_sketch(self, filepath: str) -> None:
        """Write sketch to file."""
        pass
        
    @classmethod
    @abstractmethod
    def read_sketch(cls, filepath: str) -> 'AbstractSketch':
        """Read sketch from file."""
        pass 