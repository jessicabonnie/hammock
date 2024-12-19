from abc import ABC, abstractmethod
from typing import Union, Optional

class AbstractSketch(ABC):
    """Abstract base class defining core sketching operations."""
    
    @abstractmethod
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
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
    def merge(self, other: 'AbstractSketch') -> None:
        """Merge another sketch into this one."""
        pass

class AbstractDataSketch(ABC):
    """Abstract base class for data-specific sketch wrappers."""
    
    @abstractmethod
    def __init__(self, sketch_type: str, **kwargs):
        """Initialize with specified sketch type."""
        pass
    
    @classmethod
    @abstractmethod
    def from_file(cls, filename: str, sketch_type: str, **kwargs) -> Optional['AbstractDataSketch']:
        """Create sketch from file."""
        pass 