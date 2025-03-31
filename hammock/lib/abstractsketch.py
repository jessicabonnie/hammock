from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional, Dict, List

class AbstractSketch(ABC):
    """Base class for all sketch types."""
    
    @abstractmethod
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        pass
    
    @abstractmethod
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        pass
    
    @abstractmethod
    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Estimate similarity with another sketch.
        
        Returns:
            Dictionary of similarity measures
        """
        pass
    
    # @abstractmethod
    # def merge(self, other: 'AbstractSketch') -> None:
    #     """Merge another sketch into this one."""
    #     pass
        
    @abstractmethod
    def write(self, filepath: str) -> None:
        """Write sketch to file.
        
        Args:
            filepath: Path to write the sketch to
        """
        pass
        
    @classmethod
    @abstractmethod
    def load(cls, filepath: str) -> 'AbstractSketch':
        """Load sketch from file.
        
        Args:
            filepath: Path to load the sketch from
            
        Returns:
            The loaded sketch
        """
        pass

# class AbstractDataSketch(ABC):
#     """Abstract base class for data-specific sketch wrappers."""
    
#     @abstractmethod
#     def __init__(self, sketch_type: str, **kwargs):
#         """Initialize with specified sketch type."""
#         pass
    
#     @classmethod
#     @abstractmethod
#     def from_file(cls, filename: str, sketch_type: str, **kwargs) -> Optional['AbstractDataSketch']:
#         """Create sketch from file."""
#         pass 