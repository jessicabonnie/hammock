#!/usr/bin/env python
from __future__ import annotations
import numpy as np # type: ignore
import xxhash # type: ignore
from hammock.lib.abstractsketch import AbstractSketch
from typing import List

class SetSketch(AbstractSketch):
    def __init__(self, 
                 precision: int = 8, 
                 kmer_size: int = 0, 
                 window_size: int = 0, 
                 seed: int = 0,
                 debug: bool = False):
        """Initialize SetSketch sketch.
        
        Args:
            precision: Number of registers (4 or higher)
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            debug: Whether to print debug information
        """
        if precision < 4:
            raise ValueError("Precision must be at least 4")
        
        if window_size and window_size < kmer_size:
            raise ValueError("Window size must be >= kmer size")
        
        self.precision = precision
        self.num_registers = 1 << precision
        self.seed = seed
        self.hash_size = 64  # SetSketch uses 64-bit hashing
        self.registers = np.zeros(self.num_registers, dtype=np.float64)
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.debug = debug
        
        self.item_count = 0

    # Hash functions now inherited from AbstractSketch base class

    def add_string(self, s: str) -> None:
        """Add a string to the sketch.
        
        Processes the string according to kmer_size and window_size settings.
        If kmer_size is 0, processes whole string.
        If window_size equals kmer_size, processes all k-mers.
        Otherwise uses minimizer scheme within windows.
        
        Args:
            s: String to add to sketch
        """
        if self.kmer_size == 0:
            # Whole string mode
            self.add_int(self.hash_str(s.encode()))
        else:
            # K-mer mode
            s_bytes = s.encode()
            for i in range(len(s_bytes) - self.kmer_size + 1):
                kmer = s_bytes[i:i + self.kmer_size]
                self.add_int(self.hash_str(kmer))

    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        for s in strings:
            self.add_string(s)

    def add_int(self, value: int) -> None:
        """Add an integer to the sketch."""
        raise NotImplementedError("SetSketch.add_int() not implemented")

    def estimate_cardinality(self) -> float:
        """Estimate the cardinality of the set."""
        raise NotImplementedError("SetSketch.estimate_cardinality() not implemented")

    def merge(self, other: 'SetSketch') -> None:
        """Merge another sketch into this one.
        
        Args:
            other: Another SetSketch to merge into this one
            
        Raises:
            ValueError: If sketches have different parameters
        """
        raise NotImplementedError("SetSketch.merge() not implemented")

    def estimate_intersection(self, other: 'SetSketch') -> float:
        """Estimate intersection cardinality with another sketch.
        
        Args:
            other: Another SetSketch
            
        Returns:
            Estimated size of intersection between sketches
        """
        raise NotImplementedError("SetSketch.estimate_intersection() not implemented")

    def estimate_jaccard(self, other: 'SetSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        raise NotImplementedError("SetSketch.estimate_jaccard() not implemented")

    def estimate_union(self, other: 'SetSketch') -> float:
        """Estimate union cardinality with another sketch."""
        raise NotImplementedError("SetSketch.estimate_union() not implemented")

    def write(self, filepath: str) -> None:
        """Write sketch to file in binary format.
        
        Args:
            filepath: Path to output file
        """
        np.savez_compressed(
            filepath,
            registers=self.registers,
            precision=np.array([self.precision]),
            kmer_size=np.array([self.kmer_size]),
            window_size=np.array([self.window_size]),
            seed=np.array([self.seed])
        )

    @classmethod
    def load(cls, filepath: str) -> 'SetSketch':
        """Load sketch from file in binary format.
        
        Args:
            filepath: Path to input file
            
        Returns:
            SetSketch object loaded from file
        """
        data = np.load(filepath)
        sketch = cls(
            precision=int(data['precision'][0]),
            kmer_size=int(data['kmer_size'][0]),
            window_size=int(data['window_size'][0]),
            seed=int(data['seed'][0])
        )
        sketch.registers = data['registers']
        return sketch 