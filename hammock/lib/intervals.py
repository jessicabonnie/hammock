#!/usr/bin/env python
from __future__ import annotations
from typing import Optional, List, Tuple
from hammock.lib.abstractsketch import AbstractDataSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter
from hammock.lib.exacttest import ExactTest
import pyBigWig
from numpy import floor


class IntervalSketch(AbstractDataSketch):
    """Sketch class for BED intervals."""
    
    def __init__(self, 
                 mode: str,
                 sketch_type: str = "hyperloglog",
                 precision: int = 8,
                 num_hashes: int = 128,
                 seed: int = 0,
                 expA: float = 0,
                 subsample: tuple[float, float] = (1.0, 1.0)):
        """Initialize interval sketch.
        
        Args:
            mode: A/B/C for interval/point/both comparison types
            sketch_type: "hyperloglog", "minhash", or "exact"
            precision: Precision for HyperLogLog
            num_hashes: Number of hashes for MinHash
            seed: Random seed
            expA: Exponent for interval multiplicity (mode C only)
            subsample: Tuple of (interval_rate, point_rate) between 0 and 1
        """
        if sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(precision=precision, seed=seed)
        elif sketch_type == "minhash":
            self.sketch = MinHash(num_hashes=num_hashes, seed=seed)
        elif sketch_type == "exacttest":
            self.sketch = ExactTest(seed=seed)
        else:
            raise ValueError(f"Invalid sketch type for intervals: {sketch_type}")
            
        self.mode = mode
        self.total_interval_size = 0
        self.num_intervals = 0
        self.expA = expA
        self.subsample = subsample
        
        # Validate expA usage
        if expA > 0:
            if mode != "C":
                raise ValueError("Multiplicity (expA) can only be used with mode C")
            if subsample != (1.0, 1.0):
                raise ValueError("Multiplicity (expA) cannot be used with subsampling")

    @classmethod
    def from_file(cls, filename: str, mode: str, precision: int = 8, 
                  sketch_type: str = "hyperloglog", sep: str = "-", 
                  subsample: tuple[float, float] = (1.0, 1.0),
                  expA: float = 0,
                  **kwargs) -> Optional['IntervalSketch']:
        """Create sketch from BED or BigBed file.
        
        Args:
            filename: Path to BED/BigBed file
            mode: Mode for interval processing (A/B/C)
            precision: Precision for HyperLogLog sketching
            sketch_type: Type of sketch to use
            sep: Separator for interval string representation
            subsample: Tuple of (interval_rate, point_rate) between 0 and 1,
            expA: Exponent for interval multiplicity (mode C only)
            **kwargs: Additional arguments for sketch initialization
            
        Returns:
            New IntervalSketch instance or None if error
        """
         # Validate expA usage
        if expA > 0:
            if mode != "C":
                raise ValueError("Multiplicity (expA) can only be used with mode C")
            if subsample != (1.0, 1.0):
                raise ValueError("Multiplicity (expA) cannot be used with subsampling")

        sketch = cls(mode=mode, precision=precision, sketch_type=sketch_type)
        
        try:
            if filename.endswith('.bb'):
                return cls._from_bigbed(filename, mode, sep, subsample, expA, sketch)
            else:
                return cls._from_bed(filename, mode, sep, subsample, expA, sketch)
        except Exception as e:
            print(f"Error processing file {filename}: {e}")
            return None

    @classmethod
    def _from_bigbed(cls, filename: str, mode: str, sep: str, 
                     subsample: tuple[float, float], expA: float, sketch: 'IntervalSketch') -> Optional['IntervalSketch']:
        """Process BigBed format file."""
        try:
            with pyBigWig.open(filename) as bw:
                # Get chromosomes and their sizes
                chroms = bw.chroms()
                for chrom in chroms:
                    # Fetch all intervals for this chromosome
                    intervals = bw.entries(chrom, 0, chroms[chrom])
                    if intervals:
                        for start, end, _ in intervals:
                            line = f"{chrom}\t{start}\t{end}"
                            interval, points, size = sketch.bedline(line, mode=mode, sep=sep, subsample=subsample)
                            
                            # Add interval if present and in mode A or C
                            if interval and mode in ["A", "C"]:
                                sketch.sketch.add_string(interval)
                                if expA > 0:
                                    for i in range(1, floor(10**expA)+1):
                                        sketch.sketch.add_string(interval + str(i))
                                sketch.num_intervals += 1
                                sketch.total_interval_size += size
                            
                            # Add points if present and in mode B or C
                            if points and mode in ["B", "C"]:
                                for point in points:
                                    if point is not None:
                                        sketch.sketch.add_string(point)
                                        if mode == "B":
                                            sketch.num_intervals += 1
            return sketch
        except Exception as e:
            print(f"Error processing BigBed file {filename}: {e}")
            return None

    @classmethod
    def _from_bed(cls, filename: str, mode: str, sep: str, 
                  subsample: Tuple[float, float], expA: float, sketch: 'IntervalSketch') -> Optional['IntervalSketch']:
        """Process BED format file."""
        try:
            with open(filename, 'r') as f:
                for line in f:
                    interval, points, size = sketch.bedline(line, mode=mode, sep=sep, subsample=subsample)
                    
                    # Add interval if present and in mode A or C
                    if interval and mode in ["A", "C"]:
                        print(f"\nOriginal interval: {interval}")
                        sketch.sketch.add_string(interval)
                        print(f"Cardinality after original: {sketch.sketch.estimate_cardinality()}")
                        
                        if expA > 0:
                            num_copies = int(floor(10**expA)) - 1
                            print(f"Adding {num_copies} copies...")
                            for i in range(1, num_copies + 1):
                                sketch.sketch.add_string(interval + str(i))
                            print(f"Cardinality after copies: {sketch.sketch.estimate_cardinality()}")
                        
                        sketch.num_intervals += 1
                        sketch.total_interval_size += size
                    
                    # Add points if present and in mode B or C
                    if points and mode in ["B", "C"]:
                        non_none_points = [p for p in points if p is not None]
                        print(f"\nAdding {len(non_none_points)} points")
                        print(f"First point: {non_none_points[0]}")
                        print(f"Last point: {non_none_points[-1]}")
                        for point in points:
                            if point is not None:
                                sketch.sketch.add_string(point)
                        print(f"Cardinality after points: {sketch.sketch.estimate_cardinality()}")
                        
                        if mode == "B":
                            sketch.num_intervals += 1
                    
                    print(f"\nFinal cardinality: {sketch.sketch.estimate_cardinality()}")
            return sketch
        except Exception as e:
            print(f"Error processing BED file {filename}: {e}")
            return None

    @staticmethod
    def basic_bedline(line: str) -> tuple[str, int, int]:
        """Parse a single line from a BED file into chromosome, start, and end coordinates.
        
        Args:
            line: A string containing a single line from a BED file
            
        Returns:
            Tuple of (chrom, start, end) where:
                chrom: Chromosome name (preserves 'chr' prefix if present)
                start: Integer start coordinate
                end: Integer end coordinate
                
        Raises:
            ValueError: If line has fewer than 3 tab or space-separated columns
        """
        columns = line.strip().split('\t')
        if len(columns) < 2:
            columns = line.strip().split(" ")
            if len(columns) < 2:
                raise ValueError("bedline: one of the lines in malformed")
        # Ensure chromosome has 'chr' prefix
        chrval = columns[0][3:] if columns[0].startswith('chr') else columns[0]
        return chrval, int(columns[1]), int(columns[2])

    def bedline(self, line: str, 
                mode: str, 
                sep: str, 
                subsample: tuple[float, float] = (1.0, 1.0)) -> tuple[Optional[str], list[Optional[bytes]], int]:
        """Process a single BED line.
        
        Args:
            line: BED file line
            mode: Processing mode (A/B/C)
            sep: Separator for string representation
            subsample: Tuple of (interval_rate, point_rate) between 0 and 1
            
        Returns:
            Tuple of (interval, points, interval_size) where:
                interval: String representation of interval (or None if subsampled)
                points: List of point strings as bytes (some may be None if subsampled)
                interval_size: Size of the interval
        """
        interval = None
        points = []
        chrval, start, end = self.basic_bedline(line)
        # start = min(startx, endx)
        # end = max(startx, endx)
        interval_size = end - start
        # subsample=self.subsample
        
        if mode in ["A", "C"]:
            interval = sep.join([chrval, str(start), str(end), "A"])
            # Only apply interval subsampling in mode C
            if mode == "C":
                hashv = self.sketch._hash_str(interval.encode('utf-8'), seed=777)
                if hashv % (2**32) > int(subsample[0] * (2**32)):
                    interval = None
            
        if mode in ["B", "C"]:
            # Only apply point subsampling in mode C
            subsample_rate = subsample[1] if mode == "C" else 1.0
            points = self.generate_points(chrval, start, end, sep=sep, subsample=subsample_rate)
            print(f"Generated points: {len(points)}")
            print(f"Start: {start}, End: {end}")
        
        return interval, points, interval_size

    def generate_points(self, chrval: str, 
                       start: int, 
                       end: int, 
                       sep: str = "-", 
                       subsample: float = 1, 
                       seed: int = 23) -> list[Optional[str]]:
        """Generate point strings for a given interval.
        
        Args:
            chrval: Chromosome value
            start: Start position
            end: End position
            sep: Separator for string representation
            subsample: Subsampling rate between 0 and 1
            seed: Random seed for subsampling
            
        Returns:
            List of point strings (some may be None if subsampled)
        """
        maximum = int(min(1, subsample) * (2**32))
        def gp(x):
            outstr = sep.join([str(chrval), str(x), str(x+1)])
            hashv = self.sketch._hash_str(outstr.encode('utf-8'), seed)
            if hashv % (2**32) <= maximum:
                return outstr
            return None
        # Recall that bed files are 0-indexed, so the final point should have it's start coordinate at end-1
        return [gp(x) for x in range(start, end)]

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.sketch.add_string(s)

    def estimate_jaccard(self, other: 'IntervalSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        return self.sketch.estimate_jaccard(other.sketch)