#!/usr/bin/env python
from typing import Optional, List
from hammock.lib.abstractsketch import AbstractDataSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter


class IntervalSketch(AbstractDataSketch):
    """Sketch class for BED intervals."""
    
    def __init__(self, 
                 mode: str,
                 sketch_type: str = "hyperloglog",
                 precision: int = 8,
                 num_hashes: int = 128,
                 seed: int = 0):
        """Initialize interval sketch.
        
        Args:
            mode: A/B/C for interval/point/both comparison types
            sketch_type: "hyperloglog", "minhash", or "exact"
            precision: Precision for HyperLogLog
            num_hashes: Number of hashes for MinHash
            seed: Random seed
        """
        if sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(precision=precision, seed=seed)
        elif sketch_type == "minhash":
            self.sketch = MinHash(num_hashes=num_hashes, seed=seed)
        elif sketch_type == "exact":
            self.sketch = ExactCounter(seed=seed)
        else:
            raise ValueError(f"Invalid sketch type for intervals: {sketch_type}")
            
        self.mode = mode
        self.total_interval_size = 0
        self.num_intervals = 0

    @classmethod
    def from_file(cls, filename: str, mode: str, precision: int = 8, 
                  sketch_type: str = "hyperloglog", sep: str = "-", 
                  subsample: float = 1.0, **kwargs) -> Optional['IntervalSketch']:
        """Create sketch from BED file.
        
        Args:
            filename: Path to BED file
            mode: Mode for interval processing (A/B/C)
            precision: Precision for HyperLogLog sketching
            sketch_type: Type of sketch to use
            sep: Separator for interval string representation
            subsample: Subsampling rate (-1 to 1). Negative for intervals, positive for points
            **kwargs: Additional arguments for sketch initialization
            
        Returns:
            New IntervalSketch instance or None if error
        """
        sketch = cls(mode=mode, precision=precision, sketch_type=sketch_type)
        
        try:
            with open(filename, 'r') as f:
                for line in f:
                    interval, points, size = sketch.bedline(line, mode=mode, sep=sep, subsample=subsample)
                    
                    # Add interval if present and in mode A or C
                    if interval and mode in ["A", "C"]:
                        sketch.sketch.add_string(interval.decode('utf-8'))
                        sketch.num_intervals += 1
                        sketch.total_interval_size += size
                    
                    # Add points if present and in mode B or C
                    if points and mode in ["B", "C"]:
                        for point in points:
                            if point is not None:  # Some points might be None due to subsampling
                                sketch.sketch.add_string(point.decode('utf-8'))
                                # Don't increment num_intervals for points in mode C
                                if mode == "B":
                                    sketch.num_intervals += 1
            return sketch
        except Exception as e:
            print(f"Error processing file {filename}: {e}")
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
                subsample: tuple[float, float] = (1.0, 1.0)) -> tuple[Optional[bytes], list[Optional[bytes]], int]:
        """Process a single BED line.
        
        Args:
            line: BED file line
            mode: Processing mode (A/B/C)
            sep: Separator for string representation
            subsample: Tuple of (interval_rate, point_rate) between 0 and 1
        """
        interval = None
        points = []
        chrval, startx, endx = self.basic_bedline(line)
        start = min(startx, endx)
        end = max(startx, endx)
        interval_size = end - start
        
        if mode in ["A", "C"]:
            interval = sep.join([chrval, str(start), str(end), "A"]).encode('utf-8')
            # Only apply interval subsampling in mode C
            if mode == "C":
                hashv = self.sketch._hash_str(interval, seed=777)
                if hashv % (2**32) > int(subsample[0] * (2**32)):
                    interval = None
            
        if mode in ["B", "C"]:
            # Only apply point subsampling in mode C
            subsample_rate = subsample[1] if mode == "C" else 1.0
            points = self.generate_points(chrval, start, end, sep=sep, subsample=subsample_rate)
        
        return interval, points, interval_size

    def generate_points(self, chrval: str, 
                       start: int, 
                       end: int, 
                       sep: str = "-", 
                       subsample: float = 1, 
                       seed: int = 23) -> list[Optional[bytes]]:
        maximum = int(subsample * (2**32))
        def gp(x):
            outstr = sep.join([str(chrval), str(x), str(x+1)])
            hashv = self.sketch._hash_str(outstr.encode('utf-8'), seed)
            if hashv % (2**32) <= maximum:
                return outstr.encode('utf-8')
            return None
        # Recall that bed files are 0-indexed, so the final point should have it's start coordinate at end-1
        return [gp(x) for x in range(start, end)]

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.sketch.add_string(s)

    def estimate_jaccard(self, other: 'IntervalSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        return self.sketch.estimate_jaccard(other.sketch)