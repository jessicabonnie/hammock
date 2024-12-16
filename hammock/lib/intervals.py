#!/usr/bin/env python
from typing import Optional, Tuple, List
from hammock.lib.sketchclass import Sketch
from hammock.lib.bedprocessing import basic_bedline
from itertools import islice
import gc

def bedline(line: str, 
            mode: str, 
            sep: str, 
            subsample: float = 1) -> tuple[Optional[bytes], list[Optional[bytes]], int]:
    """Process a single line describing an interval from a BED file.
    
    Args:
        line: Single line from BED file
        mode: A/B/C for interval/point/both comparison types
        sep: Separator for combining fields
        subsample: Subsampling rate. If negative, interval values are subsampled
        
    Returns:
        Tuple of (interval, points, interval_size) where:
            interval: Encoded interval string or None
            points: List of encoded point strings or empty list
            interval_size: Size of the interval
    """
    interval = None
    points = []
    chrval, startx, endx = basic_bedline(line)
    start = min(startx, endx)
    end = max(startx, endx)
    interval_size = end - start
    
    if mode in ["A","C"]:
        interval = sep.join([chrval, str(start), str(end), "A"]).encode('utf-8')
        if subsample < 0:
            hashv = HyperLogLog._hash_str(interval, seed=777)
            if hashv % (2**32) > int((1+subsample) * (2**32)):
                interval = None
    if mode in ["B","C"]:
        points = generate_points(chrval, start, end, sep=sep, subsample=abs(subsample))
        
    return interval, points, interval_size

def generate_points(chrval: str, 
                   start: int, 
                   end: int, 
                   sep: str = "-", 
                   subsample: float = 1, 
                   seed: int = 23) -> list[Optional[bytes]]:
    """Generate points from a BED interval with optional subsampling.
    
    Args:
        chrval: Chromosome value
        start: Start position
        end: End position
        sep: Separator for combining fields
        subsample: Fraction of points to sample (0-1)
        seed: Random seed for sampling
        
    Returns:
        List of encoded point strings that passed sampling
    """
    maximum = int(subsample * (2**32))
    def gp(x):
        outstr = sep.join([str(chrval), str(x), str(x+1)])
        hashv = HyperLogLog._hash_str(outstr.encode('utf-8'), seed)
        if hashv % (2**32) <= maximum:
            return outstr.encode('utf-8')
    return [gp(x) for x in range(start, end+1)]

class IntervalSketch(Sketch):
    """Sketch class for BED intervals and points."""
    
    def __init__(self, 
                 mode: str,
                 precision: int = 8,
                 num_hashes: int = 128,
                 sketch_type: str = "hyperloglog",
                 seed: int = 0):
        """Initialize IntervalSketch.
        
        Args:
            mode: A/B/C for interval/point/both comparison types
            precision: Precision for HyperLogLog sketching
            num_hashes: Number of hash functions for MinHash sketching
            sketch_type: Type of sketch to use
            seed: Random seed for hash functions
        """
        super().__init__(
            sketch_type=sketch_type,
            precision=precision,
            num_hashes=num_hashes,
            kmer_size=0,
            seed=seed
        )
        self.mode = mode
        self.total_interval_size = 0
        self.num_intervals = 0
        
    @classmethod
    def from_file(cls,
                  filename: str,
                  mode: str,
                  precision: int = 8,
                  num_hashes: int = 128,
                  sketch_type: str = "hyperloglog",
                  sep: str = "-",
                  subsample: float = 1,
                  chunk_size: int = 1000,
                  verbose: bool = False) -> Optional['IntervalSketch']:
        """Create an IntervalSketch from a BED file.
        
        Args:
            filename: Path to BED file
            mode: A/B/C for interval/point/both comparison types
            precision: Precision for HyperLogLog
            num_hashes: Number of hash functions for MinHash
            sketch_type: Type of sketch to use
            sep: Separator for combining fields
            subsample: Subsampling rate
            chunk_size: Number of lines to process at once
            verbose: Whether to print progress
            
        Returns:
            IntervalSketch object or None if file processing fails
        """
        sketch = cls(
            mode=mode,
            precision=precision,
            num_hashes=num_hashes,
            sketch_type=sketch_type
        )
        
        try:
            with open(filename, "r") as file:
                while True:
                    chunk = list(islice(file, chunk_size))
                    if not chunk:
                        break
                    
                    bedline_results = [bedline(line, mode=mode, sep=sep, subsample=subsample) 
                                     for line in chunk if not line.startswith('#')]
                    
                    sketch.total_interval_size += sum(isize for _, _, isize in bedline_results)
                    sketch.num_intervals += len(bedline_results)
                    
                    if mode in ["B", "C"]:
                        for _, b, _ in bedline_results:
                            if b is not None:
                                for point in b:
                                    if point is not None:
                                        sketch.add_string(point.decode())
                    
                    if mode in ["A", "C"]:
                        for a, _, _ in bedline_results:
                            if a is not None:
                                sketch.add_string(a.decode())
                    
                    if verbose:
                        print(f"Processed {len(chunk)} lines from {filename}")
                    gc.collect()
                    
            return sketch
            
        except Exception as e:
            print(f"Error processing file {filename}: {str(e)}")
            return None