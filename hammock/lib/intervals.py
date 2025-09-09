#!/usr/bin/env python
from __future__ import annotations
from typing import Optional, List, Tuple, TextIO, Union, Dict
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.rusthll import RustHLL
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter
from hammock.lib.exacttest import ExactTest

# Try to import FastHyperLogLog for better performance
try:
    from hammock.lib.hyperloglog_fast import FastHyperLogLog
    FAST_HLL_AVAILABLE = True
except ImportError:
    FAST_HLL_AVAILABLE = False
    FastHyperLogLog = None
import pyBigWig  # type: ignore  # no type stubs available
from numpy import floor  # type: ignore  # no type stubs available
import gzip
import pysam  # type: ignore  # no type stubs available
from multiprocessing import Pool, cpu_count
from itertools import islice
import os
import xxhash  # type: ignore
import numpy as np  # type: ignore
import mmh3  # type: ignore
import sys

def _process_chunk(chunk: List[str], mode: str, subsample: Tuple[float, float], expA: float = 0, sep: str = "-", precision: int = 12, debug: bool = False) -> Tuple[List[str], List[str], List[int]]:
    """Process a chunk of lines in parallel.
    
    Args:
        chunk: List of lines to process
        mode: Mode for comparison (A/B/C/D)
        subsample: Subsampling rate for points
        expA: Power of 10 exponent to use to multiply contribution of A-type intervals
        sep: Separator for string representation
        precision: Precision for HyperLogLog sketch
        debug: Whether to enable debug mode
        
    Returns:
        Tuple of (intervals, points, sizes)
    """
    try:
        intervals = []
        points = []
        sizes = []
        
        if debug:
            print(f"Processing chunk of {len(chunk)} lines")
        
        # Create a temporary sketch to use its bedline method
        temp_sketch = IntervalSketch(
            mode=mode,
            subsample=subsample,
            precision=precision,
            debug=debug,
            expA=expA
        )
        
        for line in chunk:
            if not line.startswith('#'):
                try:
                    # Use bedline to process the line
                    interval, point_list, size = temp_sketch.bedline(
                        line, 
                        mode=mode, 
                        sep=sep, 
                        subsample=subsample,
                        debug=debug
                    )
                    
                    if debug:
                        print(f"Processed line: {line.strip()}")
                        print(f"Got interval: {interval}")
                        print(f"Got {len(point_list)} points")
                    
                    # Process interval
                    if interval:
                        intervals.append(interval)
                        if expA > 0:
                            for i in range(1, int(10**expA)+1):
                                intervals.append(interval + str(i))
                    
                    # Process points - always extend points list in mode B or C
                    if mode in ["B", "C"] and point_list:
                        points.extend(point_list)
                        if debug:
                            print(f"Added {len(point_list)} points from interval, total points now {len(points)}")
                    
                    if size is not None:
                        sizes.append(size)
                except Exception as e:
                    if debug:
                        print(f"Error processing line: {str(e)}")
                    continue
        
        if debug:
            print(f"Chunk processing complete: {len(intervals)} intervals, {len(points)} points")
        
        return intervals, points, sizes
    except Exception as e:
        if debug:
            print(f"Error in _process_chunk: {str(e)}")
        return [], [], []

class IntervalSketch(AbstractSketch):
    """Sketch class for BED intervals."""
    
    def __init__(self, **kwargs):
        """Initialize the sketch."""
        # Extract parameters with defaults
        precision = kwargs.get('precision', 12)
        debug = kwargs.get('debug', False)
        use_cpp = kwargs.get('use_cpp', False)
        # Store mode and validation
        self.mode = kwargs.get('mode', 'A')
        if self.mode not in ['A', 'B', 'C', 'D']:
            raise ValueError(f"Invalid mode: {self.mode}")
        
        # Store and validate expA
        self.expA = kwargs.get('expA', 0)
        if not isinstance(self.expA, (int, float)) or self.expA < 0:
            raise ValueError(f"expA must be a non-negative number, got {self.expA}")
        
        # Validate expA can only be used with mode C
        if self.expA > 0 and self.mode != 'C':
            raise ValueError("Multiplicity (expA) can only be used with mode C")
        
        # Store and validate subsample
        self.subsample = kwargs.get('subsample', (1.0, 1.0))
        if not isinstance(self.subsample, tuple) or len(self.subsample) != 2:
            raise ValueError(f"subsample must be a tuple of length 2, got {self.subsample}")
        if not all(isinstance(x, (int, float)) and 0 <= x <= 1 for x in self.subsample):
            raise ValueError(f"subsample values must be between 0 and 1, got {self.subsample}")
        
        # Choose sketch implementation based on use_cpp flag and availability
        if use_cpp:
            # Use FastHyperLogLog for C++ acceleration
            if FAST_HLL_AVAILABLE:
                self.sketch = FastHyperLogLog(precision=precision, debug=debug)
                if debug:
                    print("Using FastHyperLogLog with C++ acceleration")
            else:
                # Fallback to RustHLL if FastHyperLogLog not available
                self.sketch = RustHLL(precision=precision)
                if debug:
                    print("Using RustHLL (FastHyperLogLog not available)")
        else:
            # Use pure Python HyperLogLog implementation
            self.sketch = HyperLogLog(precision=precision, debug=debug)
            if debug:
                print("Using pure Python HyperLogLog implementation")
        
        # Initialize other instance variables
        self.num_intervals = 0
        self.total_interval_size = 0

    @classmethod
    def from_file(cls, filename: str, **kwargs) -> Optional['IntervalSketch']:
        """Create a sketch from a file.
        
        Args:
            filename: Path to the file to read
            **kwargs: Additional parameters to pass to the sketch constructor
                - mode: Mode for comparison (A/B/C/D)
                - precision: Precision for HyperLogLog sketch
                - num_hashes: Number of hashes for MinHash sketch
                - kmer_size: Size of k-mers
                - window_size: Size of sliding window
                - sketch_type: Type of sketch to use
                - subsample: Subsampling rate for points
                - expA: Power of 10 exponent to use to multiply contribution of A-type intervals
                - debug: Whether to enable debug mode
                - use_cpp: Whether to use C++-accelerated implementation
        """
        debug = kwargs.get('debug', False)
        expA = kwargs.get('expA', 0)
        
        # Check if file exists
        if not os.path.exists(filename):
            if debug:
                print(f"Error: File {filename} does not exist", file=sys.stderr)
            return None
            
        # Check file extension
        supported_extensions = ['.bed', '.bed.gz']
        file_ext = os.path.splitext(filename)[1].lower()
        if file_ext.endswith('.gz'):
            base_ext = os.path.splitext(os.path.splitext(filename)[0])[1].lower()
            if base_ext not in ['.bed']:
                if debug:
                    print(f"Skipping unsupported gzipped file format: {filename}")
                return None
        elif file_ext not in supported_extensions:
            if debug:
                print(f"Skipping unsupported file format: {filename}")
            return None
            
        # Print start of file processing only in debug mode
        if debug:
            print(f"Processing file: {filename}")
            
        # Create sketch with appropriate parameters
        sketch = cls(**kwargs)
        
        # Get separator from kwargs, default to tab
        sep = kwargs.get('sep', '\t')
        
        # Read file and add intervals
        try:
            opener = gzip.open if file_ext.endswith('.gz') else open
            mode = 'rt' if file_ext.endswith('.gz') else 'r'
            
            with opener(filename, mode) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    # Skip header lines that start with common header words
                    first_word = line.split()[0].lower() if line.split() else ""
                    if first_word in ['chromosome', 'start', 'end', 'chrom', 'chromosome_name']:
                        continue
                    try:
                        interval, points, size = sketch.bedline(line, mode=kwargs.get('mode', 'A'), sep=sep, subsample=kwargs.get('subsample', (1.0, 1.0)))
                        if interval:
                            sketch.sketch.add_string(interval)
                            # Add multiple copies of the interval if expA is set
                            if expA > 0:
                                for i in range(1, int(10**expA) + 1):
                                    sketch.sketch.add_string(interval + str(i))
                            sketch.num_intervals += 1  # Increment counter for each interval
                            sketch.total_interval_size += size
                        if points:
                            for point in points:
                                if point is not None:
                                    sketch.sketch.add_string(point)
                    except ValueError as e:
                        if debug:
                            print(f"Skipping invalid line: {line.strip()[:50]}... (Error: {e})")
                        continue
            
            # Print completion of file processing only in debug mode
            if debug:
                print(f"Completed processing file: {filename}")
                print(f"  - Number of intervals: {sketch.num_intervals}")
                print(f"  - Total interval size: {sketch.total_interval_size}")
                if kwargs.get('mode', 'A') in ['B', 'C']:
                    print(f"  - Estimated cardinality: {sketch.sketch.estimate_cardinality():.2f}")
            
            return sketch
        except Exception as e:
            if debug:
                print(f"Error reading file {filename}: {e}", file=sys.stderr)
            return None

    @classmethod
    def _from_bigbed(cls, filename: str, **kwargs) -> Optional['IntervalSketch']:
        """Process BigBed format file."""
        if 'sketch' not in kwargs:
            return None
        
        try:
            # Get parameters from kwargs
            mode = kwargs.get('mode', 'A')
            sep = kwargs.get('sep', '-')
            subsample = kwargs.get('subsample', (1.0, 1.0))
            expA = kwargs.get('expA', 0)
            sketch = kwargs['sketch']
            # Use pysam to read BigBed file
            with pysam.TabixFile(filename) as bb:
                for line in bb.fetch():
                    interval, points, size = sketch.bedline(line, mode=mode, sep=sep, subsample=subsample)
                    
                    # Add interval if present and in mode A or C
                    if interval and mode in ["A", "C"]:
                        sketch.sketch.add_string(interval)
                        if expA > 0:
                            for i in range(1, int(10**expA)+1):
                                sketch.sketch.add_string(interval + str(i))
                        # sketch.num_intervals += 1
                        # sketch.total_interval_size += size
                    
                    # Add points if present and in mode B or C
                    if points and mode in ["B", "C"]:
                        for point in points:
                            if point is not None:
                                sketch.sketch.add_string(point)
                        # if mode == "B":
                        #     sketch.num_intervals += 1
            return sketch
        except Exception as e:
            print(f"Error processing BigBed file {filename}: {e}")
            return None

    @classmethod
    def _from_bed(cls, filename: str, **kwargs) -> Optional['IntervalSketch']:
        """Process BED format file."""
        if 'sketch' not in kwargs:
            return None
            
        try:
            # Get parameters from kwargs
            mode = kwargs.get('mode', 'A')
            sep = kwargs.get('sep', '-')
            subsample = kwargs.get('subsample', (1.0, 1.0))
            expA = kwargs.get('expA', 0)
            sketch = kwargs['sketch']
            debug = kwargs.get('debug', False)
            use_cpp = kwargs.get('use_cpp', False)
            
            # If using C++ implementation, be more conservative with Python-level threading
            if use_cpp:
                num_threads = 1  # Let Rust handle the threading
            else:
                num_threads = min(cpu_count(), 4)  # Limit Python threads
            
            if debug:
                print(f"Processing BED file {filename} with mode={mode}, subsample={subsample}, debug={debug}")
                print(f"Using {num_threads} Python threads (C++ threading: {use_cpp})")
            
            # Open file with gzip if it ends in .gz, otherwise normal open
            opener = gzip.open if filename.endswith('.gz') else open
            
            # Read all lines into memory for parallel processing
            with opener(filename, 'rt') as f:
                lines = []
                for line in f:
                    if not line.startswith('#'):
                        try:
                            # Validate line format
                            sketch.basic_bedline(line)
                            lines.append(line)
                        except ValueError as e:
                            if debug:
                                print(f"Error processing line: {str(e)}")
                            continue
            
            if not lines:
                print(f"No valid lines found in file: {filename}")
                return None
            
            if debug:
                print(f"Found {len(lines)} valid lines to process")
            
            # Calculate chunk size for parallel processing
            chunk_size = max(1, len(lines) // num_threads)
            
            if debug:
                print(f"Using {num_threads} threads with chunk size {chunk_size}")
            
            # If using C++, process in a single chunk to let C++ handle threading
            if use_cpp:
                chunk_args = [(lines, mode, subsample, expA, sep, sketch.sketch.precision, debug)]
            else:
                # Prepare arguments for each chunk
                chunk_args = []
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i + chunk_size]
                    args = (chunk, mode, subsample, expA, sep, sketch.sketch.precision, debug)
                    chunk_args.append(args)
                
                # Process chunks
                if num_threads > 1:
                    with Pool(processes=num_threads) as pool:
                        results = pool.starmap(_process_chunk, chunk_args)
                else:
                    results = [_process_chunk(*args) for args in chunk_args]
            
            # Combine results from all chunks
            all_intervals = []
            all_points = []
            all_sizes = []
            
            for intervals, points, sizes in results:
                all_intervals.extend(intervals)
                all_points.extend(points)
                all_sizes.extend(sizes)
            
            if debug:
                print(f"Total points collected: {len(all_points)}")
                print(f"Total intervals collected: {len(all_intervals)}")
                if len(all_points) > 0:
                    print(f"Sample points (first 5): {all_points[:5]}")
            
            # Add intervals to sketch
            if all_intervals and mode in ["A", "C"]:
                sketch.sketch.add_batch(all_intervals)
                sketch.num_intervals += len(all_intervals)
                sketch.total_interval_size += sum(all_sizes)
                if debug:
                    print(f"Added {len(all_intervals)} intervals to sketch")
                    print(f"Sketch cardinality after adding intervals: {sketch.sketch.estimate_cardinality()}")
            
            # Add points to sketch
            if all_points and mode in ["B", "C"]:
                if debug:
                    print(f"Adding {len(all_points)} points to sketch")
                    # Create a new sketch just for points to see their cardinality
                    if len(all_points) > 0:
                        from hammock.lib.hyperloglog import HyperLogLog
                        point_sketch = HyperLogLog(precision=sketch.sketch.precision)
                        point_sketch.add_batch(all_points)
                        print(f"Point-only sketch cardinality: {point_sketch.estimate_cardinality()}")
                
                sketch.sketch.add_batch(all_points)
                if debug:
                    print(f"Sketch cardinality after adding points: {sketch.sketch.estimate_cardinality()}")
                
                # if mode == "B":
                #     sketch.num_intervals = len(lines)  # Count number of intervals processed
            
            return sketch
            
        except Exception as e:
            if debug:
                print(f"Error processing BED file {filename}: {e}")
            return None

    @staticmethod
    def basic_bedline(line: str) -> tuple[str, int, int]:
        """Parse a single line from a BED file into chromosome, start, and end coordinates.
        
        Args:
            line: A string containing a single line from a BED file
            
        Returns:
            Tuple of (chrom, start, end) where:
                chrom: Chromosome name (any string)
                start: Integer start coordinate
                end: Integer end coordinate
                
        Raises:
            ValueError: If line has fewer than 3 tab or space-separated columns
            ValueError: If start or end coordinates cannot be converted to integers
        """
        # Split on either tab or space
        columns = line.strip().split()
        if len(columns) < 3:
            raise ValueError("bedline: line must have at least 3 columns")
        
        # Allow any string for chromosome name
        chrval = columns[0]
        
        try:
            start = int(columns[1])
            end = int(columns[2])
        except ValueError:
            raise ValueError("bedline: start and end coordinates must be integers")
        
        # Validate coordinates
        if start < 0 or end < 0:
            raise ValueError("bedline: coordinates must be non-negative")
        if start >= end:
            raise ValueError("bedline: start coordinate must be less than end coordinate")
        
        return chrval, start, end

    def bedline(self, line: str, 
                mode: str, 
                sep: str, 
                subsample: tuple[float, float] = (1.0, 1.0),
                debug: bool = False) -> tuple[Optional[str], list[Optional[str]], int]:
        """Process a single BED line."""
        interval = None
        points = []
        chrval, start, end = self.basic_bedline(line)
        interval_size = end - start  # Don't include end point for BED format
        
        if mode in ["A", "C"]:
            # For mode C, use integer hashing for subsampling
            if mode == "C":
                # Create the full interval string and hash it
                interval_str = sep.join([chrval, str(start), str(end), "A"])
                interval_hash = xxhash.xxh32(interval_str, seed=31337).intdigest()
                if interval_hash < 0:
                    interval_hash += 2**32
                threshold = int(subsample[0] * (2**32 - 1))
                if interval_hash > threshold:
                    interval = None
                else:
                    interval = interval_str
            else:
                interval = sep.join([chrval, str(start), str(end), "A"])
        
        if mode in ["B", "C"]:
            # Apply point subsampling in both modes B and C
            points = self.generate_points(chrval, start, end, sep=sep, subsample=subsample[1], debug=debug)
        
        return interval, points, interval_size

    def generate_points(self, chrval: str, 
                       start: int, 
                       end: int, 
                       sep: str = "-", 
                       subsample: float = 1, 
                       seed: int = 31337,
                       debug: bool = False) -> list[Optional[str]]:
        """Generate point strings for a given interval.
        
        Args:
            chrval: Chromosome name
            start: Start position (inclusive)
            end: End position (exclusive)
            sep: Separator for string representation
            subsample: Subsampling rate
            seed: Random seed for hashing
            debug: Whether to enable debug mode
            
        Returns:
            List of point strings in BED format's half-open interval [start, end)
        """
        points = []
        if subsample >= 1.0:
            # No subsampling needed, generate all points in [start, end)
            for x in range(start, end):
                points.append(sep.join([chrval, str(x)]))
        else:
            # Apply subsampling
            if debug:
                print(f"Subsampling with rate={subsample} for {end-start} points")
            
            # Generate hashes for all points in [start, end)
            if debug:
                print("Generating hashes...")
            hashes = []
            for x in range(start, end):
                point_str = sep.join([chrval, str(x)])
                # Get hash value and convert to unsigned 32-bit
                hash_val = xxhash.xxh32(point_str, seed=seed).intdigest()
                if hash_val < 0:
                    hash_val += 2**32
                hashes.append(hash_val)
            if debug:
                print(f"Generated {len(hashes)} hashes")
            
            hashes = np.array(hashes, dtype=np.uint32)
            if debug:
                print("Converted to numpy array")
            
            # Calculate threshold for subsampling
            # Use the full range of hash values [0, 2^32-1]
            threshold = int(subsample * (2**32 - 1))
            mask = hashes <= threshold
            if debug:
                print(f"Number of points passing subsampling: {np.sum(mask)}")
            
            # Generate strings only for points that pass subsampling
            for i, x in enumerate(range(start, end)):
                if mask[i]:
                    points.append(sep.join([chrval, str(x)]))
        
        return points

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.sketch.add_string(s)
        
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        # Filter out None values
        valid_strings = [s for s in strings if s is not None]
        if valid_strings:
            self.sketch.add_batch(valid_strings)
            
    def add(self, s: str) -> None:
        """Alias for add_string."""
        self.add_string(s)

    def estimate_jaccard(self, other: 'IntervalSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        return self.sketch.estimate_jaccard(other.sketch)

    def similarity_values(self, other: 'IntervalSketch') -> float:
        """Calculate similarity values between this sketch and another.
        
        Args:
            other: Another IntervalSketch to compare against
            
        Returns:
            Dictionary containing 'jaccard_similarity' and optionally 'containment'
        """
        if not isinstance(other, IntervalSketch):
            raise ValueError("Can only compare with another IntervalSketch")
        
        # Calculate Jaccard similarity
        jaccard = self.estimate_jaccard(other)
        result = {'jaccard_similarity': jaccard}
        
        # Calculate containment if expA is set
       
        containment = self.estimate_containment(other)
        result['containment'] = containment
        
        return result

    def write(self, filepath: str) -> None:
        """Write sketch to file, with gzip compression if filepath ends in .gz."""
        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'wt') as f:  # 'wt' mode for text writing to gzip
                self.sketch.write(f)
        else:
            with open(filepath, 'w') as f:
                self.sketch.write(f)

    @classmethod
    def load(cls, filepath: str) -> 'IntervalSketch':
        """Load sketch from file, handling gzipped files if filepath ends in .gz."""
        opener = gzip.open if filepath.endswith('.gz') else open
        with opener(filepath, 'rt') as f:
            # First try loading as HyperLogLog
            try:
                sketch = HyperLogLog.load(f)
                interval_sketch = cls(mode="A", sketch_type="hyperloglog")
                interval_sketch.sketch = sketch
                return interval_sketch
            except:
                pass

            # Then try as RustHLL
            try:
                sketch = RustHLL.load(f)
                interval_sketch = cls(mode="A", sketch_type="hyperloglog", use_cpp=True)
                interval_sketch.sketch = sketch
                return interval_sketch
            except:
                pass

            # Then try as MinHash
            try:
                sketch = MinHash.load(f)
                interval_sketch = cls(mode="A", sketch_type="minhash")
                interval_sketch.sketch = sketch
                return interval_sketch
            except:
                pass

            # Finally try as ExactTest
            try:
                sketch = ExactTest.load(f)
                interval_sketch = cls(mode="A", sketch_type="exacttest")
                interval_sketch.sketch = sketch
                return interval_sketch
            except:
                raise ValueError("Could not load sketch from file")

    def estimate_similarity(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Estimate similarity between interval sketches.
        
        Returns:
            Dictionary containing 'jaccard_similarity' and optionally 'containment'
        """
        if not isinstance(other, IntervalSketch):
            raise ValueError("Can only compare with another IntervalSketch")
        
        # Calculate Jaccard similarity
        jaccard = self.estimate_jaccard(other)
        result = {'jaccard_similarity': jaccard}
        
        # Calculate containment if expA is set
        if hasattr(self, 'expA') and self.expA is not None:
            containment = self.estimate_containment(other)
            result['containment'] = containment
        
        return result

    def estimate_containment(self, other: 'IntervalSketch') -> float:
        """Estimate containment of other sketch in this one."""
        if not isinstance(other, IntervalSketch):
            raise ValueError("Can only compare with another IntervalSketch")
        
        # Get intersection and self cardinality
        intersection = self.sketch.estimate_intersection(other.sketch)
        self_card = self.sketch.estimate_cardinality()
        
        # Calculate containment
        if self_card == 0:
            return 0.0
        
        # Apply exponential scaling if expA is set
        containment = intersection / self_card
        if hasattr(self, 'expA') and self.expA is not None:
            containment = containment ** self.expA
        
        return containment
            