#!/usr/bin/env python
from __future__ import annotations
from typing import Optional, List, Tuple, TextIO, Union, Dict
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter
from hammock.lib.exacttest import ExactTest
import pyBigWig  # type: ignore  # no type stubs available
from numpy import floor  # type: ignore  # no type stubs available
import gzip
import pysam  # type: ignore  # no type stubs available


class IntervalSketch(AbstractSketch):
    """Sketch class for BED intervals."""
    
    def __init__(self, mode: str = "A", precision: int = 8, 
                 sketch_type: str = "hyperloglog", expA: float = 0,
                 subsample: tuple[float, float] = (1.0, 1.0),
                 debug: bool = False,
                 **kwargs):
        """Initialize interval sketch.
        
        Args:
            mode: Mode for interval processing (A/B/C)
            precision: Precision for HyperLogLog sketching
            sketch_type: Type of sketch to use
            expA: Exponent for interval multiplicity (mode C only)
            subsample: Tuple of (interval_rate, point_rate) between 0 and 1
            debug: Whether to print debug information
        """
        super().__init__()
        # Validate mode
        if mode not in ["A", "B", "C"]:
            raise ValueError(f"Invalid mode: {mode}. Must be one of: A, B, C")
            
        # Validate expA
        if expA < 0:
            raise ValueError("expA must be non-negative")
        if expA > 0 and mode != "C":
            raise ValueError("Multiplicity (expA) can only be used with mode C")
            
        # Validate subsample rates
        if not (0 <= subsample[0] <= 1 and 0 <= subsample[1] <= 1):
            raise ValueError("Subsample rates must be between 0 and 1")
            
        self.mode = mode
        self.expA = expA
        self.subsample = subsample
        self.precision = precision
        self.sketch_type = sketch_type
        self.debug = debug
        
        # Initialize sketch based on type
        if sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(precision=precision, debug=debug, **kwargs)
        elif sketch_type == "minhash":
            self.sketch = MinHash(debug=debug, **kwargs)
        elif sketch_type == "exact":
            self.sketch = ExactCounter(**kwargs)
        else:
            raise ValueError(f"Invalid sketch type: {sketch_type}")
        
        self.total_interval_size = 0
        self.num_intervals = 0
        
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
        """Create sketch from BED/BigBed/BigWig/BAM file.
        
        Args:
            filename: Path to BED/BigBed/BigWig/BAM file
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
        sketch = cls(mode=mode, precision=precision, sketch_type=sketch_type, 
                    expA=expA, subsample=subsample)
        
        try:
            if filename.endswith('.bb'):
                return cls._from_bigbed(filename, mode, sep, subsample, expA, sketch)
            elif filename.endswith('.bw'):
                return cls._from_bigwig(filename, mode, sep, subsample, expA, sketch)
            elif filename.endswith('.bam'):
                return cls._from_bam(filename, mode, sep, subsample, expA, sketch)
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
    def _from_bigwig(cls, filename: str, mode: str, sep: str, 
                     subsample: tuple[float, float], expA: float, sketch: 'IntervalSketch') -> Optional['IntervalSketch']:
        """Process BigWig format file."""
        try:
            with pyBigWig.open(filename) as bw:
                # Get chromosomes and their sizes
                chroms = bw.chroms()
                for chrom in chroms:
                    # Get intervals for this chromosome using intervals() method
                    # This returns runs of non-zero values as intervals
                    intervals = bw.intervals(chrom, 0, chroms[chrom])
                    if intervals:
                        for start, end, value in intervals:
                            line = f"{chrom}\t{int(start)}\t{int(end)}"
                            interval, points, size = sketch.bedline(line, mode=mode, sep=sep, subsample=subsample)
                            
                            # Add interval if present and in mode A or C
                            if interval and mode in ["A", "C"]:
                                sketch.sketch.add_string(interval)
                                if expA > 0:
                                    for i in range(1, int(10**expA)+1):
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
            print(f"Error processing BigWig file {filename}: {e}")
            return None

    @classmethod
    def _from_bam(cls, filename: str, mode: str, sep: str,
                  subsample: tuple[float, float], expA: float, 
                  sketch: 'IntervalSketch') -> Optional['IntervalSketch']:
        """Process BAM format file.
        
        Args:
            filename: Path to BAM file
            mode: Processing mode (A/B/C)
            sep: Separator for string representation
            subsample: Tuple of (interval_rate, point_rate)
            expA: Exponent for interval multiplicity
            sketch: IntervalSketch instance to populate
            
        Returns:
            Populated IntervalSketch instance or None if error
        """
        try:
            with pysam.AlignmentFile(filename, "rb") as bam:
                for read in bam:
                    if read.is_unmapped:
                        continue
                        
                    # Get chromosome name without 'chr' prefix
                    chrom = read.reference_name
                    if chrom.startswith('chr'):
                        chrom = chrom[3:]
                        
                    # Get alignment coordinates
                    start = read.reference_start
                    end = read.reference_end or (start + 1)  # fallback if end is None
                    
                    # Create BED-style line and process it
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
            print(f"Error processing BAM file {filename}: {e}")
            return None

    @classmethod
    def _from_bed(cls, filename: str, mode: str, sep: str, 
                  subsample: Tuple[float, float], expA: float, sketch: 'IntervalSketch') -> Optional['IntervalSketch']:
        """Process BED format file."""
        try:
            # Open file with gzip if it ends in .gz, otherwise normal open
            opener = gzip.open if filename.endswith('.gz') else open
            with opener(filename, 'rt') as f:  # 'rt' mode for text reading from gzip
                for line in f:
                    interval, points, size = sketch.bedline(line, mode=mode, sep=sep, subsample=subsample)
                    
                    # Add interval if present and in mode A or C
                    if interval and mode in ["A", "C"]:
                        sketch.sketch.add_string(interval)
                        
                        if expA > 0:
                            num_copies = int(10**expA) - 1
                            for i in range(1, num_copies + 1):
                                sketch.sketch.add_string(interval + str(i))
                        
                        sketch.num_intervals += 1
                        sketch.total_interval_size += size
                    
                    # Add points if present and in mode B or C
                    if points and mode in ["B", "C"]:
                        for point in points:
                            if point is not None:
                                sketch.sketch.add_string(point)
            if mode in ["B"]:
                sketch.num_intervals += 1
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
            # print(f"Generated points: {len(points)}")
            # print(f"Start: {start}, End: {end}")
        
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
            