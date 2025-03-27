#!/usr/bin/env python
from __future__ import annotations
from typing import Optional, List, Tuple, TextIO, Union, Dict
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.rusthll import FastHyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter
from hammock.lib.exacttest import ExactTest
import pyBigWig  # type: ignore  # no type stubs available
from numpy import floor  # type: ignore  # no type stubs available
import gzip
import pysam  # type: ignore  # no type stubs available


class IntervalSketch(AbstractSketch):
    """Sketch class for BED intervals."""
    
    def __init__(self, **kwargs):
        """Initialize the sketch."""
        # Extract parameters with defaults
        precision = kwargs.get('precision', 12)
        debug = kwargs.get('debug', False)
        
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
        
        # Only pass precision and debug to HyperLogLog
        self.sketch = FastHyperLogLog(precision=precision, debug=debug)
        
        # Initialize other instance variables
        self.num_intervals = 0
        self.total_interval_size = 0

    @classmethod
    def from_file(cls, filename: str, **kwargs) -> Optional['IntervalSketch']:
        """Create a sketch from a file."""
        # Initialize sketch
        sketch = cls(**kwargs)
        
        # Add sketch to kwargs for helper methods
        kwargs['sketch'] = sketch
        
        # Check file extension and process accordingly
        # Split on '.' and get all extensions (handles multiple extensions like .gff.gz)
        extensions = filename.lower().split('.')
        
        # Handle GFF files (both plain and gzipped)
        if any(ext in ['gff', 'gff3'] for ext in extensions):
            return cls._from_gff(filename, **kwargs)
        
        # NOTE: rearrange this elif to start with the most common file type
        # Handle BigBed files
        elif extensions[-1] in ['bb', 'bigbed']:
            return cls._from_bigbed(filename, **kwargs)
        
        # Skip other unsupported binary formats
        elif extensions[-1] in ['cram', 'sam']:
            print(f"Skipping unsupported file format: {filename}")
            return None
        
        # Handle regular BED files (including gzipped)
        elif extensions[-1] in ["bed"] or (len(extensions) >= 2 and extensions[-2:] == ["bed", "gz"]):
            try:
                # Handle gzipped files
                if filename.endswith('.gz'):
                    opener = gzip.open
                else:
                    opener = open
                
                with opener(filename, 'rt') as f:
                    for line in f:
                        if not line.startswith('#'):
                            interval, points, size = sketch.bedline(line, mode=kwargs.get('mode', 'A'), sep=kwargs.get('sep', '-'), subsample=kwargs.get('subsample', (1.0, 1.0)))
                            
                            # Add interval if present and in mode A or C
                            if interval and kwargs.get('mode', 'A') in ["A", "C"]:
                                sketch.sketch.add_string(interval)
                                if kwargs.get('expA', 0) > 0:
                                    for i in range(1, int(10**kwargs.get('expA', 0))+1):
                                        sketch.sketch.add_string(interval + str(i))
                                sketch.num_intervals += 1
                                sketch.total_interval_size += size
                            
                            # Add points if present and in mode B or C
                            if points and kwargs.get('mode', 'A') in ["B", "C"]:
                                for point in points:
                                    if point is not None:
                                        sketch.sketch.add_string(point)
                                        if kwargs.get('mode', 'A') == "B":
                                            sketch.num_intervals += 1
                return sketch
            except Exception as e:
                print(f"Error processing file {filename}: {e}")
                return None
        else:
            raise ValueError(f"Unsupported file extension: {extensions[-1]}")

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

    @classmethod
    def _from_gff(cls, filename: str, **kwargs) -> Optional['IntervalSketch']:
        """Process GFF format file.
        
        Args:
            filename: Path to GFF file
            **kwargs: Additional arguments including:
                - mode: Sketch mode (A/B/C)
                - sep: Separator for string representation
                - subsample: Tuple of sampling rates
                - expA: Exponential scaling factor
                - feature_types: List of feature types to include (optional)
        """
        if 'sketch' not in kwargs:
            return None
        
        try:
            # Get parameters from kwargs
            mode = kwargs.get('mode', 'A')
            sep = kwargs.get('sep', '-')
            subsample = kwargs.get('subsample', (1.0, 1.0))
            expA = kwargs.get('expA', 0)
            sketch = kwargs['sketch']
            feature_types = kwargs.get('feature_types', None)  # Optional filter for specific features
            
            opener = gzip.open if filename.endswith('.gz') else open
            with opener(filename, 'rt') as f:
                for line in f:
                    # Skip comments and empty lines
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    # Parse GFF fields
                    fields = line.strip().split('\t')
                    if len(fields) != 9:
                        continue
                    
                    seqid, source, ftype, start, end, score, strand, phase, attrs = fields
                    
                    # Skip if feature type filtering is enabled and type doesn't match
                    if feature_types and ftype not in feature_types:
                        continue
                    
                    try:
                        start_pos = int(start) - 1  # Convert to 0-based coordinates
                        end_pos = int(end)
                    except ValueError:
                        continue
                    
                    # Create BED-style line and process using existing bedline method
                    bed_line = f"{seqid}\t{start_pos}\t{end_pos}"
                    interval, points, size = sketch.bedline(
                        bed_line, 
                        mode=mode, 
                        sep=sep, 
                        subsample=subsample
                    )
                    
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
            print(f"Error processing GFF file {filename}: {e}")
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
        # if not (chrval.isnumeric() or chrval in ["X","Y","M"]):
        #     raise ValueError("bedline: chromosome is not numeric")
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
        
    def add(self, s: str) -> None:
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
                sketch = FastHyperLogLog.load(f)
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
            