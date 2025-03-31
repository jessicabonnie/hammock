"""
hammock - Python Library for Cardinality Estimation and Interval Sketches
"""

from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter
from hammock.lib.sequences import SequenceSketch
from hammock.lib.intervals import IntervalSketch
from hammock.lib.minimizer import MinimizerSketch
from hammock.lib.rusthll import RustHLL
from hammock.lib.abstractsketch import AbstractSketch

# Try importing the native RustHLL implementation
try:
    import rust_hll
    RUST_HLL_AVAILABLE = True
except ImportError:
    RUST_HLL_AVAILABLE = False

if RUST_HLL_AVAILABLE:
    from hammock.lib.rusthll_compat import RustHLLWrapper

__version__ = '0.2.0'

__all__ = [
    'HyperLogLog',
    'MinHash',
    'ExactCounter',
    'SequenceSketch',
    'IntervalSketch',
    'MinimizerSketch',
    'RustHLL',
    'create_sketch',
    'load_sketch',
    'write_sketch',
    'compare_files',
    'compare_bed_files',
    'compare_sequence_files'
]
if RUST_HLL_AVAILABLE:
    __all__.append('RustHLLWrapper')

def create_sketch(sketch_type: str = "hyperloglog", **kwargs) -> AbstractSketch:
    """Create a new sketch of the specified type.
    
    Args:
        sketch_type: Type of sketch to create ('hyperloglog', 'minhash', 'minimizer', or 'exact')
        **kwargs: Additional arguments passed to the sketch constructor
        
    Returns:
        A new sketch instance of the specified type
    """
    if sketch_type == "hyperloglog":
        return HyperLogLog(**kwargs)
    elif sketch_type == "minhash":
        return MinHash(**kwargs)
    elif sketch_type == "minimizer":
        return MinimizerSketch(**kwargs)
    elif sketch_type == "exact":
        return ExactCounter(**kwargs)
    else:
        raise ValueError(f"Unknown sketch type: {sketch_type}")

def write_sketch(sketch: AbstractSketch, filepath: str) -> None:
    """Write a sketch to a file.
    
    Args:
        sketch: The sketch to write
        filepath: Path to write the sketch to
    """
    sketch.write(filepath)

def load_sketch(filepath: str) -> AbstractSketch:
    """Load a sketch from a file.
    
    Args:
        filepath: Path to load the sketch from
        
    Returns:
        The loaded sketch
    """
    # Determine sketch type from file extension or content
    if filepath.endswith('.hll'):
        return HyperLogLog.load(filepath)
    elif filepath.endswith('.mh'):
        return MinHash.load(filepath)
    elif filepath.endswith('.mnmzr'):
        return MinimizerSketch.load(filepath)
    elif filepath.endswith('.exact'):
        return ExactCounter.load(filepath)
    else:
        raise ValueError(f"Unknown sketch file type: {filepath}")

def compare_files(file1: str, file2: str, mode: str = 'A', 
                 sketch_type: str = "hyperloglog",
                 precision: int = 12,
                 num_hashes: int = 64,
                 kmer_size: int = 8,
                 window_size: int = 40,
                 subA: float = 1.0,
                 subB: float = 1.0,
                 expA: float = 0.5,
                 use_rust: bool = True) -> dict:
    """Compare two files and return similarity metrics.
    
    Args:
        file1: Path to first file
        file2: Path to second file
        mode: Comparison mode ('A', 'B', 'C', or 'D')
        sketch_type: Type of sketch to use ('hyperloglog', 'minhash', 'minimizer', or 'exact')
        precision: Precision for HyperLogLog sketching
        num_hashes: Number of hashes for MinHash sketching
        kmer_size: Size of k-mers for sequence sketching
        window_size: Size of sliding window for sequence sketching
        subA: Subsampling rate for intervals (0 to 1)
        subB: Subsampling rate for points (0 to 1)
        expA: Power of 10 exponent for A-type intervals
        use_rust: Whether to use Rust implementation for HyperLogLog
        
    Returns:
        Dictionary containing similarity metrics
    """
    if mode == 'D':
        sketch1 = SequenceSketch.from_file(
            filename=file1,
            sketch_type=sketch_type,
            kmer_size=kmer_size,
            window_size=window_size,
            precision=precision,
            num_hashes=num_hashes
        )
        sketch2 = SequenceSketch.from_file(
            filename=file2,
            sketch_type=sketch_type,
            kmer_size=kmer_size,
            window_size=window_size,
            precision=precision,
            num_hashes=num_hashes
        )
    else:
        sketch1 = IntervalSketch.from_file(
            filename=file1,
            mode=mode,
            precision=precision,
            num_hashes=num_hashes,
            kmer_size=kmer_size,
            window_size=window_size,
            sketch_type=sketch_type,
            subsample=(subA, subB),
            expA=expA,
            use_rust=use_rust
        )
        sketch2 = IntervalSketch.from_file(
            filename=file2,
            mode=mode,
            precision=precision,
            num_hashes=num_hashes,
            kmer_size=kmer_size,
            window_size=window_size,
            sketch_type=sketch_type,
            subsample=(subA, subB),
            expA=expA,
            use_rust=use_rust
        )
    
    return sketch1.similarity_values(sketch2)

def compare_bed_files(file1: str, file2: str, 
                     sketch_type: str = "hyperloglog",
                     precision: int = 12,
                     num_hashes: int = 64,
                     subA: float = 1.0,
                     subB: float = 1.0,
                     expA: float = 0.5,
                     use_rust: bool = True) -> dict:
    """Compare two BED files and return similarity metrics.
    
    This is a convenience wrapper around compare_files() specifically for BED files.
    It uses mode 'A' by default for interval comparison.
    
    Args:
        file1: Path to first BED file
        file2: Path to second BED file
        sketch_type: Type of sketch to use ('hyperloglog', 'minhash', 'minimizer', or 'exact')
        precision: Precision for HyperLogLog sketching
        num_hashes: Number of hashes for MinHash sketching
        subA: Subsampling rate for intervals (0 to 1)
        subB: Subsampling rate for points (0 to 1)
        expA: Power of 10 exponent for A-type intervals
        use_rust: Whether to use Rust implementation for HyperLogLog
        
    Returns:
        Dictionary containing similarity metrics
    """
    return compare_files(
        file1=file1,
        file2=file2,
        mode='A',
        sketch_type=sketch_type,
        precision=precision,
        num_hashes=num_hashes,
        subA=subA,
        subB=subB,
        expA=expA,
        use_rust=use_rust
    )

def compare_sequence_files(file1: str, file2: str,
                         sketch_type: str = "hyperloglog",
                         precision: int = 12,
                         num_hashes: int = 64,
                         kmer_size: int = 8,
                         window_size: int = 40) -> dict:
    """Compare two sequence files and return similarity metrics.
    
    This is a convenience wrapper around compare_files() specifically for sequence files.
    It uses mode 'D' by default for sequence comparison.
    
    Args:
        file1: Path to first sequence file
        file2: Path to second sequence file
        sketch_type: Type of sketch to use ('hyperloglog', 'minhash', 'minimizer', or 'exact')
        precision: Precision for HyperLogLog sketching
        num_hashes: Number of hashes for MinHash sketching
        kmer_size: Size of k-mers for sequence sketching
        window_size: Size of sliding window for sequence sketching
        
    Returns:
        Dictionary containing similarity metrics
    """
    return compare_files(
        file1=file1,
        file2=file2,
        mode='D',
        sketch_type=sketch_type,
        precision=precision,
        num_hashes=num_hashes,
        kmer_size=kmer_size,
        window_size=window_size
    ) 