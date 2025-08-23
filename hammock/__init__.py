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

# Try importing the native RustHLL implementation
try:
    import rust_hll
    RUST_HLL_AVAILABLE = True
except ImportError:
    RUST_HLL_AVAILABLE = False

if RUST_HLL_AVAILABLE:
    from hammock.lib.rusthll_compat import RustHLLWrapper

__version__ = '0.4.0'

__all__ = [
    'HyperLogLog',
    'MinHash',
    'ExactCounter',
    'SequenceSketch',
    'IntervalSketch',
    'MinimizerSketch',
    'RustHLL'
]
if RUST_HLL_AVAILABLE:
    __all__.append('RustHLLWrapper') 