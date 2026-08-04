"""hammock: pairwise Jaccard similarity for BED intervals (modes A/B/C) and FASTA sequences (mode D)."""

from hammock import _core
from hammock._core import HLLSketch

__all__ = ["HLLSketch", "_core"]
__version__ = "0.6.0"
