"""hammock: pairwise Jaccard similarity for BED intervals (modes A/B/C) and FASTA sequences (mode D)."""

from hammock import _core
from hammock._core import HLLSketch, BagMinHashSketch

__all__ = ["HLLSketch", "BagMinHashSketch", "_core"]
__version__ = "0.4.0"
