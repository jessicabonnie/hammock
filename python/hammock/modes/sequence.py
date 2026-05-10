"""Mode D: FASTA sequence sketching.

Implements `MinimizerSketch` matching hammock/lib/minimizer.py's algorithm
exactly. The minimizer hashes come from `digest.window_minimizer`; we then
reproduce the original's hash → str → re-hash-as-kmers idiom (orig
hyperloglog.py:160-181) so the sketch state is byte-identical.

The original uses xxh64(seed=42) per-kmer; we route every per-kmer hash
through `xxhash.xxh64` here and call `_core.HLLSketch.add_hash64`, so the
HLL register state is exactly what `_process_kmer` would produce.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import xxhash

from hammock import _core

try:
    # bioconda / VeryAmazed source build publishes as lowercase `digest`
    from digest import window_minimizer  # type: ignore
    _DIGEST_AVAILABLE = True
except ImportError:
    try:
        # local source builds may produce a capital-D module name
        from Digest import window_minimizer  # type: ignore
        _DIGEST_AVAILABLE = True
    except ImportError:
        _DIGEST_AVAILABLE = False

        def window_minimizer(seq, **kwargs):  # type: ignore
            return []


_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def canonicalize_kmer(kmer: str) -> str:
    """Lex-min of (uppercase kmer, reverse complement). Non-ACGT bases pass through.

    Matches hammock/lib/minimizer.py:21-43.
    """
    if not kmer:
        return kmer
    kmer = kmer.upper()
    rc = ''.join(_COMPLEMENT.get(b, b) for b in kmer[::-1])
    return min(kmer, rc)


def _xxh64(s: str, seed: int) -> int:
    h = xxhash.xxh64(seed=seed)
    h.update(s.encode())
    return h.intdigest()


class MinimizerSketch:
    """Two-HLL sequence sketch: minimizer-hashes + canonicalized start/end kmers.

    Faithful port of hammock/lib/minimizer.py:45-220 (the FastHyperLogLog path,
    hash_size=64, seed=42 by default).
    """

    def __init__(
        self,
        kmer_size: int = 8,
        window_size: int = 40,
        seed: int = 42,
        precision: int = 16,
    ) -> None:
        self.kmer_size = int(kmer_size)
        self.window_size = int(window_size)
        self.seed = int(seed)
        self.precision = int(precision)
        self.minimizer_hll = _core.HLLSketch(precision=self.precision)
        self.startend_hll = _core.HLLSketch(precision=self.precision)

    def _add_kmers_to(self, hll: "_core.HLLSketch", s: str) -> None:
        """Replicates HyperLogLog.add_string: skip if shorter than k; else hash
        every length-k substring via xxh64(seed) and add to `hll`. If kmer_size
        is 0, hash the whole string as one element.
        """
        k = self.kmer_size
        if k == 0:
            hll.add_hash64(_xxh64(s, self.seed))
            return
        if len(s) < k:
            return
        for i in range(len(s) - k + 1):
            hll.add_hash64(_xxh64(s[i:i + k], self.seed))

    def _process_kmer_to(self, hll: "_core.HLLSketch", kmer: str) -> None:
        """Replicates HyperLogLog._process_kmer: hash kmer once, no length check."""
        hll.add_hash64(_xxh64(kmer, self.seed))

    def add_string(self, s: str) -> None:
        """Mirror hammock/lib/minimizer.py:83-135."""
        if not s:
            return

        try:
            minimizers: List[Tuple[int, int]] = window_minimizer(
                s, k=self.kmer_size, w=self.window_size, include_hash=True,
            )
        except Exception:
            minimizers = []

        if not minimizers:
            # Empty fallback: feed whole sequence as a single kmer to BOTH HLLs
            # (orig minimizer.py:106-107 calls _process_kmer, bypassing length check).
            self._process_kmer_to(self.minimizer_hll, s)
            self._process_kmer_to(self.startend_hll, s)
            return

        # Each (_, hash_val) → str(hash_val) re-hashed as kmer-iterated string.
        # See orig minimizer.py:114-121 ("re-hashes via xxhash64 to get proper
        # bit distribution").
        for _, hash_val in minimizers:
            self._add_kmers_to(self.minimizer_hll, str(hash_val))

        # Canonicalized start+end (or whole seq if too short).
        if len(s) >= 2 * self.kmer_size:
            payload = canonicalize_kmer(s[:self.kmer_size] + s[-self.kmer_size:])
        else:
            payload = canonicalize_kmer(s)
        self._add_kmers_to(self.startend_hll, payload)

    def merged(self) -> "_core.HLLSketch":
        """Return minimizer_hll ∪ startend_hll (used for jaccard_with_ends)."""
        return self.minimizer_hll.merge_new(self.startend_hll)

    def similarity_values(self, other: "MinimizerSketch") -> Dict[str, float]:
        if (self.kmer_size != other.kmer_size
                or self.window_size != other.window_size):
            raise ValueError("Cannot compare sketches with different parameters")
        jac = self.minimizer_hll.estimate_jaccard(other.minimizer_hll)
        jac_ends = self.merged().estimate_jaccard(other.merged())
        return {
            'jaccard_similarity': jac,
            'jaccard_similarity_with_ends': jac_ends,
        }


def sketch_fasta(path: str, args) -> MinimizerSketch:
    """Read FASTA via BioPython and accumulate one MinimizerSketch over all records."""
    from Bio import SeqIO  # local import: BioPython is only needed for Mode D

    sketch = MinimizerSketch(
        kmer_size=args.kmer_size,
        window_size=args.window_size,
        seed=args.seed,
        precision=args.precision,
    )
    if path.endswith('.gz'):
        import gzip
        with gzip.open(path, 'rt') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                sketch.add_string(str(record.seq))
    else:
        for record in SeqIO.parse(path, 'fasta'):
            sketch.add_string(str(record.seq))
    return sketch
