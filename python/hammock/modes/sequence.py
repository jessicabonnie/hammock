"""Mode D: FASTA sequence sketching.

`MinimizerSketch` holds one HLL over a FASTA's (k, w) window minimizers.
Selector identities come from `digest.window_minimizer`. By default they are
domain-separated and rehashed to uniform 64-bit values before
`_core.HLLSketch.add_hash64`; `legacy-selector32` feeds the raw selector hash
for compatibility. Neither path uses orig's hash → str → re-hash-as-kmers
idiom, which silently dropped most minimizers (see divergence #6 in
CLAUDE.md). Mode D is therefore deliberately **not** byte-identical to orig.

`xxhash.xxh64(seed)` survives on exactly one path in the legacy mode: the no-minimizer fallback
for records shorter than `k + w - 1`, which hashes the whole record as a single
element. That makes such records exact-match indicators rather than graded
similarity — one substitution and they stop matching.

Through v0.5.0 this class also carried a second HLL of per-record start/end
k-mers, backing the `jaccard_similarity_with_ends` column family. That was
removed in v0.6.0; see divergence #8 and `docs/mode-d-ends-removal.md`.
"""
from __future__ import annotations

from typing import List, Tuple

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


def _xxh64(s: str, seed: int) -> int:
    h = xxhash.xxh64(seed=seed)
    h.update(s.encode())
    return h.intdigest()


# This is deliberately bytes, rather than a decimal rendering of the selector
# value: its width and byte order are part of the feature-hash contract.
_SELECTOR64_DOMAIN = b"hammock:sequence:minimizer-selector:v1\0"
_WHOLE_RECORD64_DOMAIN = b"hammock:sequence:whole-record:v1\0"


def _rehash_selector64(selector_hash: int, seed: int) -> int:
    """Domain-separate and uniformly rehash a digest 32-bit selector hash.

    The selector remains the feature identity. A
    collision already present in digest's 32-bit selector hash is therefore
    intentionally retained; this function only prevents its order-statistic
    distribution from being used as an HLL hash.
    """
    value = int(selector_hash)
    if not 0 <= value < 2**32:
        raise ValueError(f"digest selector hash is not uint32: {value}")
    h = xxhash.xxh64(seed=seed)
    h.update(_SELECTOR64_DOMAIN)
    h.update(value.to_bytes(4, byteorder="little", signed=False))
    return h.intdigest()


def _hash_whole_record64(sequence: str, seed: int, *, domain_separated: bool) -> int:
    """Hash a sub-threshold record, optionally in its own feature domain."""
    if not domain_separated:
        return _xxh64(sequence, seed)
    h = xxhash.xxh64(seed=seed)
    h.update(_WHOLE_RECORD64_DOMAIN)
    h.update(sequence.upper().encode())
    return h.intdigest()


class MinimizerSketch:
    """Single-HLL sequence sketch over a FASTA's (k, w) window minimizers.

    Port of hammock/lib/minimizer.py:45-220 (the FastHyperLogLog path,
    hash_size=64, seed=42 by default), minus the start/end ("ends") sketch —
    see divergence #8 in CLAUDE.md for why that was dropped.
    """

    def __init__(
        self,
        kmer_size: int = 8,
        window_size: int = 40,
        seed: int = 42,
        precision: int = 16,
        hash_mode: str = "rehash-selector64",
    ) -> None:
        self.kmer_size = int(kmer_size)
        self.window_size = int(window_size)
        self.seed = int(seed)
        self.precision = int(precision)
        if hash_mode not in ("legacy-selector32", "rehash-selector64"):
            raise ValueError(f"unsupported sequence HLL hash mode: {hash_mode}")
        self.hash_mode = hash_mode
        self.minimizer_hll = _core.HLLSketch(precision=self.precision)

    def _process_kmer_to(self, hll: "_core.HLLSketch", kmer: str) -> None:
        """Replicates HyperLogLog._process_kmer: hash kmer once, no length check."""
        hll.add_hash64(_hash_whole_record64(
            kmer, self.seed,
            domain_separated=self.hash_mode == "rehash-selector64",
        ))

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
            # Empty fallback: feed whole sequence as a single kmer
            # (orig minimizer.py:106-107 calls _process_kmer, bypassing length check).
            # Note this makes a sub-threshold record (len < k+w-1) an exact-match
            # indicator: it matches only a byte-identical record, so one SNP drops it.
            self._process_kmer_to(self.minimizer_hll, s)
            return

        # Feed each minimizer's selector hash directly to the HLL in legacy
        # mode. This matches
        # orig minimizer.py:118-119 (FAST_HLL path: `add_hash64(np.uint64(hash_val))`).
        # The prior refactor ported only orig's str(hash_val) slow fallback, which
        # silently dropped any minimizer whose decimal representation was shorter
        # than k (e.g. k=15 vs typical minimizer hash ~10 digits), leaving the HLL
        # empty and J=0 on random-ACGT FASTAs. See memory note
        # [[project_modeD_no_ends_zero_bug]] for full repro.
        for _, hash_val in minimizers:
            if self.hash_mode == "legacy-selector32":
                self.minimizer_hll.add_hash64(hash_val)
            else:
                self.minimizer_hll.add_hash64(
                    _rehash_selector64(hash_val, self.seed))


def sketch_fasta(path: str, args) -> MinimizerSketch:
    """Read FASTA via BioPython and accumulate one MinimizerSketch over all records."""
    from Bio import SeqIO  # local import: BioPython is only needed for Mode D

    if not _DIGEST_AVAILABLE:
        # Mode D is meaningless without `digest.window_minimizer`: every sequence
        # would fall through add_string's empty-minimizer path and be hashed whole,
        # yielding self=1.0 / cross=0.0 for ALL pairs with no error. Fail loud
        # instead of emitting silent garbage. A common cause on this cluster is a
        # stale RPATH in _core.*.so (built with an env-module gcc) shadowing
        # libstdc++ so `import digest` fails with GLIBCXX_3.4.32 not found — see
        # memory note project_modeD_zero_rpath_digest.
        raise RuntimeError(
            "Mode D (sequence) requires the `digest` module, but it failed to "
            "import. Every sketch would be whole-sequence-hashed, giving "
            "jaccard=1.0 self / 0.0 cross for all pairs. Verify `python -c "
            "'import digest'` works; if it fails with 'GLIBCXX_3.4.32 not found', "
            "your hammock _core extension has a stale RPATH pulling in an old "
            "libstdc++ — rebuild _core without a gcc env-module loaded, or "
            "patchelf --force-rpath the conda lib dir first."
        )

    sketch = MinimizerSketch(
        kmer_size=args.kmer_size,
        window_size=args.window_size,
        seed=args.seed,
        precision=args.precision,
        hash_mode=getattr(args, "sequence_hll_hash", "rehash-selector64"),
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
