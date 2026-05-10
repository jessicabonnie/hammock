"""Hammock CLI: pairwise Jaccard similarity for BED intervals (modes A/B/C) or sequences (mode D).

Output is a tab-separated file whose query/reference columns default to basenames; pass
--full-paths for normalized full paths in those columns.
"""
from __future__ import annotations

import argparse
import os
import sys

from hammock.runner import run


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="""Calculate pairwise Jaccard similarities between lists of BED or sequence files.

        Supported file formats:
        - BED format (.bed): Tab-delimited format with chromosome, start, and end positions
        - BigBed format (.bb): Binary indexed version of BED format
        - Any tab-delimited file with at least 3 columns (chr, start, end) in BED-style format
        - Sequence files (.fa, .fasta, .fna, .ffn, .faa, .frn) - automatically uses mode D

        Output: tab-separated; query and reference columns identify inputs using basenames by default.
        Pass --full-paths to use normalized full paths instead.
        """,
        epilog='Tip: use --full-paths when filenames repeat in different directories or for joining outputs to path-based manifests.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument('filepaths_file',
                   help='Text file containing paths to files to be compared')
    p.add_argument('primary_file',
                   help='Text file containing paths to primary files to compare against')

    p.add_argument('--mode', choices=['A', 'B', 'C', 'D'], default=None,
                   help='''Mode for comparison (auto-detected from inputs if omitted):
                   A: Compare intervals only (default for BED/BigBed files)
                   B: Compare points only
                   C: Compare both intervals and points
                   D: Compare sequences (default for FASTA files)''')

    p.add_argument('--outprefix', '-o', '--out', type=str, default="hammock", help='The output file prefix')
    p.add_argument('--full-paths', action='store_true',
                   help='Write normalized paths in CSV file1/file2 columns instead of basenames.')
    p.add_argument("--precision", "-p", type=int, default=18,
                   help="Precision for HyperLogLog sketching (4..24)")
    p.add_argument("--subA", type=float, default=1.0, help="Subsampling rate for intervals (0..1)")
    p.add_argument("--subB", type=float, default=1.0, help="Subsampling rate for points (0..1)")
    p.add_argument("--expA", type=float, default=0,
                   help="Power-of-10 exponent to multiply contribution of A-type intervals")
    p.add_argument("--threads", type=int, default=None,
                   help="Number of threads. Default: min(8, cpu_count()).")
    p.add_argument("--kmer_size", '-k', type=int, default=8, help="Size of k-mers for sequence sketching")
    p.add_argument("--window_size", '-w', '--window', type=int, default=40,
                   help="Size of sliding window for sequence sketching")
    p.add_argument("--seed", type=int, default=42, help="Random seed for hashing")
    p.add_argument("--verbose", action="store_true",
                   help="Report per-file sketching progress on stderr.")
    p.add_argument("--memory-limit-gb", type=float, default=0.0,
                   help="Soft memory limit in GiB. 0 disables (default).")
    p.add_argument('--mixed-stride', action='store_true',
                   help='Use mixed-stride deterministic subsampling for points (interval-independent)')

    args = p.parse_args(argv)
    # Hardcoded constants the runner still reads. Hash is always xxh64; the
    # CSV `num_hashes` column is "NA" for HLL/minimizer (only meaningful for
    # MinHash, which isn't shipped).
    args.hash_size = 64
    args.num_hashes = "NA"
    return args


_SEQ_EXTS = {".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn"}
_BED_EXTS = {".bed", ".bb", ".bigbed", ".narrowpeak", ".broadpeak", ".gappedpeak"}


def _strip_gz(path: str) -> str:
    return path[:-3] if path.lower().endswith(".gz") else path


def _first_path(list_file: str) -> str | None:
    """Return the first non-blank, non-comment path in a hammock filepaths list."""
    try:
        with open(list_file) as f:
            for line in f:
                s = line.strip()
                if s and not s.startswith('#'):
                    return s
    except OSError:
        return None
    return None


def _classify_file(path: str) -> str | None:
    """Return 'fasta', 'bed', or None if undetermined.

    Tries the file extension first (handling .gz), then falls back to peeking
    at the first non-blank line — FASTA records start with '>'.
    """
    base = _strip_gz(path)
    _, ext = os.path.splitext(base)
    ext = ext.lower()
    if ext in _SEQ_EXTS:
        return 'fasta'
    if ext in _BED_EXTS:
        return 'bed'
    # Unknown extension — peek at content. Open transparently for .gz.
    try:
        if path.lower().endswith('.gz'):
            import gzip
            opener = lambda: gzip.open(path, 'rt')  # noqa: E731
        else:
            opener = lambda: open(path, 'r')  # noqa: E731
        with opener() as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith('#') or s.startswith('track') or s.startswith('browser'):
                    continue
                if s.startswith('>') or s.startswith(';'):
                    return 'fasta'
                # First real line that looks like "chrom\tstart\tend..." → BED.
                fields = s.split('\t')
                if len(fields) >= 3 and fields[1].isdigit() and fields[2].isdigit():
                    return 'bed'
                return None
    except OSError:
        return None
    return None


def _autodetect_mode(args) -> str:
    """Pick a mode when the user didn't supply --mode.

    Strategy:
      * If the user passed --mode explicitly, honor it.
      * Otherwise, classify the first input file (extension, then content
        peek). FASTA → D, BED → A.
      * If subA/subB/expA were tweaked, prefer C over A for BED (the orig's
        Mode C is the natural home for those flags).
      * If we couldn't classify, fall back to A and warn.
    """
    if args.mode is not None:
        return args.mode

    first = _first_path(args.filepaths_file)
    kind = _classify_file(first) if first else None

    if kind == 'fasta':
        print("hammock: auto-detected FASTA input → using --mode D.", file=sys.stderr)
        return 'D'
    if kind == 'bed':
        if args.subA != 1.0 or args.subB != 1.0 or args.expA != 0:
            print("hammock: auto-detected BED input + sub/exp flags → using --mode C.",
                  file=sys.stderr)
            return 'C'
        print("hammock: auto-detected BED input → using --mode A.", file=sys.stderr)
        return 'A'

    # Couldn't tell — keep the original default but make it visible.
    print("hammock: could not auto-detect input type; defaulting to --mode A. "
          "Pass --mode explicitly to silence this warning.", file=sys.stderr)
    return 'A'


def _resolve_sketch_type(args) -> str:
    """Mode D uses the 'minimizer' label (hyperloglog backing under the hood,
    matching the orig hammock's CSV output). Everything else is HLL."""
    return "minimizer" if args.mode == "D" else "hyperloglog"


def _apply_memory_limit(gb: float) -> None:
    if gb <= 0:
        return
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        new_soft = int(gb * 1024 * 1024 * 1024)
        resource.setrlimit(resource.RLIMIT_AS, (new_soft, hard))
        print(f"Set memory limit to {gb} GiB (soft).", file=sys.stderr)
    except (ImportError, ValueError, OSError) as e:
        print(f"Warning: could not set memory limit: {e}", file=sys.stderr)


def main(argv=None) -> int:
    args = parse_args(argv)
    args.mode = _autodetect_mode(args)
    args.sketch_type = _resolve_sketch_type(args)
    if args.threads is None:
        args.threads = min(8, os.cpu_count() or 1)
    _apply_memory_limit(args.memory_limit_gb)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
