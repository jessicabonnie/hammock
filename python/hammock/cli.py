"""Hammock CLI: pairwise Jaccard similarity for BED intervals (modes A/B/C) or sequences (mode D).

Output is a tab-separated file whose query/reference columns default to basenames; pass
--full-paths for normalized full paths in those columns.
"""
from __future__ import annotations

import argparse
import os
import subprocess
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
        epilog='BED→FASTA (Mode D): pass --ref/--ref1/--ref2 to treat LIST1/LIST2 as BED '
               'files, convert them to FASTA with bedtools getfasta, and compare the '
               'sequences. A reference is a keyword (hg38, mm10, ...) or a local FASTA path; '
               'cache keywords once on a networked node with `hammock fetch-ref <keyword>`. '
               'Tip: use --full-paths when filenames repeat in different directories.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument('filepaths_file', metavar='LIST1',
                   help='Text file of paths to compare. Normally FASTA/BED files; '
                        'with a reference flag (--ref/--ref1) these are BED files '
                        'converted to FASTA via bedtools getfasta.')
    p.add_argument('primary_file', metavar='LIST2',
                   help='Text file of primary paths to compare against (see LIST1; '
                        'uses --ref/--ref2 in BED→FASTA mode).')

    p.add_argument('--mode', choices=['A', 'B', 'C', 'D'], default=None,
                   help='''Mode for comparison (auto-detected from inputs if omitted):
                   A: Compare intervals only (default for BED/BigBed files)
                   B: Compare points only
                   C: Compare both intervals and points
                   D: Compare sequences (default for FASTA files)''')

    # BED→FASTA (bed2fasta) reference flags. Presence of any of these turns the
    # two positional lists into BED lists that are converted to FASTA (Mode D).
    ref = p.add_argument_group('BED→FASTA references (Mode D)')
    ref.add_argument('--ref', default=None,
                     help='Reference for BOTH lists: a keyword (hg38, mm10, ...) '
                          'or a local FASTA path. Mutually exclusive with --ref1/--ref2.')
    ref.add_argument('--ref1', default=None,
                     help='Reference for LIST1 (keyword or local FASTA path).')
    ref.add_argument('--ref2', default=None,
                     help='Reference for LIST2 (keyword or local FASTA path).')
    ref.add_argument('--ref-cache-dir', default=None,
                     help='Directory of cached/indexed references (default: '
                          '$HAMMOCK_REF_CACHE or ~/.hammock/refs). Populate '
                          'keywords once with `hammock fetch-ref`.')
    ref.add_argument('--fasta-outdir', default=None,
                     help='Keep the generated FASTA files here instead of a '
                          'temp dir (which is auto-cleaned).')

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
    p.add_argument("--seed", type=int, default=42,
                   help="HLL ingestion seed (xxh64). Default 42.")
    p.add_argument("--gate-seed", type=int, default=31337,
                   help="Seed for the subB gate hash (xxh32) and the mixed-stride "
                        "chr->stride hash. Default 31337 matches orig hammock. "
                        "Ignored when --subB-method=single-hash (gate IS the HLL hash).")
    p.add_argument("--verbose", action="store_true",
                   help="Report per-file sketching progress on stderr.")
    p.add_argument("--memory-limit-gb", type=float, default=0.0,
                   help="Soft memory limit in GiB. 0 disables (default).")
    p.add_argument('--subB-method',
                   choices=['hash-threshold', 'mixed-stride', 'single-hash'],
                   default='mixed-stride',
                   help='subB point-sampling method (default: mixed-stride). '
                        'mixed-stride: deterministic chr-keyed stride, fastest at low subB. '
                        'hash-threshold: random gate, xxh32 seed=--gate-seed; '
                        'matches orig hammock parity. '
                        'single-hash: one xxh64 for gate+ingestion; opt-in parity divergence.')

    args = p.parse_args(argv)

    # ---- BED→FASTA reference validation -------------------------------------
    if args.ref is not None and (args.ref1 is not None or args.ref2 is not None):
        p.error("--ref is mutually exclusive with --ref1/--ref2.")
    if args.ref is not None:
        args.ref1 = args.ref2 = args.ref

    args.bed2fasta = args.ref1 is not None or args.ref2 is not None
    if args.bed2fasta:
        if args.ref1 is None or args.ref2 is None:
            p.error("BED→FASTA mode needs a reference for both lists: pass --ref "
                    "(same reference) or both --ref1 and --ref2.")
        if args.mode not in (None, "D"):
            p.error(f"reference flags imply Mode D (BED→FASTA), but --mode "
                    f"{args.mode} was given. Drop --mode (or use --mode D).")
        args.mode = "D"
    elif args.ref_cache_dir is not None or args.fasta_outdir is not None:
        p.error("--ref-cache-dir/--fasta-outdir only apply with a reference "
                "flag (--ref/--ref1/--ref2).")

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


def _fetch_ref_main(argv) -> int:
    """`hammock fetch-ref <keyword|url>`: download + index a reference into the
    cache. Meant to be run once on a networked (login) node — compute nodes are
    typically firewalled, so comparison runs never download."""
    fp = argparse.ArgumentParser(
        prog="hammock fetch-ref",
        description="Download and index a reference genome into the cache.")
    fp.add_argument("spec", help="Reference keyword (hg38, mm10, hg19, mm39, hs1) "
                                 "or an http(s) URL to a .fa.gz")
    fp.add_argument("--ref-cache-dir", default=None,
                    help="Cache directory (default: $HAMMOCK_REF_CACHE or ~/.hammock/refs)")
    fp.add_argument("--force", action="store_true",
                    help="Re-fetch even if already cached")
    a = fp.parse_args(argv)
    from hammock import refs as refs_mod
    try:
        path = refs_mod.fetch_reference(a.spec, a.ref_cache_dir, force=a.force)
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    print(path)
    return 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "fetch-ref":
        return _fetch_ref_main(argv[1:])
    args = parse_args(argv)
    args.mode = _autodetect_mode(args)
    args.sketch_type = _resolve_sketch_type(args)
    if args.threads is None:
        args.threads = min(8, os.cpu_count() or 1)
    _apply_memory_limit(args.memory_limit_gb)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
