"""Hammock CLI: pairwise Jaccard similarity in interval modes (A/B/C, BED) or sequence mode (D, FASTA).

Output is a comma-separated CSV whose file1/file2 columns default to basenames; pass
--full-paths for normalized full paths in those columns. Three output shapes, chosen by
mutually exclusive flags: bare default (jaccard_similarity_ie only, tag _ie),
--register-equality/--re (register-equality-only, tag _re), --metrics (full 8-column
block, tag _full). See --metrics/--register-equality help.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

from hammock import __version__
from hammock.runner import run


# User-facing mode names → canonical single-letter code used everywhere
# downstream (and in the CSV `mode` column / output filename, for orig parity).
# Top-level choice is `interval` (=B) vs `sequence` (=D); A/B/C are the three
# interval flavors (string / points / hybrid).
_MODE_ALIASES = {
    "a": "A", "interval-string": "A", "interval_string": "A",
    "b": "B", "interval": "B", "interval-points": "B", "interval_points": "B",
    "c": "C", "interval-hybrid": "C", "interval_hybrid": "C",
    "d": "D", "sequence": "D",
}


def _normalize_mode(value: str) -> str:
    """argparse `type` for --mode: accept a name or a letter, return the letter."""
    key = value.strip().lower()
    if key in _MODE_ALIASES:
        return _MODE_ALIASES[key]
    raise argparse.ArgumentTypeError(
        f"invalid mode '{value}'. Use 'interval' (=B, default for BED) or "
        f"'sequence' (=D, default for FASTA/--ref). Advanced interval flavors: "
        f"interval-string (A), interval-points (B), interval-hybrid (C). "
        f"The letters A/B/C/D are also accepted.")


def _precision(value: str) -> int:
    """argparse `type` for --precision: an int in 4..24.

    Must return a plain `int` -- the value reaches the output filename
    (`outprefix.get_new_prefix`) and the CSV `precision` column, so anything
    that renders differently would change filenames and break parity CSVs.

    Without the bound this reached HLLSketch's constructor unchecked: -p 64
    returned silently wrong estimates and -p 65 underflowed max-rho to
    SIZE_MAX, hanging on a ~1.8e19-iteration loop.
    """
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid int value: '{value}'")
    if not 4 <= n <= 24:
        raise argparse.ArgumentTypeError(
            f"precision must be in 4..24 (got {n})")
    return n


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="""Calculate pairwise Jaccard similarities between lists of BED or sequence files.

        LIST1 and LIST2 are text files listing one input path per line -- not inputs themselves.

        Supported inputs:
        - Interval modes (interval / interval-string / interval-hybrid, i.e. A/B/C):
          plain-text BED, or any tab-delimited file whose first three columns are
          chrom, start, end. Nothing is decompressed or decoded first: a .gz or
          binary BigBed input is detected by magic bytes and rejected (exit 2)
          rather than silently sketching to nothing. Decompress/convert first.
        - Sequence mode (D): FASTA (.fa, .fasta, .fna, .ffn, .faa, .frn), plain or .gz.
          Selected automatically for FASTA input or for any --ref/--ref1/--ref2 flag.

        Output: comma-separated CSV at <outprefix>_<sketch>_p<precision>_jacc<MODE>[..._ie|_re|_full].csv;
        the file1/file2 columns identify inputs using basenames by default.
        Pass --full-paths to use normalized full paths instead. The trailing
        _ie/_re/_full tag names the output shape (see --metrics/--register-equality).

        Sub-command: `hammock fetch-ref <keyword|url>` caches a reference genome for
        BED->FASTA runs; see `hammock fetch-ref --help`.
        """,
        epilog='BED→FASTA (sequence mode): pass --ref/--ref1/--ref2 to treat LIST1/LIST2 as BED '
               'files, convert them to FASTA with bedtools getfasta, and compare the '
               'sequences. Any reference flag forces sequence mode, so combining one with '
               '--mode A/B/C is an error. A reference is a local FASTA path, a keyword '
               '(hg38, hg19, mm10, mm39, hs1, plus aliases such as grch38/grch37/chm13/t2t), '
               'or an http(s) URL; keywords and URLs must already be cached, because a '
               'comparison run never downloads -- cache them once on a networked node with '
               '`hammock fetch-ref <keyword>`. Requires bedtools (and samtools for indexing) '
               'on PATH. '
               'Tip: use --full-paths when filenames repeat in different directories.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Same string the standalone binary reports, from the same source, so a
    # harness can tell which pair of front-ends it is holding.
    p.add_argument('--version', action='version',
                   version=f'hammock {__version__}')
    p.add_argument('filepaths_file', metavar='LIST1',
                   help='Text file listing one input path per line (NOT an input '
                        'file itself). This is the query side (side "A": the file1 '
                        'and containment_AB columns; paired with --ref1). Listed '
                        'paths are plain-text BED in interval modes (.gz and BigBed '
                        'are NOT decoded; they are detected and rejected) or FASTA '
                        '(.gz fine) in sequence mode; with '
                        'a reference flag (--ref/--ref1) they must be BED and are '
                        'converted to FASTA via bedtools getfasta.')
    p.add_argument('primary_file', metavar='LIST2',
                   help='Text file listing one input path per line, compared against '
                        'LIST1. This is the reference side (side "B": the file2 and '
                        'containment_BA columns; paired with --ref2). Same format rules '
                        'as LIST1. (May be the same file as LIST1 for all-vs-all.)')

    p.add_argument('--mode', type=_normalize_mode, default=None, metavar='MODE',
                   help='Comparison mode; auto-detected from the first path in LIST1 '
                        'when omitted. Primary choice: "interval" (= interval-points '
                        '= B) sketches every base position covered by a BED interval, '
                        'and is the default for BED input; "sequence" (= D) sketches '
                        'FASTA window minimizers, and is the default for FASTA input '
                        'or any --ref flag. Advanced interval flavors: '
                        '"interval-string" (A) hashes each whole interval record, so '
                        'only bit-identical intervals match; "interval-hybrid" (C) '
                        'sketches interval records and base positions together and is '
                        'the only mode that reads --subA/--expA. Autodetect picks C '
                        'over B when BED input is combined with --subA or --expA '
                        '(--subB alone stays B), and falls back to B with a stderr '
                        'warning when the input type cannot be determined. The bare '
                        'letters A/B/C/D are accepted anywhere a name is, and the '
                        'letter is what appears in the CSV mode column and in the '
                        'output filename.')

    # BED→FASTA (bed2fasta) reference flags. Presence of any of these turns the
    # two positional lists into BED lists that are converted to FASTA (Mode D).
    ref = p.add_argument_group('BED→FASTA references (Mode D)')
    ref.add_argument('--ref', default=None,
                     help='Reference genome for BOTH lists. Makes this a BED→FASTA '
                          'run: sequence mode is forced and LIST1/LIST2 must name BED '
                          'files. Accepts a local FASTA path, a built-in keyword '
                          '(hg38, hg19, mm10, mm39, hs1; aliases grch38, grch37, '
                          'chm13, t2t, ...), or an http(s) URL. Keywords and URLs must '
                          'already be in the cache — a comparison run never downloads, '
                          'use `hammock fetch-ref`. Mutually exclusive with '
                          '--ref1/--ref2.')
    ref.add_argument('--ref1', default=None,
                     help='Reference for LIST1 only; same accepted forms as --ref. '
                          'Pass with --ref2 for a cross-reference run '
                          '(e.g. --ref1 hg38 --ref2 mm10); both are required together.')
    ref.add_argument('--ref2', default=None,
                     help='Reference for LIST2 only; same accepted forms as --ref. '
                          'Pass with --ref1 (see above).')
    ref.add_argument('--ref-cache-dir', default=None,
                     help='Directory of cached/indexed references (default: '
                          '$HAMMOCK_REF_CACHE, else ~/.hammock/refs). Populate '
                          'keywords once with `hammock fetch-ref`. Errors out unless '
                          'a reference flag is also given.')
    ref.add_argument('--fasta-outdir', default=None,
                     help='Write the generated FASTA files here and keep them '
                          '(default: a temp dir deleted when the run ends). Errors '
                          'out unless a reference flag is also given.')

    p.add_argument('--outprefix', '-o', '--out', type=str, default="hammock",
                   help='Output filename prefix (default: hammock). The run appends '
                        'its own suffixes plus ".csv": sketch type and precision and '
                        'mode letter, then _A<subA>/_expA<expA> and _B<subB> for any '
                        'subsampling, _k<k>_w<w> in sequence mode, '
                        '_<ref1>-vs-<ref2> for a BED→FASTA run, and finally '
                        '_ie/_re/_full for the output shape (see --metrics/'
                        '--register-equality) — e.g. the bare default writes '
                        'hammock_hll_p18_jaccB_ie.csv.')
    p.add_argument('--full-paths', action='store_true',
                   help='Write normalized paths in CSV file1/file2 columns instead of basenames.')
    p.add_argument("--precision", "-p", type=_precision, default=18,
                   help="HyperLogLog precision: each sketch gets 2**p one-byte "
                        "registers, p in 4..24 (default: 18, i.e. 256 KiB/sketch). "
                        "Raising p cuts estimator error but costs memory and pairwise "
                        "time proportional to 2**p.")
    p.add_argument("--subA", type=float, default=1.0,
                   help="interval-hybrid (C) only: fraction of interval records to "
                        "sketch as whole-interval elements, 0..1 (default: 1.0, keep "
                        "all). An interval failing the gate (xxh32 of the interval "
                        "string, seeded by --gate-seed) contributes no interval "
                        "element, but its base positions are still sketched. Passing "
                        "--subA with BED input makes autodetect choose mode C.")
    p.add_argument("--subB", type=float, default=1.0,
                   help="interval / interval-hybrid (B and C): fraction of base "
                        "positions to sketch, 0..1 (default: 1.0, every base). Which "
                        "positions survive is decided by --subB-method. Note orig "
                        "hammock silently ignored --subB in mode B; here it is "
                        "honored, so such runs are not comparable to orig.")
    p.add_argument("--expA", type=float, default=0,
                   help="interval-hybrid (C) only: give each kept interval record "
                        "10**expA distinct sketch elements instead of one "
                        "(default: 0, i.e. no expansion), up-weighting whole intervals "
                        "against base positions. A nonzero --expA with BED input makes "
                        "autodetect choose mode C, and takes the place of --subA in "
                        "the output filename.")
    p.add_argument("--threads", type=int, default=None,
                   help="Threads for the sketching phase. Default: 1 in sequence mode "
                        "(D), whose per-record loop is Python and gets slower under a "
                        "pool, and min(8, cpu_count()) in interval modes, which sketch "
                        "in C++ with the GIL released. An explicit value >1 in "
                        "sequence mode is honored with a stderr warning. In a BED→FASTA "
                        "run this value also sizes the bedtools getfasta pool; left "
                        "unset, that extraction pool still uses min(8, cpu_count()) "
                        "even though sketching stays single-threaded.")
    p.add_argument("--kmer_size", '-k', type=int, default=8,
                   help="Sequence mode (D) only: minimizer k-mer size (default: 8). "
                        "Reported in the CSV kmer_size column in every mode, and in "
                        "the output filename in sequence mode. Larger k is advisable "
                        "for cross-species comparisons.")
    p.add_argument("--window_size", '-w', '--window', type=int, default=40,
                   help="Sequence mode (D) only: minimizer window size (default: 40). "
                        "A record shorter than k + w - 1 yields no minimizer and is "
                        "hashed whole, so it matches only an identical record.")
    p.add_argument("--seed", type=int, default=42,
                   help="Seed for the xxh64 hash ingested by the HLL (default: 42). "
                        "Changes every sketch but not the expected estimate; useful "
                        "for re-rolling estimator noise. With "
                        "--subB-method=single-hash this same hash also decides which "
                        "base positions are kept.")
    p.add_argument("--gate-seed", type=int, default=31337,
                   help="Seed for the xxh32 subsampling gates: the --subB point gate "
                        "under --subB-method=hash-threshold, the mixed-stride "
                        "chr->stride/offset hash, and the --subA interval gate in mode "
                        "C. Default 31337 matches orig hammock. Lets you resample "
                        "independently of --seed. Under "
                        "--subB-method=single-hash the point gate uses --seed instead, "
                        "but the mode C --subA gate still uses this one.")
    p.add_argument("--verbose", action="store_true",
                   help="Progress on stderr: per-file sketching lines, interval/point "
                        "counts from the sketcher, the output path, and a reminder "
                        "that the two jaccard columns are different estimators.")
    p.add_argument("--memory-limit-gb", type=float, default=0.0,
                   help="Soft cap on this process's address space (RLIMIT_AS) in GiB, "
                        "so an over-large run fails with MemoryError rather than being "
                        "OOM-killed. 0 disables (default).")
    p.add_argument('--subB-method',
                   choices=['hash-threshold', 'mixed-stride', 'single-hash'],
                   default='mixed-stride',
                   help='How --subB picks which base positions to keep (default: '
                        'mixed-stride; no effect at --subB 1.0, where all three keep '
                        'every base). '
                        'mixed-stride: deterministic stride of about 1/subB, with '
                        'stride choice and phase keyed by chromosome name — fastest at '
                        'low subB. '
                        'hash-threshold: keep a position when xxh32(point, --gate-seed) '
                        'falls under the rate; this is the orig-hammock-parity choice. '
                        'single-hash: one xxh64 serves as both gate and HLL input — an '
                        'opt-in parity divergence (different accepted positions, not '
                        'byte-equal to orig), announced on stderr.')

    # Output-shape flags: which similarity columns get written, mutually
    # exclusive (argparse errors if both given, same pattern as --ref vs
    # --ref1/--ref2 above). Neither flag -> the default "ie" shape.
    metrics_group = p.add_mutually_exclusive_group()
    metrics_group.add_argument(
        '--metrics', action='store_true',
        help='Emit the full 7-column metrics block: reg_eq_similarity, '
             'jaccard_similarity_ie, containment_AB, containment_BA, '
             'cosketch_geom, cosketch_arith, cosketch_max. Tags the output '
             'filename _full. Mutually exclusive with --register-equality/--re.')
    metrics_group.add_argument(
        '--register-equality', '--re', action='store_true',
        dest='register_equality',
        help='Emit only reg_eq_similarity -- the cheap register-equality-only '
             'arm; the hammock-cpp standalone binary skips the union/'
             'containment pass entirely for this shape. Tags the output '
             'filename _re. Mutually exclusive with --metrics.')

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
            p.error(f"reference flags imply sequence mode (BED→FASTA), but --mode "
                    f"{args.mode} was given. Drop --mode (or use --mode sequence).")
        args.mode = "D"
    elif args.ref_cache_dir is not None or args.fasta_outdir is not None:
        p.error("--ref-cache-dir/--fasta-outdir only apply with a reference "
                "flag (--ref/--ref1/--ref2).")

    # Hardcoded constants the runner still reads. Hash is always xxh64; the
    # CSV `num_hashes` column is "NA" for HLL/minimizer (only meaningful for
    # MinHash, which isn't shipped).
    args.hash_size = 64
    args.num_hashes = "NA"

    # Canonical output-shape field the runner reads instead of the two raw
    # flags -- "ie" (default), "re" (--register-equality/--re), or "full"
    # (--metrics). See runner._metrics_shape.
    if args.metrics:
        args.metrics_mode = "full"
    elif args.register_equality:
        args.metrics_mode = "re"
    else:
        args.metrics_mode = "ie"
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
        peek). FASTA → sequence mode (D), BED → interval mode (interval-points, B).
      * If an interval-string knob is tweaked (--subA/--expA, which only affect
        the interval-string side), prefer interval-hybrid mode (C). --subB alone
        stays in the default interval mode (B natively subsamples points).
      * If we couldn't classify, fall back to the default interval mode (B).
    """
    if args.mode is not None:
        return args.mode

    first = _first_path(args.filepaths_file)
    kind = _classify_file(first) if first else None

    if kind == 'fasta':
        print("hammock: auto-detected FASTA input → sequence mode (D).", file=sys.stderr)
        return 'D'
    if kind == 'bed':
        if args.subA != 1.0 or args.expA != 0:
            print("hammock: auto-detected BED input + --subA/--expA → "
                  "interval-hybrid mode (C).", file=sys.stderr)
            return 'C'
        print("hammock: auto-detected BED input → interval mode "
              "(interval-points, B).", file=sys.stderr)
        return 'B'

    # Couldn't tell — fall back to the default interval mode, but make it visible.
    print("hammock: could not auto-detect input type; defaulting to interval "
          "mode (B). Pass --mode explicitly to silence this warning.", file=sys.stderr)
    return 'B'


def _resolve_sketch_type(args) -> str:
    """Mode D uses the 'minimizer' label (hyperloglog backing under the hood,
    matching the orig hammock's CSV output). Everything else is HLL."""
    return "minimizer" if args.mode == "D" else "hyperloglog"


def _default_threads(mode: str) -> int:
    """Default --threads, by mode.

    Mode D's sketch loop is pure Python (`MinimizerSketch.add_string`) wrapped
    around a `digest.window_minimizer` call that does not usefully release the
    GIL, so a thread pool is a convoy: measured 2.1x slower on 4 Maurano FASTAs
    and 2.8-4.5x on 8 full-size ones. A/B/C sketch inside the C++ extension,
    which does release the GIL, so they get the pool. See
    docs/seed-mode-d-threading.md for the evidence.
    """
    return 1 if mode == "D" else min(8, os.cpu_count() or 1)


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
        description="Download, decompress and faidx-index a reference genome into "
                    "the cache, then print its path. Run this once on a networked "
                    "(login) node: comparison runs never download, they only look "
                    "in the cache. The download is atomic, so an interrupted fetch "
                    "leaves no half-written reference behind.")
    fp.add_argument("spec", help="Reference keyword (hg38, hg19, mm10, mm39, hs1; "
                                 "aliases such as grch38, grch37, chm13, t2t are "
                                 "accepted) or an http(s) URL to a .fa/.fa.gz")
    fp.add_argument("--ref-cache-dir", default=None,
                    help="Cache directory (default: $HAMMOCK_REF_CACHE, else "
                         "~/.hammock/refs). Pass the same value to the comparison run.")
    fp.add_argument("--force", action="store_true",
                    help="Re-download and re-index even if already cached")
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

    # Guard early (before autodetect): the positionals must be list-of-paths
    # files, not input files. Fails fast with a clear message.
    from hammock.runner import _guard_is_list_file
    for arg_name, list_file in (("LIST1", args.filepaths_file),
                                ("LIST2", args.primary_file)):
        msg = _guard_is_list_file(list_file, arg_name)
        if msg:
            print(f"Error: {msg}", file=sys.stderr)
            return 2

    args.mode = _autodetect_mode(args)
    args.sketch_type = _resolve_sketch_type(args)
    # Two budgets. `threads` drives sketching, where Mode D is GIL-bound;
    # `io_threads` drives the bed2fasta extraction phase, which shells out to
    # `bedtools getfasta` and genuinely parallelizes regardless of mode.
    args.io_threads = (args.threads if args.threads is not None
                       else min(8, os.cpu_count() or 1))
    # A third budget, for the C++ OpenMP pairwise phase. It must NOT inherit
    # Mode D's `threads = 1` clamp: that clamp is about the GIL convoy while
    # *sketching*, and the pairwise loop is called once from the main thread
    # with the GIL released, so clamping it to 1 would just make Mode D slower.
    # 0 means "leave OpenMP's own default alone", which keeps the no-flag path
    # exactly as it was; an explicit --threads is honored so a run inside a
    # 4-CPU cgroup stops spawning a team per core on the whole node.
    args.omp_threads = args.threads if args.threads is not None else 0
    if args.threads is None:
        args.threads = _default_threads(args.mode)
    elif args.mode == "D" and args.threads > 1:
        print(f"hammock: --threads {args.threads} in sequence mode (D) is "
              f"usually SLOWER than --threads 1 (GIL convoy; see "
              f"docs/seed-mode-d-threading.md). Proceeding as asked.",
              file=sys.stderr)
    _apply_memory_limit(args.memory_limit_gb)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
