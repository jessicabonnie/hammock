"""Sketch each input file, compute pairwise similarity, write the output CSV.

The output format (header + columns) matches the original hammock for parity.
"""
from __future__ import annotations

import csv
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List

import numpy as np

from hammock import _core
from hammock.outprefix import get_new_prefix


# Containment + co-sketch columns for any pair of HLL sketch lists.
_CONTAINMENT_COLS = ["containment_AB", "containment_BA",
                     "cosketch_geom", "cosketch_arith", "cosketch_max"]

# Three output shapes, selected by --metrics/--register-equality (see cli.py;
# `args.metrics_mode` is "ie" (default), "re", or "full"). Docs/seed-metrics-
# column-restructure.md Part 2 -- both `_metrics_shape` and `_metrics_row_values`
# must agree on column count and order, and both call sites (A/B/C `run()` and
# `_write_mode_d_csv`) must use the same `metrics_mode` for header and rows.
def _metrics_shape(metrics_mode: str) -> tuple:
    """Return (similarity_measures column names, filename tag) for
    `metrics_mode` ("ie"/"re"/"full")."""
    if metrics_mode == "full":
        return (["reg_eq_similarity", "jaccard_similarity_ie"]
                 + _CONTAINMENT_COLS,
                 "full")
    if metrics_mode == "re":
        return (["reg_eq_similarity"], "re")
    return (["jaccard_similarity_ie"], "ie")


def _metrics_row_values(metrics_mode: str, jac, jac_ie, c_ab, c_ba,
                        cs_geom, cs_arith, cs_max, i: int, j: int) -> list:
    """Per-pair values matching `_metrics_shape`'s column list, in order.

    Only reads the arrays its own branch needs: callers that skip computing
    `jac_ie`/`cs_geom`/`cs_arith`/`cs_max` for "re" (which never uses them)
    rely on that -- do not add a code path here that reads them outside the
    "full"/"ie" branches.
    """
    if metrics_mode == "full":
        return [float(jac[i, j]), float(jac_ie[i, j]), float(c_ab[i, j]), float(c_ba[i, j]),
                float(cs_geom[i, j]), float(cs_arith[i, j]), float(cs_max[i, j])]
    if metrics_mode == "re":
        return [float(jac[i, j])]
    return [float(jac_ie[i, j])]


def _cosketch_from_containments(c_ab, c_ba):
    """Return (geom_mean, arith_mean, max) of the two directional containments."""
    geom = np.sqrt(np.maximum(c_ab * c_ba, 0.0))
    arith = 0.5 * (c_ab + c_ba)
    mx = np.maximum(c_ab, c_ba)
    return geom, arith, mx


def _jaccard_ie_from_containments(c_ab, c_ba):
    """Set-Jaccard from the two directional containments.

    Both are |A n B| / |.| over the same inclusion-exclusion intersection, so
        1 / (1/C_AB + 1/C_BA - 1) == |A n B| / (|A| + |B| - |A n B|).

    Containments are clamped to 1.0 first. They can exceed it by rounding:
    `inter = cA + cB - cUnion` is a difference of three Ertl estimates, and
    while the union estimate is monotone (union registers are element-wise max)
    so `inter <= min(cA, cB)` holds mathematically, cancellation can push the
    ratio just past 1. The excess is small but **not** bounded by one ulp, as
    this said before 2026-08-06: measured max excess is ~8 ulp on nested
    subset pairs at p<=20, and ~3.4e4 ulp (1 + 7.5e-12) in the extreme
    size-ratio regime (|A| ~ 1e3 against |B| ~ 1e7, A subset of B), because
    `inter` inherits `ulp(cB)` and that dwarfs a tiny `cA`. Only the noise is
    absorbed -- the true containment cannot exceed 1 -- so the clamp is a
    correct projection, and it guarantees `denom >= 1`, making both a division
    by zero and an out-of-range result unreachable. A zero containment means
    the intersection estimate was zero -- genuinely empty, or clamped from a
    negative by the `>= 0` clamp in pairwise_metrics_hll -- and is scored 0.0.

    Note the emitted `containment_AB`/`containment_BA` columns are NOT clamped,
    so reconstructing J_ie from them can differ from this column in the last
    couple of digits -- and, in the extreme size-ratio regime above, by more
    than that. Values just over 1.0 do reach the CSV (`cosketch_max` inherits
    them, since the cosketches are also computed from the unclamped pair);
    a consumer asserting `<= 1.0` will trip. Matches
    experiments/bedtools_benchmark/estimator_compare.py, except for the clamp.

    Unlike `reg_eq_similarity` (register equality), this carries no
    chance-agreement floor, so it is comparable to `bedtools jaccard`.
    """
    c_ab = np.minimum(np.asarray(c_ab, dtype=float), 1.0)
    c_ba = np.minimum(np.asarray(c_ba, dtype=float), 1.0)
    # Gate at 1e-300 rather than 0: below that the reciprocal overflows and
    # numpy warns (or raises under np.seterr(all='raise')). Real containments
    # bottom out around 1e-19, so this only affects hand-fed input.
    ok = (c_ab > 1e-300) & (c_ba > 1e-300)
    # The inner np.where keeps the (eagerly evaluated) reciprocals off the
    # zero branch, so no divide-by-zero warning is ever emitted.
    safe_ab = np.where(ok, c_ab, 1.0)
    safe_ba = np.where(ok, c_ba, 1.0)
    denom = 1.0 / safe_ab + 1.0 / safe_ba - 1.0
    return np.where(ok, 1.0 / denom, 0.0)


def _print_estimator_note(args) -> None:
    """Remind the user that the two Jaccard columns are different estimators.

    Verbose-gated on purpose: stderr is shared with progress output, and a
    default-on banner is the kind of thing that ends up captured into someone's
    log parser. The README and docs/jaccard-definitional-gap.md are the
    canonical explanation; this is a nudge for interactive runs.

    Parameterized by `args.metrics_mode` because only "full" emits both
    columns -- the default ("ie") shape has no reg_eq_similarity column and
    "re" has no jaccard_similarity_ie column, so the two-column comparison
    below would describe an absent column in those shapes.
    """
    if not args.verbose:
        return
    mode = getattr(args, "metrics_mode", "full")
    if mode == "ie":
        print("note: jaccard_similarity_ie is set-Jaccard (inclusion-exclusion), "
              "comparable to bedtools.", file=sys.stderr)
        return
    if mode == "re":
        print("note: reg_eq_similarity is "
              "register-equality -- biased high, and the bias depends on both "
              "sketch load and |A|/|B|, so rank only within comparable pairs.",
              file=sys.stderr)
        return
    print("note: reg_eq_similarity is register-equality -- "
          "biased high,\n"
          "      and the bias depends on both sketch load and |A|/|B|, so rank "
          "only within\n"
          "      comparable pairs. The _ie column is set-Jaccard, comparable "
          "to bedtools.",
          file=sys.stderr)


def _read_paths(list_file: str) -> List[str]:
    paths = []
    with open(list_file) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            paths.append(s)
    return paths


def _looks_like_input_data(line: str) -> str | None:
    """Return 'FASTA'/'BED' if `line` is sequence/interval DATA rather than a
    filesystem path — i.e. the user passed an input file where a list-of-paths
    was expected. Returns None for a normal path line."""
    if line.startswith('>') or line.startswith(';'):
        return 'FASTA'
    fields = line.split('\t')
    if len(fields) >= 3 and fields[1].isdigit() and fields[2].isdigit():
        return 'BED'
    return None


def _guard_is_list_file(list_file: str, arg_name: str) -> str | None:
    """If `list_file` looks like an actual BED/FASTA input file (not a text file
    listing one path per line), return an error message; else None."""
    try:
        with open(list_file) as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith('#') or s.startswith(('track', 'browser')):
                    continue
                kind = _looks_like_input_data(s)
                if kind:
                    return (
                        f"{arg_name} ({list_file}) looks like a {kind} file, not a "
                        f"list of paths. Pass a text file that lists one input path "
                        f"per line, e.g.  ls *.bed > list.txt  (then: hammock list.txt "
                        f"other_list.txt ...).")
                return None  # first real line looks like a path — OK
    except OSError:
        return None  # let the normal open/read error surface downstream
    return None


# Interval modes hand the path straight to a bare `std::ifstream` in C++
# (cpp/src/processing_modes.cpp), which neither decompresses nor decodes
# anything. A gzipped or BigBed input therefore parses to zero intervals and
# scores 0.0 against everything -- including itself -- while exiting 0. Reject
# it up front rather than emit a plausible-looking all-zero CSV.
#
# Entries are (magic bytes, description, conversion template). bgzip/BAM/.tbi
# are gzip-framed and so are caught by the gzip entry. Prefer a LONGER magic
# over a shorter one: these are matched against the head of a file that is
# usually legitimate text, so a short magic risks rejecting real data. bzip2 is
# spelled out per compression level ("BZh1".."BZh9") for exactly that reason --
# a bare b"BZh" also matches a BED whose first chromosome is named "BZh...".
_BINARY_MAGIC = (
    (b"\x1f\x8b", "gzip-compressed", "gunzip -c {p} > {out}"),
    (b"\x28\xb5\x2f\xfd", "zstd-compressed", "zstd -dc {p} > {out}"),
    (b"\xfd7zXZ", "xz-compressed", "xz -dc {p} > {out}"),
    (b"\xeb\xf2\x89\x87", "a binary BigBed", "bigBedToBed {p} {out}"),
    (b"\x87\x89\xf2\xeb", "a binary BigBed", "bigBedToBed {p} {out}"),
    (b"\x26\xfc\x8f\x88", "a binary BigWig", "bigWigToBedGraph {p} {out}"),
    (b"\x88\x8f\xfc\x26", "a binary BigWig", "bigWigToBedGraph {p} {out}"),
) + tuple(
    (b"BZh" + bytes([level]), "bzip2-compressed", "bzip2 -dc {p} > {out}")
    for level in b"123456789"
)


def _guard_plain_text_bed(paths: List[str], which: str) -> str | None:
    """Return an error message if any interval-mode input is compressed or
    binary rather than plain text; else None.

    Magic-byte matching only, and deliberately conservative: it catches the
    formats people actually hand us by accident, not every possible one. A
    binary format not listed here still reaches the C++ parser and still
    sketches to nothing -- the general fix is a post-sketch "parsed 0
    intervals" check, which this does not replace.
    """
    for p in paths:
        # Only stat/read regular files. Reading 4 bytes from a FIFO consumes
        # them and lets the writer finish, after which the C++ reopen blocks
        # forever -- and a FIFO fed by `zcat foo.bed.gz > pipe` is exactly the
        # workaround this error message invites.
        try:
            if not os.path.isfile(p):
                continue
            with open(p, 'rb') as fh:
                head = fh.read(8)
        except OSError:
            continue  # let the normal open error surface downstream
        for magic, what, fix in _BINARY_MAGIC:
            if head.startswith(magic):
                return (
                    f"{which} input {p} is {what}. Interval modes (A/B/C) read "
                    f"plain-text BED only -- nothing is decompressed or decoded "
                    f"first, so this file would sketch to nothing and score 0.0 "
                    f"against every file including itself. Convert it first, "
                    f"e.g.  " + fix.format(p=p, out=_decompressed_name(p)))
    return None


def _decompressed_name(path: str) -> str:
    """A conversion target that cannot collide with an existing file -- the
    obvious `foo.bed.gz` -> `foo.bed` would overwrite a real `foo.bed`."""
    base = path[:-3] if path.lower().endswith('.gz') else path
    base = os.path.splitext(base)[0] or base
    cand = base + '.plain.bed'
    n = 2
    while os.path.exists(cand):
        cand = f"{base}.plain{n}.bed"
        n += 1
    return cand


def _label(path: str, full_paths: bool) -> str:
    return os.path.normpath(path) if full_paths else os.path.basename(path)


def _sketch_one_file(path: str, args, sketch_threads: int = 0):
    if args.mode == "D":
        # Mode D's sketch_type is labelled "minimizer" (orig parity) but the
        # backing storage is still our HLL — see modes/sequence.py. Mode D
        # doesn't sketch via sketch_bed_file_hll, so sketch_threads (A/B/C
        # only) doesn't apply here.
        if args.sketch_type not in ("minimizer", "hyperloglog"):
            raise NotImplementedError(
                f"Mode D doesn't support sketch_type '{args.sketch_type}' yet."
            )
        from hammock.modes.sequence import sketch_fasta
        return sketch_fasta(path, args)
    if args.sketch_type != "hyperloglog":
        raise NotImplementedError(
            f"Sketch type '{args.sketch_type}' is not implemented yet (Phase 1: HLL only)."
        )
    return _core.sketch_bed_file_hll(
        path=path,
        mode=args.mode,
        precision=args.precision,
        sub_a=args.subA,
        sub_b=args.subB,
        exp_a=args.expA,
        subB_method=args.subB_method,
        seed=args.seed,
        gate_seed=args.gate_seed,
        verbose=args.verbose,
        threads=sketch_threads,
    )


def _parallel_map(n: int, fn, threads: int) -> list:
    """Run ``fn(i)`` for ``i`` in ``range(n)``, returning results in index
    order. Uses a thread pool when ``threads > 1`` and ``n > 1``; otherwise
    runs sequentially. Any exception raised by ``fn`` propagates (aborting the
    whole map) via ``future.result()`` — results stay index-aligned, never
    completion-ordered.
    """
    results = [None] * n
    if threads <= 1 or n <= 1:
        for i in range(n):
            results[i] = fn(i)
        return results
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futures = {ex.submit(fn, i): i for i in range(n)}
        for fut in as_completed(futures):
            results[futures[fut]] = fut.result()
    return results


def _split_thread_budget(total_threads: int, n: int) -> tuple:
    """Split a `--threads` budget between `_sketch_many`'s outer file-dispatch
    pool and the inner OMP team each file's sketch phase gets.

    Returns `(outer_workers, sketch_threads)`. `outer_workers` mirrors
    `_parallel_map`'s own branch condition (`threads <= 1 or n <= 1`, else a
    pool of `max_workers=threads`) — it must be called with the same
    `total_threads` that gets passed into `_parallel_map`, and the two are
    kept in sync by construction, not an enforced contract, so update both if
    either changes.

    `sketch_threads` is `total_threads` floor-divided by `outer_workers`, so
    the product is always `<= total_threads` — it can only *under*-use the
    budget when `n` doesn't divide evenly (e.g. threads=8, n=3 gives 3*2=6,
    not 8), never over-use it — floored at 1 so a large `n` never starves an
    individual file's sketch phase down to 0 threads (a `num_threads(0)`
    clause is invalid). See `_sketch_many`'s docstring for why this split
    needs to exist at all.
    """
    outer_workers = min(total_threads, n) if (total_threads > 1 and n > 1) else 1
    sketch_threads = max(1, total_threads // outer_workers)
    return outer_workers, sketch_threads


def _sketch_many(paths: List[str], args, label: str) -> list:
    """Sketch each path; parallelize with threads if --threads > 1.

    Threads are sound here: `_core.sketch_bed_file_hll` already releases the
    GIL (A/B/C), and `digest.window_minimizer` is a C++ extension that does
    too (Mode D), so we get real parallelism for the sketch phase.

    When `args.verbose` is set, prints `[i/N] <basename> (<elapsed>s)` to
    stderr as each file finishes (results stay in original order regardless).

    **Divides `--threads` between this outer file-dispatch pool and the inner
    OMP team** (A/B/C only; Mode D's per-file cost is single-threaded Python +
    digest, not OMP) via `_split_thread_budget` — see its docstring for the
    split itself. Before this, `sketch_bed_file_hll` had no thread parameter
    at all, so every pool worker spawned an *unclamped* OMP team at OpenMP's
    ambient default (all cores) — `--threads 4` on a 4-CPU cgroup could use
    ~12 cores concurrently. See docs/seed-hammock-cpp-file-dispatch.md Part 1:
    confirmed by code inspection (no code path carried `--threads` into the
    sketching region) and corroborated empirically (single-file isolation
    measured up to 1459% CPU at requested `--threads` values of 2-16).
    """
    n = len(paths)
    total_threads = args.threads or 1
    _outer_workers, sketch_threads = _split_thread_budget(total_threads, n)
    progress_lock = threading.Lock()
    done = [0]
    t0 = time.monotonic()

    def _one(i: int):
        s = _sketch_one_file(paths[i], args, sketch_threads)
        if args.verbose:
            with progress_lock:
                done[0] += 1
                print(f"  [{done[0]}/{n}] {label}: {os.path.basename(paths[i])} "
                      f"({time.monotonic() - t0:.1f}s)", file=sys.stderr, flush=True)
        return s

    return _parallel_map(n, _one, args.threads or 1)


def _build_header(args, similarity_measures: List[str]) -> List[str]:
    header = ["file1", "file2", "sketch_type", "mode"]
    if args.sketch_type in ("hyperloglog", "minhash", "minimizer"):
        header.extend(["precision", "num_hashes", "kmer_size", "window_size"])
    if args.mode == "C":
        header.extend(["subA", "subB", "expA"])
    header.extend(similarity_measures)
    return header


def _row_prefix(args, qlabel: str, rlabel: str) -> List:
    """Build the per-row metadata prefix matching the original CSV layout."""
    row: List = [qlabel, rlabel, args.sketch_type, args.mode]
    if args.sketch_type in ("hyperloglog", "minhash", "minimizer"):
        precision = args.precision if args.sketch_type in ("hyperloglog", "minimizer") else "NA"
        num_hashes = args.num_hashes if args.sketch_type == "minhash" else "NA"
        if args.mode == "D" or args.kmer_size > 0:
            row.extend([precision, num_hashes, args.kmer_size, args.window_size])
        else:
            row.extend([precision, num_hashes, "NA", "NA"])
    if args.mode == "C":
        row.extend([args.subA, args.subB, args.expA])
    return row


def _write_mode_d_csv(args, queries, refs, query_sketches, ref_sketches,
                      query_labels=None, ref_labels=None,
                      ref1="NA", ref2="NA") -> int:
    """Mode D output: Jaccard + containment + cosketch on the minimizer sketch.

    `query_labels`/`ref_labels` override the CSV file1/file2 text (used by the
    bed2fasta path so rows are labelled by the original BED files, not the
    generated FASTA temp names). `ref1`/`ref2` populate the always-present
    trailing reference columns (`"NA"` for plain FASTA runs)."""
    if query_labels is None:
        query_labels = [_label(q, args.full_paths) for q in queries]
    if ref_labels is None:
        ref_labels = [_label(r, args.full_paths) for r in refs]
    if args.verbose:
        print("Computing pairwise minimizer Jaccard + containment (Mode D)...",
              file=sys.stderr)
    minimizer_query = [s.minimizer_hll for s in query_sketches]
    minimizer_ref = [s.minimizer_hll for s in ref_sketches]

    # `pairwise_metrics_hll` is called unconditionally regardless of shape --
    # the Python binding always computes the fused union pass as a side
    # effect (see docs/seed-metrics-column-restructure.md Part 2's cost
    # note), so there is no cheaper call for the "re" shape on this front-end.
    # Unlike hammock-cpp (Part 4), which branches per-shape and skips the
    # union pass entirely for --register-equality, the Python CLI cannot skip
    # it without a binding change -- not a regression, since it already paid
    # this cost for every shape before this restructure.
    jac, c_ab, c_ba = _core.pairwise_metrics_hll(
        minimizer_query, minimizer_ref, threads=getattr(args, "omp_threads", 0) or 0)
    metrics_mode = getattr(args, "metrics_mode", "full")
    cs_geom = cs_arith = cs_max = jac_ie = None
    if metrics_mode in ("full", "ie"):
        jac_ie = _jaccard_ie_from_containments(c_ab, c_ba)
    if metrics_mode == "full":
        cs_geom, cs_arith, cs_max = _cosketch_from_containments(c_ab, c_ba)

    similarity_measures, tag = _metrics_shape(metrics_mode)

    out_path = get_new_prefix(
        outprefix=args.outprefix,
        sketch_type=args.sketch_type,
        mode=args.mode,
        num_hashes=args.num_hashes,
        precision=args.precision,
        subA=args.subA,
        subB=args.subB,
        expA=args.expA,
        kmer_size=args.kmer_size,
        window_size=args.window_size,
        sequence_hll_hash=getattr(args, "sequence_hll_hash", "rehash-selector64"),
        metrics_tag=tag,
    ) + ".csv"

    # ref1/ref2 are always emitted (trailing, "NA" outside bed2fasta mode) so
    # the Mode D header stays fixed and cross-reference provenance is recorded.
    header = _build_header(args, similarity_measures) + ["ref1", "ref2"]

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, qlabel in enumerate(query_labels):
            for j, rlabel in enumerate(ref_labels):
                row = _row_prefix(args, qlabel, rlabel)
                row.extend(_metrics_row_values(
                    metrics_mode, jac, jac_ie, c_ab, c_ba,
                    cs_geom, cs_arith, cs_max, i, j))
                # ref1/ref2 trailing columns are Mode-D-specific (see
                # _build_header above) and independent of metrics_mode --
                # always appended after the similarity columns, for every
                # shape, matching the header built above.
                row.extend([ref1, ref2])
                w.writerow(row)

    if args.verbose:
        print(f"Wrote {out_path}", file=sys.stderr)
    _print_estimator_note(args)
    return 0


def _ref_tag(spec: str) -> str:
    """Short filesystem-safe tag for a reference spec, for the output filename."""
    from hammock import refs as refs_mod
    kw = refs_mod.canonical_keyword(spec)
    if kw:
        return kw
    base = os.path.basename(spec.rstrip("/")) or "ref"
    base = base.split("?")[0]
    for ext in (".gz", ".fa", ".fasta", ".fna"):
        if base.lower().endswith(ext):
            base = base[: -len(ext)]
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in base) or "ref"


def _warn_duplicate_labels(labels: List[str], which: str, full_paths: bool) -> None:
    if full_paths:
        return
    seen = set()
    dups = {x for x in labels if x in seen or seen.add(x)}
    if dups:
        print(f"[hammock] warning: duplicate {which} basenames "
              f"({', '.join(sorted(dups))}); rows will be ambiguous. "
              f"Pass --full-paths to disambiguate.", file=sys.stderr)


def _run_bed2fasta(args, queries: List[str], refs: List[str]) -> int:
    """BED lists → per-list FASTA via bedtools getfasta → existing Mode D path.

    Reference resolution never downloads (see hammock.refs); conversion,
    sketching, and the CSV write all happen inside one TemporaryDirectory so
    the generated FASTAs outlive sketching and are cleaned up afterwards.
    """
    import contextlib
    import tempfile

    from hammock import bed2fasta as b2f
    from hammock import refs as refs_mod

    # Labels come from the ORIGINAL BED paths, captured before we swap in the
    # generated FASTA paths — otherwise CSV rows would be named after temp files.
    query_labels = [_label(q, args.full_paths) for q in queries]
    ref_labels = [_label(r, args.full_paths) for r in refs]
    _warn_duplicate_labels(query_labels, "query", args.full_paths)
    _warn_duplicate_labels(ref_labels, "primary", args.full_paths)

    cache = args.ref_cache_dir or refs_mod.default_cache_dir()
    try:
        ref1_fasta = refs_mod.resolve_reference(args.ref1, cache)
        ref2_fasta = refs_mod.resolve_reference(args.ref2, cache)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2

    print(f"[hammock] bed2fasta: bedtools getfasta ({b2f.bedtools_version()}); "
          f"ref1={args.ref1}, ref2={args.ref2}", file=sys.stderr)

    # Extraction shells out to `bedtools getfasta`, so it parallelizes for real
    # even though this is Mode D — it is exempt from the Mode D thread clamp and
    # uses the I/O budget (`cli._default_threads` / `args.io_threads`). getattr
    # keeps hand-built Namespaces (tests, library callers) working.
    n_threads = getattr(args, "io_threads", None) or args.threads or 1
    with contextlib.ExitStack() as stack:
        if args.fasta_outdir:
            os.makedirs(args.fasta_outdir, exist_ok=True)
            base_dir = args.fasta_outdir
        else:
            base_dir = stack.enter_context(
                tempfile.TemporaryDirectory(prefix="hammock_b2f_"))
        out1 = os.path.join(base_dir, "list1")
        out2 = os.path.join(base_dir, "list2")

        try:
            query_fastas = b2f.convert_list(queries, ref1_fasta, out1,
                                            n_threads, args.verbose)
            ref_fastas = b2f.convert_list(refs, ref2_fasta, out2,
                                          n_threads, args.verbose)
        except b2f.ConversionError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 2

        if args.verbose:
            print(f"Sketching {len(query_fastas)} query + {len(ref_fastas)} "
                  f"primary FASTA files (mode D)...", file=sys.stderr)
        query_sketches = _sketch_many(query_fastas, args, label="query")
        ref_sketches = _sketch_many(ref_fastas, args, label="primary")

        # Tag the output filename with both refs so cross-reference runs to the
        # same -o don't silently overwrite each other.
        args.outprefix = f"{args.outprefix}_{_ref_tag(args.ref1)}-vs-{_ref_tag(args.ref2)}"
        return _write_mode_d_csv(
            args, query_fastas, ref_fastas, query_sketches, ref_sketches,
            query_labels=query_labels, ref_labels=ref_labels,
            ref1=args.ref1, ref2=args.ref2,
        )


def run(args) -> int:
    queries = _read_paths(args.filepaths_file)
    refs = _read_paths(args.primary_file)

    if not queries:
        print(f"Error: no paths found in {args.filepaths_file}", file=sys.stderr)
        return 2
    if not refs:
        print(f"Error: no paths found in {args.primary_file}", file=sys.stderr)
        return 2

    if getattr(args, "bed2fasta", False):
        return _run_bed2fasta(args, queries, refs)

    if args.mode in ("A", "B", "C"):
        for paths, which in ((queries, "Query"), (refs, "Primary")):
            msg = _guard_plain_text_bed(paths, which)
            if msg:
                print(f"Error: {msg}", file=sys.stderr)
                return 2

    n_threads = args.threads or 1
    if args.subB_method == 'single-hash' and args.subB < 1.0:
        print("[hammock] note: --subB-method=single-hash diverges from orig-parity "
              "(one xxh64 does both gate and HLL ingestion; accepted-position set "
              "differs from hash-threshold). CSV results are not byte-equal to orig.",
              file=sys.stderr)
    if args.verbose:
        print(f"Sketching {len(queries)} query files "
              f"(mode {args.mode}, p={args.precision}, threads={n_threads})...",
              file=sys.stderr)
    query_sketches = _sketch_many(queries, args, label="query")

    if args.verbose:
        print(f"Sketching {len(refs)} primary files...", file=sys.stderr)
    ref_sketches = _sketch_many(refs, args, label="primary")

    if args.mode == "D":
        return _write_mode_d_csv(args, queries, refs, query_sketches, ref_sketches)

    if args.verbose:
        print("Computing pairwise Jaccard + containment...", file=sys.stderr)
    # See the matching comment in _write_mode_d_csv: this call is
    # unconditional regardless of shape -- the Python binding always
    # computes the fused union pass, so "re" is not cheaper than "full" on
    # this front-end (docs/seed-metrics-column-restructure.md Part 2).
    jaccard, c_ab, c_ba = _core.pairwise_metrics_hll(
        query_sketches, ref_sketches, threads=getattr(args, "omp_threads", 0) or 0)
    metrics_mode = getattr(args, "metrics_mode", "full")
    cs_geom = cs_arith = cs_max = jac_ie = None
    if metrics_mode in ("full", "ie"):
        jac_ie = _jaccard_ie_from_containments(c_ab, c_ba)
    if metrics_mode == "full":
        cs_geom, cs_arith, cs_max = _cosketch_from_containments(c_ab, c_ba)

    similarity_measures, tag = _metrics_shape(metrics_mode)

    out_path = get_new_prefix(
        outprefix=args.outprefix,
        sketch_type=args.sketch_type,
        mode=args.mode,
        num_hashes=args.num_hashes,
        precision=args.precision,
        subA=args.subA,
        subB=args.subB,
        expA=args.expA,
        kmer_size=args.kmer_size,
        window_size=args.window_size,
        sequence_hll_hash=getattr(args, "sequence_hll_hash", "rehash-selector64"),
        metrics_tag=tag,
    ) + ".csv"

    header = _build_header(args, similarity_measures)

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, q in enumerate(queries):
            qlabel = _label(q, args.full_paths)
            for j, r in enumerate(refs):
                rlabel = _label(r, args.full_paths)
                row = _row_prefix(args, qlabel, rlabel)
                row.extend(_metrics_row_values(
                    metrics_mode, jaccard, jac_ie, c_ab, c_ba,
                    cs_geom, cs_arith, cs_max, i, j))
                w.writerow(row)

    if args.verbose:
        print(f"Wrote {out_path}", file=sys.stderr)
    _print_estimator_note(args)
    return 0
