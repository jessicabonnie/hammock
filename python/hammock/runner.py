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

    Containments are clamped to 1.0 first. They can exceed it by a few ulp:
    `inter = cA + cB - cUnion` is a difference of three Ertl estimates, and
    while the union estimate is monotone (union registers are element-wise max)
    so `inter <= min(cA, cB)` holds mathematically, rounding can push the ratio
    just past 1 -- measured 0.2% of random pairs across p=4..20, max excess one
    ulp. Clamping guarantees `denom >= 1`, so neither a division by zero nor an
    out-of-range result is reachable. A zero containment means the intersection
    estimate was zero -- genuinely empty, or clamped from a negative in
    HLLSketch::intersection_size -- and is scored 0.0.

    Note the emitted `containment_AB`/`containment_BA` columns are NOT clamped,
    so reconstructing J_ie from them can differ from this column in the last
    couple of digits. Matches
    experiments/bedtools_benchmark/estimator_compare.py, except for the clamp.

    Unlike `jaccard_similarity` (register equality), this carries no
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
    """
    if not args.verbose:
        return
    print("note: jaccard_similarity is register-equality -- "
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


def _label(path: str, full_paths: bool) -> str:
    return os.path.normpath(path) if full_paths else os.path.basename(path)


def _sketch_one_file(path: str, args):
    if args.mode == "D":
        # Mode D's sketch_type is labelled "minimizer" (orig parity) but the
        # backing storage is still our HLL — see modes/sequence.py.
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


def _sketch_many(paths: List[str], args, label: str) -> list:
    """Sketch each path; parallelize with threads if --threads > 1.

    Threads are sound here: `_core.sketch_bed_file_hll` already releases the
    GIL (A/B/C), and `digest.window_minimizer` is a C++ extension that does
    too (Mode D), so we get real parallelism for the sketch phase.

    When `args.verbose` is set, prints `[i/N] <basename> (<elapsed>s)` to
    stderr as each file finishes (results stay in original order regardless).
    """
    n = len(paths)
    progress_lock = threading.Lock()
    done = [0]
    t0 = time.monotonic()

    def _one(i: int):
        s = _sketch_one_file(paths[i], args)
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

    jac, c_ab, c_ba = _core.pairwise_metrics_hll(minimizer_query, minimizer_ref)
    cs_geom, cs_arith, cs_max = _cosketch_from_containments(c_ab, c_ba)
    jac_ie = _jaccard_ie_from_containments(c_ab, c_ba)

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
    ) + ".csv"

    similarity_measures = (
        ["jaccard_similarity", "jaccard_similarity_ie"] + _CONTAINMENT_COLS
    )
    # ref1/ref2 are always emitted (trailing, "NA" outside bed2fasta mode) so
    # the Mode D header stays fixed and cross-reference provenance is recorded.
    header = _build_header(args, similarity_measures) + ["ref1", "ref2"]

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, qlabel in enumerate(query_labels):
            for j, rlabel in enumerate(ref_labels):
                row = _row_prefix(args, qlabel, rlabel)
                row.extend([
                    float(jac[i, j]), float(jac_ie[i, j]),
                    float(c_ab[i, j]), float(c_ba[i, j]),
                    float(cs_geom[i, j]), float(cs_arith[i, j]), float(cs_max[i, j]),
                    ref1, ref2,
                ])
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
    jaccard, c_ab, c_ba = _core.pairwise_metrics_hll(query_sketches, ref_sketches)
    cs_geom, cs_arith, cs_max = _cosketch_from_containments(c_ab, c_ba)
    jac_ie = _jaccard_ie_from_containments(c_ab, c_ba)

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
    ) + ".csv"

    similarity_measures = (["jaccard_similarity", "jaccard_similarity_ie"]
                           + _CONTAINMENT_COLS)
    header = _build_header(args, similarity_measures)

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, q in enumerate(queries):
            qlabel = _label(q, args.full_paths)
            for j, r in enumerate(refs):
                rlabel = _label(r, args.full_paths)
                row = _row_prefix(args, qlabel, rlabel)
                row.extend([
                    float(jaccard[i, j]), float(jac_ie[i, j]),
                    float(c_ab[i, j]), float(c_ba[i, j]),
                    float(cs_geom[i, j]), float(cs_arith[i, j]), float(cs_max[i, j]),
                ])
                w.writerow(row)

    if args.verbose:
        print(f"Wrote {out_path}", file=sys.stderr)
    _print_estimator_note(args)
    return 0
