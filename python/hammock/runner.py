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


def _read_paths(list_file: str) -> List[str]:
    paths = []
    with open(list_file) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            paths.append(s)
    return paths


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


def _sketch_many(paths: List[str], args, label: str) -> list:
    """Sketch each path; parallelize with threads if --threads > 1.

    Threads are sound here: `_core.sketch_bed_file_hll` already releases the
    GIL (A/B/C), and `digest.window_minimizer` is a C++ extension that does
    too (Mode D), so we get real parallelism for the sketch phase.

    When `args.verbose` is set, prints `[i/N] <basename> (<elapsed>s)` to
    stderr as each file finishes (results stay in original order regardless).
    """
    n_threads = args.threads or 1
    n = len(paths)
    sketches = [None] * n
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
        return i, s

    if n_threads <= 1 or n <= 1:
        for i in range(n):
            _, s = _one(i)
            sketches[i] = s
        return sketches

    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        for fut in as_completed([ex.submit(_one, i) for i in range(n)]):
            i, s = fut.result()
            sketches[i] = s
    return sketches


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


def _write_mode_d_csv(args, queries, refs, query_sketches, ref_sketches) -> int:
    """Mode D output: Jaccard + containment + cosketch on both the minimizer
    and merged (`_with_ends`) sketches."""
    if args.verbose:
        print("Computing pairwise minimizer Jaccard + containment (Mode D)...",
              file=sys.stderr)
    minimizer_query = [s.minimizer_hll for s in query_sketches]
    minimizer_ref = [s.minimizer_hll for s in ref_sketches]
    merged_query = [s.merged() for s in query_sketches]
    merged_ref = [s.merged() for s in ref_sketches]

    jac, c_ab, c_ba = _core.pairwise_metrics_hll(minimizer_query, minimizer_ref)
    cs_geom, cs_arith, cs_max = _cosketch_from_containments(c_ab, c_ba)

    jac_e, c_ab_e, c_ba_e = _core.pairwise_metrics_hll(merged_query, merged_ref)
    cs_geom_e, cs_arith_e, cs_max_e = _cosketch_from_containments(c_ab_e, c_ba_e)

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
        ["jaccard_similarity"] + _CONTAINMENT_COLS
        + ["jaccard_similarity_with_ends"]
        + [c + "_with_ends" for c in _CONTAINMENT_COLS]
    )
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
                    float(jac[i, j]),
                    float(c_ab[i, j]), float(c_ba[i, j]),
                    float(cs_geom[i, j]), float(cs_arith[i, j]), float(cs_max[i, j]),
                    float(jac_e[i, j]),
                    float(c_ab_e[i, j]), float(c_ba_e[i, j]),
                    float(cs_geom_e[i, j]), float(cs_arith_e[i, j]), float(cs_max_e[i, j]),
                ])
                w.writerow(row)

    if args.verbose:
        print(f"Wrote {out_path}", file=sys.stderr)
    return 0


def run(args) -> int:
    queries = _read_paths(args.filepaths_file)
    refs = _read_paths(args.primary_file)

    if not queries:
        print(f"Error: no paths found in {args.filepaths_file}", file=sys.stderr)
        return 2
    if not refs:
        print(f"Error: no paths found in {args.primary_file}", file=sys.stderr)
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
    jaccard, c_ab, c_ba = _core.pairwise_metrics_hll(query_sketches, ref_sketches)
    cs_geom, cs_arith, cs_max = _cosketch_from_containments(c_ab, c_ba)

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

    similarity_measures = ["jaccard_similarity"] + _CONTAINMENT_COLS
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
                    float(jaccard[i, j]),
                    float(c_ab[i, j]), float(c_ba[i, j]),
                    float(cs_geom[i, j]), float(cs_arith[i, j]), float(cs_max[i, j]),
                ])
                w.writerow(row)

    if args.verbose:
        print(f"Wrote {out_path}", file=sys.stderr)
    return 0
