#!/usr/bin/env python3
"""Sweep subB at multiple size classes for hammock-cpp Mode B with --mixed-stride.

Iterates (size_class × subB × replicate), runs hammock-cpp Mode B on the 6 BEDs
of each size class (all 15 pairs), captures wall time, CPU time, and max RSS
via /usr/bin/time, and parses per-pair Jaccards from the CSV output.

Writes a long-form CSV to results/. Each row is one (size_class, subB, rep,
file_a, file_b). Timing fields are attached identically to every pair within a
run (a hammock-cpp invocation produces all 15 pairs at once).
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "bedtools_benchmark"))
from benchmark_cpp_vs_bedtools import (  # noqa: E402
    find_hammock_cpp,
    get_system_info,
    run_with_time,
    SKETCH_RE,
    PAIR_RE,
)

SUBB_VALUES = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01]
SIZE_CLASSES = ["10k", "100k", "1M"]
FILES_PER_CLASS = 6
DEFAULT_REPS = 5
DEFAULT_PRECISION = 18
DEFAULT_THREADS = 8

DATA_DIR = SCRIPT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR / "results"

ROW_COLS = [
    "size_class", "num_intervals", "subB", "rep", "run_id",
    "file_a", "file_b", "jaccard",
    "wall_time", "cpu_time", "max_rss_mb",
    "sketch_creation_time", "comparison_time",
    "precision", "threads", "mixed_stride",
]

SIZE_TO_N = {"10k": 10_000, "100k": 100_000, "1M": 1_000_000}


def bed_paths_for(size_class: str) -> List[Path]:
    return [DATA_DIR / f"bed_{size_class}_{i}.bed" for i in range(FILES_PER_CLASS)]


def parse_hammock_csv(path: str) -> List[Dict[str, Any]]:
    """Hammock-cpp Mode B output is TAB-separated: query\\treference\\tjaccard_similarity."""
    rows: List[Dict[str, Any]] = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            j = r.get("jaccard_similarity") or r.get("jaccard")
            if j is None:
                continue
            rows.append({
                "file_a": os.path.basename(r.get("query") or r.get("file_a") or ""),
                "file_b": os.path.basename(r.get("reference") or r.get("file_b") or ""),
                "jaccard": float(j),
            })
    return rows


def run_one(
    binary: str,
    size_class: str,
    subB: float,
    precision: int,
    threads: int,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run hammock-cpp Mode B on the 6 BEDs of size_class; return timing + per-pair rows."""
    beds = bed_paths_for(size_class)
    for p in beds:
        if not p.exists():
            raise FileNotFoundError(f"missing input: {p}; run generate_data.py first")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as fh:
        listfile = fh.name
        for p in beds:
            fh.write(str(p) + "\n")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".out", delete=False) as fh:
        out_prefix = fh.name

    try:
        # mixed-stride is the default for --subB-method since v.X; this
        # sweep was written when it was an opt-in flag.
        cmd = [
            binary, listfile, listfile,
            "--mode", "B",
            "--subB", str(subB),
            "--subB-method", "mixed-stride",
            "-p", str(precision),
            "-t", str(threads),
            "-o", out_prefix,
            "--verbose",
        ]
        if verbose:
            print("  +", " ".join(cmd), file=sys.stderr)
        r = run_with_time(cmd)

        sketch_s: Optional[float] = None
        pair_s: Optional[float] = None
        for line in r["stderr"].splitlines():
            m = SKETCH_RE.match(line)
            if m:
                sketch_s = int(m.group(1)) / 1000.0
                continue
            m = PAIR_RE.match(line)
            if m:
                pair_s = int(m.group(1)) / 1000.0
        r["sketch_creation_time"] = sketch_s
        r["comparison_time"] = pair_s

        csvs = [f for f in glob.glob(out_prefix + "*") if f.endswith(".csv")]
        if not csvs:
            raise RuntimeError(f"no CSV produced; stderr=\n{r['stderr'][-500:]}")
        pairs = parse_hammock_csv(csvs[0])
        # keep only unordered self-pairs once: file_a < file_b (lex), drop equal
        seen = set()
        kept = []
        for row in pairs:
            a, b = row["file_a"], row["file_b"]
            if a == b:
                continue
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            kept.append({"file_a": key[0], "file_b": key[1], "jaccard": row["jaccard"]})

        r["pairs"] = kept
        return r
    finally:
        os.unlink(listfile)
        for f in glob.glob(out_prefix + "*"):
            try:
                os.remove(f)
            except OSError:
                pass


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size-classes", nargs="+", default=SIZE_CLASSES,
                        choices=SIZE_CLASSES, help="Subset of size classes to run")
    parser.add_argument("--subB-values", nargs="+", type=float, default=SUBB_VALUES,
                        help="subB values to sweep")
    parser.add_argument("--reps", type=int, default=DEFAULT_REPS,
                        help="Replicate runs per (size_class, subB) cell")
    parser.add_argument("--precision", "-p", type=int, default=DEFAULT_PRECISION,
                        help="HLL precision")
    parser.add_argument("--threads", "-t", type=int, default=DEFAULT_THREADS,
                        help="Threads passed to hammock-cpp")
    parser.add_argument("--binary", default=None, help="Override hammock-cpp path")
    parser.add_argument("--output", default=None,
                        help="Output CSV path (default: results/sweep_<timestamp>.csv)")
    parser.add_argument("--verbose", action="store_true", help="Print every command")
    parser.add_argument("--smoke", action="store_true",
                        help="Smoke test: one rep at one subB on smallest class")
    args = parser.parse_args()

    if args.smoke:
        args.size_classes = ["10k"]
        args.subB_values = [0.1]
        args.reps = 1

    binary = args.binary or find_hammock_cpp()
    print(f"hammock-cpp: {binary}")
    print(f"size_classes: {args.size_classes}")
    print(f"subB values:  {args.subB_values}")
    print(f"reps:         {args.reps}")
    print(f"precision:    {args.precision}")
    print(f"threads:      {args.threads}")
    print(f"system:       {get_system_info()}")
    print()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = Path(args.output) if args.output else RESULTS_DIR / f"sweep_{stamp}.csv"

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ROW_COLS)
        w.writeheader()

        run_id = 0
        total = len(args.size_classes) * len(args.subB_values) * args.reps
        t_start = time.time()

        for size_class in args.size_classes:
            n_int = SIZE_TO_N[size_class]
            for subB in args.subB_values:
                for rep in range(args.reps):
                    run_id += 1
                    elapsed = time.time() - t_start
                    print(f"[{run_id:>3}/{total}] size={size_class} subB={subB} rep={rep+1}/{args.reps}  (elapsed {elapsed:.1f}s)",
                          flush=True)
                    r = run_one(binary, size_class, subB,
                                args.precision, args.threads, verbose=args.verbose)
                    for pair in r["pairs"]:
                        w.writerow({
                            "size_class": size_class,
                            "num_intervals": n_int,
                            "subB": subB,
                            "rep": rep,
                            "run_id": run_id,
                            "file_a": pair["file_a"],
                            "file_b": pair["file_b"],
                            "jaccard": pair["jaccard"],
                            "wall_time": r["wall_time"],
                            "cpu_time": r["cpu_time"],
                            "max_rss_mb": r["max_rss_mb"],
                            "sketch_creation_time": r["sketch_creation_time"],
                            "comparison_time": r["comparison_time"],
                            "precision": args.precision,
                            "threads": args.threads,
                            "mixed_stride": 1,
                        })
                    fh.flush()
                    print(f"      wall={r['wall_time']:.2f}s sketch={r['sketch_creation_time']}s pair={r['comparison_time']}s pairs={len(r['pairs'])}",
                          flush=True)

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
