#!/usr/bin/env python3
"""Time pairwise bedtools jaccard on the Maurano corpus.

Reference baseline for the subB-method sweep: how fast is bedtools vs
hammock at any subB on real data?

Pre-sorts the 20 Maurano BEDs once into data/maurano_sorted/ (sort time is
not charged to bedtools, matching bedtools_benchmark conventions). Then
times N replicate batches of pairwise `bedtools jaccard` across all 190
unique pairs, using GNU parallel for the 8-way fan-out.

Writes results/sweep_maurano_bedtools_<stamp>.csv with the same schema as
run_sweep.py output, so analyze.R can consume it as a peer dataset.
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "bedtools_benchmark"))
from benchmark_cpp_vs_bedtools import run_with_time  # noqa: E402

DATA_DIR = SCRIPT_DIR / "data"
SORTED_DIR = DATA_DIR / "maurano_sorted"
RESULTS_DIR = SCRIPT_DIR / "results"
BEDTOOLS_SCRIPT = SCRIPT_DIR.parent / "bedtools_benchmark" / "bedtools.sh"

ROW_COLS = [
    "corpus", "method", "size_class", "num_intervals", "subB", "rep", "run_id",
    "file_a", "file_b", "jaccard",
    "wall_time", "cpu_time", "max_rss_mb",
    "sketch_creation_time", "comparison_time",
    "precision", "threads",
]


def sort_maurano() -> List[Path]:
    """Sort the symlinked Maurano BEDs into data/maurano_sorted/ once."""
    src_dir = DATA_DIR / "maurano"
    srcs = sorted(src_dir.glob("*.bed"))
    if not srcs:
        raise FileNotFoundError(
            f"no BEDs in {src_dir}; run prepare_maurano.sh first"
        )
    SORTED_DIR.mkdir(parents=True, exist_ok=True)
    sorted_paths: List[Path] = []
    for src in srcs:
        dst = SORTED_DIR / src.name
        if dst.exists() and dst.stat().st_mtime >= src.resolve().stat().st_mtime:
            sorted_paths.append(dst)
            continue
        print(f"sort  {src.name}", flush=True)
        with open(dst, "w") as fh:
            subprocess.run(["sort", "-k1,1", "-k2,2n", str(src)],
                           stdout=fh, check=True)
        sorted_paths.append(dst)
    return sorted_paths


def parse_bedtools_output(stdout: str) -> Dict[tuple, float]:
    """bedtools.sh emits: <path_a>\\t<path_b>\\t<jaccard>. Dedupe to unordered pairs."""
    out: Dict[tuple, float] = {}
    for line in stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        a, b, j = parts
        a, b = os.path.basename(a), os.path.basename(b)
        if a == b:
            continue
        key = tuple(sorted((a, b)))
        if key in out:
            continue
        try:
            out[key] = float(j)
        except ValueError:
            continue
    return out


def run_one_batch(sorted_beds: List[Path], threads: int,
                  verbose: bool = False) -> Dict[str, Any]:
    """Time one full 190-pair bedtools jaccard batch."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as fh:
        listfile = fh.name
        for p in sorted_beds:
            fh.write(str(p) + "\n")
    try:
        cmd = ["bash", str(BEDTOOLS_SCRIPT), listfile, listfile, str(threads)]
        if verbose:
            print("  +", " ".join(cmd), file=sys.stderr)
        r = run_with_time(cmd)
        pairs_dict = parse_bedtools_output(r["stdout"])
        r["pairs"] = [{"file_a": k[0], "file_b": k[1], "jaccard": v}
                      for k, v in pairs_dict.items()]
        return r
    finally:
        os.unlink(listfile)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=5,
                        help="Replicate batches (default 5)")
    parser.add_argument("--threads", "-t", type=int, default=8,
                        help="GNU parallel jobs (default 8)")
    parser.add_argument("--output", default=None,
                        help="Output CSV path (default: results/sweep_maurano_bedtools_<stamp>.csv)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    sorted_beds = sort_maurano()
    print(f"\nbedtools on {len(sorted_beds)} pre-sorted BEDs, "
          f"{len(sorted_beds)*(len(sorted_beds)-1)//2} unique pairs, "
          f"{args.threads} threads, {args.reps} reps\n", flush=True)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = Path(args.output) if args.output \
        else RESULTS_DIR / f"sweep_maurano_bedtools_{stamp}.csv"

    t_start = time.time()
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ROW_COLS)
        w.writeheader()
        for rep in range(args.reps):
            elapsed = time.time() - t_start
            print(f"[{rep+1}/{args.reps}] (elapsed {elapsed:.1f}s)", flush=True)
            r = run_one_batch(sorted_beds, args.threads, verbose=args.verbose)
            print(f"      wall={r['wall_time']:.2f}s pairs={len(r['pairs'])}",
                  flush=True)
            for pair in r["pairs"]:
                w.writerow({
                    "corpus": "maurano",
                    "method": "bedtools",
                    "size_class": "maurano",
                    "num_intervals": "",
                    "subB": 1.0,           # placeholder; bedtools has no subsample
                    "rep": rep,
                    "run_id": rep + 1,
                    "file_a": pair["file_a"],
                    "file_b": pair["file_b"],
                    "jaccard": pair["jaccard"],
                    "wall_time": r["wall_time"],
                    "cpu_time": r["cpu_time"],
                    "max_rss_mb": r["max_rss_mb"],
                    "sketch_creation_time": "",
                    "comparison_time": "",
                    "precision": "",
                    "threads": args.threads,
                })
            fh.flush()

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
