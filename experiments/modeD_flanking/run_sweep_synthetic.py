#!/usr/bin/env python3
"""[Part 2 skeleton] Mode D sweep on the synthetic FASTA pairs.

For each pair listed in data/synthetic/manifest.tsv, run hammock Mode D
across the (k, w, p) grid and capture both Jaccard columns into a tidy
long-form CSV at results/part2_sweep.csv.

Skeleton — fill in the grid you actually want before launching, or just
let it inherit from workflow/config.yaml when we wire Snakemake up.

Each row in the output CSV:
  label, n_intervals, mean_len, dist, mutation, k, w, p,
  j_no_ends, j_with_ends, wall_s
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import List, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_HAMMOCK = "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock"
MANIFEST = SCRIPT_DIR / "data" / "synthetic" / "manifest.tsv"
OUT      = SCRIPT_DIR / "results" / "part2_sweep.csv"

DEFAULT_KS = (8, 10, 15, 20)
DEFAULT_WS = (10, 20, 50, 100, 200, 500)
DEFAULT_PS = (16, 20, 24)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hammock", type=str, default=DEFAULT_HAMMOCK)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--ks", type=int, nargs="+", default=list(DEFAULT_KS))
    ap.add_argument("--ws", type=int, nargs="+", default=list(DEFAULT_WS))
    ap.add_argument("--ps", type=int, nargs="+", default=list(DEFAULT_PS))
    ap.add_argument("--dry-run", action="store_true")
    return ap.parse_args()


def parse_mode_d_csv(path: Path) -> List[Tuple[str, str, float, float]]:
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            if "jaccard_similarity_with_ends" not in r:
                raise SystemExit(
                    "This experiment requires hammock <= 0.5.0: it compares "
                    "jaccard_similarity against jaccard_similarity_with_ends, "
                    "and v0.6.0 removed the latter (CLAUDE.md divergence #8, "
                    "docs/mode-d-ends-removal.md). The archived results under "
                    "results/ were produced with the old schema and stand; this "
                    "driver cannot be re-run against current hammock."
                )
            rows.append((
                r["file1"], r["file2"],
                float(r["jaccard_similarity"]),
                float(r["jaccard_similarity_with_ends"]),
            ))
    return rows


def main() -> int:
    args = parse_args()
    if not MANIFEST.exists():
        sys.exit(f"manifest missing: {MANIFEST}; run generate_synthetic.py first")
    if not Path(args.hammock).exists():
        sys.exit(f"hammock not found: {args.hammock}")

    with MANIFEST.open() as fh:
        pairs = list(csv.DictReader(fh, delimiter="\t"))

    plan = [(p, k, w) for p in args.ps for k in args.ks for w in args.ws if w >= k]
    print(f"Plan: {len(pairs)} pairs × {len(plan)} (k,w,p) configs = "
          f"{len(pairs)*len(plan)} hammock invocations", file=sys.stderr)
    if args.dry_run:
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "n_intervals", "mean_len", "dist", "mutation",
                    "k", "w", "p", "j_no_ends", "j_with_ends", "wall_s"])

        for r in pairs:
            label = r["label"]
            fa_A = Path(r["fasta_A"])
            fa_B = Path(r["fasta_B"])
            # Hammock wants list-files. Build a single-line listfile per FASTA.
            with tempfile.TemporaryDirectory() as td:
                td = Path(td)
                listA = td / "A.txt"; listA.write_text(str(fa_A) + "\n")
                listB = td / "B.txt"; listB.write_text(str(fa_B) + "\n")
                for (p, k, ws) in plan:
                    outprefix = td / f"out_{label}_p{p}_k{k}_w{ws}"
                    t0 = time.perf_counter()
                    rc = subprocess.run([
                        args.hammock,
                        str(listA), str(listB),
                        "--mode", "D",
                        "--precision", str(p),
                        "-k", str(k), "-w", str(ws),
                        "--threads", str(args.threads),
                        "--outprefix", str(outprefix),
                    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
                    wall = time.perf_counter() - t0
                    if rc != 0:
                        print(f"  FAILED {label} p={p} k={k} w={ws}", file=sys.stderr)
                        continue
                    csv_path = Path(str(outprefix) + f"_mnmzr_p{p}_jaccD_k{k}_w{ws}.csv")
                    rows = parse_mode_d_csv(csv_path)
                    for (f1, f2, j_n, j_e) in rows:
                        # Only the A-vs-B row is interesting (self pairs trivially 1.0).
                        if f1 == fa_A.name and f2 == fa_B.name:
                            w.writerow([label, r["n_intervals"], r["mean_len"],
                                        r["dist"], r["mutation"], k, ws, p,
                                        f"{j_n:.6f}", f"{j_e:.6f}",
                                        f"{wall:.3f}"])
            print(f"  done {label}", file=sys.stderr)

    print(f"Wrote {OUT}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
