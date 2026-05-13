#!/usr/bin/env python3
"""[Part 2 skeleton] Generate synthetic FASTA pairs for the flanking study.

For each (n_intervals, mean_length, length_dist, mutation_rate) configuration
we emit a *pair* of FASTAs (A, B) such that B is a mutated copy of A at the
specified per-base mutation rate. n_intervals records per FASTA, with
lengths drawn from `length_dist` parameterised by `mean_length`.

Why pairs? The flanking question is about how the two columns track
similarity. We need a controlled-divergence pair where the true k-mer
Jaccard is computable exactly.

Outputs (under data/synthetic/):
  config_<n>_<mean_len>_<dist>_<mut>_A.fa
  config_<n>_<mean_len>_<dist>_<mut>_B.fa
  manifest.tsv   one row per generated pair with axis values

This is a skeleton. TODO before running:
  - Decide on the size of each axis (defaults inside are conservative).
  - Wire to exact_kmer_jaccard.py for the ground truth pass.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from pathlib import Path
from typing import List

SCRIPT_DIR = Path(__file__).resolve().parent
SYN_DIR = SCRIPT_DIR / "data" / "synthetic"

BASES = "ACGT"
RNG_BASE_SEED = 20260513   # bump if you want a fresh corpus

# Default axes. Override with --axes-csv to load from a YAML/CSV instead.
DEFAULT_N_INTERVALS = [10, 100, 1000, 10000]
DEFAULT_MEAN_LEN    = [50, 150, 500, 5000]
DEFAULT_DISTS       = ["constant", "lognormal_0.5", "lognormal_1.5"]
DEFAULT_MUTATIONS   = [0.01, 0.05, 0.15, 0.30]


def random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(BASES, k=length))


def mutate(rng: random.Random, seq: str, rate: float) -> str:
    if rate <= 0:
        return seq
    out = []
    for b in seq:
        if rng.random() < rate:
            out.append(rng.choice([c for c in BASES if c != b]))
        else:
            out.append(b)
    return "".join(out)


def sample_lengths(rng: random.Random, n: int, mean: int, dist: str) -> List[int]:
    if dist == "constant":
        return [mean] * n
    if dist.startswith("lognormal_"):
        sigma = float(dist.split("_", 1)[1])
        mu = math.log(mean) - sigma * sigma / 2
        out = []
        for _ in range(n):
            v = int(round(rng.lognormvariate(mu, sigma)))
            out.append(max(1, v))
        return out
    raise ValueError(f"unknown dist: {dist}")


def write_fasta(path: Path, seqs: List[str], header_prefix: str) -> None:
    with path.open("w") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">{header_prefix}_{i}\n{s}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=SYN_DIR)
    ap.add_argument("--n-intervals", type=int, nargs="+",
                    default=DEFAULT_N_INTERVALS)
    ap.add_argument("--mean-len", type=int, nargs="+",
                    default=DEFAULT_MEAN_LEN)
    ap.add_argument("--dists", type=str, nargs="+", default=DEFAULT_DISTS)
    ap.add_argument("--mutations", type=float, nargs="+",
                    default=DEFAULT_MUTATIONS)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out_dir / "manifest.tsv"

    rows = []
    for n in args.n_intervals:
        for mlen in args.mean_len:
            for dist in args.dists:
                for mut in args.mutations:
                    label = f"n{n}_L{mlen}_{dist}_m{mut:.2f}"
                    rows.append((n, mlen, dist, mut, label))

    print(f"Will generate {len(rows)} synthetic FASTA pairs in {args.out_dir}",
          file=sys.stderr)
    if args.dry_run:
        for r in rows:
            print("  ", r, file=sys.stderr)
        return 0

    with manifest_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["n_intervals", "mean_len", "dist", "mutation",
                    "label", "fasta_A", "fasta_B",
                    "total_len_A", "total_len_B"])
        for (n, mlen, dist, mut, label) in rows:
            seed = hash((label, RNG_BASE_SEED)) & 0xFFFFFFFF
            rng = random.Random(seed)
            lengths = sample_lengths(rng, n, mlen, dist)
            seqs_A = [random_seq(rng, L) for L in lengths]
            seqs_B = [mutate(rng, s, mut) for s in seqs_A]
            total_len_A = sum(len(s) for s in seqs_A)
            total_len_B = sum(len(s) for s in seqs_B)
            fa_A = args.out_dir / f"{label}_A.fa"
            fa_B = args.out_dir / f"{label}_B.fa"
            write_fasta(fa_A, seqs_A, header_prefix=label + "_A")
            write_fasta(fa_B, seqs_B, header_prefix=label + "_B")
            w.writerow([n, mlen, dist, f"{mut:.2f}", label,
                        str(fa_A), str(fa_B),
                        total_len_A, total_len_B])
            print(f"  wrote {label}  total_len A={total_len_A} B={total_len_B}",
                  file=sys.stderr)

    print(f"Wrote {manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
