#!/usr/bin/env python3
"""[Part 2 skeleton] Compute exact k-mer set Jaccard for synthetic FASTA pairs.

Reads data/synthetic/manifest.tsv (built by generate_synthetic.py),
extracts canonical k-mers from each FASTA, computes |A∩B| / |A∪B| exactly,
and writes results/part2_synthetic_truth.tsv.

This is the *correct* ground truth for the flanking-k-mer question, where
bedtools is the wrong yardstick (boundary k-mers come from sequence content,
not interval coordinates).

Skeleton — TODO before running:
  - Decide whether to use canonical (min(kmer, revcomp(kmer))) or strand-
    sensitive k-mers. Mode D uses canonical, so default canonical = True.
  - For large k or many intervals, the set sizes explode; consider switching
    to xxhash + Python sets, or to numpy / Rust if memory becomes an issue.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterator, Set

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MANIFEST = SCRIPT_DIR / "data" / "synthetic" / "manifest.tsv"
DEFAULT_OUT      = SCRIPT_DIR / "results" / "part2_synthetic_truth.tsv"

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def iter_fasta_records(path: Path) -> Iterator[str]:
    seq_parts: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_parts:
                    yield "".join(seq_parts)
                    seq_parts = []
            else:
                seq_parts.append(line.strip())
        if seq_parts:
            yield "".join(seq_parts)


def kmer_set(path: Path, k: int, canonical: bool = True) -> Set[str]:
    """Return the set of (canonical) k-mers in the FASTA. Ns are skipped."""
    s: Set[str] = set()
    for rec in iter_fasta_records(path):
        rec = rec.upper()
        for i in range(len(rec) - k + 1):
            kmer = rec[i:i + k]
            if "N" in kmer:
                continue
            if canonical:
                rc = revcomp(kmer)
                if rc < kmer:
                    kmer = rc
            s.add(kmer)
    return s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--out",      type=Path, default=DEFAULT_OUT)
    ap.add_argument("--ks",       type=int, nargs="+", default=[8, 10, 15, 20],
                    help="k values to compute exact Jaccards at.")
    ap.add_argument("--no-canonical", action="store_true",
                    help="Use strand-sensitive k-mers (Mode D uses canonical).")
    args = ap.parse_args()

    if not args.manifest.exists():
        sys.exit(f"manifest missing: {args.manifest}; run generate_synthetic.py first")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    print(f"Computing exact k-mer Jaccards for {len(rows)} pairs × "
          f"{len(args.ks)} k values ({len(rows)*len(args.ks)} pairs total)",
          file=sys.stderr)

    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "k", "canonical", "n_intervals", "mean_len",
                    "dist", "mutation", "n_kmers_A", "n_kmers_B",
                    "n_intersection", "n_union", "jaccard"])
        for r in rows:
            fa_A = Path(r["fasta_A"])
            fa_B = Path(r["fasta_B"])
            for k in args.ks:
                set_A = kmer_set(fa_A, k, canonical=not args.no_canonical)
                set_B = kmer_set(fa_B, k, canonical=not args.no_canonical)
                inter = len(set_A & set_B)
                union = len(set_A | set_B)
                j = inter / union if union > 0 else 0.0
                w.writerow([r["label"], k, not args.no_canonical,
                            r["n_intervals"], r["mean_len"], r["dist"],
                            r["mutation"], len(set_A), len(set_B),
                            inter, union, f"{j:.6f}"])
                print(f"  {r['label']} k={k}  J={j:.4f}", file=sys.stderr)

    print(f"Wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
