#!/usr/bin/env python3
"""Single clean experiment: demonstrate when `_with_ends` adds value.

Premise (user-stated): `_with_ends` was designed to capture similarity
between strings that is *not* immediately apparent via minimizers — i.e.
shared content that the randomly-sampled minimizer set happens to miss.

Cleanest demonstration: construct pairs where the *interior* of every
record is mutated at rate `r_int` but the *first and last (k-1) bases*
of every record are kept *identical* between A and B. Interior minimizers
diverge with r_int; boundary k-mers stay perfectly shared. The exact
canonical-k-mer Jaccard between A and B reflects both effects.

Expected behaviour at (k=10, w=20, p=24):
  - Ĵ_no_ends  ↘  with r_int (minimizers see only the diverging interior)
  - Ĵ_with_ends ↘ much more slowly (boundary k-mers anchor the sketch)
  - Exact J = ground truth: sits between the two, closer to with_ends

The "value" of `with_ends` is the gap between its line and `no_ends` —
that's the shared-boundary signal it captures and `no_ends` misses.

Stays inside the known-working hammock range (k=10, w ≤ 30); see
[[project_modeD_no_ends_zero_bug]] for why we avoid k ≥ 12 / large w.

Outputs:
  results/value_demo.csv
  figures/value_demo.png (drawn by value_demo.R)
"""

from __future__ import annotations

import csv
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
HAMMOCK   = "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock"
CORPUS    = SCRIPT_DIR / "data" / "value_demo"
OUT_CSV   = SCRIPT_DIR / "results" / "value_demo.csv"

# Fixed design.
N_REPS       = 30
N_INTERVALS  = 100
MEAN_LEN     = 1000          # interior k-mers per record well above the noise floor
K            = 10            # working range of hammock Mode D
W            = 20
PRECISION    = 24
R_INTERIORS  = (0.0, 0.05, 0.10, 0.20, 0.30, 0.50)
BASES = "ACGT"

# Bring exact canonical k-mer Jaccard into this module rather than re-importing.
sys.path.insert(0, str(SCRIPT_DIR))
import exact_kmer_jaccard as ekj   # noqa: E402


def random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(BASES, k=length))


def mutate_interior(rng: random.Random, seq: str, rate: float, k: int) -> str:
    """Mutate `seq` at per-base rate `rate`, BUT preserve the first and last
    (k-1) bases exactly. Substitution-only (no indels)."""
    out = list(seq)
    L = len(seq)
    if L <= 2 * (k - 1):
        # Record too short to have an "interior"; leave it alone.
        return seq
    for i in range(k - 1, L - (k - 1)):
        if rng.random() < rate:
            b = out[i]
            out[i] = rng.choice([c for c in BASES if c != b])
    return "".join(out)


def make_pair(rep: int, r_int: float, out_dir: Path) -> tuple[Path, Path]:
    rng = random.Random(20260514 + rep * 1000 + int(r_int * 1000))
    seqs_A = [random_seq(rng, MEAN_LEN) for _ in range(N_INTERVALS)]
    seqs_B = [mutate_interior(rng, s, r_int, K) for s in seqs_A]
    label = f"rep{rep:02d}_r{int(r_int*100):03d}"
    fa_A = out_dir / f"{label}_A.fa"
    fa_B = out_dir / f"{label}_B.fa"
    with fa_A.open("w") as fh:
        for i, s in enumerate(seqs_A):
            fh.write(f">{label}_A_{i}\n{s}\n")
    with fa_B.open("w") as fh:
        for i, s in enumerate(seqs_B):
            fh.write(f">{label}_B_{i}\n{s}\n")
    return fa_A, fa_B


def run_hammock(fa_A: Path, fa_B: Path) -> tuple[float, float, float]:
    """Returns (j_no_ends, j_with_ends, wall_sec)."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        listA = td / "A.txt"; listA.write_text(str(fa_A) + "\n")
        listB = td / "B.txt"; listB.write_text(str(fa_B) + "\n")
        outprefix = td / "out"
        t0 = time.perf_counter()
        rc = subprocess.run([
            HAMMOCK,
            str(listA), str(listB),
            "--mode", "D",
            "--precision", str(PRECISION),
            "-k", str(K), "-w", str(W),
            "--threads", "4",
            "--outprefix", str(outprefix),
        ], capture_output=True, text=True)
        wall = time.perf_counter() - t0
        if rc.returncode != 0:
            raise RuntimeError(f"hammock failed: {rc.stderr[:500]}")
        csv_path = td / f"out_mnmzr_p{PRECISION}_jaccD_k{K}_w{W}.csv"
        with csv_path.open() as fh:
            for r in csv.DictReader(fh):
                if r["file1"] == fa_A.name and r["file2"] == fa_B.name:
                    return (float(r["jaccard_similarity"]),
                            float(r["jaccard_similarity_with_ends"]),
                            wall)
    raise RuntimeError("hammock produced no A-vs-B row")


def main() -> int:
    CORPUS.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    total = N_REPS * len(R_INTERIORS)
    print(f"Plan: {total} pairs at k={K}, w={W}, p={PRECISION}", file=sys.stderr)

    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rep", "r_int", "k", "w", "p",
                         "j_no_ends", "j_with_ends", "j_truth",
                         "n_kmers_A", "n_kmers_B", "wall_s"])
        fh.flush()

        for r_int in R_INTERIORS:
            for rep in range(N_REPS):
                fa_A, fa_B = make_pair(rep, r_int, CORPUS)
                # Exact canonical k-mer Jaccard as truth.
                set_A = ekj.canonical_kmer_set(fa_A, K)
                set_B = ekj.canonical_kmer_set(fa_B, K)
                inter = len(set_A & set_B)
                union = len(set_A | set_B)
                j_truth = inter / union if union > 0 else 0.0
                j_no, j_with, wall = run_hammock(fa_A, fa_B)
                writer.writerow([rep, f"{r_int:.2f}", K, W, PRECISION,
                                 f"{j_no:.6f}", f"{j_with:.6f}", f"{j_truth:.6f}",
                                 len(set_A), len(set_B), f"{wall:.3f}"])
                fh.flush()
            print(f"  r_int={r_int:.2f} done ({N_REPS} reps)", file=sys.stderr,
                  flush=True)

    print(f"Wrote {OUT_CSV}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
