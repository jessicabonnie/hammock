#!/usr/bin/env python3
"""subB x subB-method x precision grid on the 20 Maurano BEDs, reading the IE column.

Answers: does --subB subsampling degrade `jaccard_similarity_ie`?

No production code is touched. hammock-cpp has emitted jaccard_similarity_ie
natively and by default since v0.7.0, and exact per-pair truth already exists in
docs/data/maurano_bedtools_ref.tsv, so this is just an invocation harness.

The output name deliberately does NOT begin with `sweep_`: analyze.R,
headline_figures.R and pareto_variants.R all pick_latest on ^sweep_ globs and
would adopt any such file as the newest sweep.
"""
from __future__ import annotations

import csv
import glob
import itertools
import os
import subprocess
import sys
import tempfile
import time

REPO = "/home/jbonnie1/interval_sketch/hammock_claude"
BINARY = os.path.join(REPO, "build/cp310-cp310-linux_x86_64/hammock-cpp")
BEDS = sorted(glob.glob(os.path.join(REPO, "experiments/subB_mixed_stride/data/maurano/*.bed")))

PRECISIONS = [18, 23]
SUBBS = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01]
METHODS = ["mixed-stride", "hash-threshold", "single-hash"]
THREADS = 8

OUT = sys.argv[1] if len(sys.argv) > 1 else "ie_maurano.csv"

COLS = ["jaccard_similarity", "jaccard_similarity_ie",
        "containment_AB", "containment_BA"]


def run_cell(p: int, sub_b: float, method: str, tmp: str):
    listfile = os.path.join(tmp, "files.txt")
    with open(listfile, "w") as fh:
        fh.write("\n".join(BEDS) + "\n")
    prefix = os.path.join(tmp, "out")
    cmd = [BINARY, listfile, listfile, "--mode", "B", "-p", str(p),
           "-t", str(THREADS), "-o", prefix,
           "--subB", repr(sub_b), "--subB-method", method,
           # Needs the full block: reads jaccard_similarity, jaccard_similarity_ie,
           # containment_AB, containment_BA (COLS below), none of which the
           # bare/no-flag default (jaccard_similarity_ie alone) or
           # --register-equality (jaccard_similarity alone) would fully supply.
           "--metrics"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"hammock-cpp failed (p={p} subB={sub_b} {method}):\n{r.stderr[-2000:]}")
    wall = time.time() - t0
    produced = [f for f in glob.glob(prefix + "*") if f.endswith(".csv")]
    if len(produced) != 1:
        raise SystemExit(f"expected 1 CSV, got {produced}")
    rows = []
    with open(produced[0]) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in COLS if c not in (rd.fieldnames or [])]
        if missing:
            raise SystemExit(f"output lacks {missing}; header={rd.fieldnames}")
        for row in rd:
            a, b = os.path.basename(row["query"]), os.path.basename(row["reference"])
            if a == b:
                continue
            rec = {"precision": p, "subB": sub_b, "method": method,
                   "file_a": a, "file_b": b}
            rec.update({c: row[c] for c in COLS})
            rows.append(rec)
    for f in glob.glob(prefix + "*"):
        os.remove(f)
    return rows, wall


def main() -> int:
    if len(BEDS) != 20:
        raise SystemExit(f"expected 20 Maurano BEDs, found {len(BEDS)}")
    # subB=1.0 short-circuits the gate, so all three methods are identical there
    # (hammock-cpp --help says so explicitly). Run it once and label it.
    cells = [(p, 1.0, "none") for p in PRECISIONS]
    cells += [(p, s, m) for p, s, m in itertools.product(
        PRECISIONS, [s for s in SUBBS if s < 1.0], METHODS)]

    all_rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for i, (p, s, m) in enumerate(cells, 1):
            meth = "mixed-stride" if m == "none" else m
            rows, wall = run_cell(p, s, meth, tmp)
            for r in rows:
                r["method"] = m
            all_rows.extend(rows)
            print(f"[{i:>2}/{len(cells)}] p={p} subB={s:<5g} {m:<14} "
                  f"{len(rows)} pairs  {wall:6.1f}s", flush=True)

    fields = ["precision", "subB", "method", "file_a", "file_b"] + COLS
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nWrote {OUT}  ({len(all_rows)} rows, {len(cells)} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
