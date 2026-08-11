#!/usr/bin/env python3
"""Summarize a --metrics (jaccard_similarity_ie) synthetic subB sweep.

Counterpart to the Maurano summary behind docs/data/maurano_subB_ie_summary.csv,
but wall-time only: the synthetic corpus has no exact bedtools per-pair truth
file the way Maurano does via docs/data/maurano_bedtools_ref.tsv (see
docs/seed-subsampling-synthetic-supplement.md), so there is no MAE column here
-- this measures the subB speedup on hammock's own wall time, which is what the
outline/draft paragraph actually claims.

Usage: summarize_synthetic_ie.py <sweep_synthetic_ie_*.csv> [output.csv]
"""
import csv
import statistics
import sys
from collections import defaultdict

IN = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else None

rows_by_cell = defaultdict(list)  # (method, size_class, subB) -> [wall_time, ...]
n_pairs_by_cell = defaultdict(set)
threads_by_cell = {}

with open(IN, newline="") as fh:
    for r in csv.DictReader(fh):
        if not r["metrics"] or r["metrics"] != "True":
            continue
        key = (r["method"], r["size_class"], float(r["subB"]))
        # wall_time is attached identically to every pair row within a run;
        # dedupe by run_id so each replicate contributes once.
        n_pairs_by_cell[key].add((r["file_a"], r["file_b"]))
        threads_by_cell[key] = r["threads"]

# Re-walk to pull one wall_time per (cell, run_id) -- avoids double counting
# across the per-pair rows within a run.
seen_run = defaultdict(set)
with open(IN, newline="") as fh:
    for r in csv.DictReader(fh):
        if not r["metrics"] or r["metrics"] != "True":
            continue
        key = (r["method"], r["size_class"], float(r["subB"]))
        run_key = (key, r["run_id"])
        if run_key in seen_run:
            continue
        seen_run[run_key] = True
        rows_by_cell[key].append(float(r["wall_time"]))

out_rows = []
nosub_wall = {}
for key, walls in rows_by_cell.items():
    method, size_class, subB = key
    if subB == 1.0:
        nosub_wall[(method, size_class)] = statistics.median(walls)

for key, walls in sorted(rows_by_cell.items(), key=lambda kv: (kv[0][1], -kv[0][2])):
    method, size_class, subB = key
    med = statistics.median(walls)
    mean = statistics.mean(walls)
    sd = statistics.stdev(walls) if len(walls) > 1 else 0.0
    base = nosub_wall.get((method, size_class))
    out_rows.append({
        "method": method,
        "size_class": size_class,
        "subB": subB,
        "n_reps": len(walls),
        "threads": threads_by_cell[key],
        "wall_median": med,
        "wall_mean": mean,
        "wall_sd": sd,
        "n_pairs": len(n_pairs_by_cell[key]),
        "wall_nosub": base,
        "speedup_vs_nosub": (base / med) if base else "",
    })

fieldnames = ["method", "size_class", "subB", "n_reps", "threads",
              "wall_median", "wall_mean", "wall_sd", "n_pairs",
              "wall_nosub", "speedup_vs_nosub"]

out_fh = open(OUT, "w", newline="") if OUT else sys.stdout
w = csv.DictWriter(out_fh, fieldnames=fieldnames)
w.writeheader()
for row in out_rows:
    w.writerow(row)
    print(row, file=sys.stderr)
if OUT:
    out_fh.close()
    print(f"wrote {OUT}", file=sys.stderr)
