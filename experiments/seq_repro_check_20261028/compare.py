#!/usr/bin/env python3
"""Compare new vs archived Mode D CSVs on estimate columns (read-only on archives)."""
import csv
import glob
import math
import os
import sys

CHECK_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_ROOT = "/vast/blangme2/jbonnie/hammock-seq-repro-20261028/results"
KEY = ["file1", "file2"]
ESTIMATES = ["reg_eq_similarity", "jaccard_similarity_ie", "containment_AB", "containment_BA",
             "cosketch_geom", "cosketch_arith", "cosketch_max"]
META = ["sketch_type", "mode", "precision", "num_hashes", "kmer_size", "window_size"]


def load(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {(r["file1"], r["file2"]): r for r in rows}, rows


def main():
    jobs = list(csv.DictReader(open(os.path.join(CHECK_DIR, "jobs.tsv")), delimiter="\t"))
    out_rows = []
    for job in jobs:
        pattern = os.path.join(RESULT_ROOT, f"task{int(job['task']):02d}_{job['label']}_*", "*.csv")
        found = glob.glob(pattern)
        rec = {k: job[k] for k in ("task", "label", "runtime", "precision", "seed", "k", "w")}
        if len(found) != 1:
            rec.update(status=f"missing({len(found)})", rows="", n_estimate_cells="",
                       exact_text_equal="", max_abs_diff="")
            out_rows.append(rec)
            continue
        new, new_list = load(found[0])
        old, old_list = load(job["archived_csv"])
        same_keys = set(new) == set(old) and len(new_list) == len(old_list)
        same_order = [tuple(r[c] for c in KEY) for r in new_list] == [tuple(r[c] for c in KEY) for r in old_list]
        text_equal = True
        max_diff = 0.0
        cells = 0
        meta_equal = True
        for key in old:
            if key not in new:
                continue
            for col in META:
                meta_equal &= new[key][col] == old[key][col]
            for col in ESTIMATES:
                a, b = new[key][col], old[key][col]
                cells += 1
                if a != b:
                    text_equal = False
                    fa, fb = float(a), float(b)
                    if math.isnan(fa) and math.isnan(fb):
                        continue
                    max_diff = max(max_diff, abs(fa - fb))
        status = "IDENTICAL" if (same_keys and same_order and text_equal and meta_equal) else "DIFFERS"
        rec.update(status=status, rows=len(new_list), n_estimate_cells=cells,
                   exact_text_equal=text_equal and meta_equal and same_keys and same_order,
                   max_abs_diff=f"{max_diff:.3g}", new_csv=found[0], archived_csv=job["archived_csv"])
        out_rows.append(rec)
    fields = ["task", "label", "runtime", "precision", "seed", "k", "w", "status", "rows",
              "n_estimate_cells", "exact_text_equal", "max_abs_diff", "new_csv", "archived_csv"]
    dest = os.path.join(CHECK_DIR, "comparison.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    for r in out_rows:
        print("\t".join(str(r.get(c, "")) for c in fields[:12]))
    return 0 if all(r["status"] == "IDENTICAL" for r in out_rows) else 1


if __name__ == "__main__":
    sys.exit(main())
