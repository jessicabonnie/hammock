#!/usr/bin/env python3
"""Validate and summarize the corrected cross-reference Mode D rerun."""
from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


SEEDS = (1, 42, 99)
PEAK_TYPES = ("broad", "narrow")
EXPECTED_COLUMNS = {
    "file1", "file2", "mode", "precision", "kmer_size", "window_size",
    "jaccard_similarity_ie", "containment_AB", "containment_BA",
}


def parse_label(path: str) -> tuple[str, str]:
    p = Path(path)
    return p.stem, p.parent.name


def read_cells(path: Path) -> list[tuple[int, int]]:
    with path.open() as handle:
        cells = [(int(r["k"]), int(r["w"])) for r in csv.DictReader(handle, delimiter="\t")]
    if len(cells) != 63 or len(set(cells)) != 63:
        raise ValueError(f"expected 63 unique cells, found {len(cells)} rows/{len(set(cells))} unique")
    return cells


def read_tissues(path: Path) -> dict[str, str]:
    with path.open() as handle:
        return {r["sample_id"]: r["tissue"] for r in csv.DictReader(handle, delimiter="\t")}


def validate_matrix(path: Path, k: int, w: int) -> list[dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle)
        missing = EXPECTED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows = list(reader)
    if len(rows) != 81:
        raise ValueError(f"{path}: expected 81 rows, found {len(rows)}")
    values: dict[tuple[str, str], float] = {}
    for row in rows:
        if row["mode"] != "D" or int(row["precision"]) != 24:
            raise ValueError(f"{path}: wrong mode/precision")
        if int(row["kmer_size"]) != k or int(row["window_size"]) != w:
            raise ValueError(f"{path}: wrong k/w metadata")
        value = float(row["jaccard_similarity_ie"])
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{path}: invalid Jaccard value {value}")
        values[(row["file1"], row["file2"])] = value
    if len(values) != 81:
        raise ValueError(f"{path}: duplicate file pairs")
    for (a, b), value in values.items():
        other = values.get((b, a))
        if other is None or not math.isclose(value, other, abs_tol=1e-12):
            raise ValueError(f"{path}: asymmetric pair {a}, {b}")
        if a == b and not math.isclose(value, 1.0, abs_tol=1e-12):
            raise ValueError(f"{path}: non-unit diagonal for {a}")
    return rows


def cell_stats(rows: list[dict[str, str]], tissues: dict[str, str]) -> dict[str, float | int | bool]:
    xref: dict[str, list[float]] = {"reg": [], "ie": []}
    diff: dict[str, list[float]] = {"reg": [], "ie": []}
    for row in rows:
        sample1, ref1 = parse_label(row["file1"])
        sample2, ref2 = parse_label(row["file2"])
        if sample1 == sample2 and ref1 == ref2:
            continue
        values = {
            "reg": float(row.get("reg_eq_similarity", row.get("jaccard_similarity", "nan"))),
            "ie": float(row["jaccard_similarity_ie"]),
        }
        if sample1 == sample2 and ref1 != ref2:
            for tag, value in values.items():
                xref[tag].append(value)
        elif tissues[sample1] != tissues[sample2]:
            for tag, value in values.items():
                diff[tag].append(value)
    if (len(xref["ie"]), len(diff["ie"])) != (18, 54):
        raise ValueError(
            f"expected 18/54 ordered comparisons, found {len(xref['ie'])}/{len(diff['ie'])}")
    out: dict[str, float | int | bool] = {
        "n_xref": len(xref["ie"]), "n_diff": len(diff["ie"]),
    }
    for tag in ("reg", "ie"):
        out[f"median_xref_{tag}"] = statistics.median(xref[tag])
        out[f"median_diff_{tag}"] = statistics.median(diff[tag])
        out[f"delta_{tag}"] = statistics.median(xref[tag]) - statistics.median(diff[tag])
        out[f"separated_{tag}"] = min(xref[tag]) > max(diff[tag])
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    cells = read_cells(args.cells)
    tissues = read_tissues(args.metadata)
    rows_out: list[dict[str, object]] = []
    for seed in SEEDS:
        for peak_type in PEAK_TYPES:
            for k, w in cells:
                path = (args.results / f"seed{seed}" / peak_type / f"k{k}_w{w}" /
                        f"exp_a_mnmzr_p24_jaccD_k{k}_w{w}_rehash-selector64_full.csv")
                if not path.is_file():
                    raise FileNotFoundError(path)
                stats = cell_stats(validate_matrix(path, k, w), tissues)
                rows_out.append({"seed": seed, "peak_type": peak_type, "cell": f"k{k}_w{w}",
                                 "k": k, "w": w, **stats, "matrix": str(path)})

    ranks: dict[tuple[int, str, str], int] = {}
    grouped: dict[tuple[int, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows_out:
        grouped[(int(row["seed"]), str(row["peak_type"]))].append(row)
    for (seed, peak_type), group in grouped.items():
        for rank, row in enumerate(sorted(group, key=lambda r: (-float(r["delta_ie"]), str(r["cell"]))), 1):
            ranks[(seed, peak_type, str(row["cell"]))] = rank
            row["rank_ie"] = rank

    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_seed = args.output_dir / "per_seed_cells.csv"
    with per_seed.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows_out[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)

    sensitivity_rows: list[dict[str, object]] = []
    by_cell: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows_out:
        by_cell[(str(row["peak_type"]), str(row["cell"]))].append(row)
    for (peak_type, cell), group in sorted(by_cell.items()):
        deltas = [float(r["delta_ie"]) for r in group]
        cell_ranks = [int(r["rank_ie"]) for r in group]
        sensitivity_rows.append({
            "peak_type": peak_type, "cell": cell, "k": group[0]["k"], "w": group[0]["w"],
            "delta_ie_min": min(deltas), "delta_ie_median": statistics.median(deltas),
            "delta_ie_max": max(deltas), "delta_ie_range": max(deltas) - min(deltas),
            "rank_best": min(cell_ranks), "rank_worst": max(cell_ranks),
            "all_seeds_separated": all(bool(r["separated_ie"]) for r in group),
        })
    sensitivity = args.output_dir / "seed_sensitivity.csv"
    with sensitivity.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(sensitivity_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(sensitivity_rows)

    seed42 = args.output_dir / "exp_a_estimator_delta_expanded_rehash_seed42.csv"
    seed42_rows = [r for r in rows_out if r["seed"] == 42]
    with seed42.open("w", newline="") as handle:
        fields = ["peak_type", "cell", "k", "w", "n_xref", "n_diff",
                  "median_xref_reg", "median_diff_reg", "delta_reg", "separated_reg",
                  "median_xref_ie", "median_diff_ie", "delta_ie", "separated_ie"]
        writer = csv.DictWriter(
            handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(seed42_rows)

    for peak_type in PEAK_TYPES:
        leaders = sorted((r for r in seed42_rows if r["peak_type"] == peak_type),
                         key=lambda r: int(r["rank_ie"]))[:3]
        print(f"{peak_type} seed42 leaders: " + ", ".join(
            f"{r['cell']} delta={float(r['delta_ie']):.6f}" for r in leaders))
    print(f"validated {len(rows_out)} matrices; wrote {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
