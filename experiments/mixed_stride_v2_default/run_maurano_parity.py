#!/usr/bin/env python3
"""Maurano parity, timing, and accuracy gate for the mixed-stride v2 default."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import time
from datetime import datetime
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DEFAULT_DATA = REPO_ROOT / "experiments" / "subB_mixed_stride" / "data" / "maurano"
RATES = (0.3, 0.1, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001)
INTEGRAL_PARITY_RATES = frozenset((0.1, 0.01, 0.008, 0.005, 0.001))
METHODS = ("mixed-stride", "mixed-stride-v1")


def find_binary() -> Path:
    hits = list((REPO_ROOT / "build").glob("*/hammock-cpp"))
    if not hits:
        raise FileNotFoundError("hammock-cpp not found under build/*/")
    return max(hits, key=lambda path: path.stat().st_mtime)


def read_pairs(path: Path) -> dict[tuple[str, str], float]:
    pairs: dict[tuple[str, str], float] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "jaccard_similarity_ie" not in (reader.fieldnames or []):
            raise RuntimeError(f"missing jaccard_similarity_ie in {path}")
        for row in reader:
            left = Path(row["query"]).name
            right = Path(row["reference"]).name
            if left == right:
                continue
            key = tuple(sorted((left, right)))
            value = float(row["jaccard_similarity_ie"])
            previous = pairs.setdefault(key, value)
            if previous != value:
                raise RuntimeError(f"asymmetric pairwise result for {key}: {previous} != {value}")
    return pairs


def run_cell(binary: Path, file_list: Path, run_dir: Path, rate: float,
             method: str, precision: int, threads: int) -> tuple[Path, float]:
    label = f"{method}_subB_{rate:g}".replace(".", "p")
    prefix = run_dir / "raw" / label
    cmd = [
        str(binary), str(file_list), str(file_list), "--mode", "B",
        "--precision", str(precision), "--threads", str(threads),
        "--subB", repr(rate), "--subB-method", method,
        "--output", str(prefix), "--verbose",
    ]
    started = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    wall_seconds = time.perf_counter() - started
    (run_dir / "raw" / f"{label}.stderr.txt").write_text(proc.stderr)
    if proc.returncode:
        raise RuntimeError(f"hammock-cpp failed for {method} subB={rate}:\n{proc.stderr}")
    written = [line[6:].strip() for line in proc.stderr.splitlines()
               if line.startswith("Wrote ")]
    if len(written) != 1:
        raise RuntimeError(f"expected one output path for {label}, got {written}")
    return Path(written[0]), wall_seconds


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=None)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--precision", type=int, default=18)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--results-dir", type=Path, default=SCRIPT_DIR / "results")
    args = parser.parse_args()

    binary = (args.binary or find_binary()).resolve()
    beds = sorted(args.data_dir.resolve().glob("*.bed"))
    if len(beds) != 20:
        raise SystemExit(f"expected 20 Maurano BEDs in {args.data_dir}, found {len(beds)}")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = args.results_dir / f"run_{stamp}"
    (run_dir / "raw").mkdir(parents=True)
    file_list = run_dir / "maurano_files.txt"
    file_list.write_text("".join(f"{path.resolve()}\n" for path in beds))

    version = subprocess.run([str(binary), "--version"], check=True,
                             capture_output=True, text=True).stdout.strip()
    print(f"binary: {binary} ({version})", flush=True)
    print(f"output: {run_dir}", flush=True)

    baseline_path, baseline_wall = run_cell(
        binary, file_list, run_dir, 1.0, "mixed-stride",
        args.precision, args.threads)
    baseline = read_pairs(baseline_path)
    if len(baseline) != 190:
        raise RuntimeError(f"expected 190 unique Maurano pairs, got {len(baseline)}")

    rows = [{
        "subB": 1.0, "method": "mixed-stride", "wall_seconds": baseline_wall,
        "pair_count": len(baseline), "mae_vs_subB_1": 0.0,
        "max_abs_error_vs_subB_1": 0.0, "exactly_equals_other_method": "",
        "unequal_pairs_vs_other_method": "", "mae_vs_other_method": "",
        "max_abs_error_vs_other_method": "",
    }]
    results: dict[tuple[float, str], dict[tuple[str, str], float]] = {}

    for rate in RATES:
        for method in METHODS:
            output, wall = run_cell(binary, file_list, run_dir, rate, method,
                                    args.precision, args.threads)
            pairs = read_pairs(output)
            if pairs.keys() != baseline.keys():
                raise RuntimeError(f"pair set differs from baseline for {method} subB={rate}")
            errors = [abs(pairs[key] - baseline[key]) for key in baseline]
            results[(rate, method)] = pairs
            rows.append({
                "subB": rate, "method": method, "wall_seconds": wall,
                "pair_count": len(pairs), "mae_vs_subB_1": math.fsum(errors) / len(errors),
                "max_abs_error_vs_subB_1": max(errors),
                "exactly_equals_other_method": "pending",
                "unequal_pairs_vs_other_method": "pending",
                "mae_vs_other_method": "pending",
                "max_abs_error_vs_other_method": "pending",
            })
            print(f"subB={rate:<6g} {method:<16} wall={wall:8.3f}s "
                  f"MAE={rows[-1]['mae_vs_subB_1']:.6g}", flush=True)

        public = results[(rate, "mixed-stride")]
        legacy = results[(rate, "mixed-stride-v1")]
        unequal = [key for key in public if public[key] != legacy[key]]
        method_errors = [abs(public[key] - legacy[key]) for key in public]
        equal = not unequal
        for row in rows[-2:]:
            row["exactly_equals_other_method"] = str(equal).lower()
            row["unequal_pairs_vs_other_method"] = len(unequal)
            row["mae_vs_other_method"] = math.fsum(method_errors) / len(method_errors)
            row["max_abs_error_vs_other_method"] = max(method_errors)
        print(f"  v1/v2: unequal={len(unequal)}/{len(public)} "
              f"MAE={rows[-1]['mae_vs_other_method']:.6g}", flush=True)
        if rate in INTEGRAL_PARITY_RATES and unequal:
            raise RuntimeError(
                f"v1/v2 parity failed at subB={rate}: {len(unequal)} unequal pairs; "
                f"first={unequal[0]} ({public[unequal[0]]} != {legacy[unequal[0]]})")
        if rate not in INTEGRAL_PARITY_RATES and not unequal:
            raise RuntimeError(
                f"v1/v2 unexpectedly produced identical pairwise output at "
                f"non-integral reciprocal rate subB={rate}")

    summary = run_dir / "summary.csv"
    fields = list(rows[0])
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"all integral-rate parity checks passed; wrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
