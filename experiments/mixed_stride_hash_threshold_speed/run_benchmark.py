#!/usr/bin/env python3
"""Paired, randomized timing of mixed-stride versus hash-threshold."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "bedtools_benchmark"))
from benchmark_cpp_vs_bedtools import (  # noqa: E402
    PAIRONLY_RE,
    PAIR_RE,
    SKETCH_RE,
    WRITE_RE,
    check_binary_version,
    find_hammock_cpp,
    get_system_info,
    run_with_time,
)


METHODS = ("mixed-stride", "hash-threshold")
RATES = (0.1, 0.06, 0.03, 0.01, 0.006, 0.003, 0.001)
COUNT_RE = re.compile(
    r"^Processed\s+\d+\s+intervals,\s+(\d+)/(\d+)\s+points\s+"
    r"\(subB=[^,]+,\s+method=([^)]+)\)"
)
FIELDS = (
    "round", "order_in_round", "method", "subB", "wall_seconds",
    "cpu_seconds", "max_rss_mb", "sketch_seconds", "pairwise_write_seconds",
    "pairwise_seconds", "write_seconds", "sampled_points", "total_points",
    "realized_fraction", "output_sha256",
)


def write_input_manifest(path: Path, beds: list[Path]) -> None:
    """Record portable input identities independently of machine-local paths."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("basename", "size_bytes", "sha256"),
            lineterminator="\n",
        )
        writer.writeheader()
        for bed in beds:
            writer.writerow({
                "basename": bed.name,
                "size_bytes": bed.stat().st_size,
                "sha256": sha256_path(bed),
            })


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_stderr(stderr: str, expected_method: str) -> dict[str, float | int]:
    sketch = pair_write = pair = write = None
    sampled = total = 0
    count_lines = 0
    for line in stderr.splitlines():
        match = SKETCH_RE.match(line)
        if match:
            sketch = int(match.group(1)) / 1e6
            continue
        match = PAIR_RE.match(line)
        if match:
            pair_write = int(match.group(1)) / 1e6
            continue
        match = PAIRONLY_RE.match(line)
        if match:
            pair = int(match.group(1)) / 1e6
            continue
        match = WRITE_RE.match(line)
        if match:
            write = int(match.group(1)) / 1e6
            continue
        match = COUNT_RE.match(line)
        if match:
            if match.group(3) != expected_method:
                raise RuntimeError(
                    f"method provenance mismatch: expected {expected_method}, "
                    f"stderr reported {match.group(3)}"
                )
            sampled += int(match.group(1))
            total += int(match.group(2))
            count_lines += 1
    if sketch is None or pair_write is None or pair is None or write is None:
        raise RuntimeError(f"could not parse all phase timings from stderr:\n{stderr}")
    if count_lines == 0 or total <= 0:
        raise RuntimeError(f"could not parse sampled-point counts from stderr:\n{stderr}")
    return {
        "sketch_seconds": sketch,
        "pairwise_write_seconds": pair_write,
        "pairwise_seconds": pair,
        "write_seconds": write,
        "sampled_points": sampled,
        "total_points": total,
    }


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def run_cell(binary: Path, list_file: Path, method: str, rate: float,
             precision: int, threads: int, seed: int, gate_seed: int,
             env: dict[str, str],
             archive_path: Path | None) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="hammock-ms-ht-") as temporary:
        prefix = Path(temporary) / "result"
        command = [
            str(binary), str(list_file), str(list_file),
            "--mode", "B", "--subB", str(rate), "--subB-method", method,
            "--precision", str(precision), "--threads", str(threads),
            "--seed", str(seed), "--gate-seed", str(gate_seed),
            "--output", str(prefix), "--verbose",
        ]
        measured = run_with_time(command, env=env)
        outputs = sorted(Path(temporary).glob("result*.csv"))
        if len(outputs) != 1:
            raise RuntimeError(f"expected one hammock output, got {outputs}")
        output = outputs[0]
        if archive_path is not None:
            archive_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(output, archive_path)
        parsed = parse_stderr(measured["stderr"], method)
        return {
            "wall_seconds": measured["wall_time"],
            "cpu_seconds": measured["cpu_time"],
            "max_rss_mb": measured["max_rss_mb"],
            **parsed,
            "realized_fraction": parsed["sampled_points"] / parsed["total_points"],
            "output_sha256": sha256_path(output),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--data-dir", type=Path,
                        default=REPO_ROOT / "experiments/subB_mixed_stride/data/maurano")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--rounds", type=int, default=10)
    parser.add_argument("--warmup-rounds", type=int, default=1)
    parser.add_argument("--precision", type=int, default=18)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gate-seed", type=int, default=31337)
    parser.add_argument("--shuffle-seed", type=int, default=20260901)
    parser.add_argument("--allow-shared", action="store_true",
                        help="Allow a non-Slurm/non-exclusive plumbing run; timings are not publishable")
    parser.add_argument("--smoke", action="store_true",
                        help="One rate, one measured round, no warm-up; implies --allow-shared")
    args = parser.parse_args()

    if args.rounds < 1 or args.warmup_rounds < 0:
        parser.error("--rounds must be positive and --warmup-rounds nonnegative")
    if args.smoke:
        args.rounds = 1
        args.warmup_rounds = 0
        args.allow_shared = True
    if not args.allow_shared and not os.environ.get("SLURM_JOB_ID"):
        parser.error("publishable timings require a Slurm allocation; use --allow-shared only for plumbing checks")

    binary = (args.binary or Path(find_hammock_cpp())).resolve()
    version = check_binary_version(str(binary))
    beds = sorted(args.data_dir.resolve().glob("*.bed"))
    if len(beds) != 20:
        raise RuntimeError(f"expected 20 Maurano BED files in {args.data_dir}, got {len(beds)}")

    stamp = datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    output_dir = (args.output_dir or SCRIPT_DIR / "results" / stamp).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    list_file = output_dir / "files.txt"
    list_file.write_text("".join(f"{path.resolve()}\n" for path in beds))
    write_input_manifest(output_dir / "input_manifest.csv", beds)

    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=REPO_ROOT,
        check=True, capture_output=True, text=True,
    ).stdout.strip())
    rates = (RATES[0],) if args.smoke else RATES
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_sha": git_sha,
        "git_dirty": dirty,
        "binary": str(binary),
        "binary_version": version,
        "binary_sha256": sha256_path(binary),
        "methods": METHODS,
        "rates": rates,
        "rounds": args.rounds,
        "warmup_rounds": args.warmup_rounds,
        "precision": args.precision,
        "threads": args.threads,
        "seed": args.seed,
        "gate_seed": args.gate_seed,
        "shuffle_seed": args.shuffle_seed,
        "publishable_allocation_required": not args.allow_shared,
        "system": get_system_info(),
    }
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    env = dict(os.environ)
    env.update({
        "OMP_NUM_THREADS": str(args.threads),
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": "close",
        "OMP_PLACES": "cores",
    })
    cells = [(method, rate) for method in METHODS for rate in rates]

    for warmup in range(args.warmup_rounds):
        order = list(cells)
        random.Random(args.shuffle_seed - 1 - warmup).shuffle(order)
        print(f"warm-up round {warmup + 1}/{args.warmup_rounds}", flush=True)
        for method, rate in order:
            run_cell(binary, list_file, method, rate, args.precision,
                     args.threads, args.seed, args.gate_seed, env,
                     archive_path=None)

    rows: list[dict[str, object]] = []
    raw_path = output_dir / "raw_runs.csv"
    for round_index in range(args.rounds):
        order = list(cells)
        random.Random(args.shuffle_seed + round_index).shuffle(order)
        print(f"measured round {round_index + 1}/{args.rounds}", flush=True)
        for order_index, (method, rate) in enumerate(order):
            rate_label = format(rate, ".6g").replace(".", "p")
            archive = None
            if round_index == 0:
                archive = output_dir / "estimates" / f"{method}_subB_{rate_label}.tsv"
            measured = run_cell(binary, list_file, method, rate, args.precision,
                                args.threads, args.seed, args.gate_seed, env,
                                archive)
            row = {
                "round": round_index,
                "order_in_round": order_index,
                "method": method,
                "subB": rate,
                **measured,
            }
            rows.append(row)
            write_rows(raw_path, rows)
            print(
                f"  {method:14s} subB={rate:g} "
                f"sketch={measured['sketch_seconds']:.3f}s "
                f"wall={measured['wall_seconds']:.3f}s",
                flush=True,
            )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
