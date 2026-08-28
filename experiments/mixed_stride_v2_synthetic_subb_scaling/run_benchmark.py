#!/usr/bin/env python3
"""Matched-corpus BEDTools, register-equality, and +IE synthetic benchmark."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import glob
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "bedtools_benchmark"))

from benchmark_cpp_vs_bedtools import (  # noqa: E402
    _rotate,
    bedtools_serial_ms,
    check_binary_version,
    derive_seed,
    generate_bed_file,
    get_system_info,
    run_bedtools,
    run_hammock,
)


N_VALUES = (2, 4, 8, 16, 32, 64, 128, 256, 512)
SUBB_VALUES = (1.0, 0.1, 0.01, 0.001)
FIELDS = (
    "num_files", "rep", "tool", "similarity", "subB", "num_pairs",
    "num_intervals_per_file", "precision", "threads", "corpus_seed",
    "wall_time", "cpu_time", "max_rss_mb", "sort_time",
    "sketch_creation_time", "comparison_time", "pair_time", "write_time",
    "bedtools_serial_ms", "bedtools_parallel_eff",
    "mae_vs_bedtools", "max_err_vs_bedtools",
    "hostname", "cpu_model", "git_sha", "slurm_job_id",
)


def sort_one(path: Path) -> None:
    sorted_path = path.with_suffix(path.suffix + ".sorted")
    with sorted_path.open("w") as output:
        subprocess.run(["sort", "-k1,1", "-k2,2n", str(path)],
                       stdout=output, check=True)
    sorted_path.replace(path)


def parse_bedtools(stdout: str) -> dict[tuple[str, str], float]:
    values = {}
    for line in stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        a, b, value = parts
        values[(Path(a).name, Path(b).name)] = float(value)
    return values


def parse_hammock(path: Path, column: str) -> dict[tuple[str, str], float]:
    values = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            values[(Path(row["query"]).name, Path(row["reference"]).name)] = float(row[column])
    return values


def clean_hammock_output(result: dict) -> None:
    prefix = result.get("_out_prefix")
    if prefix:
        for path in glob.glob(prefix + "*"):
            try:
                os.remove(path)
            except OSError:
                pass


def expected_reps(n: int, small_runs: int, large_runs: int) -> int:
    return small_runs if n <= 32 else large_runs


def completed(path: Path, n: int, rep: int, args: argparse.Namespace) -> bool:
    if not path.exists():
        return False
    try:
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        return (
            len(rows) == 9
            and {int(row["num_files"]) for row in rows} == {n}
            and {int(row["rep"]) for row in rows} == {rep}
            and {int(row["num_intervals_per_file"]) for row in rows} == {args.num_intervals}
            and {int(row["precision"]) for row in rows} == {args.precision}
            and {int(row["threads"]) for row in rows} == {args.threads}
            and {int(row["corpus_seed"]) for row in rows} == {args.corpus_seed}
        )
    except (OSError, ValueError, KeyError):
        return False


def base_row(n: int, rep: int, args: argparse.Namespace, system: dict) -> dict:
    return {
        "num_files": n,
        "rep": rep,
        "num_pairs": n * n,
        "num_intervals_per_file": args.num_intervals,
        "precision": args.precision,
        "threads": args.threads,
        "corpus_seed": args.corpus_seed,
        "hostname": system["hostname"],
        "cpu_model": system["cpu_model"],
        "git_sha": system["git_sha"],
        "slurm_job_id": system["slurm_job_id"],
    }


def run_replicate(n: int, rep: int, args: argparse.Namespace, system: dict) -> list[dict]:
    with tempfile.TemporaryDirectory() as tmp_name:
        tmp = Path(tmp_name)
        queries, references = [], []
        for index in range(n):
            query = tmp / f"set1_{index}.bed"
            reference = tmp / f"set2_{index}.bed"
            generate_bed_file(args.num_intervals, str(query),
                              derive_seed(args.corpus_seed, rep, index, 0))
            generate_bed_file(args.num_intervals, str(reference),
                              derive_seed(args.corpus_seed, rep, index, 1))
            queries.append(query)
            references.append(reference)

        sort_start = time.time()
        paths = queries + references
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=min(args.threads, len(paths))) as executor:
            list(executor.map(sort_one, paths))
        sort_time = time.time() - sort_start

        query_list = tmp / "queries.txt"
        reference_list = tmp / "references.txt"
        query_list.write_text("\n".join(map(str, queries)) + "\n")
        reference_list.write_text("\n".join(map(str, references)) + "\n")

        arms = [("register-equality", rate) for rate in SUBB_VALUES]
        arms += [("inclusion-exclusion", rate) for rate in SUBB_VALUES]
        execution = [("bedtools", None)] + arms
        results = {}
        estimates = {}
        truth = None

        for similarity, rate in _rotate(execution, rep):
            if similarity == "bedtools":
                print("    BEDTools", flush=True)
                result = run_bedtools(str(query_list), str(reference_list), args.threads)
                truth = parse_bedtools(result["stdout"])
                serial_ms = bedtools_serial_ms(str(query_list), str(reference_list))
                result["bedtools_serial_ms"] = serial_ms
                result["bedtools_parallel_eff"] = (
                    (n * n * serial_ms / 1000.0) / (result["wall_time"] * args.threads)
                )
                results[(similarity, rate)] = result
            else:
                print(f"    {similarity} subB={rate:g}", flush=True)
                result = run_hammock(
                    args.binary,
                    str(query_list),
                    str(reference_list),
                    args.precision,
                    args.threads,
                    sub_b=rate,
                    ie_only=(similarity == "inclusion-exclusion"),
                    keep_output=True,
                )
                try:
                    column = ("reg_eq_similarity" if similarity == "register-equality"
                              else "jaccard_similarity_ie")
                    estimates[(similarity, rate)] = parse_hammock(
                        Path(result["output_csv"]), column)
                finally:
                    clean_hammock_output(result)
                results[(similarity, rate)] = result

        if truth is None or len(truth) != n * n:
            raise RuntimeError(f"BEDTools returned {0 if truth is None else len(truth)} of {n*n} pairs")

        rows = []
        common = base_row(n, rep, args, system)
        bedtools_result = results[("bedtools", None)]
        rows.append({
            **common,
            "tool": "bedtools",
            "similarity": "exact",
            "subB": "",
            "wall_time": bedtools_result["wall_time"],
            "cpu_time": bedtools_result["cpu_time"],
            "max_rss_mb": bedtools_result["max_rss_mb"],
            "sort_time": sort_time,
            "bedtools_serial_ms": bedtools_result["bedtools_serial_ms"],
            "bedtools_parallel_eff": bedtools_result["bedtools_parallel_eff"],
            "mae_vs_bedtools": 0.0,
            "max_err_vs_bedtools": 0.0,
        })

        for similarity, rate in arms:
            result = results[(similarity, rate)]
            values = estimates[(similarity, rate)]
            if values.keys() != truth.keys():
                missing = len(truth.keys() - values.keys())
                extra = len(values.keys() - truth.keys())
                raise RuntimeError(
                    f"pair-key mismatch for {similarity} subB={rate}: missing={missing}, extra={extra}"
                )
            errors = [abs(values[key] - truth[key]) for key in truth]
            rows.append({
                **common,
                "tool": "hammock",
                "similarity": similarity,
                "subB": rate,
                "wall_time": result["wall_time"],
                "cpu_time": result["cpu_time"],
                "max_rss_mb": result["max_rss_mb"],
                "sort_time": sort_time,
                "sketch_creation_time": result["sketch_creation_time"],
                "comparison_time": result["comparison_time"],
                "pair_time": result["pair_time"],
                "write_time": result["write_time"],
                "mae_vs_bedtools": sum(errors) / len(errors),
                "max_err_vs_bedtools": max(errors),
            })
        return rows


def write_atomic(path: Path, rows: list[dict]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--num-files", default=",".join(map(str, N_VALUES)))
    parser.add_argument("--small-runs", type=int, default=20)
    parser.add_argument("--large-runs", type=int, default=3)
    parser.add_argument("--num-intervals", type=int, default=10_000)
    parser.add_argument("--precision", type=int, default=18)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--corpus-seed", type=int, default=20_260_810)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    check_binary_version(args.binary)
    system = get_system_info()
    n_values = tuple(int(value) for value in args.num_files.split(","))
    if any(value not in N_VALUES for value in n_values):
        raise ValueError(f"N must be selected from {N_VALUES}")

    print(f"system: {system}")
    print(f"N values: {n_values}; subB values: {SUBB_VALUES}")
    print(f"p={args.precision}; threads={args.threads}; seed={args.corpus_seed}")
    for n in n_values:
        reps = expected_reps(n, args.small_runs, args.large_runs)
        for rep in range(reps):
            output = args.output_dir / f"N{n}_rep{rep:02d}.csv"
            if completed(output, n, rep, args):
                print(f"[{n=:>5} rep={rep+1:>2}/{reps}] complete; skipping", flush=True)
                continue
            print(f"[{n=:>5} rep={rep+1:>2}/{reps}] generating matched corpus", flush=True)
            rows = run_replicate(n, rep, args, system)
            write_atomic(output, rows)
            print(f"    wrote {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
