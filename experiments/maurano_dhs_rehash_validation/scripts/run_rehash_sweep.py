#!/usr/bin/env python3
"""Plan or run isolated, atomic, resumable Mode D rehash jobs."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from common import EXPERIMENT, atomic_json, fasta_paths, grid, load_config, sha256, strip_fasta, write_csv
from validation import validate_hll_csv


def parse_time_report(path: Path) -> dict[str, float]:
    text = path.read_text()
    patterns = {
        "user_cpu_s": r"User time \(seconds\): ([0-9.]+)",
        "system_cpu_s": r"System time \(seconds\): ([0-9.]+)",
        "maximum_rss_kb": r"Maximum resident set size \(kbytes\): ([0-9.]+)",
    }
    parsed = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            raise ValueError(f"GNU time report lacks {name}: {path}")
        parsed[name] = float(match.group(1))
    parsed["cpu_s"] = parsed["user_cpu_s"] + parsed["system_cpu_s"]
    return parsed


def normalize_csv_lf(path: Path) -> bool:
    """Normalize Python csv.writer CRLF for tracked, diff-checkable artifacts."""
    data = path.read_bytes()
    normalized = data.replace(b"\r\n", b"\n")
    if normalized == data:
        return False
    path.write_bytes(normalized)
    return True


def source_commit(hammock: Path) -> str:
    executable = hammock.resolve()
    # The isolated console script lives in ENV/bin; source is supplied explicitly
    # so callers cannot silently use the historical shared environment.
    result = subprocess.run([str(executable), "--version"], check=True,
                            text=True, capture_output=True)
    return result.stdout.strip() or result.stderr.strip()


def followup_cells(config) -> list[tuple[int, int]]:
    leaders_path = EXPERIMENT / config["followup"]["leaders_file"]
    if not leaders_path.is_file():
        raise SystemExit(f"follow-up leaders file missing: {leaders_path}")
    leaders = json.loads(leaders_path.read_text())
    if leaders.get("status") != "frozen":
        raise SystemExit("follow-up is gated until p=18 leaders have status=frozen")
    records = [*config["followup"]["historical_cells"],
               *leaders["exact_leaders"], *leaders["robustness_leaders"]]
    cells = sorted({(int(record["k"]), int(record["w"])) for record in records})
    invalid = set(cells) - set(grid(config))
    if invalid:
        raise SystemExit(f"follow-up leaders are outside the frozen grid: {sorted(invalid)}")
    return cells


def requested_seeds(spec) -> list[int]:
    requested = set(range(int(spec["seed_start"]), int(spec["seed_stop"]) + 1))
    requested.update(map(int, spec.get("additional_seeds", [])))
    return sorted(requested)


def extension_seeds(config) -> list[int]:
    """Return only the new seeds, leaving the frozen eight-seed phases intact."""
    return sorted(set(requested_seeds(config["extension"])) - set(map(int, config["seeds"])))


def jobs(config, base: Path, hammock: Path, commit: str, phase: str):
    list_path = (base / config["inputs"]["fasta_list"]).resolve()
    short_commit = commit[:12]
    if phase == "primary":
        precisions = [int(config["primary_precision"])]
        cells = grid(config)
        seeds = list(map(int, config["seeds"]))
    elif phase == "followup":
        # p=18 is already complete for every possible leader in the primary sweep.
        precisions = [int(value) for value in config["precisions"]
                      if int(value) != int(config["primary_precision"])]
        cells = followup_cells(config)
        seeds = list(map(int, config["seeds"]))
    elif phase == "extension":
        precisions = list(map(int, config["extension"]["precisions"]))
        cells = [(int(cell["k"]), int(cell["w"])) for cell in config["extension"]["cells"]]
        seeds = extension_seeds(config)
    else:
        precisions = list(map(int, config["interpolation"]["precisions"]))
        cells = [(int(cell["k"]), int(cell["w"])) for cell in config["interpolation"]["cells"]]
        seeds = requested_seeds(config["interpolation"])
    resources = (config[phase].get("resources", config["resources"])
                 if phase in {"extension", "interpolation"} else config["resources"])
    for precision in precisions:
        for seed in seeds:
            for k, w in cells:
                job_id = f"rehash_p{precision}_seed{seed:05d}_k{k}_w{w}_{short_commit}"
                result = (EXPERIMENT / "results" / "rehash" / f"p{precision}" /
                          f"seed{seed:05d}" / f"k{k}_w{w}_{short_commit}.csv")
                outprefix = "hammock"
                command = [str(hammock.resolve()), str(list_path), str(list_path),
                           "--mode", "D", "--sequence-hll-hash", config["hash_scheme"],
                           "--seed", str(seed), "--precision", str(precision),
                           "-k", str(k), "-w", str(w), "--threads", "1", "--metrics",
                           "--outprefix", outprefix]
                yield {"job_id": job_id, "phase": phase, "hash_scheme": config["hash_scheme"],
                       "source_commit": commit, "seed": seed, "precision": precision,
                       "k": k, "w": w, "cpus": resources["cpus_per_job"],
                       "memory_mb": resources["memory_mb"],
                       "time_minutes": resources["time_minutes"],
                       "output": str(result), "command": command}


def plan_fingerprint(planned: list[dict]) -> str:
    payload = [{key: value for key, value in job.items() if key != "command"} |
               {"command": job["command"]} for job in planned]
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def quarantine(path: Path) -> Path:
    destination = (EXPERIMENT / "results" / "rehash" / "quarantine" /
                   f"{path.name}.{dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(path, destination)
    return destination


def run_job(job, fastas: list[Path], samples: list[str]) -> None:
    output = Path(job["output"])
    completion = output.with_suffix(".complete.json")
    if output.exists() and completion.exists():
        validate_hll_csv(output, samples, k=job["k"], w=job["w"], precision=job["precision"])
        recorded = json.loads(completion.read_text())
        identity_keys = ("phase", "hash_scheme", "source_commit", "seed", "precision", "k", "w")
        if (recorded.get("sha256") == sha256(output) and
                all(recorded.get(key) == job[key] for key in identity_keys)):
            print(f"[skip validated] {job['job_id']}")
            return
    quarantined = []
    for stale in (output, completion):
        if stale.exists():
            quarantined.append(str(quarantine(stale)))
    output.parent.mkdir(parents=True, exist_ok=True)
    log_dir = EXPERIMENT / "logs" / job["job_id"]
    log_dir.mkdir(parents=True, exist_ok=True)
    started = dt.datetime.now(dt.timezone.utc)
    with tempfile.TemporaryDirectory(prefix=f"{job['job_id']}.", dir=output.parent) as raw_stage:
        stage = Path(raw_stage)
        command = [*job["command"]]
        command[-1] = str(stage / "hammock")
        timed = ["/usr/bin/time", "-v", "-o", str(log_dir / "time.txt"), *command]
        before = time.monotonic()
        with (log_dir / "stdout.txt").open("wb") as stdout, (log_dir / "stderr.txt").open("wb") as stderr:
            process = subprocess.run(timed, stdout=stdout, stderr=stderr)
        elapsed = time.monotonic() - before
        timing = parse_time_report(log_dir / "time.txt")
        status = {**{key: job[key] for key in ("job_id", "phase", "hash_scheme", "source_commit", "seed", "precision", "k", "w")},
                  "command": command, "started_utc": started.isoformat(), "elapsed_s": elapsed,
                  "returncode": process.returncode, "quarantined": quarantined,
                  "timing": timing, "stdout_log": str(log_dir / "stdout.txt"),
                  "stderr_log": str(log_dir / "stderr.txt"), "time_log": str(log_dir / "time.txt")}
        if process.returncode != 0:
            atomic_json(log_dir / "failed.json", status)
            raise RuntimeError(f"{job['job_id']} failed with rc={process.returncode}")
        candidates = list(stage.glob("*.csv"))
        if len(candidates) != 1:
            status["csv_candidates"] = [str(path) for path in candidates]
            atomic_json(log_dir / "failed.json", status)
            raise RuntimeError(f"expected one staged CSV, found {len(candidates)}")
        validation = validate_hll_csv(candidates[0], samples, k=job["k"], w=job["w"], precision=job["precision"])
        raw_cli_sha256 = sha256(candidates[0])
        normalized_line_endings = normalize_csv_lf(candidates[0])
        os.replace(candidates[0], output)
    status.update({"status": "complete", "output": str(output), "sha256": sha256(output),
                   "raw_cli_sha256": raw_cli_sha256,
                   "tracked_normalization": "CRLF converted to LF" if normalized_line_endings else "none",
                   "validation": validation, "input_fastas": [str(path) for path in fastas]})
    atomic_json(completion, status)
    print(f"[complete] {job['job_id']} ({elapsed:.1f}s)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXPERIMENT / "config.yaml")
    parser.add_argument("--hammock", type=Path, required=True,
                        help="console script installed from this disposable clone")
    parser.add_argument("--source-commit", required=True,
                        help="full git commit of the disposable source build")
    parser.add_argument("--phase", choices=["primary", "followup", "extension", "interpolation"],
                        default="primary")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--plan-sha256",
                        help="required for execution; fingerprint printed by the reviewed dry run")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--job-id", help="run exactly one planned job")
    selection.add_argument("--job-index", type=int,
                           help="run one zero-based row from the frozen job plan")
    args = parser.parse_args()
    config, base = load_config(args.config)
    if not args.hammock.is_file():
        raise SystemExit(f"isolated hammock executable missing: {args.hammock}")
    repo = EXPERIMENT.parents[1]
    actual_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    if args.source_commit != actual_commit:
        source_paths = ["python", "cpp", "CMakeLists.txt", "pyproject.toml"]
        comparison = subprocess.run(
            ["git", "diff", "--quiet", args.source_commit, actual_commit, "--", *source_paths],
            cwd=repo,
        )
        if comparison.returncode != 0:
            raise SystemExit(
                f"source commit mismatch in production code: requested {args.source_commit}, "
                f"checkout {actual_commit}"
            )
    planned = list(jobs(config, base, args.hammock, args.source_commit, args.phase))
    if args.phase == "primary" and len(planned) != 296:
        raise SystemExit("primary job-plan gate failed: expected 296 jobs")
    if args.phase == "extension":
        expected = (len(config["extension"]["precisions"]) *
                    len(config["extension"]["cells"]) * len(extension_seeds(config)))
        if len(planned) != expected:
            raise SystemExit(f"extension job-plan gate failed: expected {expected} jobs")
    if args.phase == "interpolation":
        expected = (len(config["interpolation"]["precisions"]) *
                    len(config["interpolation"]["cells"]) *
                    len(requested_seeds(config["interpolation"])))
        if len(planned) != expected:
            raise SystemExit(f"interpolation job-plan gate failed: expected {expected} jobs")
    if not planned or len({job["output"] for job in planned}) != len(planned):
        raise SystemExit(f"{args.phase} job-plan gate failed: outputs are empty or collide")
    fingerprint = plan_fingerprint(planned)
    rows = [{**{key: value for key, value in job.items() if key != "command"},
             "command": shlex.join(job["command"])} for job in planned]
    if args.dry_run:
        write_csv(EXPERIMENT / "results" / "manifests" / f"rehash_{args.phase}_jobs.csv",
                  ["job_id", "phase", "hash_scheme", "source_commit", "seed", "precision", "k", "w",
                   "cpus", "memory_mb", "time_minutes", "output", "command"], rows)
        print(f"planned {len(planned)} unique {args.phase} jobs; plan_sha256={fingerprint}; "
              "scheduler submission not performed")
        return 0
    if args.plan_sha256 != fingerprint:
        raise SystemExit(f"reviewed plan fingerprint mismatch: expected {fingerprint}")
    if args.job_index is not None:
        if not 0 <= args.job_index < len(planned):
            raise SystemExit(f"job index out of range: {args.job_index}")
        selected = [planned[args.job_index]]
    elif args.job_id:
        selected = [job for job in planned if job["job_id"] == args.job_id]
        if len(selected) != 1:
            raise SystemExit(f"unknown or ambiguous job id: {args.job_id}")
    else:
        raise SystemExit("refusing bulk local run: pass --job-id or --job-index after reviewing the manifest")
    fastas = fasta_paths(config, base)
    samples = [strip_fasta(str(path)) for path in fastas]
    run_job(selected[0], fastas, samples)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
