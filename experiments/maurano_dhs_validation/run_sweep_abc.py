#!/usr/bin/env python3
"""A/B/C sweep on the Maurano DHS corpus, vs bedtools ground truth.

Sweep grid (default):
  - Mode A, precisions {18, 21, 23}
  - Mode B, precisions {18, 21, 23}                    (mixed-stride, subB=1.0)
  - Mode C, precision  21, expA in {0.1, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0}
  - Mode C, precision  21, subB in {0.01, 0.05, 0.1, 0.25, 0.5, 0.8}

Each hammock invocation emits its own CSV into ``results/raw_abc/``; this
driver only writes a tiny index CSV summarising (config, wall, cpu, RSS,
output_path). analyze.R is responsible for joining the raw CSVs against
the bedtools reference.

Outputs (relative to this experiment dir):
  results/raw_abc/hammock_hll_pXX_jaccM[_expA*][_BXX].csv  -- one per config
  results/sweep_abc_<stamp>.csv                            -- index

Use --dry-run to print the full plan without launching hammock.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR / "results"
RAW_DIR = RESULTS_DIR / "raw_abc"
FILES_LIST = DATA_DIR / "maurano_files.txt"

# Refactored hammock lives in the claude-ref-comparison conda env (see
# memory/reference_new_hammock_env.md). Override with --hammock if needed.
DEFAULT_HAMMOCK = "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock"

# (precision, expA, subB) — None means "leave at default 0 / 1.0"
DEFAULT_GRID: List[Tuple[str, int, float, float]] = []
for p in (18, 21, 23):
    DEFAULT_GRID.append(("A", p, 0.0, 1.0))
    DEFAULT_GRID.append(("B", p, 0.0, 1.0))
for expA in (0.1, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0):
    DEFAULT_GRID.append(("C", 21, expA, 1.0))
for subB in (0.01, 0.05, 0.1, 0.25, 0.5, 0.8):
    DEFAULT_GRID.append(("C", 21, 0.0, subB))


def find_hammock(override: str = "") -> str:
    if override:
        if not Path(override).exists():
            sys.exit(f"--hammock={override} does not exist")
        return override
    if Path(DEFAULT_HAMMOCK).exists():
        return DEFAULT_HAMMOCK
    exe = shutil.which("hammock")
    if not exe:
        sys.exit(f"hammock not found at {DEFAULT_HAMMOCK} or on PATH; "
                 f"activate claude-ref-comparison or pass --hammock")
    return exe


def output_csv_name(mode: str, precision: int, expA: float, subB: float) -> str:
    name = f"hammock_hll_p{precision}_jacc{mode}"
    if expA > 0:
        name += f"_expA{expA:.2f}"
    if subB != 1.0:
        name += f"_B{subB:.2f}"
    # "_full" matches the --metrics flag added below (python/hammock/
    # outprefix.py now always tags Python-CLI output _ie/_re/_full;
    # docs/seed-metrics-column-restructure.md) -- analyze.R's parse_abc_name
    # requires this exact suffix.
    return name + "_full.csv"


def run_with_time(cmd: List[str], log_path: Path) -> Tuple[float, float, float, int]:
    """Run cmd; return (wall_s, cpu_s, max_rss_mb, rc). Uses /usr/bin/time -v if available."""
    time_bin = "/usr/bin/time"
    use_time = Path(time_bin).exists()
    with log_path.open("w") as logf:
        t0 = time.perf_counter()
        if use_time:
            proc = subprocess.run(
                [time_bin, "-v", *cmd], stdout=logf, stderr=subprocess.STDOUT,
            )
        else:
            proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
        wall = time.perf_counter() - t0
    cpu = -1.0
    rss = -1.0
    if use_time:
        with log_path.open() as fh:
            for line in fh:
                m = re.search(r"User time .*: ([\d.]+)", line)
                if m:
                    cpu = float(m.group(1))
                m = re.search(r"Maximum resident set size .*: (\d+)", line)
                if m:
                    rss = float(m.group(1)) / 1024.0  # KB -> MB
    return wall, cpu, rss, proc.returncode


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--hammock", type=str, default="",
                   help=f"Path to hammock binary (default: {DEFAULT_HAMMOCK} "
                        f"if it exists, else first on PATH).")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the planned invocations; don't launch hammock.")
    p.add_argument("--subset", type=str, default="",
                   help="Comma-separated mode filter (e.g. 'A,B' or 'C'). Default all.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not FILES_LIST.exists():
        sys.exit(f"{FILES_LIST} missing — run ./prepare_data.sh first")

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    logs_dir = SCRIPT_DIR / "logs"
    logs_dir.mkdir(exist_ok=True)

    hammock = find_hammock(args.hammock)
    print(f"Using hammock: {hammock}", file=sys.stderr)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_path = RESULTS_DIR / f"sweep_abc_{stamp}.csv"

    wanted_modes = set(args.subset.split(",")) if args.subset else None
    plan = [g for g in DEFAULT_GRID if (wanted_modes is None or g[0] in wanted_modes)]

    print(f"Planning {len(plan)} hammock invocations on {FILES_LIST}.", file=sys.stderr)
    if args.dry_run:
        for mode, prec, expA, subB in plan:
            print(f"  mode={mode} p={prec} expA={expA} subB={subB} -> "
                  f"{output_csv_name(mode, prec, expA, subB)}", file=sys.stderr)
        return 0

    with summary_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mode", "precision", "expA", "subB", "wall_s", "cpu_s",
                    "max_rss_mb", "csv_path", "rc"])
        for mode, prec, expA, subB in plan:
            outname = output_csv_name(mode, prec, expA, subB)
            outpath = RAW_DIR / outname
            outprefix = RAW_DIR / "hammock"
            cmd = [
                hammock,
                str(FILES_LIST), str(FILES_LIST),
                "--mode", mode,
                "--precision", str(prec),
                "--threads", str(args.threads),
                "--outprefix", str(outprefix),
                # analyze.R reads the full containment/cosketch block, not
                # just jaccard_similarity -- needs the full shape, not the
                # new bare default or --register-equality.
                "--metrics",
            ]
            if mode == "C":
                cmd += ["--expA", str(expA), "--subB", str(subB)]
            log_path = logs_dir / f"abc_{outname}.log"
            print(f"[run] {' '.join(shlex.quote(c) for c in cmd)}", file=sys.stderr)
            wall, cpu, rss, rc = run_with_time(cmd, log_path)
            if rc != 0:
                print(f"  FAILED rc={rc} (see {log_path})", file=sys.stderr)
            else:
                print(f"  wall={wall:.1f}s cpu={cpu:.1f}s rss={rss:.0f}MB -> {outpath}",
                      file=sys.stderr)
            w.writerow([mode, prec, expA, subB, f"{wall:.3f}", f"{cpu:.3f}",
                        f"{rss:.1f}", str(outpath), rc])

    print(f"Wrote {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
