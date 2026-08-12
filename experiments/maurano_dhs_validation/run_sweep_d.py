#!/usr/bin/env python3
"""Mode D (k, w, precision) sweep on the Maurano DHS FASTAs.

Default grid mirrors the original dnase1-hypersensitivity experiment:
  k         in {8, 10, 15, 20, 25}
  w         in {8, 10, 20, 30, 50, 100, 200}     (w >= k enforced)
  precision in {20, 22, 23, 24}

Each hammock invocation writes one CSV
``hammock_mnmzr_pXX_jaccD_kKK_wWW.csv`` into ``results/raw_d/``. The driver
records wall/cpu/RSS into an index CSV. analyze.R joins all of these against
``data/maurano_bedtools_ref.tsv`` for accuracy + clustering metrics.

The headline finding to look for: ``jaccard_similarity`` at k=20, w=20,
precision=24 hits Spearman 0.9997 vs bedtools on the ``_ie`` flavor.
(This note used to point at ``jaccard_similarity_with_ends``; that column was
removed in hammock v0.6.0 -- see CLAUDE.md divergence #8.)
"""

from __future__ import annotations

import argparse
import csv
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
RAW_DIR = RESULTS_DIR / "raw_d"
FASTAS = DATA_DIR / "maurano_fastas.txt"

# Refactored hammock lives in the claude-ref-comparison conda env (see
# memory/reference_new_hammock_env.md). Override with --hammock if needed.
DEFAULT_HAMMOCK = "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock"

DEFAULT_KS = (8, 10, 15, 20, 25)
DEFAULT_WS = (8, 10, 20, 30, 50, 100, 200)
DEFAULT_PS = (20, 22, 23, 24)


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


def run_with_time(cmd: List[str], log_path: Path) -> Tuple[float, float, float, int]:
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
                    rss = float(m.group(1)) / 1024.0
    return wall, cpu, rss, proc.returncode


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--hammock", type=str, default="",
                   help=f"Path to hammock binary (default: {DEFAULT_HAMMOCK} "
                        f"if it exists, else first on PATH).")
    p.add_argument("--ks", type=int, nargs="+", default=list(DEFAULT_KS))
    p.add_argument("--ws", type=int, nargs="+", default=list(DEFAULT_WS))
    p.add_argument("--ps", type=int, nargs="+", default=list(DEFAULT_PS))
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not FASTAS.exists():
        sys.exit(f"{FASTAS} missing — run ./prepare_data.sh first (need FASTAs).")

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    logs_dir = SCRIPT_DIR / "logs"
    logs_dir.mkdir(exist_ok=True)

    hammock = find_hammock(args.hammock)
    print(f"Using hammock: {hammock}", file=sys.stderr)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_path = RESULTS_DIR / f"sweep_d_{stamp}.csv"

    plan: List[Tuple[int, int, int]] = []
    for k in args.ks:
        for w in args.ws:
            if w < k:
                continue
            for p in args.ps:
                plan.append((k, w, p))

    print(f"Planning {len(plan)} Mode D invocations on {FASTAS}.", file=sys.stderr)
    if args.dry_run:
        for k, w, p in plan:
            print(f"  k={k} w={w} p={p} -> hammock_mnmzr_p{p}_jaccD_k{k}_w{w}_full.csv",
                  file=sys.stderr)
        return 0

    outprefix = RAW_DIR / "hammock"
    with summary_path.open("w", newline="") as fh:
        wcsv = csv.writer(fh)
        wcsv.writerow(["k", "w", "precision", "wall_s", "cpu_s",
                       "max_rss_mb", "csv_path", "rc"])
        for k, w, p in plan:
            # "_full" matches the --metrics flag below (python/hammock/
            # outprefix.py now always tags Python-CLI output _ie/_re/_full;
            # docs/seed-metrics-column-restructure.md) -- analyze.R's
            # parse_d_name requires this exact suffix.
            outname = f"hammock_mnmzr_p{p}_jaccD_k{k}_w{w}_full.csv"
            outpath = RAW_DIR / outname
            cmd = [
                hammock,
                str(FASTAS), str(FASTAS),
                "--mode", "D",
                "--precision", str(p),
                "-k", str(k),
                "-w", str(w),
                "--threads", str(args.threads),
                "--outprefix", str(outprefix),
                # analyze.R's Mode D path reads the full containment/cosketch
                # block (with a fallback that reconstructs jaccard_similarity_ie
                # from containments if absent) -- needs the full shape.
                "--metrics",
            ]
            log_path = logs_dir / f"d_{outname}.log"
            print(f"[run] {' '.join(shlex.quote(c) for c in cmd)}", file=sys.stderr)
            wall, cpu, rss, rc = run_with_time(cmd, log_path)
            if rc != 0:
                print(f"  FAILED rc={rc} (see {log_path})", file=sys.stderr)
            else:
                print(f"  wall={wall:.1f}s cpu={cpu:.1f}s rss={rss:.0f}MB -> {outpath}",
                      file=sys.stderr)
            wcsv.writerow([k, w, p, f"{wall:.3f}", f"{cpu:.3f}", f"{rss:.1f}",
                           str(outpath), rc])

    print(f"Wrote {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
