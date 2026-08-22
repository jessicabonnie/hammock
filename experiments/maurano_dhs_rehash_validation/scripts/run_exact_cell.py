#!/usr/bin/env python3
"""Dispatch one zero-based exact-grid cell for a bounded scheduler array."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from common import EXPERIMENT, grid, load_config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXPERIMENT / "config.yaml")
    parser.add_argument("--cell-index", type=int, required=True)
    args = parser.parse_args()
    config, _ = load_config(args.config)
    cells = grid(config)
    if not 0 <= args.cell_index < len(cells):
        raise SystemExit(f"cell index out of range: {args.cell_index}")
    k, w = cells[args.cell_index]
    command = [sys.executable, str(EXPERIMENT / "scripts" / "extract_exact_features.py"),
               "--config", str(args.config.resolve()), "--cell", str(k), str(w)]
    print(f"exact cell index={args.cell_index} k={k} w={w}", flush=True)
    return subprocess.run(command).returncode


if __name__ == "__main__":
    raise SystemExit(main())
