#!/usr/bin/env python3
"""Validate exact and HLL products without trusting filenames alone."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import EXPERIMENT, fasta_paths, grid, load_config, strip_fasta
from validation import validate_hll_csv


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXPERIMENT / "config.yaml")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("paths", type=Path, nargs="*")
    args = parser.parse_args()
    config, base = load_config(args.config)
    samples = [strip_fasta(str(path)) for path in fasta_paths(config, base)]
    checked = []
    paths = list(args.paths)
    if args.all:
        paths.extend((EXPERIMENT / "results" / "rehash").glob("p*/seed*/k*_w*_*.csv"))
    for path in sorted(set(paths)):
        completion_path = path.with_suffix(".complete.json")
        if not completion_path.is_file():
            raise SystemExit(f"completion manifest missing: {completion_path}")
        metadata = json.loads(completion_path.read_text())
        result = validate_hll_csv(path, samples, k=int(metadata["k"]), w=int(metadata["w"]),
                                  precision=int(metadata["precision"]))
        checked.append({"path": str(path), **result})
    print(json.dumps({"validated": len(checked), "outputs": checked}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
