#!/usr/bin/env python3
"""Generate the curated BED set for the subB / mixed-stride experiment.

18 BED files: 6 files per size class, 3 size classes (10k, 100k, 1M intervals).
Seeded so the corpus is reproducible across machines and reruns.
"""

import argparse
import os
import random
from pathlib import Path

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
SIZE_CLASSES = {"10k": 10_000, "100k": 100_000, "1M": 1_000_000}
FILES_PER_CLASS = 6
BASE_SEED = 20260511


def write_bed(path: Path, num_intervals: int, seed: int) -> None:
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for _ in range(num_intervals):
            chrom = rng.choice(CHROMS)
            start = rng.randint(0, 10_000_000)
            end = rng.randint(start + 100, start + 10_000)
            fh.write(f"{chrom}\t{start}\t{end}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default=str(Path(__file__).resolve().parent / "data"),
        help="Output directory for BED files (default: ./data)",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing files")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    seed_counter = 0
    for label, n in SIZE_CLASSES.items():
        for idx in range(FILES_PER_CLASS):
            path = data_dir / f"bed_{label}_{idx}.bed"
            seed = BASE_SEED + seed_counter
            seed_counter += 1
            if path.exists() and not args.force:
                print(f"skip   {path.name} (exists)")
                continue
            write_bed(path, n, seed)
            print(f"wrote  {path.name}  intervals={n:>8}  seed={seed}")

    print(f"\nData dir: {data_dir}")
    print(f"Total files: {len(SIZE_CLASSES) * FILES_PER_CLASS}")


if __name__ == "__main__":
    main()
