#!/usr/bin/env python3
"""
Build broadPeak BED path lists grouped by reference build.

Reads an ENCODE file manifest TSV (typically ``data/histone-marks/bed/manifest.tsv``),
selects broadPeak BED rows, maps assembly/build, and writes:

  path_lists/histone_broadpeak_bed_paths_<build>.txt

Each list contains one absolute BED path per line for files present under ``--bed-dir``.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from histone_encode_manifest_common import (  # noqa: E402
    REFERENCE_BUILD_TO_ENCODE_ASSEMBLY,
    encode_assembly_to_bed_list_slug,
    normalize_organism_label,
    write_bed_path_lists_by_reference_build,
)


def _effective_build(assembly: str, organism: str) -> str:
    a = (assembly or "").strip()
    if a:
        return a
    org = normalize_organism_label(organism)
    if org == "Mus musculus":
        return "mm10"
    if org == "Homo sapiens":
        return "GRCh38"
    return ""


def _is_broadpeak(row: Dict[str, str]) -> bool:
    ftype = (row.get("File format type") or "").strip().lower()
    file_type = (row.get("File type") or "").strip().lower()
    out_type = (row.get("Output type") or "").strip().lower()
    if "broadpeak" in ftype:
        return True
    if "broadpeak" in file_type:
        return True
    if "broad" in out_type and "peak" in out_type:
        return True
    return False


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Build histone broadPeak BED path lists by reference build."
    )
    p.add_argument(
        "--bed-manifest",
        type=Path,
        required=True,
        help="ENCODE file manifest TSV (e.g. data/histone-marks/bed/manifest.tsv)",
    )
    p.add_argument(
        "--bed-dir",
        type=Path,
        required=True,
        help="Directory containing local ENCFF*.bed.gz files",
    )
    p.add_argument(
        "--path-lists-dir",
        type=Path,
        required=True,
        help="Directory to write histone_broadpeak_bed_paths_<build>.txt",
    )
    args = p.parse_args(argv)

    bed_manifest = args.bed_manifest.resolve()
    bed_dir = args.bed_dir.resolve()
    path_lists_dir = args.path_lists_dir.resolve()
    if not bed_manifest.is_file():
        print("error: bed manifest not found: {}".format(bed_manifest), file=sys.stderr)
        return 2
    if not bed_dir.is_dir():
        print("error: bed dir not found: {}".format(bed_dir), file=sys.stderr)
        return 2

    rows_out = []  # type: List[Dict[str, str]]
    missing = 0
    with bed_manifest.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if not _is_broadpeak(row):
                continue
            acc = (row.get("File accession") or "").strip()
            if not acc:
                continue
            bed_path = bed_dir / "{}.bed.gz".format(acc)
            if not bed_path.is_file():
                missing += 1
                continue
            raw_assembly = (row.get("File assembly") or "").strip()
            assembly = REFERENCE_BUILD_TO_ENCODE_ASSEMBLY.get(raw_assembly, raw_assembly)
            organism = (row.get("Biosample organism") or "").strip()
            build = _effective_build(assembly, organism)
            rows_out.append(
                {
                    "bed_path": str(bed_path.resolve()),
                    "assembly": build,
                    "reference_build": build,
                    "organism": organism,
                }
            )

    list_paths, n_unknown = write_bed_path_lists_by_reference_build(
        rows_out,
        output_dir=path_lists_dir,
        basename_prefix="histone_broadpeak_bed_paths",
    )
    for lp in sorted(list_paths):
        print("Wrote {}".format(lp))
    print("BroadPeak rows with local BEDs: {}".format(len(rows_out)))
    if missing:
        print(
            "warning: {} broadPeak accessions missing local .bed.gz".format(missing),
            file=sys.stderr,
        )
    if n_unknown:
        print(
            "warning: {} rows had no resolvable build".format(n_unknown), file=sys.stderr
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
