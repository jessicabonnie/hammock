#!/usr/bin/env python3
"""
Build ``histone_fastq_manifest.tsv`` from an ENCODE cart experiment report and a FASTQ directory.

The report must contain a header row (line 2 after the timestamp line) with a ``Files`` column
listing ``/files/ENCFF.../`` paths. Each ``ENCFFxxxxxx.fastq.gz`` under ``--fastq-dir`` is joined
to the experiment row that references that accession. Includes ``lab``, ``project``, and library
construction columns from the cart for batch / provenance tracking.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from histone_encode_manifest_common import batch_fields_from_experiment_row, parse_files_cell


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Link local FASTQs to ENCODE experiment report rows.")
    p.add_argument(
        "--experiment-report",
        type=Path,
        required=True,
        help="ENCODE cart TSV (first line timestamp/URL, second line header)",
    )
    p.add_argument("--fastq-dir", type=Path, required=True, help="Directory with ENCFFxxxxxx.fastq.gz")
    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV (default: <fastq-dir>/../histone_fastq_manifest.tsv)",
    )
    args = p.parse_args(argv)

    report = args.experiment_report.resolve()
    fastq_dir = args.fastq_dir.resolve()
    if not report.is_file():
        print(f"error: report not found: {report}", file=sys.stderr)
        return 2
    if not fastq_dir.is_dir():
        print(f"error: fastq dir not found: {fastq_dir}", file=sys.stderr)
        return 2

    out = args.output
    if out is None:
        out = (fastq_dir.parent / "histone_fastq_manifest.tsv").resolve()
    else:
        out = out.resolve()

    encff_to_row = {}  # type: Dict[str, Dict[str, str]]
    with report.open(newline="", encoding="utf-8") as f:
        f.readline()  # skip timestamp / URL
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for encff in parse_files_cell(row.get("Files", "")):
                encff_to_row[encff] = row

    cols = [
        "fastq_path",
        "file_accession",
        "experiment_accession",
        "lab",
        "project",
        "library_construction_platform",
        "library_construction_method",
        "biosample_accession",
        "life_stage",
        "biosample_age",
        "tissue_type",
        "histone_target",
        "organism",
        "biological_replicate",
        "technical_replicate",
        "replicate_uris",
        "linked_antibody",
    ]

    rows_out = []  # type: List[Dict[str, str]]
    missing = []  # type: List[str]
    for pth in sorted(fastq_dir.glob("*.fastq.gz")):
        name = pth.name
        if not name.endswith(".fastq.gz"):
            continue
        encff = name[: -len(".fastq.gz")]
        meta = encff_to_row.get(encff)
        if not meta:
            missing.append(encff)
            continue
        batch = batch_fields_from_experiment_row(meta)
        rows_out.append(
            {
                "fastq_path": str(pth.resolve()),
                "file_accession": encff,
                "experiment_accession": (meta.get("Accession") or "").strip(),
                **batch,
                "biosample_accession": (meta.get("Biosample accession") or "").strip(),
                "life_stage": (meta.get("Life stage") or "").strip(),
                "biosample_age": (meta.get("Biosample age") or "").strip(),
                "tissue_type": (meta.get("Biosample term name") or "").strip(),
                "histone_target": (meta.get("Target of assay") or "").strip(),
                "organism": (meta.get("Organism") or "").strip(),
                "biological_replicate": (meta.get("Biological replicate") or "").strip(),
                "technical_replicate": (meta.get("Technical replicate") or "").strip(),
                "replicate_uris": (meta.get("Replicates") or "").strip(),
                "linked_antibody": (meta.get("Linked antibody") or "").strip(),
            }
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows_out)

    print(f"Wrote {out} ({len(rows_out)} rows)")
    if missing:
        print(f"warning: {len(missing)} FASTQs not found in report Files column", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
