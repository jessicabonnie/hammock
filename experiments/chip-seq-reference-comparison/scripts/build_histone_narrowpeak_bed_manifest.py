#!/usr/bin/env python3
"""
Build ``histone_narrowpeak_bed_manifest.tsv`` from an ENCODE cart experiment report and a BED directory.

The report must contain a header row (line 2 after the timestamp line) with a ``Files`` column
listing ``/files/ENCFF.../`` paths. Each ``ENCFFxxxxxx.bed.gz`` under ``--bed-dir`` is joined to the
experiment row that references that accession (ENCODE cart narrowPeak files typically use the
``.bed.gz`` suffix). Includes ``lab``, ``project``, and library construction columns from the cart
for batch / provenance tracking.

Also writes ``histone_narrowpeak_bed_paths_<build>.txt`` under ``path_lists/`` (one absolute
``bed_path`` per line per resolved reference build: ``mm10``, ``GRCh38``, ``hg19``, ``mm9``,
``mm10_minimal``, …).

The manifest includes ``assembly`` (empty unless ``--fetch-assembly``) and ``reference_build``
(per-file when fetched, otherwise mm10/GRCh38 inferred from ``organism``). Use
``--fetch-assembly`` before ``histone_bedpaths_to_fasta.py --manifest`` when the cart mixes
hg19/GRCh38 or mm9/mm10/mm10-minimal.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from histone_encode_manifest_common import (
    batch_fields_from_experiment_row,
    fetch_encode_file_assembly,
    parse_files_cell,
    reference_build_for_manifest,
    write_bed_path_lists_by_reference_build,
)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Link local narrowPeak .bed.gz files to ENCODE experiment report rows."
    )
    p.add_argument(
        "--experiment-report",
        type=Path,
        required=True,
        help="ENCODE cart TSV (first line timestamp/URL, second line header)",
    )
    p.add_argument(
        "--bed-dir",
        type=Path,
        required=True,
        help="Directory with ENCFFxxxxxx.bed.gz narrowPeak files",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV (default: <bed-dir>/../histone_narrowpeak_bed_manifest.tsv)",
    )
    p.add_argument(
        "--fetch-assembly",
        action="store_true",
        help="Query encodeproject.org per file for assembly (mm10, hg19, GRCh38, …); slow but "
        "needed for correct bedtools getfasta references on mixed carts",
    )
    p.add_argument(
        "--path-lists-dir",
        type=Path,
        default=None,
        help="Directory for bed path lists (default: <output-dir>/path_lists)",
    )
    args = p.parse_args(argv)

    report = args.experiment_report.resolve()
    bed_dir = args.bed_dir.resolve()
    if not report.is_file():
        print("error: report not found: {}".format(report), file=sys.stderr)
        return 2
    if not bed_dir.is_dir():
        print("error: bed dir not found: {}".format(bed_dir), file=sys.stderr)
        return 2

    out = args.output
    if out is None:
        out = (bed_dir.parent / "histone_narrowpeak_bed_manifest.tsv").resolve()
    else:
        out = out.resolve()

    encff_to_row = {}  # type: Dict[str, Dict[str, str]]
    with report.open(newline="", encoding="utf-8") as f:
        f.readline()
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for encff in parse_files_cell(row.get("Files", "")):
                encff_to_row[encff] = row

    cols = [
        "bed_path",
        "file_accession",
        "assembly",
        "reference_build",
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
    suffix = ".bed.gz"
    bed_paths = sorted(bed_dir.glob("*.bed.gz"))
    for i, pth in enumerate(bed_paths):
        name = pth.name
        if not name.endswith(suffix):
            continue
        encff = name[: -len(suffix)]
        meta = encff_to_row.get(encff)
        if not meta:
            missing.append(encff)
            continue
        batch = batch_fields_from_experiment_row(meta)
        assembly = ""
        if args.fetch_assembly:
            try:
                assembly = fetch_encode_file_assembly(encff)
            except Exception as e:  # noqa: BLE001
                print(
                    "warning: could not fetch assembly for {}: {}".format(encff, e),
                    file=sys.stderr,
                )
            if (i + 1) % 25 == 0 or (i + 1) == len(bed_paths):
                print("  ... assembly {}/{}".format(i + 1, len(bed_paths)), file=sys.stderr)
        organism = (meta.get("Organism") or "").strip()
        ref_build = reference_build_for_manifest(assembly, organism)
        rows_out.append(
            {
                "bed_path": str(pth.resolve()),
                "file_accession": encff,
                "assembly": assembly,
                "reference_build": ref_build,
                "experiment_accession": (meta.get("Accession") or "").strip(),
                **batch,
                "biosample_accession": (meta.get("Biosample accession") or "").strip(),
                "life_stage": (meta.get("Life stage") or "").strip(),
                "biosample_age": (meta.get("Biosample age") or "").strip(),
                "tissue_type": (meta.get("Biosample term name") or "").strip(),
                "histone_target": (meta.get("Target of assay") or "").strip(),
                "organism": organism,
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

    print("Wrote {} ({} rows)".format(out, len(rows_out)))

    path_lists_dir = (args.path_lists_dir or (out.parent / "path_lists")).resolve()
    list_paths, n_unknown = write_bed_path_lists_by_reference_build(
        rows_out,
        output_dir=path_lists_dir,
        basename_prefix="histone_narrowpeak_bed_paths",
    )
    for lp in sorted(list_paths):
        print("Wrote {}".format(lp))
    if n_unknown:
        print(
            "warning: {} manifest rows had no resolvable reference build (omitted from per-build .txt lists)".format(
                n_unknown
            ),
            file=sys.stderr,
        )

    if missing:
        print(
            "warning: {} .bed.gz files not found in report Files column".format(len(missing)),
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
