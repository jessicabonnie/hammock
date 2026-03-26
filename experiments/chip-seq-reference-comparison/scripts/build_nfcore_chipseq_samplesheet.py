#!/usr/bin/env python3
"""
Create nf-core/chipseq samplesheet and intermediate ENCODE artifacts from local FASTQs.

Typical flow
------------
1. Place IP (histone ChIP) FASTQs under ``--fastq-dir`` as ``ENCFFxxxxxx.fastq.gz``.
2. Optionally pass ``--manifest`` (e.g. ``histone_fastq_manifest.tsv``) to restrict rows.
3. Run this script (network required unless metadata cache is complete).
4. Run the emitted ``download_control_fastqs_<species>.sh`` scripts to fetch input controls.
5. Run nf-core/chipseq separately per species with ``--input`` pointing at
   ``nfcore_chipseq_samplesheet_Mus_musculus.csv`` or ``nfcore_chipseq_samplesheet_Homo_sapiens.csv``.

Example::

    python3 build_nfcore_chipseq_samplesheet.py \\
        --fastq-dir /path/to/histone-marks/fastq \\
        --output-dir /path/to/histone-marks \\
        --manifest /path/to/histone_fastq_manifest.tsv
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from chipseq_nfcore_samplesheet import (
    build_samplesheet,
    collect_local_ip_accessions,
    discover_fastqs_from_dir,
    load_metadata_jsonl,
    read_manifest_fastqs,
    refresh_file_metadata,
    save_metadata_jsonl,
    split_samplesheet_by_organism,
    write_control_manifest_tsv,
    write_download_script,
    write_samplesheet_csv,
)


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build nf-core/chipseq samplesheet from ENCODE ChIP FASTQs on disk."
    )
    p.add_argument(
        "--fastq-dir",
        type=Path,
        required=True,
        help="Directory containing IP FASTQs named ENCFFxxxxxx.fastq.gz",
    )
    p.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Optional TSV with columns file_accession and fastq_path (e.g. histone_fastq_manifest.tsv)",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for cache, TSV manifest, shell downloader, and samplesheet CSV",
    )
    p.add_argument(
        "--control-dir",
        type=Path,
        default=None,
        help="Directory for input/control FASTQs (default: <output-dir>/fastq_controls)",
    )
    p.add_argument(
        "--encode-cache",
        type=Path,
        default=None,
        help="JSONL cache of ENCODE /files/{ENCFF}/ metadata (default: <output-dir>/encode_file_metadata_cache.jsonl)",
    )
    p.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Do not query ENCODE; require complete metadata in --encode-cache",
    )
    p.add_argument(
        "--refresh-metadata",
        action="store_true",
        help="Re-fetch metadata for all IP FASTQs even if present in cache",
    )
    p.add_argument(
        "--sleep",
        type=float,
        default=0.1,
        help="Seconds to sleep between ENCODE requests (default: 0.1)",
    )
    p.add_argument(
        "--write-combined",
        action="store_true",
        help="Also write a single nfcore_chipseq_samplesheet.csv (all species; default: species-split only)",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = _parse_args(argv)
    fastq_dir = args.fastq_dir.resolve()
    out_dir = args.output_dir.resolve()
    control_dir = (args.control_dir or (out_dir / "fastq_controls")).resolve()
    cache_path = (args.encode_cache or (out_dir / "encode_file_metadata_cache.jsonl")).resolve()

    if not fastq_dir.is_dir():
        print(f"error: --fastq-dir is not a directory: {fastq_dir}", file=sys.stderr)
        return 2

    organism_by_acc = None
    if args.manifest:
        manifest_path = args.manifest.resolve()
        if not manifest_path.is_file():
            print(f"error: --manifest not found: {manifest_path}", file=sys.stderr)
            return 2
        pairs, organism_by_acc = read_manifest_fastqs(manifest_path)
        missing = [str(p) for _, p in pairs if not p.is_file()]
        if missing:
            print("warning: manifest paths missing on disk:", file=sys.stderr)
            for m in missing[:10]:
                print(f"  {m}", file=sys.stderr)
            if len(missing) > 10:
                print(f"  ... and {len(missing) - 10} more", file=sys.stderr)
        pairs = [(a, p.resolve()) for a, p in pairs if p.is_file()]
    else:
        pairs = discover_fastqs_from_dir(fastq_dir)
        organism_by_acc = None

    accessions = collect_local_ip_accessions(pairs)
    if not accessions:
        print("error: no *.fastq.gz files found (or manifest had no existing paths)", file=sys.stderr)
        return 2

    if args.skip_fetch:
        by_acc = load_metadata_jsonl(cache_path)
        missing = [a for a in accessions if a not in by_acc or by_acc[a].get("_error")]
        if missing:
            print(
                "error: --skip-fetch but cache missing or errored for: " + ", ".join(missing[:20]),
                file=sys.stderr,
            )
            return 2
    else:
        by_acc = refresh_file_metadata(
            accessions,
            cache_path=cache_path,
            sleep_s=args.sleep,
            only_missing=not args.refresh_metadata,
        )
        bad = [a for a in accessions if by_acc.get(a, {}).get("_error")]
        if bad:
            print("error: ENCODE metadata fetch failed for:", file=sys.stderr)
            for a in bad[:15]:
                print(f"  {a}: {by_acc[a].get('_error')}", file=sys.stderr)
            return 2

    try:
        control_rows, ip_rows, _ctrl_info = build_samplesheet(
            pairs,
            encode_by_acc=by_acc,
            control_dir=control_dir,
            sleep_s=args.sleep,
            organism_from_manifest=organism_by_acc,
        )
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    by_species = split_samplesheet_by_organism(control_rows, ip_rows)
    if not by_species:
        print("error: no IP rows after organism split (unexpected)", file=sys.stderr)
        return 2

    save_metadata_jsonl(cache_path, by_acc)

    print("Wrote:")
    print(f"  {cache_path}")

    for slug, (ctrl_sub, ip_sub) in sorted(by_species.items()):
        suf = "_{}".format(slug)
        ss_path = out_dir / "nfcore_chipseq_samplesheet{}.csv".format(suf)
        tsv_path = out_dir / "nfcore_control_fastqs_encode{}.tsv".format(suf)
        sh_path = out_dir / "download_control_fastqs{}.sh".format(suf)
        write_samplesheet_csv(ss_path, ctrl_sub, ip_sub)
        write_control_manifest_tsv(tsv_path, ctrl_sub)
        write_download_script(sh_path, ctrl_sub)
        print(f"  {ss_path}  ({len(ip_sub)} IP rows, {len(ctrl_sub)} control FASTQs)")
        print(f"  {tsv_path}")
        print(f"  {sh_path}")

    if args.write_combined:
        samplesheet_path = out_dir / "nfcore_chipseq_samplesheet.csv"
        control_tsv = out_dir / "nfcore_control_fastqs_encode.tsv"
        dl_script = out_dir / "download_control_fastqs.sh"
        write_samplesheet_csv(samplesheet_path, control_rows, ip_rows)
        write_control_manifest_tsv(control_tsv, control_rows)
        write_download_script(dl_script, control_rows)
        print(f"  {samplesheet_path}  (combined, --write-combined)")
        print(f"  {control_tsv}")
        print(f"  {dl_script}")

    print(
        "Total IP FASTQ rows: {}; unique control FASTQs: {}".format(
            len(ip_rows), len(control_rows)
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
