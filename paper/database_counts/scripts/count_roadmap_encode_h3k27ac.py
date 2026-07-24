#!/usr/bin/env python3
"""Count Roadmap hg19 and ENCODE GRCh38 H3K27ac peak BED files.

The output separates directly comparable within-resource pairs from cross-resource
pairs that cannot be compared with ordinary coordinate overlap because the two
collections use different human reference assemblies.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import ssl
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Iterable

ROADMAP_NARROWPEAK_URL = (
    "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"
)
ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/"
USER_AGENT = "hammock-database-counts/1.1"


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_ssl_context(ca_bundle: str | None, insecure: bool) -> ssl.SSLContext:
    if insecure:
        return ssl._create_unverified_context()
    candidates = [
        ca_bundle,
        os.environ.get("SSL_CERT_FILE"),
        os.environ.get("REQUESTS_CA_BUNDLE"),
        os.environ.get("CURL_CA_BUNDLE"),
    ]
    try:
        import certifi  # type: ignore
        candidates.append(certifi.where())
    except ImportError:
        pass
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return ssl.create_default_context(cafile=candidate)
    return ssl.create_default_context()


def open_url(request: urllib.request.Request, ssl_context: ssl.SSLContext):
    try:
        return urllib.request.urlopen(request, timeout=180, context=ssl_context)
    except urllib.error.URLError as exc:
        if isinstance(exc.reason, ssl.SSLCertVerificationError):
            raise RuntimeError(
                "TLS certificate verification failed. Pass --ca-bundle /path/to/ca.pem, "
                "set SSL_CERT_FILE, or use --insecure only as a last resort."
            ) from exc
        raise


def fetch_text(url: str, ssl_context: ssl.SSLContext) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with open_url(request, ssl_context) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_json(url: str, ssl_context: ssl.SSLContext) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    with open_url(request, ssl_context) as response:
        return json.load(response)


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)


def count_roadmap(url: str, target: str, ssl_context: ssl.SSLContext) -> list[dict[str, Any]]:
    parser = LinkParser()
    parser.feed(fetch_text(url, ssl_context))
    pattern = re.compile(rf"^(E\d+)-{re.escape(target)}\.narrowPeak\.gz$", re.IGNORECASE)
    rows = []
    for href in sorted(set(parser.links)):
        name = urllib.parse.unquote(Path(href).name)
        match = pattern.match(name)
        if match:
            rows.append({
                "repository": "Roadmap Epigenomics",
                "dataset_id": match.group(1),
                "file_accession": "",
                "target": target,
                "assembly": "hg19",
                "file_format": "narrowPeak",
                "output_type": "consolidated narrow peaks",
                "output_category": "annotation",
                "download_url": urllib.parse.urljoin(url, href),
            })
    return rows


def encode_search_url(target: str, assembly: str) -> str:
    params = [
        ("type", "File"),
        ("status", "released"),
        ("file_format", "bed"),
        ("assembly", assembly),
        ("assay_title", "Histone ChIP-seq"),
        ("target.label", target),
        ("limit", "all"),
        ("format", "json"),
        ("frame", "embedded"),
    ]
    return ENCODE_SEARCH_URL + "?" + urllib.parse.urlencode(params)


def extract_accession(value: Any, collection: str) -> str:
    """Extract an accession from an embedded object or /collection/ACCESSION/ path."""
    if isinstance(value, dict):
        accession = value.get("accession")
        if isinstance(accession, str) and accession:
            return accession
        value = value.get("@id", "")
    if isinstance(value, str):
        match = re.search(rf"/{re.escape(collection)}/([^/]+)/?", value)
        if match:
            return match.group(1)
    return ""


def count_encode(
    target: str,
    assembly: str,
    ssl_context: ssl.SSLContext,
) -> tuple[list[dict[str, Any]], dict[str, int], list[dict[str, Any]]]:
    payload = fetch_json(encode_search_url(target, assembly), ssl_context)
    candidates: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []

    for record in payload.get("@graph", []):
        output_type = str(record.get("output_type", ""))
        file_type = str(record.get("file_type", ""))
        output_category = str(record.get("output_category", ""))
        dataset_id = extract_accession(record.get("dataset"), "experiments")
        accession = str(record.get("accession", ""))
        href = str(record.get("href", ""))

        reason = ""
        if "peak" not in output_type.casefold():
            reason = "output_type_not_peak"
        elif "bed" not in file_type.casefold() and record.get("file_format") != "bed":
            reason = "not_bed"
        elif not dataset_id:
            reason = "missing_experiment_accession"
        elif output_category and output_category.casefold() != "annotation":
            reason = f"unexpected_output_category:{output_category}"
        elif not accession:
            reason = "missing_file_accession"
        elif not href:
            reason = "missing_download_href"

        row = {
            "repository": "ENCODE",
            "dataset_id": dataset_id,
            "file_accession": accession,
            "target": target,
            "assembly": assembly,
            "file_format": file_type or "bed",
            "output_type": output_type,
            "output_category": output_category,
            "download_url": urllib.parse.urljoin("https://www.encodeproject.org", href),
            "preferred_default": bool(record.get("preferred_default", False)),
            "date_created": str(record.get("date_created", "")),
            "rejection_reason": reason,
        }
        (rejected if reason else candidates).append(row)

    rank_terms = (
        "optimal IDR thresholded peaks",
        "conservative IDR thresholded peaks",
        "IDR thresholded peaks",
        "replicated peaks",
        "pseudoreplicated peaks",
        "peaks",
    )
    rank = {name.casefold(): index for index, name in enumerate(rank_terms)}
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in candidates:
        grouped.setdefault(row["dataset_id"], []).append(row)

    selected = []
    for rows in grouped.values():
        rows.sort(key=lambda row: (
            not row["preferred_default"],
            rank.get(row["output_type"].casefold(), len(rank)),
            row["date_created"],
            row["file_accession"],
        ))
        selected.append(rows[0])

    diagnostics = {
        "api_file_records": len(payload.get("@graph", [])),
        "eligible_peak_bed_records": len(candidates),
        "rejected_records": len(rejected),
        "unique_experiments": len(grouped),
        "experiments_with_multiple_candidate_files": sum(
            1 for rows in grouped.values() if len(rows) > 1
        ),
        "duplicate_candidate_files_removed": len(candidates) - len(grouped),
        "selected_files": len(selected),
        "selected_missing_dataset_id": sum(1 for row in selected if not row["dataset_id"]),
    }
    return sorted(selected, key=lambda row: (row["dataset_id"], row["file_accession"])), diagnostics, rejected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="H3K27ac")
    parser.add_argument("--roadmap-url", default=ROADMAP_NARROWPEAK_URL)
    parser.add_argument("--encode-assembly", default="GRCh38")
    parser.add_argument("--ca-bundle")
    parser.add_argument("--insecure", action="store_true")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results",
    )
    args = parser.parse_args()

    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    ssl_context = build_ssl_context(args.ca_bundle, args.insecure)
    roadmap = count_roadmap(args.roadmap_url, args.target, ssl_context)
    encode, encode_diagnostics, encode_rejected = count_encode(
        args.target, args.encode_assembly, ssl_context
    )
    if not roadmap:
        raise RuntimeError("No Roadmap narrowPeak files matched the requested target")
    if not encode:
        raise RuntimeError("No ENCODE experiment-level peak BED files matched the requested target")

    stem = f"roadmap_encode_{args.target.lower()}"
    manifest_path = args.output_dir / f"{stem}_manifest.tsv"
    summary_path = args.output_dir / f"{stem}_summary.tsv"
    provenance_path = args.output_dir / f"{stem}_provenance.json"
    audit_path = args.output_dir / f"{stem}_encode_audit.tsv"

    manifest = roadmap + encode
    for row in manifest:
        row["retrieved_at_utc"] = retrieved_at
    write_tsv(
        manifest_path,
        manifest,
        [
            "repository", "dataset_id", "file_accession", "target", "assembly",
            "file_format", "output_type", "output_category", "download_url",
            "retrieved_at_utc",
        ],
    )

    for row in encode_rejected:
        row["retrieved_at_utc"] = retrieved_at
    write_tsv(
        audit_path,
        encode_rejected,
        [
            "dataset_id", "file_accession", "target", "assembly", "file_format",
            "output_type", "output_category", "preferred_default", "rejection_reason",
            "download_url", "retrieved_at_utc",
        ],
    )

    n_roadmap = len(roadmap)
    n_encode = len(encode)
    summary = [
        {
            "row_type": "repository",
            "repository": "Roadmap Epigenomics",
            "target": args.target,
            "assembly": "hg19",
            "bed_file_count": n_roadmap,
            "within_reference_pairwise_comparisons": pair_count(n_roadmap),
            "cross_repository_incompatible_pairs": "",
            "retrieved_at_utc": retrieved_at,
        },
        {
            "row_type": "repository",
            "repository": "ENCODE",
            "target": args.target,
            "assembly": args.encode_assembly,
            "bed_file_count": n_encode,
            "within_reference_pairwise_comparisons": pair_count(n_encode),
            "cross_repository_incompatible_pairs": "",
            "retrieved_at_utc": retrieved_at,
        },
        {
            "row_type": "cross_repository",
            "repository": "Roadmap Epigenomics × ENCODE",
            "target": args.target,
            "assembly": f"hg19 × {args.encode_assembly}",
            "bed_file_count": n_roadmap + n_encode,
            "within_reference_pairwise_comparisons": pair_count(n_roadmap) + pair_count(n_encode),
            "cross_repository_incompatible_pairs": n_roadmap * n_encode,
            "retrieved_at_utc": retrieved_at,
        },
    ]
    write_tsv(summary_path, summary, list(summary[-1].keys()))

    env_ca_bundle = next(
        (
            os.environ.get(name, "")
            for name in ("SSL_CERT_FILE", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE")
            if os.environ.get(name)
        ),
        "",
    )
    provenance = {
        "retrieved_at_utc": retrieved_at,
        "target": args.target,
        "roadmap": {
            "source_url": args.roadmap_url,
            "assembly": "hg19",
            "selection": "consolidated narrowPeak filenames matching E###-TARGET.narrowPeak.gz",
            "file_count": n_roadmap,
        },
        "encode": {
            "search_url": encode_search_url(args.target, args.encode_assembly),
            "assembly": args.encode_assembly,
            "selection": (
                "one released experiment-linked peak BED per experiment; prefer "
                "preferred_default, then IDR/replicated output rank"
            ),
            "file_count": n_encode,
            "diagnostics": encode_diagnostics,
        },
        "cross_repository_pair_policy": (
            "Roadmap_count multiplied by ENCODE unique-experiment count; pairs share target "
            "and species but use incompatible assemblies"
        ),
        "tls": {
            "certificate_verification": not args.insecure,
            "ca_bundle_argument": args.ca_bundle or "",
            "environment_ca_bundle": env_ca_bundle,
        },
        "outputs": [str(manifest_path), str(summary_path), str(audit_path)],
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"Roadmap hg19 {args.target} narrowPeak files: {n_roadmap:,}")
    print(f"ENCODE {args.encode_assembly} unique experiments/files selected: {n_encode:,}")
    print("ENCODE audit: " + ", ".join(f"{k}={v:,}" for k, v in encode_diagnostics.items()))
    print(f"Cross-resource incompatible pairs: {n_roadmap * n_encode:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
