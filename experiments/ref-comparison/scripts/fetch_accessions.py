#!/usr/bin/env python3
"""
scripts/fetch_accessions.py
===========================
Queries the ENCODE portal and NCBI SRA to confirm and resolve accession numbers
for the samples listed in config/config.yaml.

Outputs an updated accession table: config/confirmed_accessions.tsv

Run this BEFORE executing the Snakemake workflow to populate the sra_map
in config.yaml with confirmed SRA run IDs.

Usage:
    python scripts/fetch_accessions.py

Dependencies:
    pip install requests pandas pysradb
"""

import json
import sys
import time
import requests
import pandas as pd

ENCODE_API   = "https://www.encodeproject.org"
ENCODE_HEADS = {"Accept": "application/json"}

# =============================================================================
# Roadmap Epigenomics human samples
# GEO Superseries: GSE16256
# Roadmap Epigenomics Consortium (2015) Nature 518:317-330
# Epigenome IDs and corresponding GEO sample IDs from:
#   https://egg2.wustl.edu/roadmap/web_portal/meta.php
# =============================================================================

ROADMAP_SAMPLES = {
    # EID   : (tissue,          mark,     GEO_sample_accession)
    "E095"  : ("heart",         "H3K27ac", "GSM1127115"),
    "E066"  : ("liver",         "H3K27ac", "GSM1127068"),
    "E096"  : ("lung",          "H3K27ac", "GSM1127116"),
    "E067"  : ("brain_angular", "H3K27ac", "GSM1127069"),
    "E113"  : ("spleen",        "H3K27ac", "GSM1127132"),
    "E109"  : ("small_int",     "H3K27ac", "GSM1127128"),
    # H3K4me3 equivalents for multi-mark expansion
    "E095_K4me3": ("heart",    "H3K4me3", "GSM1127114"),
    "E066_K4me3": ("liver",    "H3K4me3", "GSM1127067"),
}

# =============================================================================
# Mouse ENCODE LICR ChIP-seq samples
# GEO Series: GSE49847
# Yue et al. (2014) Nature 515:355-364
# =============================================================================

MOUSE_ENCODE_SAMPLES = {
    # sample_id        : (tissue,   mark,     GEO_sample_accession — confirm via GEO)
    "mm_heart_H3K27ac" : ("heart",  "H3K27ac", None),  # TODO: confirm GSM
    "mm_liver_H3K27ac" : ("liver",  "H3K27ac", None),  # TODO: confirm GSM
    "mm_lung_H3K27ac"  : ("lung",   "H3K27ac", None),  # TODO: confirm GSM
    "mm_brain_H3K27ac" : ("brain",  "H3K27ac", None),  # TODO: confirm GSM
    "mm_spleen_H3K27ac": ("spleen", "H3K27ac", None),  # TODO: confirm GSM
}


def query_encode_experiment(accession: str) -> dict:
    """Retrieve ENCODE experiment metadata."""
    url = f"{ENCODE_API}/experiments/{accession}/?format=json"
    r = requests.get(url, headers=ENCODE_HEADS, timeout=30)
    r.raise_for_status()
    return r.json()


def geo_to_sra(gsm: str) -> list[str]:
    """Use pysradb to convert GEO sample (GSM) to SRA run IDs (SRR)."""
    try:
        from pysradb import SRAweb
        db = SRAweb()
        df = db.gsm_to_srr(gsm)
        return list(df["run_accession"].values)
    except Exception as e:
        print(f"  [WARN] pysradb failed for {gsm}: {e}", file=sys.stderr)
        return []


def main():
    rows = []

    print("=== Resolving Roadmap Epigenomics human samples ===")
    for eid, (tissue, mark, gsm) in ROADMAP_SAMPLES.items():
        sample_id = f"{eid}_{mark}"
        print(f"  {sample_id} | tissue={tissue} | GSM={gsm}")
        srr_list = []
        if gsm:
            srr_list = geo_to_sra(gsm)
            print(f"    SRR runs: {srr_list}")
        rows.append({
            "sample_id": sample_id,
            "species":   "human",
            "tissue":    tissue,
            "mark":      mark,
            "geo_series": "GSE16256",
            "geo_sample": gsm,
            "sra_runs":   ",".join(srr_list),
            "reference":  "GRCh38",
            "source":     "Roadmap2015",
        })
        time.sleep(0.5)

    print("\n=== Resolving Mouse ENCODE LICR samples ===")
    for sample_id, (tissue, mark, gsm) in MOUSE_ENCODE_SAMPLES.items():
        print(f"  {sample_id} | tissue={tissue} | GSM={gsm}")
        srr_list = []
        if gsm:
            srr_list = geo_to_sra(gsm)
            print(f"    SRR runs: {srr_list}")
        rows.append({
            "sample_id": sample_id,
            "species":   "mouse",
            "tissue":    tissue,
            "mark":      mark,
            "geo_series": "GSE49847",
            "geo_sample": gsm or "TODO",
            "sra_runs":   ",".join(srr_list),
            "reference":  "GRCm39",
            "source":     "MouseENCODE_Yue2014",
        })
        time.sleep(0.5)

    df = pd.DataFrame(rows)
    out = "config/confirmed_accessions.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nSaved confirmed accessions to {out}")
    print("\nTODO: rows with 'TODO' in geo_sample need manual confirmation.")
    print("      Check: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49847")


if __name__ == "__main__":
    main()

