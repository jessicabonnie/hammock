from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "source_data" / "sra_annual_experiment_counts.tsv"
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def fetch_count(term: str, email: str | None, api_key: str | None) -> tuple[int, str]:
    params = {
        "db": "sra",
        "term": term,
        "rettype": "count",
        "retmode": "json",
        "tool": "hammock_database_counts",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    url = f"{ESEARCH_URL}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "hammock-database-counts/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.load(response)
    return int(payload["esearchresult"]["count"]), url


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Query NCBI Entrez for cumulative SRA experiment counts at the end "
            "of each calendar year. Entrez SRA search results are experiment-level records."
        )
    )
    parser.add_argument("--start-year", type=int, default=2009)
    parser.add_argument("--end-year", type=int, default=2025)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--email", help="Email passed to NCBI E-utilities.")
    parser.add_argument("--api-key", help="Optional NCBI API key.")
    args = parser.parse_args()

    if args.start_year > args.end_year:
        raise ValueError("start-year must not exceed end-year")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for year in range(args.start_year, args.end_year + 1):
        # Publication date is used because it reflects when an experiment became
        # publicly searchable. The cumulative query includes all records through
        # December 31 of the specified year.
        term = f"1900/01/01:{year}/12/31[PDAT]"
        count, query_url = fetch_count(term, args.email, args.api_key)
        rows.append(
            {
                "year": year,
                "cumulative_experiments": count,
                "query_term": term,
                "retrieval_date": date.today().isoformat(),
                "query_url": query_url,
            }
        )
        time.sleep(0.34 if args.api_key else 0.67)

    counts = [int(row["cumulative_experiments"]) for row in rows]
    if any(current < previous for previous, current in zip(counts, counts[1:])):
        raise ValueError("Cumulative SRA experiment count decreased across years.")

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "year",
                "cumulative_experiments",
                "query_term",
                "retrieval_date",
                "query_url",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
