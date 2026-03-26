"""
Build nf-core/chipseq samplesheets from local ENCODE FASTQs.

Resolves input controls via file ``controlled_by`` when set, otherwise via the
ChIP experiment's ``possible_controls`` (preferring input-library controls) and
matching biological replicate in the control experiment's FASTQs.

See: https://nf-co.re/chipseq/usage#samplesheet-input
"""

import csv
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from histone_encode_manifest_common import (  # noqa: E402
    ORGANISM_HUMAN,
    ORGANISM_MOUSE,
    ORGANISM_SLUG,
    normalize_organism_label,
)

ENCODE_HOST = "https://www.encodeproject.org"
DEFAULT_USER_AGENT = "hammock-chipseq-samplesheet/1.0"


def encode_file_download_url(accession: str) -> str:
    return f"{ENCODE_HOST}/files/{accession}/@@download/{accession}.fastq.gz"


def _request_json(url: str, timeout: int = 120) -> Any:
    req = urllib.request.Request(url, headers={"User-Agent": DEFAULT_USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def experiment_accession_from_dataset(dataset: Optional[str]) -> str:
    if not dataset:
        return ""
    m = re.search(r"/experiments/(ENCSR[A-Z0-9]+)/", dataset)
    return m.group(1) if m else ""


def parse_controlled_by(controlled_by: Optional[List[str]]) -> List[str]:
    if not controlled_by:
        return []
    return re.findall(r"/files/(ENCFF\w+)/", ",".join(controlled_by))


def antibody_label(file_doc: Dict[str, Any]) -> str:
    t = file_doc.get("target")
    if isinstance(t, dict):
        lab = t.get("label") or t.get("name") or "UNKNOWN"
        return lab.replace(" ", "_").replace(",", "")
    return "UNKNOWN"


def biological_replicate_int(file_doc: Dict[str, Any]) -> int:
    reps = file_doc.get("biological_replicates") or [1]
    try:
        return int(reps[0])
    except (TypeError, ValueError, IndexError):
        return 1


def organism_from_experiment_biosample_summary(summary: Optional[str]) -> str:
    if not summary:
        return ""
    if "Mus musculus" in summary:
        return ORGANISM_MOUSE
    if "Homo sapiens" in summary:
        return ORGANISM_HUMAN
    return ""


class EncodeCaches(object):
    """In-memory cache for ENCODE experiment search results and file JSON."""

    def __init__(self):
        self.experiments = {}  # type: Dict[str, Dict[str, Any]]
        self.control_fastqs = {}  # type: Dict[str, List[Dict[str, Any]]]
        self.files = {}  # type: Dict[str, Dict[str, Any]]


def resolve_ip_organism(
    ip_acc: str,
    ip_doc: Dict[str, Any],
    *,
    organism_from_manifest: Optional[Dict[str, str]],
    sleep_s: float,
    caches: EncodeCaches,
    fetch_json: Callable[[str], Any],
) -> str:
    if organism_from_manifest:
        lab = normalize_organism_label(organism_from_manifest.get(ip_acc, ""))
        if lab:
            return lab
    chip_exp = experiment_accession_from_dataset(ip_doc.get("dataset"))
    if not chip_exp:
        return ""
    exp_doc = fetch_experiment_metadata(
        chip_exp, sleep_s=sleep_s, caches=caches, fetch_json=fetch_json
    )
    summary = exp_doc.get("biosample_summary") or ""
    return organism_from_experiment_biosample_summary(summary)


def load_metadata_jsonl(path: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    if not path.is_file():
        return out
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            acc = obj.get("accession")
            data = obj.get("data")
            if isinstance(acc, str) and isinstance(data, dict):
                out[acc] = data
    return out


def save_metadata_jsonl(path: Path, by_acc: Dict[str, Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for acc in sorted(by_acc.keys()):
            f.write(json.dumps({"accession": acc, "data": by_acc[acc]}, default=str) + "\n")


def fetch_file_metadata(
    accession: str,
    *,
    sleep_s: float = 0.0,
    fetch_json: Callable[[str], Any] = _request_json,
) -> Dict[str, Any]:
    if sleep_s > 0:
        time.sleep(sleep_s)
    return fetch_json(f"{ENCODE_HOST}/files/{accession}/?format=json")


def fetch_experiment_metadata(
    experiment_accession: str,
    *,
    sleep_s: float = 0.0,
    caches: EncodeCaches,
    fetch_json: Callable[[str], Any] = _request_json,
) -> Dict[str, Any]:
    if experiment_accession in caches.experiments:
        return caches.experiments[experiment_accession]
    if sleep_s > 0:
        time.sleep(sleep_s)
    doc = fetch_json(f"{ENCODE_HOST}/experiments/{experiment_accession}/?format=json")
    caches.experiments[experiment_accession] = doc
    return doc


def search_control_fastqs(
    control_experiment: str,
    *,
    sleep_s: float = 0.0,
    caches: EncodeCaches,
    fetch_json: Callable[[str], Any] = _request_json,
) -> List[Dict[str, Any]]:
    if control_experiment in caches.control_fastqs:
        return caches.control_fastqs[control_experiment]
    if sleep_s > 0:
        time.sleep(sleep_s)
    url = (
        f"{ENCODE_HOST}/search/?type=File&status=released"
        f"&dataset=/experiments/{control_experiment}/&file_type=fastq&format=json&limit=all"
    )
    data = fetch_json(url)
    rows = list(data.get("@graph", []))
    caches.control_fastqs[control_experiment] = rows
    return rows


def pick_control_fastq_accession(
    control_experiment: str,
    bio_rep: int,
    *,
    sleep_s: float = 0.0,
    caches: EncodeCaches,
    fetch_json: Callable[[str], Any] = _request_json,
) -> Optional[str]:
    files = search_control_fastqs(
        control_experiment, sleep_s=sleep_s, caches=caches, fetch_json=fetch_json
    )
    cand = [f for f in files if biological_replicate_int(f) == bio_rep]
    if not cand:
        cand = files
    if not cand:
        return None

    def sort_key(f: Dict[str, Any]) -> Tuple[str, str]:
        tr = (f.get("technical_replicates") or ["99_99"])[0]
        acc = f.get("accession") or ""
        return (str(tr), str(acc))

    cand.sort(key=sort_key)
    acc = cand[0].get("accession")
    return str(acc) if acc else None


def choose_input_control_experiment(chip_experiment_doc: Dict[str, Any]) -> Optional[str]:
    pcs = chip_experiment_doc.get("possible_controls") or []
    ctl: Optional[str] = None
    for pc in pcs:
        if pc.get("control_type") == "input library":
            ctl = pc.get("accession")
            break
        if "Control" in (pc.get("assay_title") or ""):
            ctl = pc.get("accession")
            break
    if not ctl and pcs:
        ctl = pcs[0].get("accession")
    return str(ctl) if ctl else None


def resolve_control_accession(
    ip_accession: str,
    ip_doc: Dict[str, Any],
    *,
    sleep_s: float,
    caches: EncodeCaches,
    fetch_json: Callable[[str], Any],
) -> Optional[str]:
    cb = parse_controlled_by(ip_doc.get("controlled_by"))
    if cb:
        return cb[0]
    chip_exp = experiment_accession_from_dataset(ip_doc.get("dataset"))
    if not chip_exp:
        return None
    exp_doc = fetch_experiment_metadata(
        chip_exp, sleep_s=sleep_s, caches=caches, fetch_json=fetch_json
    )
    ctl_exp = choose_input_control_experiment(exp_doc)
    if not ctl_exp:
        return None
    bio_rep = biological_replicate_int(ip_doc)
    return pick_control_fastq_accession(
        ctl_exp, bio_rep, sleep_s=sleep_s, caches=caches, fetch_json=fetch_json
    )


def refresh_file_metadata(
    accessions: Iterable[str],
    *,
    cache_path: Path,
    sleep_s: float = 0.1,
    only_missing: bool = True,
    fetch_json: Callable[[str], Any] = _request_json,
) -> Dict[str, Dict[str, Any]]:
    by_acc = load_metadata_jsonl(cache_path) if only_missing else {}
    todo = [a for a in sorted(set(accessions)) if a not in by_acc]
    for i, acc in enumerate(todo):
        try:
            by_acc[acc] = fetch_file_metadata(acc, sleep_s=sleep_s, fetch_json=fetch_json)
        except urllib.error.HTTPError as e:
            by_acc[acc] = {"_error": f"HTTPError {e.code}: {e.reason}"}
        except Exception as e:  # noqa: BLE001 — surface network/parse issues in cache
            by_acc[acc] = {"_error": str(e)}
        if (i + 1) % 25 == 0:
            save_metadata_jsonl(cache_path, by_acc)
    save_metadata_jsonl(cache_path, by_acc)
    return by_acc


def read_manifest_fastqs(
    manifest_path: Path,
) -> Tuple[List[Tuple[str, Path]], Dict[str, str]]:
    """
    Parse histone_fastq_manifest-style TSV.

    Returns (pairs, organism_by_file_accession) where organism values are raw strings
    from the ``organism`` column (normalized later via ``normalize_organism_label``).
    """
    rows: List[Tuple[str, Path]] = []
    organism_by_acc: Dict[str, str] = {}
    with manifest_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = (row.get("file_accession") or "").strip()
            p = (row.get("fastq_path") or "").strip()
            if acc and p:
                rows.append((acc, Path(p)))
                org = (row.get("organism") or "").strip()
                if org:
                    organism_by_acc[acc] = org
    return rows, organism_by_acc


def discover_fastqs_from_dir(fastq_dir: Path) -> List[Tuple[str, Path]]:
    out: List[Tuple[str, Path]] = []
    for p in sorted(fastq_dir.glob("*.fastq.gz")):
        name = p.name
        if not name.endswith(".fastq.gz"):
            continue
        acc = name[: -len(".fastq.gz")]
        out.append((acc, p.resolve()))
    return out


def collect_local_ip_accessions(
    pairs: List[Tuple[str, Path]],
) -> List[str]:
    return [acc for acc, path in pairs if path.is_file()]


def build_samplesheet(
    ip_pairs: List[Tuple[str, Path]],
    *,
    encode_by_acc: Dict[str, Dict[str, Any]],
    control_dir: Path,
    sleep_s: float = 0.08,
    fetch_json: Callable[[str], Any] = _request_json,
    organism_from_manifest: Optional[Dict[str, str]] = None,
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], Dict[str, Dict[str, Any]]]:
    """
    Returns (control_rows, ip_rows, ctrl_info) where ctrl_info maps control ENCFF
    to {control_exp, bio_rep}.

    Each row includes ``encode_organism`` (``Mus musculus`` or ``Homo sapiens``) for splitting.
    Organism is taken from ``organism_from_manifest`` when provided, else from the IP file's
    experiment ``biosample_summary`` via ENCODE.
    """
    caches = EncodeCaches()

    ip_rows_meta: List[Dict[str, Any]] = []
    for acc, path in sorted(ip_pairs, key=lambda x: x[1].as_posix()):
        doc = encode_by_acc.get(acc) or {}
        if doc.get("_error"):
            raise ValueError(f"Missing or errored ENCODE metadata for {acc}: {doc.get('_error')}")
        chip_exp = experiment_accession_from_dataset(doc.get("dataset"))
        ab = antibody_label(doc)
        bio_rep = biological_replicate_int(doc)
        ctrl = resolve_control_accession(acc, doc, sleep_s=sleep_s, caches=caches, fetch_json=fetch_json)
        if not ctrl:
            raise ValueError(f"Could not resolve control FASTQ for IP {acc} (experiment {chip_exp})")
        org = resolve_ip_organism(
            acc,
            doc,
            organism_from_manifest=organism_from_manifest,
            sleep_s=sleep_s,
            caches=caches,
            fetch_json=fetch_json,
        )
        if not org:
            raise ValueError(
                f"Could not determine organism for IP {acc} (add 'organism' to manifest or ensure "
                f"experiment biosample_summary mentions Mus musculus / Homo sapiens)"
            )
        ip_rows_meta.append(
            {
                "ip_acc": acc,
                "ip_path": str(path.resolve()),
                "chip_exp": chip_exp,
                "bio_rep": bio_rep,
                "antibody": ab,
                "control_file": ctrl,
                "organism": org,
            }
        )

    ctrl_info: Dict[str, Dict[str, Any]] = {}
    control_organism: Dict[str, str] = {}
    for r in ip_rows_meta:
        cf = r["control_file"]
        org_ip = str(r["organism"])
        if cf in control_organism and control_organism[cf] != org_ip:
            raise ValueError(
                f"Control {cf} is paired with IPs from different organisms "
                f"({control_organism[cf]} vs {org_ip})"
            )
        control_organism[cf] = org_ip
        if cf not in ctrl_info:
            if cf in encode_by_acc and not encode_by_acc[cf].get("_error"):
                fm = encode_by_acc[cf]
            elif cf in caches.files:
                fm = caches.files[cf]
            else:
                if sleep_s > 0:
                    time.sleep(sleep_s)
                fm = fetch_file_metadata(cf, sleep_s=0.0, fetch_json=fetch_json)
                caches.files[cf] = fm
                encode_by_acc[cf] = fm
            ctrl_info[cf] = {
                "control_exp": experiment_accession_from_dataset(fm.get("dataset")),
                "bio_rep": biological_replicate_int(fm),
            }

    control_rows: List[Dict[str, str]] = []
    for cf in sorted(ctrl_info.keys()):
        info = ctrl_info[cf]
        cexp = str(info["control_exp"])
        br = int(info["bio_rep"])
        sample = f"INPUT_{cexp}"
        row = {
            "sample": sample,
            "fastq_1": str((control_dir / f"{cf}.fastq.gz").resolve()),
            "fastq_2": "",
            "replicate": str(br),
            "antibody": "",
            "control": "",
            "control_replicate": "",
            "encode_ip_file_accession": "",
            "encode_control_file_accession": cf,
            "encode_ip_download_url": "",
            "encode_control_download_url": encode_file_download_url(cf),
            "chip_experiment_accession": "",
            "control_experiment_accession": cexp,
            "encode_organism": control_organism.get(cf, ""),
        }
        control_rows.append(row)

    ip_rows: List[Dict[str, str]] = []
    for r in ip_rows_meta:
        cf = r["control_file"]
        cexp = str(ctrl_info[cf]["control_exp"])
        cbr = int(ctrl_info[cf]["bio_rep"])
        sample = f"{r['chip_exp']}_{r['antibody']}"
        acc = r["ip_acc"]
        org = str(r["organism"])
        ip_rows.append(
            {
                "sample": sample,
                "fastq_1": r["ip_path"],
                "fastq_2": "",
                "replicate": str(int(r["bio_rep"])),
                "antibody": r["antibody"],
                "control": f"INPUT_{cexp}",
                "control_replicate": str(cbr),
                "encode_ip_file_accession": acc,
                "encode_control_file_accession": cf,
                "encode_ip_download_url": encode_file_download_url(acc),
                "encode_control_download_url": encode_file_download_url(cf),
                "chip_experiment_accession": r["chip_exp"],
                "control_experiment_accession": cexp,
                "encode_organism": org,
            }
        )

    return control_rows, ip_rows, ctrl_info


def split_samplesheet_by_organism(
    control_rows: List[Dict[str, str]],
    ip_rows: List[Dict[str, str]],
) -> Dict[str, Tuple[List[Dict[str, str]], List[Dict[str, str]]]]:
    """
    Partition INPUT + IP rows by ``encode_organism``.

    Returns slug -> (control_rows_subset, ip_rows_subset). Slugs are ``Mus_musculus`` and
    ``Homo_sapiens``. Omits species with no IP rows.
    """
    out: Dict[str, Tuple[List[Dict[str, str]], List[Dict[str, str]]]] = {}
    for canon, slug in ORGANISM_SLUG.items():
        ip_sub = [r for r in ip_rows if r.get("encode_organism") == canon]
        if not ip_sub:
            continue
        ctrl_needed = {
            r.get("encode_control_file_accession", "")
            for r in ip_sub
            if r.get("encode_control_file_accession")
        }
        ctrl_sub = [
            r
            for r in control_rows
            if r.get("encode_control_file_accession") in ctrl_needed
        ]
        out[slug] = (ctrl_sub, ip_sub)
    return out


def write_samplesheet_csv(path: Path, control_rows: List[Dict[str, str]], ip_rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not control_rows and not ip_rows:
        raise ValueError("write_samplesheet_csv: no rows")
    fieldnames = list(control_rows[0].keys()) if control_rows else list(ip_rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in sorted(control_rows, key=lambda r: (r["sample"], int(r["replicate"]), r["fastq_1"])):
            w.writerow(row)
        for row in sorted(ip_rows, key=lambda r: (r["sample"], int(r["replicate"]), r["fastq_1"])):
            w.writerow(row)


def write_control_manifest_tsv(
    path: Path,
    control_rows: List[Dict[str, str]],
    ip_rows: Optional[List[Dict[str, str]]] = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    target_by_control: Dict[str, List[str]] = {}
    if ip_rows:
        tmp: Dict[str, set] = {}
        for row in ip_rows:
            cf = (row.get("encode_control_file_accession") or "").strip()
            tg = (row.get("antibody") or "").strip()
            if not cf or not tg:
                continue
            tmp.setdefault(cf, set()).add(tg)
        target_by_control = {k: sorted(v) for k, v in tmp.items()}
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "file_accession",
                "control_experiment_accession",
                "biological_replicate",
                "histone_target",
                "encode_download_url",
                "local_path_after_download",
            ]
        )
        for row in sorted(control_rows, key=lambda r: (r["control_experiment_accession"], int(r["replicate"]), r["fastq_1"])):
            m = re.search(r"(ENCFF\w+)\.fastq\.gz$", row["fastq_1"] or "")
            if not m:
                continue
            cf = m.group(1)
            w.writerow(
                [
                    cf,
                    row["control_experiment_accession"],
                    row["replicate"],
                    ",".join(target_by_control.get(cf, [])),
                    row["encode_control_download_url"],
                    row["fastq_1"],
                ]
            )


def write_download_script(path: Path, control_rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
    ]
    dirs = sorted({str(Path(r["fastq_1"]).parent) for r in control_rows})
    for d in dirs:
        lines.append(f'mkdir -p "{d}"')
    for row in sorted(control_rows, key=lambda r: r["fastq_1"]):
        url = row["encode_control_download_url"]
        dest = row["fastq_1"]
        lines.append(f'curl -fL --retry 3 -o "{dest}" "{url}"')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | 0o111)
