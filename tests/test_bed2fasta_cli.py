"""End-to-end + validation tests for the bed2fasta CLI (--ref/--ref1/--ref2).

The end-to-end conversion test needs both `hammock` and `bedtools` on PATH; the
argparse-validation tests only need `hammock`.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

OURS = shutil.which("hammock")
HAS_BEDTOOLS = shutil.which("bedtools") is not None

pytestmark = pytest.mark.skipif(OURS is None, reason="hammock not on PATH")

_REF = ">chr1\n" + "ACGTACGTAC" * 20 + "\n>chr2\n" + "TTTTGGGGCC" * 20 + "\n"


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def _setup(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(_REF)
    b1 = tmp_path / "a.bed"
    b1.write_text("chr1\t0\t80\nchr2\t0\t80\n")
    b2 = tmp_path / "b.bed"
    b2.write_text("chr1\t10\t90\n")
    l1 = tmp_path / "list1.txt"
    l1.write_text(f"{b1}\n")
    l2 = tmp_path / "list2.txt"
    l2.write_text(f"{b2}\n")
    return ref, l1, l2


@pytest.mark.skipif(not HAS_BEDTOOLS, reason="bedtools not on PATH")
def test_bed2fasta_end_to_end(tmp_path: Path) -> None:
    ref, l1, l2 = _setup(tmp_path)
    r = _run([OURS, str(l1), str(l2), "--ref", str(ref),
              "-p", "14", "-k", "8", "-w", "40",
              "-o", str(tmp_path / "out")], tmp_path)
    assert r.returncode == 0, r.stderr
    csv = next(tmp_path.glob("out*.csv")).read_text().splitlines()
    header = csv[0].split(",")
    row = csv[1].split(",")
    # Row is labelled by the original BED basenames, not the temp FASTA names.
    assert row[0] == "a.bed"
    assert row[1] == "b.bed"
    # ref1/ref2 provenance columns carry the reference spec.
    assert row[header.index("ref1")] == str(ref)
    assert row[header.index("ref2")] == str(ref)
    # Mode was forced to D.
    assert row[header.index("mode")] == "D"


@pytest.mark.skipif(not HAS_BEDTOOLS, reason="bedtools not on PATH")
def test_bed2fasta_self_comparison_is_one(tmp_path: Path) -> None:
    ref, l1, _ = _setup(tmp_path)
    r = _run([OURS, str(l1), str(l1), "--ref", str(ref),
              "-p", "14", "-o", str(tmp_path / "self")], tmp_path)
    assert r.returncode == 0, r.stderr
    csv = next(tmp_path.glob("self*.csv")).read_text().splitlines()
    header = csv[0].split(",")
    row = csv[1].split(",")
    assert float(row[header.index("jaccard_similarity")]) == pytest.approx(1.0)


def test_ref_mutually_exclusive_with_ref1(tmp_path: Path) -> None:
    _, l1, l2 = _setup(tmp_path)
    r = _run([OURS, str(l1), str(l2), "--ref", "hg38", "--ref1", "mm10"], tmp_path)
    assert r.returncode != 0
    assert "mutually exclusive" in r.stderr


def test_ref_plus_explicit_non_d_mode_errors(tmp_path: Path) -> None:
    _, l1, l2 = _setup(tmp_path)
    r = _run([OURS, str(l1), str(l2), "--ref", "hg38", "--mode", "interval"], tmp_path)
    assert r.returncode != 0
    assert "sequence mode" in r.stderr


def test_only_one_of_ref1_ref2_errors(tmp_path: Path) -> None:
    _, l1, l2 = _setup(tmp_path)
    r = _run([OURS, str(l1), str(l2), "--ref1", "hg38"], tmp_path)
    assert r.returncode != 0
    assert "both" in r.stderr.lower()


def test_missing_keyword_gives_fetch_hint(tmp_path: Path) -> None:
    _, l1, l2 = _setup(tmp_path)
    cache = tmp_path / "emptycache"
    r = _run([OURS, str(l1), str(l2), "--ref", "mm10",
              "--ref-cache-dir", str(cache)], tmp_path)
    assert r.returncode != 0
    assert "fetch-ref mm10" in r.stderr


def test_cache_dir_without_ref_flag_errors(tmp_path: Path) -> None:
    _, l1, l2 = _setup(tmp_path)
    r = _run([OURS, str(l1), str(l2), "--ref-cache-dir", str(tmp_path)], tmp_path)
    assert r.returncode != 0
    assert "only apply with a reference" in r.stderr
