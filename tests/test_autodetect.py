"""Unit tests for CLI input-type auto-detection (cli._classify_file / _autodetect_mode).

These run pure-Python and don't need the C++ extension or the installed binary,
so they're safe to run while a rebuild is pending.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from hammock import cli

DATA = Path(__file__).parent / "data"


def _ns(mode=None, list_path=None, **overrides) -> argparse.Namespace:
    base = dict(mode=mode, filepaths_file=str(list_path) if list_path else "",
                subA=1.0, subB=1.0, expA=0)
    base.update(overrides)
    return argparse.Namespace(**base)


def test_classify_fasta_by_extension(tmp_path: Path) -> None:
    assert cli._classify_file(str(DATA / "tiny.fa")) == "fasta"
    p = tmp_path / "x.fasta"
    p.write_text(">id\nACGT\n")
    assert cli._classify_file(str(p)) == "fasta"


def test_classify_bed_by_extension() -> None:
    assert cli._classify_file(str(DATA / "tiny_a.bed")) == "bed"


def test_classify_unknown_ext_uses_content_sniff(tmp_path: Path) -> None:
    fa = tmp_path / "weird.txt"
    fa.write_text(">id\nACGT\n")
    assert cli._classify_file(str(fa)) == "fasta"

    bed = tmp_path / "regions.tsv"
    bed.write_text("chr1\t100\t200\n")
    assert cli._classify_file(str(bed)) == "bed"


def test_classify_handles_gz(tmp_path: Path) -> None:
    import gzip
    p = tmp_path / "regions.bed.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("chr1\t100\t200\n")
    assert cli._classify_file(str(p)) == "bed"

    p2 = tmp_path / "seq.fa.gz"
    with gzip.open(p2, "wt") as fh:
        fh.write(">id\nACGT\n")
    assert cli._classify_file(str(p2)) == "fasta"


def test_autodetect_fasta_picks_d(tmp_path: Path) -> None:
    lst = tmp_path / "files.txt"
    lst.write_text(str(DATA / "tiny.fa") + "\n")
    assert cli._autodetect_mode(_ns(list_path=lst)) == "D"


def test_autodetect_bed_picks_a(tmp_path: Path) -> None:
    lst = tmp_path / "files.txt"
    lst.write_text(str(DATA / "tiny_a.bed") + "\n")
    assert cli._autodetect_mode(_ns(list_path=lst)) == "A"


def test_autodetect_bed_with_subsampling_picks_c(tmp_path: Path) -> None:
    lst = tmp_path / "files.txt"
    lst.write_text(str(DATA / "tiny_a.bed") + "\n")
    assert cli._autodetect_mode(_ns(list_path=lst, subA=0.5)) == "C"


def test_explicit_mode_wins(tmp_path: Path) -> None:
    lst = tmp_path / "files.txt"
    lst.write_text(str(DATA / "tiny.fa") + "\n")
    # FASTA input but --mode B explicitly: honor user's choice.
    assert cli._autodetect_mode(_ns(mode="B", list_path=lst)) == "B"


def test_autodetect_skips_comments_and_blanks(tmp_path: Path) -> None:
    lst = tmp_path / "files.txt"
    lst.write_text("# header\n\n   \n" + str(DATA / "tiny.fa") + "\n")
    assert cli._autodetect_mode(_ns(list_path=lst)) == "D"


def test_autodetect_missing_list_falls_back_to_a(tmp_path: Path) -> None:
    assert cli._autodetect_mode(_ns(list_path=tmp_path / "nope.txt")) == "A"
