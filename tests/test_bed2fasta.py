"""Tests for hammock.bed2fasta (bedtools getfasta extraction).

Extraction tests are skipped when bedtools isn't on PATH (`ml bedtools`); the
pure-Python helpers (_fasta_stats, _out_name) run unconditionally.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from hammock import bed2fasta as b2f

HAS_BEDTOOLS = shutil.which("bedtools") is not None
bedtools_only = pytest.mark.skipif(not HAS_BEDTOOLS, reason="bedtools not on PATH")

_REF = (
    ">chr1\n" + "ACGTACGTAC" * 5 + "\n"
    ">chr2\n" + "TTTTGGGGCC" * 5 + "\n"
)


def _write_ref(tmp_path: Path) -> str:
    ref = tmp_path / "ref.fa"
    ref.write_text(_REF)
    return str(ref)


# ---- pure-Python helpers (no bedtools) ------------------------------------

def test_out_name_is_unique_per_index() -> None:
    a = b2f._out_name(0, "/x/peaks.bed")
    b = b2f._out_name(1, "/y/peaks.bed")  # same basename, different dir
    assert a != b
    assert a == "0000_peaks.fa"
    assert b == "0001_peaks.fa"


def test_fasta_stats_counts_records_and_n(tmp_path: Path) -> None:
    fa = tmp_path / "s.fa"
    fa.write_text(">a\nACGTNNNN\n>b\nACGT\n")
    records, n_frac = b2f._fasta_stats(str(fa))
    assert records == 2
    assert n_frac == pytest.approx(4 / 12)


# ---- extraction (needs bedtools) ------------------------------------------

@bedtools_only
def test_bed_to_fasta_extracts_expected(tmp_path: Path) -> None:
    ref = _write_ref(tmp_path)
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t0\t8\nchr2\t0\t4\n")
    out = tmp_path / "out.fa"
    b2f.bed_to_fasta(str(bed), ref, str(out))
    text = out.read_text()
    assert text.count(">") == 2
    assert "ACGTACGT" in text


@bedtools_only
def test_bed_to_fasta_chrom_mismatch_raises(tmp_path: Path) -> None:
    ref = _write_ref(tmp_path)
    bed = tmp_path / "p.bed"
    bed.write_text("1\t0\t8\n")  # '1' not 'chr1' → bedtools skips, warns
    out = tmp_path / "out.fa"
    with pytest.raises(b2f.ConversionError):
        b2f.bed_to_fasta(str(bed), ref, str(out))


@bedtools_only
def test_convert_list_preserves_order_and_no_collision(tmp_path: Path) -> None:
    ref = _write_ref(tmp_path)
    d1 = tmp_path / "expA"
    d2 = tmp_path / "expB"
    d1.mkdir(); d2.mkdir()
    # Two BEDs with the SAME basename in different dirs, different content.
    (d1 / "peaks.bed").write_text("chr1\t0\t8\n")
    (d2 / "peaks.bed").write_text("chr2\t0\t8\n")
    beds = [str(d1 / "peaks.bed"), str(d2 / "peaks.bed")]
    out = b2f.convert_list(beds, ref, str(tmp_path / "fa"), threads=2)
    assert len(out) == 2
    assert out[0] != out[1]  # distinct output paths, no overwrite
    assert Path(out[0]).read_text() != Path(out[1]).read_text()
