"""Unit tests for hammock.refs (reference resolution). No network access."""
from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from hammock import refs


def test_canonical_keyword_and_aliases() -> None:
    assert refs.canonical_keyword("hg38") == "hg38"
    assert refs.canonical_keyword("HG38") == "hg38"
    assert refs.canonical_keyword("GRCh38") == "hg38"
    assert refs.canonical_keyword("chm13") == "hs1"
    assert refs.canonical_keyword("t2t") == "hs1"
    assert refs.canonical_keyword("not-a-genome") is None


def _fake_cached_ref(cache: Path, keyword: str) -> Path:
    """Create a fake complete cache entry (fasta + .fai + .done)."""
    fa = cache / f"{keyword}.fa"
    fa.write_text(">chr1\nACGTACGTACGT\n")
    (cache / f"{keyword}.fa.fai").write_text("chr1\t12\t6\t12\t13\n")
    (cache / f"{keyword}.fa.done").write_text("source=test\nsize=1\n")
    return fa


def test_resolve_keyword_hits_cache(tmp_path: Path) -> None:
    fa = _fake_cached_ref(tmp_path, "hg38")
    assert refs.resolve_reference("hg38", str(tmp_path)) == str(fa)
    # Alias resolves to the same canonical cache entry.
    assert refs.resolve_reference("GRCh38", str(tmp_path)) == str(fa)


def test_resolve_keyword_missing_raises_with_fetch_hint(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as e:
        refs.resolve_reference("mm10", str(tmp_path))
    msg = str(e.value)
    assert "fetch-ref mm10" in msg
    assert str(tmp_path) in msg


def test_resolve_incomplete_cache_missing_done_raises(tmp_path: Path) -> None:
    # fasta + fai present but no .done sentinel → treated as incomplete.
    (tmp_path / "hg38.fa").write_text(">chr1\nACGT\n")
    (tmp_path / "hg38.fa.fai").write_text("chr1\t4\t6\t4\t5\n")
    with pytest.raises(ValueError):
        refs.resolve_reference("hg38", str(tmp_path))


def test_resolve_unknown_spec_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as e:
        refs.resolve_reference("wat", str(tmp_path))
    assert "Unrecognised reference" in str(e.value)


def test_resolve_url_in_run_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as e:
        refs.resolve_reference("https://example.com/hg38.fa.gz", str(tmp_path))
    assert "fetch-ref" in str(e.value)


def test_resolve_local_non_fasta_path_raises(tmp_path: Path) -> None:
    p = tmp_path / "notfasta.txt"
    p.write_text("hello\n")
    with pytest.raises(ValueError) as e:
        refs.resolve_reference(str(p), str(tmp_path / "cache"))
    assert "does not look like a FASTA" in str(e.value)


@pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools to faidx")
def test_resolve_local_fasta_builds_index(tmp_path: Path) -> None:
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "ACGT" * 10 + "\n")
    out = refs.resolve_reference(str(fa), str(tmp_path / "cache"))
    assert out == str(fa)
    assert (tmp_path / "ref.fa.fai").exists()


def test_fetch_reference_http_scheme_only(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as e:
        refs.fetch_reference("ftp://example.com/genome.fa.gz", str(tmp_path))
    assert "http(s)" in str(e.value)


def test_safe_key_is_traversal_free() -> None:
    key = refs._safe_key("../../etc/passwd")
    assert "/" not in key and ".." not in key
    assert key.startswith("cache_")
