"""Interval modes read plain-text BED only: the C++ parser is a bare
`std::ifstream`, so a compressed or binary input sketches to nothing and scores
0.0 against every file *including itself*, with exit 0. `runner._guard_plain_text_bed`
rejects the common offenders up front.

Pure-Python (the guard never touches the extension), so these run while a
rebuild is pending.
"""
from __future__ import annotations

import gzip
import os

import pytest

from hammock import runner


# (filename, leading bytes, substring expected in the message)
BINARY_CASES = [
    ("a.bed.gz", b"\x1f\x8b\x08\x00rest", "gzip-compressed"),
    ("a.bed.zst", b"\x28\xb5\x2f\xfdrest", "zstd-compressed"),
    ("a.bed.bz2", b"BZh9rest", "bzip2-compressed"),
    ("a.bed.xz", b"\xfd7zXZ\x00rest", "xz-compressed"),
    ("a.bb", b"\xeb\xf2\x89\x87rest", "BigBed"),
    ("a.bigbed", b"\x87\x89\xf2\xebrest", "BigBed"),
    ("a.bw", b"\x26\xfc\x8f\x88rest", "BigWig"),
]


@pytest.mark.parametrize("name,head,expected", BINARY_CASES)
def test_binary_inputs_are_rejected(tmp_path, name, head, expected):
    p = tmp_path / name
    p.write_bytes(head)
    msg = runner._guard_plain_text_bed([str(p)], "Query")
    assert msg is not None, f"{name} should have been rejected"
    assert expected in msg
    assert "Query" in msg


# Plain text must pass, including the degenerate lengths that a naive
# `head[:4]` comparison could mishandle.
@pytest.mark.parametrize("content", [
    b"chr1\t100\t200\n",
    b"",                       # empty file
    b"ch",                     # shorter than the longest magic
    b"\xef\xbb\xbfchr1\t1\t2\n",   # UTF-8 BOM
    b"chr1\t100\t200\r\n",     # CRLF
    b"track name=x\nchr1\t1\t2\n",
])
def test_plain_text_passes(tmp_path, content):
    p = tmp_path / "a.bed"
    p.write_bytes(content)
    assert runner._guard_plain_text_bed([str(p)], "Query") is None


def test_non_regular_files_are_skipped(tmp_path):
    """Reading from a FIFO consumes bytes the C++ reopen then blocks forever
    waiting for -- and `zcat foo.bed.gz > pipe` is exactly the workaround the
    rejection message invites. The guard must not touch it."""
    fifo = tmp_path / "pipe"
    os.mkfifo(fifo)
    # No writer attached: if the guard opened this it would block and the test
    # would hang rather than fail.
    assert runner._guard_plain_text_bed([str(fifo)], "Query") is None


def test_missing_file_defers_to_downstream_error(tmp_path):
    assert runner._guard_plain_text_bed([str(tmp_path / "nope.bed")], "Query") is None


def test_suggested_command_does_not_clobber_an_existing_file(tmp_path):
    """`foo.bed.gz` -> `foo.bed` would overwrite a real neighbouring foo.bed."""
    real = tmp_path / "a.bed"
    real.write_text("chr1\t1\t2\n")
    gz = tmp_path / "a.bed.gz"
    gz.write_bytes(gzip.compress(b"chr1\t1\t2\n"))

    msg = runner._guard_plain_text_bed([str(gz)], "Query")
    assert msg is not None
    # Check the redirect TARGET specifically: a naive substring test would match
    # the source path, since ".../a.bed" is a prefix of ".../a.bed.gz".
    target = msg.split(">")[-1].strip()
    assert target != str(real)
    assert not os.path.exists(target)


def test_message_names_the_side(tmp_path):
    p = tmp_path / "a.bed.gz"
    p.write_bytes(b"\x1f\x8b\x08\x00rest")
    assert "Primary" in runner._guard_plain_text_bed([str(p)], "Primary")
