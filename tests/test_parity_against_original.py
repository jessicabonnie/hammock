"""Parity tests: our `hammock` CLI must match `hammock-orig` (the upstream
hammock 0.4.0 installed via pipx) byte-for-byte on the same fixtures.

Skipped automatically if `hammock-orig` isn't on PATH.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"

ORIG = shutil.which("hammock-orig")
OURS = shutil.which("hammock")

pytestmark = pytest.mark.skipif(
    ORIG is None or OURS is None,
    reason="Both hammock-orig and hammock must be on PATH for parity tests.",
)


def _run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def _files_list(tmp_path: Path) -> Path:
    f = tmp_path / "files.txt"
    f.write_text(f"{DATA / 'tiny_a.bed'}\n{DATA / 'tiny_b.bed'}\n")
    return f


_PROJECTED_OUT = {
    # Orig 0.4.0 emitted an unreliable `containment` column; we no longer
    # emit that name. Our well-defined replacements are containment_AB /
    # containment_BA and the three cosketch summaries, none of which the
    # orig has.
    "containment",
    "containment_AB", "containment_BA",
    "cosketch_geom", "cosketch_arith", "cosketch_max",
    # A *Jaccard* column is projected out here, which looks wrong at a glance.
    # It isn't: jaccard_similarity_ie is a second, inclusion-exclusion
    # estimator emitted alongside the register-equality one. Orig has no
    # counterpart, so there is nothing to be unfaithful to. The column that
    # must stay byte-equal -- jaccard_similarity -- is still compared.
    "jaccard_similarity_ie",
}


def _projected_rows(csv_text: str) -> list[tuple]:
    """Drop containment/cosketch and inclusion-exclusion columns before comparing.

    Parity is required for `jaccard_similarity` (the register-equality column
    orig also emits); orig 0.4.0 has no counterpart for our containment,
    cosketch, or inclusion-exclusion surface.
    """
    lines = csv_text.strip().split("\n")
    header = lines[0].split(",")
    keep = [i for i, name in enumerate(header) if name not in _PROJECTED_OUT]
    return [tuple(line.split(",")[i] for i in keep) for line in lines]


# NOTE on parity scope. Three ways our hammock differs from hammock-orig:
#   - Mode B + --subB: orig's intervals.py:542 silently ignores subsample[1]
#     in Mode B. Our Mode B honors --subB; we accept the divergence.
#   - --subB-method=mixed-stride: hammock_claude extension. Default since
#     we dropped orig-parity as the default. To compare against orig we
#     must explicitly pass --subB-method=hash-threshold.
#   - --subB-method=single-hash: opt-in parity divergence (one xxh64 for
#     gate+ingestion). Not tested against orig.
@pytest.mark.parametrize("mode,extra", [
    ("A", []),
    ("B", []),
    ("C", []),
    # subB-using Mode C tests must override the default mixed-stride back to
    # hash-threshold to match orig's contract.
    ("C", ["--subB", "0.5", "--subB-method", "hash-threshold"]),
    ("C", ["--subA", "0.5"]),
    ("C", ["--expA", "1.0"]),
    ("C", ["--subB", "0.3", "--expA", "0.5", "--subB-method", "hash-threshold"]),
])
def test_jaccard_byte_equal(tmp_path: Path, mode: str, extra: list[str]) -> None:
    files = _files_list(tmp_path)
    # Drop our-only flags (and their args) before invoking orig.
    our_only_flags = {"--subB-method"}
    orig_extra: list[str] = []
    skip_next = False
    for a in extra:
        if skip_next:
            skip_next = False
            continue
        if a in our_only_flags:
            skip_next = True
            continue
        orig_extra.append(a)
    common = [str(files), str(files), "--mode", mode, "-p", "14"]
    _run([ORIG, *common, *orig_extra, "-o", str(tmp_path / "orig")], tmp_path)
    _run([OURS, *common, *extra, "-o", str(tmp_path / "ours")], tmp_path)

    orig_csvs = list(tmp_path.glob("orig*.csv"))
    ours_csvs = list(tmp_path.glob("ours*.csv"))
    assert len(orig_csvs) == 1, f"expected 1 orig CSV, got {orig_csvs}"
    assert len(ours_csvs) == 1, f"expected 1 ours CSV, got {ours_csvs}"

    orig_rows = _projected_rows(orig_csvs[0].read_text())
    ours_rows = _projected_rows(ours_csvs[0].read_text())
    assert orig_rows == ours_rows, (
        f"Jaccard-column mismatch for mode={mode} extra={extra}\n"
        f"orig (first 3): {orig_rows[:3]}\n"
        f"ours (first 3): {ours_rows[:3]}"
    )


# Functional regression tests: these verify behaviors that intentionally
# diverge from hammock-orig (broken there, fixed here).

@pytest.mark.skipif(OURS is None, reason="hammock not on PATH")
def test_mode_b_subB_actually_subsamples(tmp_path: Path) -> None:
    """Our Mode B honors --subB (orig silently ignores it). Verify the
    Jaccard with subB=0.25 differs meaningfully from subB=1.0 — i.e. the
    flag is doing something."""
    files = _files_list(tmp_path)
    common = [str(files), str(files), "--mode", "B", "-p", "14"]
    _run([OURS, *common, "-o", str(tmp_path / "full")], tmp_path)
    _run([OURS, *common, "--subB", "0.25", "-o", str(tmp_path / "sub")], tmp_path)

    full = next(tmp_path.glob("full*.csv")).read_text().splitlines()
    sub = next(tmp_path.glob("sub*.csv")).read_text().splitlines()
    # cross-pair Jaccard column should differ between full sampling and 25% sampling
    jac_col = full[0].split(",").index("jaccard_similarity")
    full_jac = full[2].split(",")[jac_col]
    sub_jac = sub[2].split(",")[jac_col]
    assert full_jac != sub_jac, f"--subB 0.25 produced same Jaccard as no subsampling ({full_jac})"


@pytest.mark.skipif(OURS is None, reason="hammock not on PATH")
@pytest.mark.parametrize("mode", ["A", "B", "C"])
def test_threads_match_sequential(tmp_path: Path, mode: str) -> None:
    """--threads N must produce byte-equal output to single-threaded run for
    every interval mode (sketch_bed_file_hll releases the GIL, so the parallel
    path must produce the same registers)."""
    files = _files_list(tmp_path)
    common = [str(files), str(files), "--mode", mode, "-p", "14"]
    _run([OURS, *common, "-o", str(tmp_path / "seq")], tmp_path)
    _run([OURS, *common, "--threads", "4", "-o", str(tmp_path / "par")], tmp_path)
    seq = next(tmp_path.glob("seq*.csv")).read_text()
    par = next(tmp_path.glob("par*.csv")).read_text()
    assert seq == par


@pytest.mark.skipif(OURS is None, reason="hammock not on PATH")
def test_mode_b_mixed_stride_is_deterministic(tmp_path: Path) -> None:
    """Re-running --mode B with mixed-stride (the default) must produce
    byte-equal output across invocations."""
    files = _files_list(tmp_path)
    common = [str(files), str(files), "--mode", "B", "-p", "14",
              "--subB", "0.25"]
    _run([OURS, *common, "-o", str(tmp_path / "run1")], tmp_path)
    _run([OURS, *common, "-o", str(tmp_path / "run2")], tmp_path)
    a = next(tmp_path.glob("run1*.csv")).read_text()
    b = next(tmp_path.glob("run2*.csv")).read_text()
    assert a == b
