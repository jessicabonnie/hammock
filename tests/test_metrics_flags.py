"""Tests for the three-shape output contract: bare default (_ie),
--register-equality/--re (_re), --metrics (_full).

See docs/seed-metrics-column-restructure.md Part 2. Exercises Mode B
(interval) and Mode D (sequence) on both front-ends' shared CSV writer
(runner.py). Follows the subprocess + shutil.which("hammock") + skipif
convention already used by tests/test_bed2fasta_cli.py and tests/test_mode_d.py,
so the tests exercise the real console-script entrypoint the same way those do.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"

OURS = shutil.which("hammock")

pytestmark = pytest.mark.skipif(OURS is None, reason="hammock not on PATH")

_FULL_TAIL = ["reg_eq_similarity", "jaccard_similarity_ie",
              "containment_AB", "containment_BA",
              "cosketch_geom", "cosketch_arith", "cosketch_max",
              "register_equality_similarity"]
_RE_TAIL = ["reg_eq_similarity", "register_equality_similarity"]
_IE_TAIL = ["jaccard_similarity_ie"]

# (flag args, filename tag, expected trailing similarity columns)
_SHAPES = [
    ([], "ie", _IE_TAIL),
    (["--register-equality"], "re", _RE_TAIL),
    (["--metrics"], "full", _FULL_TAIL),
]


def _files_list(tmp_path: Path, *names: str) -> Path:
    f = tmp_path / (f"files_{'_'.join(names)}.txt")
    f.write_text("\n".join(str(DATA / n) for n in names) + "\n")
    return f


def _run(cmd: list[str], cwd: Path):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def _read_csv(path: Path):
    lines = path.read_text().splitlines()
    header = lines[0].split(",")
    rows = [line.split(",") for line in lines[1:]]
    return header, rows


def _run_shape(tmp_path: Path, extra_flags: list, mode: str, *fasta_or_bed: str,
              precision: str = "12") -> Path:
    files = _files_list(tmp_path, *fasta_or_bed)
    out = tmp_path / f"out_{mode}_{'_'.join(extra_flags).replace('-', '')}"
    cmd = [OURS, str(files), str(files), "--mode", mode, "-p", precision]
    if mode == "D":
        cmd += ["-k", "8", "-w", "40"]
    cmd += extra_flags + ["-o", str(out)]
    r = _run(cmd, tmp_path)
    assert r.returncode == 0, r.stderr
    matches = list(tmp_path.glob(f"{out.name}*.csv"))
    assert len(matches) == 1, (matches, r.stderr)
    return matches[0]


@pytest.mark.parametrize("flags,tag,tail", _SHAPES)
def test_mode_b_shape_header_and_filename(tmp_path: Path, flags, tag, tail) -> None:
    csv_path = _run_shape(tmp_path, flags, "B", "tiny_a.bed", "tiny_b.bed")
    assert csv_path.name.endswith(f"_{tag}.csv"), csv_path.name
    header, rows = _read_csv(csv_path)
    assert header[-len(tail):] == tail
    for row in rows:
        assert len(row) == len(header), (header, row)


@pytest.mark.parametrize("flags,tag,tail", _SHAPES)
def test_mode_d_shape_header_filename_and_ref_columns(tmp_path: Path, flags, tag, tail) -> None:
    csv_path = _run_shape(tmp_path, flags, "D", "tiny.fa", "tiny2.fa")
    # Tag comes after the Mode D k/w component in the filename.
    assert csv_path.name.endswith(f"_k8_w40_{tag}.csv"), csv_path.name
    header, rows = _read_csv(csv_path)
    # ref1/ref2 are Mode-D-specific trailing columns, present in every shape,
    # appended after the similarity block -- not dropped by any shape.
    assert header[-2:] == ["ref1", "ref2"]
    assert header[-2 - len(tail):-2] == tail
    for row in rows:
        # Catches a shape silently dropping/adding columns (row/header length
        # mismatch) even if the specific similarity columns checked above
        # happen to still be found correctly.
        assert len(row) == len(header), (header, row)
        assert row[-2:] == ["NA", "NA"]  # plain FASTA run, no --ref


def test_register_equality_matches_metrics_reg_eq_similarity(tmp_path: Path) -> None:
    """reg_eq_similarity must be byte-identical between --re and --metrics."""
    re_csv = _run_shape(tmp_path, ["--register-equality"], "B", "tiny_a.bed", "tiny_b.bed")
    full_csv = _run_shape(tmp_path, ["--metrics"], "B", "tiny_a.bed", "tiny_b.bed")
    re_header, re_rows = _read_csv(re_csv)
    full_header, full_rows = _read_csv(full_csv)
    re_val = re_rows[0][re_header.index("reg_eq_similarity")]
    full_val = full_rows[0][full_header.index("reg_eq_similarity")]
    assert re_val == full_val


def test_default_matches_metrics_jaccard_similarity_ie(tmp_path: Path) -> None:
    """jaccard_similarity_ie must be byte-identical between the bare default
    and --metrics."""
    ie_csv = _run_shape(tmp_path, [], "B", "tiny_a.bed", "tiny_b.bed")
    full_csv = _run_shape(tmp_path, ["--metrics"], "B", "tiny_a.bed", "tiny_b.bed")
    ie_header, ie_rows = _read_csv(ie_csv)
    full_header, full_rows = _read_csv(full_csv)
    ie_val = ie_rows[0][ie_header.index("jaccard_similarity_ie")]
    full_val = full_rows[0][full_header.index("jaccard_similarity_ie")]
    assert ie_val == full_val


def test_register_equality_similarity_duplicates_reg_eq_similarity(tmp_path: Path) -> None:
    for flags in (["--register-equality"], ["--metrics"]):
        csv_path = _run_shape(tmp_path, flags, "B", "tiny_a.bed", "tiny_b.bed")
        header, rows = _read_csv(csv_path)
        j = rows[0][header.index("reg_eq_similarity")]
        re_sim = rows[0][header.index("register_equality_similarity")]
        assert j == re_sim, (flags, header, rows)


def test_metrics_and_register_equality_are_mutually_exclusive(tmp_path: Path) -> None:
    files = _files_list(tmp_path, "tiny_a.bed", "tiny_b.bed")
    r = _run([OURS, str(files), str(files), "--mode", "B", "-p", "12",
              "--metrics", "--register-equality"], tmp_path)
    assert r.returncode != 0
    assert "not allowed" in r.stderr or "mutually exclusive" in r.stderr


def test_re_is_an_alias_for_register_equality(tmp_path: Path) -> None:
    csv_path = _run_shape(tmp_path, ["--re"], "B", "tiny_a.bed", "tiny_b.bed")
    assert csv_path.name.endswith("_re.csv")
    header, _ = _read_csv(csv_path)
    assert header[-2:] == _RE_TAIL
