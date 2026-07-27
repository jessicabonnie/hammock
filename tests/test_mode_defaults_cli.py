"""End-to-end checks that the default/aliased --mode resolves correctly and
lands the right letter in the CSV `mode` column."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"
OURS = shutil.which("hammock")

pytestmark = pytest.mark.skipif(OURS is None, reason="hammock not on PATH")


def _bed_list(tmp_path: Path) -> Path:
    f = tmp_path / "beds.txt"
    f.write_text(str(DATA / "tiny_a.bed") + "\n")
    return f


def _mode_col(csv_path: Path) -> str:
    lines = csv_path.read_text().splitlines()
    header = lines[0].split(",")
    return lines[1].split(",")[header.index("mode")]


def test_default_bed_mode_is_interval_points_b(tmp_path: Path) -> None:
    beds = _bed_list(tmp_path)
    r = subprocess.run([OURS, str(beds), str(beds), "-p", "12", "-o", str(tmp_path / "d")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    csv = next(tmp_path.glob("d*.csv"))
    assert _mode_col(csv) == "B"                 # interval-points is the default
    assert "jaccB" in csv.name                    # filename keeps the letter


def test_mode_name_alias_interval_string_maps_to_a(tmp_path: Path) -> None:
    beds = _bed_list(tmp_path)
    r = subprocess.run([OURS, str(beds), str(beds), "--mode", "interval-string",
                        "-p", "12", "-o", str(tmp_path / "s")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert _mode_col(next(tmp_path.glob("s*.csv"))) == "A"


def test_mode_name_alias_sequence_rejects_bad_value(tmp_path: Path) -> None:
    beds = _bed_list(tmp_path)
    r = subprocess.run([OURS, str(beds), str(beds), "--mode", "bogus"],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "invalid mode" in r.stderr
