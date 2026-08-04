"""End-to-end checks that the default/aliased --mode resolves correctly and
lands the right letter in the CSV `mode` column."""
from __future__ import annotations

import os
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


def _fasta_list(tmp_path: Path) -> Path:
    f = tmp_path / "fastas.txt"
    f.write_text(str(DATA / "tiny.fa") + "\n")
    return f


_CONVOY_WARNING = "GIL convoy"


def test_mode_d_default_threads_is_silent(tmp_path: Path) -> None:
    """No --threads in mode D: default resolves to 1, no warning."""
    fas = _fasta_list(tmp_path)
    r = subprocess.run([OURS, str(fas), str(fas), "--mode", "D",
                        "-p", "12", "-k", "8", "-w", "40", "-o", str(tmp_path / "d")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert _CONVOY_WARNING not in r.stderr
    assert _mode_col(next(tmp_path.glob("d*.csv"))) == "D"


def test_mode_d_explicit_threads_warns_but_runs(tmp_path: Path) -> None:
    """--threads > 1 in mode D is still honored, with a one-line stderr note."""
    fas = _fasta_list(tmp_path)
    r = subprocess.run([OURS, str(fas), str(fas), "--mode", "D", "--threads", "4",
                        "-p", "12", "-k", "8", "-w", "40", "-o", str(tmp_path / "t")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert _CONVOY_WARNING in r.stderr


def test_interval_mode_threads_never_warn(tmp_path: Path) -> None:
    """The convoy note is Mode D only — A/B/C threading is real parallelism."""
    beds = _bed_list(tmp_path)
    r = subprocess.run([OURS, str(beds), str(beds), "--mode", "B", "--threads", "4",
                        "-p", "12", "-o", str(tmp_path / "b")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert _CONVOY_WARNING not in r.stderr


def _resolved_args(argv, monkeypatch):
    """Run cli.main far enough to resolve defaults, capturing the args it built."""
    from hammock import cli
    captured = {}

    def _fake_run(args):
        captured["args"] = args
        return 0

    # cli.py does `from hammock.runner import run`, so patch the name in cli.
    monkeypatch.setattr(cli, "run", _fake_run)
    assert cli.main(argv) == 0
    return captured["args"]


def test_mode_d_resolves_sketch_threads_1_but_keeps_io_threads(tmp_path, monkeypatch) -> None:
    """bed2fasta extraction shells out to bedtools and must NOT be clamped to 1
    by the Mode D convoy fix — it reads args.io_threads (runner._run_bed2fasta)."""
    fas = _fasta_list(tmp_path)
    args = _resolved_args([str(fas), str(fas), "--mode", "D", "-o", str(tmp_path / "x")],
                          monkeypatch)
    assert args.threads == 1
    assert args.io_threads == min(8, os.cpu_count() or 1)


def test_explicit_threads_drives_both_budgets(tmp_path, monkeypatch) -> None:
    fas = _fasta_list(tmp_path)
    args = _resolved_args([str(fas), str(fas), "--mode", "D", "--threads", "4",
                           "-o", str(tmp_path / "x")], monkeypatch)
    assert args.threads == 4 and args.io_threads == 4


def test_mode_name_alias_sequence_rejects_bad_value(tmp_path: Path) -> None:
    beds = _bed_list(tmp_path)
    r = subprocess.run([OURS, str(beds), str(beds), "--mode", "bogus"],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "invalid mode" in r.stderr
