"""Mode D process dispatch must work under spawn and preserve exact HLL state."""
from __future__ import annotations

import gzip
import multiprocessing as mp
from pathlib import Path
from types import SimpleNamespace

import pytest

from hammock import _core, bed2fasta, cli, runner


DATA = Path(__file__).parent / "data"


def _block_forever() -> None:
    """Picklable spawn target for shutdown escalation coverage."""
    import time

    time.sleep(60)


def _args(threads: int, sequence_hll_hash: str = "rehash-selector64") -> SimpleNamespace:
    return SimpleNamespace(
        mode="D",
        threads=threads,
        kmer_size=8,
        window_size=40,
        seed=42,
        precision=12,
        sequence_hll_hash=sequence_hll_hash,
        sketch_type="minimizer",
        verbose=False,
    )


@pytest.mark.parametrize(
    "sequence_hll_hash", ["rehash-selector64", "legacy-selector32"]
)
def test_mode_d_spawn_processes_match_serial_and_keep_order(
    sequence_hll_hash: str,
) -> None:
    # Repetition makes completion order potentially differ from input order.
    paths = [str(DATA / "tiny2.fa"), str(DATA / "tiny.fa"), str(DATA / "tiny2.fa")]
    serial = runner._sketch_many(paths, _args(1, sequence_hll_hash), "query")
    spawned = runner._sketch_many(paths, _args(2, sequence_hll_hash), "query")

    assert [
        _core._hll_transport_state(sketch.minimizer_hll)
        for sketch in spawned
    ] == [
        _core._hll_transport_state(sketch.minimizer_hll)
        for sketch in serial
    ]
    assert all(sketch.hash_mode == sequence_hll_hash for sketch in spawned)


def test_mode_d_worker_failure_identifies_path_and_leaves_no_child(tmp_path: Path) -> None:
    missing = str(tmp_path / "missing.fa")
    with pytest.raises(runner.ModeDSketchWorkerError, match="missing.fa"):
        runner._sketch_many([missing, str(DATA / "tiny.fa")], _args(2), "query")

    assert not mp.active_children()


def test_cli_worker_failure_is_concise_and_publishes_no_csv(tmp_path: Path, capsys) -> None:
    missing = tmp_path / "missing.fa"
    paths = tmp_path / "paths.txt"
    paths.write_text(f"{missing}\n{DATA / 'tiny.fa'}\n")

    rc = cli.main([
        str(paths), str(paths), "--mode", "D", "--threads", "2",
        "-p", "12", "-k", "8", "-w", "40", "-o", str(tmp_path / "out"),
    ])

    captured = capsys.readouterr()
    assert rc == 1
    assert "missing.fa" in captured.err
    assert "Traceback" not in captured.err
    assert not list(tmp_path.glob("out*.csv"))
    assert not list(tmp_path.glob(".*.tmp"))


def test_mode_d_spawn_processes_support_gzipped_fasta(tmp_path: Path) -> None:
    gz_path = tmp_path / "tiny.fa.gz"
    with open(DATA / "tiny.fa", "rb") as source, gzip.open(gz_path, "wb") as target:
        target.write(source.read())

    serial = runner._sketch_many([str(gz_path)], _args(1), "query")
    # Two entries exercise both spawn dispatch and duplicate-path ordering.
    spawned = runner._sketch_many([str(gz_path), str(gz_path)], _args(2), "query")
    expected = _core._hll_transport_state(serial[0].minimizer_hll)
    assert [_core._hll_transport_state(s.minimizer_hll) for s in spawned] == [expected, expected]


def test_cli_process_csv_is_byte_identical_and_atomic(tmp_path: Path) -> None:
    paths = tmp_path / "paths.txt"
    paths.write_text(f"{DATA / 'tiny.fa'}\n{DATA / 'tiny2.fa'}\n")
    common = [str(paths), str(paths), "--mode", "D", "-p", "12", "-k", "8", "-w", "40"]

    assert cli.main([*common, "--threads", "1", "-o", str(tmp_path / "serial")]) == 0
    assert cli.main([*common, "--threads", "2", "-o", str(tmp_path / "process")]) == 0

    serial = tmp_path / "serial_mnmzr_p12_jaccD_k8_w40_rehash-selector64_ie.csv"
    process = tmp_path / "process_mnmzr_p12_jaccD_k8_w40_rehash-selector64_ie.csv"
    assert serial.read_bytes() == process.read_bytes()
    assert not list(tmp_path.glob(".*.tmp"))


def test_bounded_shutdown_terminates_a_blocked_spawn_worker() -> None:
    context = mp.get_context("spawn")
    tasks = context.Queue(maxsize=1)
    results = context.Queue(maxsize=1)
    worker = context.Process(target=_block_forever)
    worker.start()

    runner._stop_mode_d_workers([worker], tasks, results, graceful=False)

    assert not worker.is_alive()
    assert worker.exitcode is not None


def test_cli_interrupt_exits_130_without_traceback(tmp_path: Path, monkeypatch, capsys) -> None:
    paths = tmp_path / "paths.txt"
    paths.write_text(f"{DATA / 'tiny.fa'}\n")

    def interrupt(_args):
        raise KeyboardInterrupt

    monkeypatch.setattr(cli, "run", interrupt)
    rc = cli.main([str(paths), str(paths), "--mode", "D"])

    captured = capsys.readouterr()
    assert rc == 130
    assert captured.err == "Interrupted.\n"


def test_bed2fasta_then_spawn_sketching_keeps_io_budget(tmp_path: Path, monkeypatch) -> None:
    """Extraction can use its own pool before Mode D starts spawn workers."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTACGTACGTACGT\n")
    list1 = tmp_path / "list1.txt"
    list2 = tmp_path / "list2.txt"
    list1.write_text("first.bed\nsecond.bed\n")
    list2.write_text("third.bed\nfourth.bed\n")
    seen_threads = []

    def fake_convert_list(_beds, _ref, _outdir, threads, _verbose):
        seen_threads.append(threads)
        return [str(DATA / "tiny.fa"), str(DATA / "tiny2.fa")]

    monkeypatch.setattr(bed2fasta, "convert_list", fake_convert_list)
    rc = cli.main([
        str(list1), str(list2), "--ref", str(ref), "--threads", "2",
        "-p", "12", "-k", "8", "-w", "40", "-o", str(tmp_path / "out"),
    ])

    assert rc == 0
    assert seen_threads == [2, 2]
    assert list(tmp_path.glob("out*.csv"))
