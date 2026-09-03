"""Public mixed-stride spellings resolve consistently in the Python front end."""

from pathlib import Path

import pytest

from hammock import _core
from hammock.cli import parse_args


def _sketch(path: Path, rate: float, method: str | None, threads: int = 1):
    kwargs = dict(
        path=str(path), mode="B", precision=18, separator="\t",
        sub_a=1.0, sub_b=rate, exp_a=0.0, seed=42, gate_seed=31337,
        verbose=False, threads=threads,
    )
    if method is not None:
        kwargs["subB_method"] = method
    return _core.sketch_bed_file_hll(**kwargs)


def _same_registers(left, right) -> bool:
    return left.estimate_reg_eq_similarity(right) == 1.0


def test_binding_default_and_v2_alias_are_v2_at_non_integral_rate(tmp_path: Path):
    bed = tmp_path / "span.bed"
    bed.write_text("chrX\t0\t36000\n")

    default = _sketch(bed, 0.3, None)
    public = _sketch(bed, 0.3, "mixed-stride")
    v2_alias = _sketch(bed, 0.3, "mixed-stride-v2")
    legacy = _sketch(bed, 0.3, "mixed-stride-v1")

    assert _same_registers(default, public)
    assert _same_registers(public, v2_alias)
    assert not _same_registers(public, legacy)


def test_integral_rates_preserve_legacy_registers(tmp_path: Path):
    bed = tmp_path / "span.bed"
    bed.write_text("chr7\t135\t30170\n")

    for rate in (0.1, 0.01):
        public = _sketch(bed, rate, "mixed-stride")
        v2_alias = _sketch(bed, rate, "mixed-stride-v2")
        legacy = _sketch(bed, rate, "mixed-stride-v1")
        assert _same_registers(public, v2_alias)
        assert _same_registers(public, legacy)


def test_v2_default_is_thread_invariant(tmp_path: Path):
    bed = tmp_path / "partitioned.bed"
    bed.write_text("chrX\t0\t317\nchrX\t317\t811\nchrX\t811\t36000\n")
    assert _same_registers(_sketch(bed, 0.3, None, threads=1),
                           _sketch(bed, 0.3, None, threads=4))


def test_python_cli_accepts_all_public_method_names():
    for method in ("hash-threshold", "mixed-stride", "mixed-stride-v1",
                   "mixed-stride-v2", "single-hash"):
        args = parse_args(["queries.txt", "refs.txt", "--subB-method", method])
        assert args.subB_method == method


@pytest.mark.parametrize("option", ["--subA", "--subB"])
@pytest.mark.parametrize("value", ["-0.01", "1.01", "nan", "inf", "-inf"])
def test_python_cli_rejects_invalid_sampling_fractions(option, value):
    with pytest.raises(SystemExit):
        parse_args(["queries.txt", "refs.txt", option, value])


@pytest.mark.parametrize("argument", ["sub_a", "sub_b"])
@pytest.mark.parametrize("value", [-0.01, 1.01, float("nan"), float("inf"), -float("inf")])
def test_binding_rejects_invalid_sampling_fractions(tmp_path: Path, argument, value):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    kwargs = dict(path=str(bed), mode="B", precision=12)
    kwargs[argument] = value
    with pytest.raises(ValueError, match=argument + r" must be finite and in \[0, 1\]"):
        _core.sketch_bed_file_hll(**kwargs)


@pytest.mark.parametrize("method", ["mixed-stride", "mixed-stride-v1", "mixed-stride-v2"])
def test_python_cli_rejects_unrepresentable_mixed_stride_reciprocal(method):
    with pytest.raises(SystemExit):
        parse_args(["queries.txt", "refs.txt", "--subB", "1e-20",
                    "--subB-method", method])


@pytest.mark.parametrize("method", ["hash-threshold", "single-hash"])
def test_python_cli_allows_tiny_rates_without_integer_strides(method):
    args = parse_args(["queries.txt", "refs.txt", "--subB", "1e-20",
                       "--subB-method", method])
    assert args.subB == 1e-20


@pytest.mark.parametrize("method", ["mixed-stride", "mixed-stride-v1", "mixed-stride-v2"])
def test_binding_rejects_unrepresentable_mixed_stride_reciprocal(tmp_path: Path, method):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    with pytest.raises(ValueError, match="reciprocal exceeds INT64_MAX"):
        _core.sketch_bed_file_hll(
            path=str(bed), mode="B", precision=12, sub_b=1e-20,
            subB_method=method)
