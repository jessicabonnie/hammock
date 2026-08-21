"""Validated HLL transport used by Mode D spawn workers."""
from __future__ import annotations

import pytest

from hammock import _core


def _populated(precision: int) -> "_core.HLLSketch":
    sketch = _core.HLLSketch(precision)
    for value in range(10_000):
        sketch.add_hash64((value * 0x9E3779B185EBCA87) & ((1 << 64) - 1))
    return sketch


@pytest.mark.parametrize("precision", [4, 12, 18, 24])
def test_transport_round_trip_is_byte_exact(precision: int) -> None:
    original = _populated(precision)
    state = _core._hll_transport_state(original)
    restored = _core._hll_from_transport_state(state)

    assert _core._hll_transport_state(restored) == state
    assert restored.estimate_cardinality() == original.estimate_cardinality()
    assert restored.estimate_reg_eq_similarity(original) == 1.0
    assert restored.estimate_intersection(original) == original.estimate_cardinality()


def test_transport_round_trips_valid_32_bit_state() -> None:
    state = ("hammock.HLLSketch.transport", 1, 4, 32, bytes(range(16)))
    restored = _core._hll_from_transport_state(state)

    assert _core._hll_transport_state(restored) == state


@pytest.mark.parametrize(
    "state, message",
    [
        (("hammock.HLLSketch.transport", 1, 4, 64), "length"),
        (("wrong", 1, 4, 64, bytes(16)), "magic"),
        (("hammock.HLLSketch.transport", 2, 4, 64, bytes(16)), "version"),
        (("hammock.HLLSketch.transport", True, 4, 64, bytes(16)), "types"),
        (("hammock.HLLSketch.transport", 1, 3, 64, bytes(8)), "precision"),
        (("hammock.HLLSketch.transport", 1, 4, 33, bytes(16)), "hash_size"),
        (("hammock.HLLSketch.transport", 1, 4, 64, bytes(15)), "length"),
        (("hammock.HLLSketch.transport", 1, 4, 64, bytes([62]) + bytes(15)),
         "out-of-range"),
    ],
)
def test_transport_rejects_malformed_state(state, message: str) -> None:
    with pytest.raises((TypeError, ValueError, OverflowError), match=message):
        _core._hll_from_transport_state(state)
