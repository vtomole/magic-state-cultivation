import collections

import sinter

import gen
import cultiv
from ._stats_util import compute_expected_injection_growth_volume, split_by_custom_count


def test_compute_expected_injection_growth_volume():
    circuit = cultiv.make_inject_and_cultivate_circuit(dcolor=3, basis="Y", inject_style="unitary")
    circuit = gen.NoiseModel.uniform_depolarizing(1e-3).noisy_circuit_skipping_mpp_boundaries(circuit)
    v = compute_expected_injection_growth_volume(circuit)
    assert 70 <= v <= 72

    circuit = cultiv.make_end2end_cultivation_circuit(
        dcolor=3,
        dsurface=11,
        r_growing=2,
        r_end=2,
        basis="Y",
        inject_style="unitary",
    )
    circuit = gen.NoiseModel.uniform_depolarizing(1e-3).noisy_circuit_skipping_mpp_boundaries(circuit)
    v = compute_expected_injection_growth_volume(circuit)
    assert 1550 <= v <= 1650


def test_split_by_custom_count_preserves_non_numeric_values():
    stats = [
        sinter.TaskStats(
            strong_id="source",
            decoder="decoder",
            json_metadata={},
            shots=1,
            errors=0,
            discards=0,
            seconds=0,
            custom_counts=collections.Counter({"basis=unknown": 1}),
        )
    ]

    result = split_by_custom_count(stats)

    assert len(result) == 1
    assert result[0].json_metadata["basis"] == "unknown"
    assert result[0].json_metadata["hits"] == 1
