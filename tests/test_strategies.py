from itertools import islice
from typing import Dict

import pytest
from mablibs.strategies import *


@pytest.fixture()
def position_residue_mapping() -> Dict[int, str]:
    return dict(zip(range(1, 5), ["A C D E".split()] * 4))


def test_randomization_strategy(position_residue_mapping) -> None:
    randomization = RandomizationStrategy(position_residue_mapping)


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            2,
            [
                ((1, "A"), (2, "A")),
                ((1, "A"), (2, "C")),
                ((1, "A"), (2, "D")),
                ((1, "A"), (2, "E")),
            ],
        ),
        (
            3,
            [
                ((1, "A"), (2, "A"), (3, "A")),
                ((1, "A"), (2, "A"), (3, "C")),
                ((1, "A"), (2, "A"), (3, "D")),
                ((1, "A"), (2, "A"), (3, "E")),
            ],
        ),
        (
            4,
            [
                ((1, "A"), (2, "A"), (3, "A"), (4, "A")),
                ((1, "A"), (2, "A"), (3, "A"), (4, "C")),
                ((1, "A"), (2, "A"), (3, "A"), (4, "D")),
                ((1, "A"), (2, "A"), (3, "A"), (4, "E")),
            ],
        ),
    ],
)
def test_randomization_strategy_get_mutations(
    position_residue_mapping, test_input, expected
) -> None:
    randomization = RandomizationStrategy(position_residue_mapping, test_input)
    assert list(it.islice(randomization.get_mutations(), 4)) == expected


@pytest.mark.parametrize("test_input,expected", [(2, 96), (3, 256), (4, 256)])
def test_randomization_strategy_len(
    position_residue_mapping, test_input, expected
) -> None:
    randomization = RandomizationStrategy(position_residue_mapping, test_input)
    assert len(randomization) == expected


def test_randomization_strategy_zip_keys_to_vals(
    position_residue_mapping,
) -> None:
    r = RandomizationStrategy(position_residue_mapping)
    assert r.zip_keys_to_vals(position_residue_mapping) == [
        [(1, "A"), (1, "C"), (1, "D"), (1, "E")],
        [(2, "A"), (2, "C"), (2, "D"), (2, "E")],
        [(3, "A"), (3, "C"), (3, "D"), (3, "E")],
        [(4, "A"), (4, "C"), (4, "D"), (4, "E")],
    ]


def test_randomization_strategy_sample(position_residue_mapping) -> None:
    r = RandomizationStrategy(position_residue_mapping, n=2)
    assert list(r.sample(3, seed=42)) == [
        ((1, "A"), (2, "D")),
        ((1, "D"), (3, "D")),
        ((2, "E"), (3, "C")),
    ]


def test_randomization_strategymul_lens(position_residue_mapping) -> None:
    r = RandomizationStrategy(position_residue_mapping, n=2)
    assert r.mul_lens([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) == 27
