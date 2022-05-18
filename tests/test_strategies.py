from itertools import islice
from mablibs.strategies import *
import pytest


@pytest.fixture()
def position_residue_mapping():
    return dict(zip(range(1, 5), ["A B C D".split()] * 4))


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            2,
            [
                ((1, "A"), (2, "A")),
                ((1, "A"), (2, "B")),
                ((1, "A"), (2, "C")),
                ((1, "A"), (2, "D")),
            ],
        ),
        (
            3,
            [
                ((1, "A"), (2, "A"), (3, "A")),
                ((1, "A"), (2, "A"), (3, "B")),
                ((1, "A"), (2, "A"), (3, "C")),
                ((1, "A"), (2, "A"), (3, "D")),
            ],
        ),
        (
            4,
            [
                ((1, "A"), (2, "A"), (3, "A"), (4, "A")),
                ((1, "A"), (2, "A"), (3, "A"), (4, "B")),
                ((1, "A"), (2, "A"), (3, "A"), (4, "C")),
                ((1, "A"), (2, "A"), (3, "A"), (4, "D")),
            ],
        ),
    ],
)
def test_randomization_strategy_get_mutations(
    position_residue_mapping, test_input, expected
):
    randomization = RandomizationStrategy(position_residue_mapping, test_input)
    assert list(it.islice(randomization.get_mutations(), 4)) == expected


@pytest.mark.parametrize("test_input,expected", [(2, 96), (3, 256), (4, 256)])
def test_randomization_strategy_library_size(
    position_residue_mapping, test_input, expected
):
    randomization = RandomizationStrategy(position_residue_mapping, test_input)
    assert randomization.library_size == expected


def test_randomization_strategy_zip_keys_to_vals(position_residue_mapping):
    r = RandomizationStrategy(position_residue_mapping)
    assert r._zip_keys_to_vals(position_residue_mapping) == [
        [(1, "A"), (1, "B"), (1, "C"), (1, "D")],
        [(2, "A"), (2, "B"), (2, "C"), (2, "D")],
        [(3, "A"), (3, "B"), (3, "C"), (3, "D")],
        [(4, "A"), (4, "B"), (4, "C"), (4, "D")],
    ]


def test_randomization_strategy_sample(position_residue_mapping):
    r = RandomizationStrategy(position_residue_mapping, n=2)
    assert list(r.sample(3, seed=42)) == [
        ((1, "A"), (2, "C")),
        ((1, "C"), (3, "C")),
        ((2, "D"), (3, "B")),
    ]


def test_randomization_strategy_mul_lens(position_residue_mapping):
    r = RandomizationStrategy(position_residue_mapping, n=2)
    assert r._mul_lens([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) == 27
