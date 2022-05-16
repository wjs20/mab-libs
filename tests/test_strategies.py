from mablibs.strategies import *
import pytest


@pytest.fixture()
def position_residue_mapping():
    return {1: ["A", "C", "D", "E"], 3: ["A", "C", "D", "E"], 6: ["A", "C", "D", "E"]}


@pytest.fixture
def twomerrandomization(position_residue_mapping):
    return NmerRandomizationStrategy(2, position_residue_mapping)


@pytest.fixture
def randomtrimerrandomization(position_residue_mapping):
    return RandomTrimerRandomizationStrategy(position_residue_mapping)


def test_randomization_zip_keys_to_vals(position_residue_mapping):
    r = RandomizationStrategy(position_residue_mapping)
    assert r.zip_keys_to_vals(position_residue_mapping) == [
        [(1, "A"), (1, "C"), (1, "D"), (1, "E")],
        [(3, "A"), (3, "C"), (3, "D"), (3, "E")],
        [(6, "A"), (6, "C"), (6, "D"), (6, "E")],
    ]


def test_twomerrandomization_get_mutations(twomerrandomization):
    assert twomerrandomization.get_mutations(generator=False, first=10) == [
        ((1, "A"), (3, "A")),
        ((1, "A"), (3, "C")),
        ((1, "A"), (3, "D")),
        ((1, "A"), (3, "E")),
        ((1, "A"), (6, "A")),
        ((1, "A"), (6, "C")),
        ((1, "A"), (6, "D")),
        ((1, "A"), (6, "E")),
        ((1, "C"), (3, "A")),
        ((1, "C"), (3, "C")),
    ]


def test_randomtrimerrandomization_get_mutations(randomtrimerrandomization):
    assert randomtrimerrandomization.get_mutations(generator=False, first=10) == [
        ((1, "A"), (3, "A"), (6, "A")),
        ((1, "A"), (3, "A"), (6, "C")),
        ((1, "A"), (3, "A"), (6, "D")),
        ((1, "A"), (3, "A"), (6, "E")),
        ((1, "A"), (3, "C"), (6, "A")),
        ((1, "A"), (3, "C"), (6, "C")),
        ((1, "A"), (3, "C"), (6, "D")),
        ((1, "A"), (3, "C"), (6, "E")),
        ((1, "A"), (3, "D"), (6, "A")),
        ((1, "A"), (3, "D"), (6, "C")),
    ]


def test_randomtrimerrandomization_library_size(randomtrimerrandomization):
    assert randomtrimerrandomization.library_size == 64


def test_twomerrandomization_library_size(twomerrandomization):
    assert twomerrandomization.library_size == 48


def test_sample(twomerrandomization):
    assert [mutations for mutations in twomerrandomization.sample(10, seed=42)] == [
        ((1, "A"), (3, "C")),
        ((1, "A"), (6, "A")),
        ((1, "C"), (3, "D")),
        ((1, "C"), (6, "C")),
        ((1, "D"), (6, "A")),
        ((1, "E"), (6, "D")),
        ((3, "A"), (6, "A")),
        ((3, "A"), (6, "E")),
        ((3, "D"), (6, "D")),
    ]


def test_mul_lens(twomerrandomization):
    assert twomerrandomization._mul_lens([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) == 27
