import pytest
from mablibs.mutagenesis import *


@pytest.mark.parametrize(
    "test_input,expected",
    [("human", "TTC"), ("yeast", "TTA"), ("e_coli", "TTT"), ("hamster", "TTC")],
)
def test_get_preferred_codons(test_input, expected):
    get_preferred_codons(test_input)["F"] == expected
