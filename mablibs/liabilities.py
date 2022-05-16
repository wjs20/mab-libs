import re
from collections import namedtuple

amino_acid_liability_specification = (
    ("GLYCOSYLATION_MOTIF", r"N[^P][ST]"),
    ("DEAMIDATION_MOTIF", r"N[GSA]"),
    ("ISOMERIZATION_MOTIF", r"D[GS]"),
    ("CLEAVAGE_MOTIF", r"DP"),
    ("OXIDATION_MOTIF", r"M"),
)

amino_acid_liability_regex = re.compile(
    "|".join("(?P<%s>%s)" % pair for pair in amino_acid_liability_specification)
)

Liability = namedtuple("Liability", "kind match start end")


def find_liabilities(seq, report=False):
    matches = list(amino_acid_liability_regex.finditer(seq))
    if not matches:
        return
    if not report:
        return True
    return [
        Liability(mo.lastgroup, mo.group(0), mo.start(), mo.end()) for mo in matches
    ]


if __name__ == "__main__":
    test_seqs = ["ATG"]

    print(find_liabilities(test_seqs[0]))
