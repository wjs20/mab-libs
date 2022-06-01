import re
from collections import namedtuple

amino_acid_ptm_specification = (
    ("GLYCOSYLATION_MOTIF", r"N[^P][ST]"),
    ("DEAMIDATION_MOTIF", r"N[GSA]"),
    ("ISOMERIZATION_MOTIF", r"D[GS]"),
    ("CLEAVAGE_MOTIF", r"DP"),
    ("OXIDATION_MOTIF", r"M"),
)

amino_acid_ptm_regex = re.compile(
    "|".join(
        f"(?P<{name}>{pattern})"
        for name, pattern in amino_acid_ptm_specification
    )
)

PTMmotif = namedtuple("PTMmotif", "kind match start end")


def find_ptm_motifs(seq, report=False):
    if not (matches := list(amino_acid_ptm_regex.finditer(seq))):
        return
    return [
        PTMmotif(mo.lastgroup, mo.group(0), mo.start(), mo.end())
        for mo in matches
    ]
