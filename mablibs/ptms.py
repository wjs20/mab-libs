import re
from collections import namedtuple
from typing import Union

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


def find_ptm_motifs(seq: str, verbose: bool = True) -> Union[PTMmotif, bool]:
    if not (matches := tuple(amino_acid_ptm_regex.finditer(seq))):
        return False
    if verbose:
        return [
            PTMmotif(mo.lastgroup, mo.group(0), mo.start(), mo.end())
            for mo in matches
        ]
    return True
