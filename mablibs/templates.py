import json
import operator
import pkgutil
from typing import List, Tuple

AA2CODON = json.loads(
    pkgutil.get_data("mablibs", "amino_acid2codon.json").decode("utf-8")
)

CODON2AA = json.loads(
    pkgutil.get_data("mablibs", "codon2amino_acid.json").decode("utf-8")
)


def translate(sequence_nt: str) -> str:
    """Method converts nucleotide sequences to amino acid sequences"""
    seq_len = len(sequence_nt)
    assert seq_len % 3 == 0, "Ensure sequence contains complete reading frame"
    return "".join(
        [
            CODON2AA[codon]
            for codon in [sequence_nt[i : i + 3] for i in range(0, seq_len, 3)]
        ]
    )


class Template:
    def __init__(self, sequence_nt: str) -> None:
        self.sequence_nt = sequence_nt
        self.sequence_aa = translate(sequence_nt)

    def get_mutant(self, mutations: List[Tuple[str]]) -> str:
        """Method returns the template sequence with the residues at the
        positions specified at initialization replaced by residues listed in
        'mutations.'
        """
        codons_copy = list(self.codons())
        for position, residue in mutations:
            replacement_codon = AA2CODON[residue][0]
            codons_copy[position] = replacement_codon
        return "".join(codons_copy)

    def codons(self):
        return [self.sequence_nt[i : i + 3] for i in range(0, len(self.sequence_nt), 3)]

    def __repr__(self):
        return f"Amino acids: {self.sequence_aa}\nNucleotides: {self.sequence_nt}"


if __name__ == "__main__":
    t = Template("TTTCTTATTCCTCCCATTAATGGTCATGAACTTCTC")
    assert t.sequence_aa == "FLIPPINGHELL"
    assert translate(t.get_mutant([(1, "T")])) == "FTIPPINGHELL"
