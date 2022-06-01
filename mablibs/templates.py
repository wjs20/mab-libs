from typing import List

from mablibs import codons


class Template:
    def __init__(self, nucleotides) -> None:
        self.nucleotides = nucleotides

    @property
    def amino_acids(self) -> str:
        try:
            return "".join([codons.CODON2AA[codon] for codon in self.codons()])
        except KeyError as e:
            raise NotImplementedError(
                f"Amino acid could not be found, error: {e}"
            )

    def codons(self) -> List[str]:
        return [
            self.nucleotides[i : i + 3]
            for i in range(0, len(self.nucleotides), 3)
        ]

    def __repr__(self):
        return f"Template({self.nucleotides})"
