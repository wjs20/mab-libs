import json
import operator
import pkgutil
from typing import List, Tuple
from dataclasses import dataclass
from mablibs import codons

@dataclass
class Template:
    nucleotides: str

    def codons(self) -> List[str]:
        return [self.nucleotides[i : i + 3] for i in range(0, len(self.nucleotides), 3)]

    @property
    def amino_acids(self) -> str:
        try:
            return "".join([codons.CODON2AA[codon] for codon in self.codons()])
        except KeyError as e:
            raise NotImplementedError(f'Amino acid could not be found, error: {e}')
