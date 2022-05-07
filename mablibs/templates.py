import json
import operator
import pkgutil
from typing import List, Tuple

CODON_TABLE = json.loads(pkgutil.get_data('mablibs', 'codon_table.json').decode('utf-8'))

class Template:
    def __init__(self, sequence_nt: str) -> None:
        self.sequence_nt = sequence_nt
        self.sequence_aa = self.translate(sequence_nt)

    def get_mutant(self, mutations: List[Tuple[str]]) -> str:
        """Method returns the template sequence with the residues at the 
        positions specified at initialization replaced by residues listed in 
        'mutations.'
        """
        codons_copy = list(self.codons()) 
        for position, residue in mutations:
            replacement_codon = CODON_TABLE[residue][0] 
            codons_copy[position] = replacement_codon
        return ''.join(codons_copy)

    def codons(self):
        return [self.sequence_nt[i:i+3] for i in range(0, len(self.sequence_nt), 3)]

    @staticmethod
    def translate(sequence_nt: str) -> str:
        """Method converts nucleotide sequences to amino acid sequences"""
        seq_len = len(sequence_nt)
        assert seq_len%3 == 0, 'Ensure sequence contains complete reading frame'
        sequence_aa = []

        for codon_start in range(0, seq_len, 3):
            current_codon = sequence_nt[codon_start:codon_start + 3]
            for aa, codon_list in CODON_TABLE.items():
                if current_codon in codon_list:
                    sequence_aa.append(aa)
                    break
            else:
                raise ValueError('Could not find aa, check sequence is valid')
        return ''.join(sequence_aa)

if __name__ == '__main__':
    t = Template('TTTCTTATTCCTCCCATTAATGGTCATGAACTTCTC', [0, 1, 2])
    assert t.sequence_aa == 'FLIPPINGHELL'
    assert t.translate(t.mutate([(1, 'T')])) == 'FTIPPINGHELL'
