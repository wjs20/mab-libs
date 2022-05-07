import re
from collections import namedtuple 

# amino acid - codon regexprs 
N = 'AA[TC]'
G = 'GG[TCAG]'
S = 'TC[TCAG]|AG[TC]'
A = 'GC[TCAG]'
D = 'GA[TC]'
P = 'CC[TCAG]'
M = 'ATG'
T = 'AC[TCAG]'

# antibody liabilities identified by their codon regexrps
liability_specification = (
        ('GLYCOSYLATION_MOTIF', r'AA[TC](?!CC[TCAG])[ATCG]{3}(TC[TCAG]|AG[TC]|AC[TCAG])'),
        ('DEAMIDATION_MOTIF', r'AA[TC](GG[TCAG]|TC[TCAG]|AG[TC]|AC[TCAG])'),
        ('ISOMERIZATION_MOTIF', r'GA[TC](GG[TCAG]|TC[TCAG]|AG[TC])'),
        ('CLEAVAGE_MOTIF', r'GA[TC]CC[TCAG]'),
        ('OXIDATION_MOTIF', r'ATG')
    )

liability_regex = re.compile('|'.join('(?P<%s>%s)' % pair for pair in liability_specification))

Liability = namedtuple('Liability', 'kind start end')

def find_liabilities(seq, report=False):
    matches = list(liability_regex.finditer(seq))
    if not matches: return 
    if not report: return True
    return [
        Liability(mo.lastgroup, mo.start(), mo.end())
        for mo in matches
    ]

if __name__ == '__main__':
    test_seqs = [
            'ATG'
        ]

    print(find_liabilities(test_seqs[0]))
