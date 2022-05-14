import randomization
from templates import Template
import string
import itertools as it
from time import time
import random

def main():
    r = randomization.RandomTrimerRandomization({0: 'A C D E'.split(), 3: 'A C D E'.split(), 6: 'A C D E'.split()})
    t = Template('CCCTTTTATAAAGGTGGGAGA')
    for m in r.get_mutations():
        print(t.translate(t.get_mutant(m)))
 
if __name__ == '__main__':
