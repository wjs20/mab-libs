import functools
import itertools as it
import random
import re
from typing import Callable, Dict, List

from mablibs import codons, enzymes
from mablibs.templates import Template


def compile_restriction_site_regex(*restriction_sites):
    return re.compile(
        "|".join(
            f"(?P<{name}>{pattern})"
            for name, pattern in enzymes.ENZYMES.items()
            if name in restriction_sites
        )
    )


def compile_nucleotide_repeat_regex():
    single_repeat = r"([ATCG])\1{3,}"
    double_repeat = "|".join(
        [
            f'({"".join(dinucleotide)})' + "{4,}"
            for dinucleotide in it.product("ATCG", repeat=2)
        ]
    )
    return re.compile("|".join((single_repeat, double_repeat)))


def pattern_not_found(seq, patterns=None):
    if patterns is None:
        raise NotImplementedError('"patterns" must me provided')

    if patterns.search(seq) is None:
        return True
    else:
        return False


def is_not_palindromic(seq):
    half_seq_len = len(seq) // 2
    front, back = seq[:half_seq_len], seq[-half_seq_len:]
    return front != "".join(
        codons.COMPLEMENTARY_BASES[base] for base in back[::-1]
    )


def is_below_gc_content_threshold(seq, threshold=0.65):
    return (seq.count("G") + seq.count("C")) / len(seq) <= threshold


@functools.lru_cache(maxsize=None)
def get_synonymous_codons(
    codon: str, synonymous_codons: List[List[str]]
) -> List[str]:
    for codons in synonymous_codons:
        if codon in codons:
            return codons
    else:
        raise NotImplementedError("codon not in synonymous codons")


class DNAOptimizer:
    def __init__(
        self,
        constraints: List[Callable],
        species: str,
        seed=None,
    ) -> None:
        self.constraints = constraints
        self.codon_ref = tuple(codons.AA2CODON.values())
        self.codon_frequencies = codons.CODON_FREQUENCIES[species.upper()]
        self.seed = seed

    def change_codon(self, i: int, template: Template) -> Template:
        codons = template.codons()
        codon_to_replace = codons[i]
        synonymous_codons = get_synonymous_codons(
            codon_to_replace, self.codon_ref
        )
        weights = [self.codon_frequencies[codon] for codon in synonymous_codons]
        if self.seed:
            random.seed(self.seed)
        replacement_codon = random.choices(
            synonymous_codons, weights=weights, k=1
        )
        codons[i] = replacement_codon[0]  # codon comes wrapped in list
        nucleotides = "".join(codons)
        return Template(nucleotides)

    def optimize_template(self, template: Template) -> str:
        if all(
            constraint(template.nucleotides) for constraint in self.constraints
        ):
            return template

        n_codons = len(template.codons())
        codon_indices = range(n_codons)
        for i in it.cycle(codon_indices):
            template = self.change_codon(i, template)
            if all(
                constraint(template.nucleotides)
                for constraint in self.constraints
            ):
                return template
