from mablibs import codons
from mablibs.ptms import find_ptm_motifs
from mablibs.templates import Template


def get_preferred_codons(species):
    key = codons.CODON_FREQUENCIES[species.upper()].get
    return {aa: max(co, key=key) for aa, co in codons.AA2CODON.items()}


class Mutagenesis:
    def __init__(
        self,
        randomization_strategy,
        template,
        species,
        ptms_to_exclude=None,
        optimizer=None,
        sample_size=None,
    ):
        self.randomization_strategy = randomization_strategy
        self.template = template
        self.aa2codon = get_preferred_codons(species)
        self.ptms_to_exclude = ptms_to_exclude
        self.optimizer = optimizer
        self.sample_size = sample_size

    def mutate(self, template, mutations) -> str:
        codons = template.codons()
        for position, substitution in mutations:
            codons[position] = self.aa2codon[substitution]

        nucleotides = "".join(codons)
        return Template(nucleotides)

    def generate_library(self):
        if self.sample_size is None:
            randomization = self.randomization_strategy.get_mutations()
        else:
            randomization = self.randomization_strategy.sample(self.sample_size)

        for mutations in randomization:
            new_template = self.mutate(self.template, mutations)

            if self.ptms_to_exclude is not None:
                if ptm_motifs := find_ptm_motifs(
                    "".join(new_template.amino_acids)
                ):
                    # skip mutation if it generates a sequence with a ptm site
                    if any(
                        ptm.kind in self.ptms_to_exclude for ptm in ptm_motifs
                    ):
                        continue

            if self.optimizer is None:
                yield new_template

            else:
                yield self.optimizer.optimize_template(new_template)


if __name__ == "__main__":
    pass
