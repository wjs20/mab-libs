import itertools as it
import string
import operator
import re
import functools
import warnings

MAX_SEQUENCES_IN_MEM = 10_000

class Randomization:
    def __init__(self, position_residue_mapping):
        self.position_residue_mapping = position_residue_mapping
    
    def _get_position_residue_pairs(self):
        return it.chain.from_iterable(
            self.zip_keys_to_vals(self.position_residue_mapping)
        )
    
    def _create_combinations(self):
        return 
    
    def get_mutations(self, generator=True, first=None):
        assert self._create_combinations() is not None, 'Check you have overloaded _create_combinations method'
        mutations = self._create_combinations()
        
        if not generator:
            if self.library_size <= MAX_SEQUENCES_IN_MEM:
                return list(it.islice(mutations, 0, first))
            elif self.library_size >= MAX_SEQUENCES_IN_MEM:
                warnings.warn(f"""
                Library size of {self.library_size} may exhuast RAM.
                first 10,000 sequences will be returned.
                """)
                return list(it.islice(mutations, 0, first or MAX_SEQUENCES_IN_MEM))
        else:
            return mutations
    
    @staticmethod
    def zip_keys_to_vals(position_residue_mapping):
        return [
            list(it.zip_longest([k], v, fillvalue=k))
            for k, v in position_residue_mapping.items()
        ]

class TwomerRandomization(Randomization):
    
    def _non_redundant_pair(self, ms):
        if ms[0][0] != ms[1][0]: 
            return True

    def _filter_redundant_mutations(self, mutatations):
        return filter(self._non_redundant_pair, mutatations)
    
    def _create_combinations(self):
        return self._filter_redundant_mutations(
            it.combinations(self._get_position_residue_pairs(), 2)
        )
    
    @functools.cached_property
    def library_size(self):
        return sum(
            len(position_one_residues) * len(position_two_residues) 
            for (_, position_one_residues), (_, position_two_residues) 
            in it.combinations(self.position_residue_mapping.items(), 2)
        )

class TrimerRandomization(Randomization):
    def _group_by_position(self, mutations):
        position_getter = operator.itemgetter(0)
        return [
            list(g) 
            for _,g in it.groupby(mutations, key=position_getter)
        ]
    
    def _create_combinations(self):
        return it.product(
            *self._group_by_position(self._get_position_residue_pairs())
        )
    
    @functools.cached_property
    def library_size(self):
        n_res_by_position = [
            len(res) 
            for res in self.position_residue_mapping.values()
        ]
        
        return functools.reduce(
            operator.mul, n_res_by_position
        )

if __name__ == '__main__':

    n_positions = 7
    position_residue_mapping = dict(zip(range(n_positions), [string.ascii_letters[:18]]*n_positions))

    twomer = TwomerRandomization(position_residue_mapping)
    assert twomer.library_size == 6804
    assert twomer.get_mutations(generator=False, first=1) == [((0, 'a'), (1, 'a'))]

    trimer = TrimerRandomization(position_residue_mapping)
    assert trimer.library_size == 612220032
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert trimer.get_mutations(generator=False, first=1) == [((0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'a'), (5, 'a'), (6, 'a'))]
