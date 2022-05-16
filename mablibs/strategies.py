import itertools as it
import operator
import re
import functools
import random
import string
import warnings
from typing import Dict, Iterator, Union, Tuple, List, Sequence, Any

Mutation = Tuple[int, str]


class RandomizationStrategy:
    """Top level class containing methods for generating the mutations
    associated with a particular randomization strategy.

    This class is meant to be subclassed to create different Randomization strategies.
    Any subclass of RandomizationStrategy must overload the _create_combinations
    method and the library_size property.
    """

    position_getter = operator.itemgetter(0)

    def __init__(self, position_residue_mapping: Dict) -> None:
        """
        Args:
            position_residue_mapping (dict): A dictionary mapping positions to
            amino acid substitutions. e.g.

                {
                    1: ['A, 'D' ..., 'V', 'W'],
                    3: ['A, 'D' ..., 'V', 'W']
                }
        """
        self.position_residue_mapping = position_residue_mapping

    def _get_position_residue_pairs(self) -> Iterator[Mutation]:
        """
        Method to convert position residue mapping dict into iterator of mutation
        tuples (position, residue) e.g.

            {
                1: ['A, 'D' ..., 'V', 'W'],
                3: ['A, 'D' ..., 'V', 'W']
            }

            ->

            (1,'A'), (1,'D'), ..., (3, 'V'), (3, 'W')
        """
        return it.chain.from_iterable(
            self.zip_keys_to_vals(self.position_residue_mapping)
        )

    def _create_combinations(self) -> None:
        """
        Overwrite method to specify how mutations should be combined to create
        the library.
        """
        return

    def get_mutations(
        self, generator: bool = True, first: Union[int, None] = None
    ) -> Iterator[Tuple[Mutation]]:
        """
        Method overrides_create_combinations method to generate an iterator
        of mutations.

        Args:
            generator (bool): if True, result will be returned as a generator.
            This prevents large lists of mutations being generated in memory.

            first (Union[int, None]): if None generator or list containing all mutations will be returned.
            if an int is passed to first, the mutations up to 'first' will be returned.

        Returns:
            Union[Iterator[Tuple[Mutation]], list[Tuple[Mutation]]: an iterator
            or list containing a mutations e.g.

            -> ((1, 'A'), (3, 'A')), ((1, 'A'), (3, 'D')) ... ((1, 'W'), (3, 'W'))

            each tuple represents all mutations that will be made to a sequence,
            mapping a set of positions to a set of amino acids.

        """
        assert (
            self._create_combinations() is not None
        ), "Check you have overloaded _create_combinations method"

        mutations = self._create_combinations()

        if first is not None:
            mutations = it.islice(mutations, 0, first)

        # Use for debugging/inspecting, may cause memory issues if library size is large
        if not generator:
            return list(mutations)

        else:
            return mutations

    def sample(
        self, k: int, seed: Union[int, None] = None
    ) -> Iterator[Tuple[Mutation]]:
        """
        Method takes a random sample of size k mutations from the total set of
        mutations. Useful for investigating properties of very large libraries.

        Args:
            k (int): sample size

        Returns:
            Iterator[Tuple[Mutation]]: a sample of mutations
        """
        mutations = self.get_mutations()
        library_idxs = range(0, self.library_size)

        if seed:
            random.seed(42)

        sample_indices = set(random.choices(library_idxs, k=k))
        for i, mutation in enumerate(mutations):
            if i in sample_indices:
                yield mutation

    @functools.cached_property
    def library_size() -> None:
        """
        Override method to determine size of the library
        """
        return

    @staticmethod
    def zip_keys_to_vals(
        position_residue_mapping: Dict[int, List[str]]
    ) -> List[List[Mutation]]:
        """
        Method converts a dict to a list such that items in each value list are
        each paired with their key e.g.


        {
            1: ['A, 'D' ..., 'V', 'W'],
            3: ['A, 'D' ..., 'V', 'W']
        }

        ->

        [
            [(1, 'A'), (1, 'D'), ..., (1, 'V'), (1, 'W')],
            [(3, 'A'), (3, 'D'), ..., (3, 'V'), (3, 'W')]
        ]

        Used to create the input for _create_combinations class method.
        """
        return [
            list(it.zip_longest([position], residues, fillvalue=position))
            for position, residues in position_residue_mapping.items()
        ]

    @staticmethod
    def _mul_lens(iterable: Sequence[Sequence[Any]]) -> int:
        """
        Helper method to multiply together the lengths of a collection of
        iterables

        Args:
            iterator (Sequence[Sequence[Any]]): an iterable of iterables e.g.
            [(1, 2, 3), (1, 2, 3)]

        Returns:
            int: the product of the lengths of all the iterables

        """
        return functools.reduce(operator.mul, [len(i) for i in iterable])


class NmerRandomizationStrategy(RandomizationStrategy):
    """
    Nmer randomization strategy will only mutatate n positions and any one time
    out of all positions to be mutated.
    """

    def __init__(self, n: int, *args) -> None:
        """
        Args:
            n (int): number of mutations to include in each combination e.g.
            n = 2 -> every unique combination of 2 mutations.
        """
        self.n = n
        super().__init__(*args)

    def _non_redundant_pair(self, mutations: List[Tuple[Mutation]]) -> bool:
        """Method checks to see if any mutations in a set of mutations to be
        applied to a single sequence will mutate the same position i.e. are redundant.

        Args:
            ms (List[Tuple[Mutation]]): a set of mutations to be applied to a single
            sequence e.g.

        Returns:
            bool: wether the mutation set is not redundant e.g.
            ((1, 'A'), (1, 'D')) -> False (mutate the same position)
            ((1, 'A'), (3, 'D')) -> True (mutate the different position)

        """
        positions = list(map(self.position_getter, mutations))

        if len(positions) == len(set(positions)):
            return True
        else:
            return False

    def _filter_redundant_mutations(
        self, mutatations: Iterator[Tuple[Mutation]]
    ) -> Iterator[Tuple[Mutation]]:
        """
        Method applies non redudant filter to a collection of mutations to ensure
        the same position is not being mutated twice in the same sequence.

        Args:
            mutations (Iterator[Tuple[Mutation]]): an iterator generating
            a collection of mutations e.g. ((1, 'A'), (3, 'A')), ...

        Returns:
            (Iterator[Tuple[Mutation]]): filtered iterator of mutations
        """
        return filter(self._non_redundant_pair, mutatations)

    def _create_combinations(self) -> Iterator[Tuple[Mutation]]:
        """
        Method returns an Iterator of all unique n length combinations of mutations.
        """
        return self._filter_redundant_mutations(
            it.combinations(self._get_position_residue_pairs(), self.n)
        )

    @functools.cached_property
    def library_size(self) -> int:
        """
        Property calculates the total library size for the specified randomization strategy.

        Returns (int):
            library size (number of unique sequences)
        """
        return sum(
            self._mul_lens(residues)
            for residues in it.combinations(
                self.position_residue_mapping.values(), self.n
            )
        )


class RandomTrimerRandomizationStrategy(RandomizationStrategy):
    """
    Randomization strategy in which all positions are mutated simultaneously.
    """

    def _group_by_position(self, mutations: Iterator[Mutation]) -> List[List[Mutation]]:
        """
        Method groups mutations by position.

        Args:
            mutations (Iterator[Mutation]): an iterator containing mutations e.g.

            (1, 'A'), (1, 'D'), (2, 'A'), (2, 'D')

        Returns:
            List[List[Mutation]]: [(1, 'A'), (1, 'D')], [(2, 'A'), (2, 'D')]
        """

        return [list(g) for _, g in it.groupby(mutations, key=self.position_getter)]

    def _create_combinations(self):
        """
        Method returns an Iterator of all unique n length combinations of mutations.
        """
        return it.product(*self._group_by_position(self._get_position_residue_pairs()))

    @functools.cached_property
    def library_size(self):
        return self._mul_lens(self.position_residue_mapping.values())
