import itertools
from typing import List
from plum import dispatch

from rna_library.structure.motif import Motif
from rna_library.core.enums import *
from rna_library.core.util import *


class Helix(Motif):
    """
    Represents a helix or stack in an RNA structure. Inherits from :class:`Motif()`.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.HELIX
        self.size_ = (len(self.sequence()) - 1) // 2
        self.structure_ = f'{"("*self.size_}&{")"*self.size_}'
        self.token_ = f"Helix{self.size_}"
        self.pairs_ = []

        for index in range(self.size_):
            self.pairs_.append(self.sequence()[index] + self.sequence()[-index - 1])

    @dispatch
    def size(self) -> int:
        """
        Returns the size of the :class:`Helix()` which is just the number of 
        pairs in the stack.

        :return: size
        :rtype: int
        """
        return self.size_

    @dispatch
    def size(self, val):
        """
        Sets the current size for the :class:`Helix()`. 

        :param int val: the new size of the helix.

        """
        # TODO some kind of validation
        self.size_ = val

    def buffer(self):
        """
        Returns the buffer of the :class:`Helix()` which is just the number of 
        pairs in the stack.

        :return: buffer
        :rtype: int
        """
        return self.size_

    def pairs(self) -> List[str]:
        """
        Returns the basepairs in the stack as a list of strings of length 2.
        Pairs are returned in order of lowest 3 prime starting index. 
        
        :return: pairs
        :rtype: List[str]
        """
        return self.pairs_

    def is_helix(self):
        """
        Indicates that the :class:`Motif()` is of type :class:`Helix()`.

        :return: is_helix
        :rtype: bool
        """
        return True

    def recursive_structure(self):
        """
        Builds and returns the continguous sequence of the structure viewing the current
        :class:`Motif()` as the root of the structure. The returned sequence will be part of 
        the main sequence.

        :return: sequence
        :rtype: str
        """
        result = self.structure().split("&")

        result = result[0] + self.children()[0].recursive_structure() + result[1]
        if len(self.children()) == 2:
            result += self.children()[1].recursive_structure()

        return result

    def recursive_sequence(self):
        """
        Builds and returns the continguous structure of the structure viewing the current
        :class:`Motif()` as the root of the structure. The returned structure will be part of 
        the main structure.

        :return: structure
        :rtype: str
        """
        result = self.sequence().split("&")

        result = result[0] + self.children()[0].recursive_sequence() + result[1]
        if len(self.children()) == 2:
            result += self.children()[1].recursive_sequence()

        return result

    def has_non_canonical(self):
        """
        Checks if any of the basepairs are non-canonical (i.e. non- AU/UA, GU/UG, GC/CG).

        :return: has_non_canonical
        :rtype: bool
        """
        # TODO should reference ALLOWED_PAIRS
        for pair in self.pairs():
            if pair not in ALLOWED_PAIRS:
                return True
        return False

    def generate_sequences(self):
        """
        Generates all possible sequences for the :class:`Helix()` that are compatible with
        the constraints for the motif.
        """

        bp_codes = []
        for idx, bp in enumerate(self.pairs()):
            n_count = bp.count("N")
            if n_count == 0:
                bp_codes.append([int(BASEPAIR_MAPPER[bp])])
            elif n_count == 2:
                bp_codes.append(BP_VALS)
            elif n_count == 1:
                (left, right) = bp
                allowed = []
                if left != "N":
                    allowed = list(filter(lambda pr: pr[0] == left, BPS))
                elif right != "N":
                    allowed = list(filter(lambda pr: pr[1] == right, BPS))
                bp_codes.append(
                    list(map(lambda code: int(BASEPAIR_MAPPER[code]), allowed))
                )
        nt_combos = list(itertools.product(*bp_codes))
        sequences = list(map(bp_codes_to_sequence, nt_combos))
        self.sequences(sequences)

    def change_outer_flanking( self, new_cp : str ) -> None:
        """Changes the outer closing pair to the supplied new_cp. new_cp must be a string of length 2.
		Consider the below example:
		>>> h = Helix(sequence='12&34',structure='((&))'))
		>>> h.sequence()
		'12&34'
		>>> h.change_outer_flanking('AB')
		>>> h.sequence()
		'A2&3B'
		"""
        assert len(new_cp) == 2
        tks : List[str] = list(self.sequence())
        tks[0], tks[-1] = new_cp[0], new_cp[1] 
        self.sequence_ = ''.join( tks )

    def change_inner_flanking( self, new_cp : str ) -> None:
        """Changes the inner closing pair to the supplied new_cp. new_cp must be a string of length 2.
		Consider the below example:
		>>> h = Helix(sequence='12&34',structure='((&))'))
		>>> h.sequence()
		'12&34'
		>>> h.change_inner_flanking('AB')
		>>> h.sequence()
		'1A&B4'
		"""
        assert len(new_cp) == 2
        it : int  = self.sequence().find('&')
        tks : List[str] = list(self.sequence())
        tks[it-1], tks[it+1] = new_cp[0], new_cp[1] 
        self.sequence_ = ''.join( tks )
