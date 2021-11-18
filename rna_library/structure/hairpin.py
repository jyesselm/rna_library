import itertools

from rna_library.core.util import *
from rna_library.core.enums import *
from rna_library.structure.motif import Motif


class Hairpin(Motif):
    """Represents a hairpin loop in an RNA structure. Inherits from :class:`Motif()`."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.HAIRPIN

        size = len(self.sequence()) - 2
        self.structure_ = "(" + "." * size + ")"
        self.token_ = f"Hairpin{size}"

    def buffer(self) -> int:
        """
        For the :class:`Hairpin()` type, this is simply the size of the 
        closing helix meaning the number of closing pairs.

        :return: buffer
        :rtype: int
        """
        return self.parent().buffer()

    def is_hairpin(self) -> bool:
        """
        Indicates that the :class:`Motif()` is of type :class:`Hairpin()`.

        :return: is_hairpin
        :rtype: bool
        """
        return True

    def recursive_structure(self) -> str:
        """
        Returns the owned portion of the structure. In this coding of structure 
        it is just the loop portion and does not include the closing pair.

        :return: recursive_structure
        :rtype: str
        """
        return self.sequence()[1:-1]

    def recursive_sequence(self) -> str:
        """
        Returns the owned portion of the sequence. In this coding of sequence 
        it is just the loop portion and does not include the closing pair.

        :return: recursive_sequence
        :rtype: str
        """
        return self.sequence()[1:-1]

    def has_non_canonical(self) -> bool:
        """
        Returns whether or not the closing pair is canonical (i.e. is AU/UA, CG/GC, GU/UG).

        :return: has_non_canonical
        :rtype: bool
        """
        seq = self.sequence()
        pair = seq[0] + seq[-1]
        return pair not in ALLOWED_PAIRS

    def generate_sequences(self):
        """
        Generates all possible sequences for the :class:`Hairpin()` that are compatible with
        the constraints for the motif.
        """
        nts = []
        for n in self.recursive_sequence():
            if n != "N":
                nts.append([int(NUCLEOTIDE_MAPPER[n])])
            else:
                nts.append(NT_VALS)
        nt_combos = list(itertools.product(*nts))
        sequences = list(map(nt_codes_to_sequences, nt_combos))
        self.sequences(["N" + seq + "N" for seq in sequences])
