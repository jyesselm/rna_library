from rna_library.core import *
from .motif import Motif


class SingleStrand(Motif):
    """
    Represents a single stranded region in an RNA structure. Does not include unpaired regions that are part of a :class:`Junction()` or :class:`Helix()`.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.SINGLESTRAND

        self.structure_ = "." * len(self.sequence())
        self.token_ = f"SingleStrand{len(self.sequence())}"

    def buffer(self) -> int:
        """
        For the :class:`SingleStrand()` type, this does not have any meaning and is always the value 
        ``-1``.
        
        :return: buffer
        :rtype: int
        """
        return -1

    def is_singlestrand(self) -> bool:
        """
        Indicates that the :class:`Motif()` is of type :class:`SingleStrand()`.

        :return: is_singlestrand 
        :rtype: bool
        """
        return True

    def recursive_structure(self) -> str:
        """
        Returns the owned portion of the structure. In this coding of structure 
        it is just the nucleotides in the single strand plus its child if it exists.

        :return: recursive_structure
        :rtype: str
        """
        result = self.structure()
        for child in self.children():
            result += child.recursive_structure()
        return result

    def recursive_sequence(self) -> str:
        """
        Returns the owned portion of the sequence. In this coding of sequence 
        it is just the nucleotides in the single strand plus its child if it exists.

        :return: recursive_sequence
        :rtype: str
        """
        result = self.sequence()
        for child in self.children():
            result += child.recursive_sequence()
        return result

    def has_non_canonical(self):
        """
        Because there are no pairs "owned" by :class:`SingleStrand()`'s, it always returns ``False``.
        
        :return: has_non_canonical
        :rtype: bool
        """
        return False

    def generate_sequences(self):
        """
        Generates all possible sequences for the :class:`SingleStrand()` that are compatible with
        the constraints for the motif.
        """
        nts = []
        for nt in self.sequence():
            if nt == "N":
                nts.append(NT_VALS)
            else:
                nts.append([int(NUCLEOTIDE_MAPPER[nt])])
        nt_combos = list(itertools.product(*nts))
        sequences = list(map(nt_codes_to_sequences, nt_combos))
        self.sequences(sequences)
