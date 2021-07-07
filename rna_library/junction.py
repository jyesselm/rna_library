from .motif import Motif

from .enums import *

from typing import List

class Junction(Motif):
    """
    Represents a junction of any size in an RNA structure including bulges and multi-loops.    
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._Motif__type = MotifType.JUNCTION
        self.__closing_pairs = []
        self.__gaps = [len(strand) - 2 for strand in self.strands()]

        secstruct = ["."] * len(self.sequence())
        self.__closing_pairs.append(self.sequence()[0] + self.sequence()[-1])

        for index, nt in enumerate(self.sequence()):
            if nt == "&":
                secstruct[index - 1] = "("
                secstruct[index] = "&"
                secstruct[index + 1] = ")"
                self.__closing_pairs.append(
                    self.sequence()[index - 1] + self.sequence()[index + 1]
                )

        secstruct[0] = "("
        secstruct[-1] = ")"

        self.structure( "".join( secstruct ) )
        self.__token = f"Junction{len(self.strands())}_" + "|".join(
            [str(len(strand) - 2) for strand in self.strands()]
        ) 
        self.__num_branches = len(self.strands())
        self.__symmetric = len(set([len(strand) - 2 for strand in self.strands()])) == 1

    def buffer(self) -> List[int]:
        """
        For the :class:`Junction()` type this is a :class:`list()` of :class:`int()`'s where the first
        is the size of the parent :class:`Helix()` and then they are arranged in 3' to 5' order.
        Will have the same size as number of branches in the :class:`Jucntion()`.

        :return: buffers
        :rtype: List[int]
        """
        buffers = [self.parent.buffer()]
        for child in self.children():
            buffers.append(child.buffer())
        return buffers

    def gaps(self) -> List[int]:
        """
        Returns a :class:`list()` of :class:`int()`'s of gap sizes in 3' to 5' order.
        Will have the same size as number of branches in the :class:`Jucntion()`.
        
        :return: gaps
        :rtype: List[int]
        """
        return self.__gaps

    def is_junction(self) -> bool:
        """
        Indicates that the :class:`Motif()` is of type :class:`Junction()`.

        :return: is_hairpin
        :rtype: bool
        """
        return True

    def recursive_structure(self) -> str:
        """
        Returns the owned portion of the structure. In this coding of structure 
        it is the closing pairs as well as the child :class:`Helix()`'s and their children.

        :return: recursive_structure
        :rtype: str
        """
        result = str()
        secstruct_chunks = [subseq[1:-1] for subseq in self.structure().split("&")]
        for secstruct_chunk, child in zip(secstruct_chunks, self.children()):
            result += secstruct_chunk
            result += child.recursive_structure()
        result += secstruct_chunks[-1]
        return result

    def recursive_sequence(self) -> str:
        """
        Returns the owned portion of the sequence. In this coding of structure 
        it is the closing pairs as well as the child :class:`Helix()`'s and their children.

        :return: recursive_sequence
        :rtype: str
        """
        result = str()
        seq_chunks = [subseq[1:-1] for subseq in self.sequence().split("&")]
        for seq_chunk, child in zip(seq_chunks, self.children()):
            result += seq_chunk
            result += child.recursive_sequence()
        result += seq_chunks[-1]
        return result

    def closing_pairs(self) -> List[str]:
        """
        Returns a :class:`list()` of :class:`str()`'s that correspond to the closing pairs in 
        the :class:`Junction()` Motif. 

        :return: closing_pairs
        :rtype: List[str]
        """
        return self.__closing_pairs

    def has_non_canonical(self) -> bool:
        """
        Returns whether or not any of the closing pairs are non-canonical (i.e. not AU/UA, CG/GC, GU/UG).

        :return: has_non_canonical
        :rtype: bool
        """
        for pair in self.closing_pairs_:
            if pair not in ALLOWED_PAIRS:
                return True
        return False

    def number_branches(self) -> int:
        """
        Returns the number of branches in the current :class:`Junction()`.

        :return: number_branches
        :rtype: int
        """
        return self.__num_branches

    def symmetric(self) -> bool:
        """
        Indicates if the current :class:`Junction()` is symmetric, that is the sizes of all of the 
        gaps are the same.

        :return: is_symmetric
        :rtype: bool
        """
        return self.__symmetric
    
    def generate_sequences( self ):
        """
        Would generate all possible sequences for the :class:`Junction()` that are compatible with
        the constraints for the motif. **Not currently implemented.**

        :raises: TypeError
        """
        raise TypeError(f"The method generate_sequences() is not supported for the Junction type")

