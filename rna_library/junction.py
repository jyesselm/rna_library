from .motif import Motif

from .enums import *

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
        self.num_branches_ = len(self.strands())
        self.symmetric_ = len(set([len(strand) - 2 for strand in self.strands()])) == 1

    def buffer(self):
        buffers = [self.parent.buffer()]
        for child in self.children():
            buffers.append(child.buffer())
        return buffers

    def gaps(self):
        return self.gaps_

    def is_junction(self):
        return True

    def recursive_structure(self):
        result = str()
        secstruct_chunks = [subseq[1:-1] for subseq in self.structure().split("&")]
        for secstruct_chunk, child in zip(secstruct_chunks, self.children()):
            result += secstruct_chunk
            result += child.recursive_structure()
        result += secstruct_chunks[-1]
        return result

    def recursive_sequence(self):
        result = str()
        seq_chunks = [subseq[1:-1] for subseq in self.sequence().split("&")]
        for seq_chunk, child in zip(seq_chunks, self.children()):
            result += seq_chunk
            result += child.recursive_sequence()
        result += seq_chunks[-1]
        return result

    def closing_pairs(self):
        return self.__closing_pairs

    def has_non_canonical(self):
        for pair in self.closing_pairs_:
            if pair not in ALLOWED_PAIRS:
                return True
        return False

    def number_branches(self):
        return self.num_branches_

    def symmetric(self):
        return self.symmetric_
    
    def generate_sequences( self ):
        raise TypeError(f"The method generate_sequences() is not supported for the Junction type")

