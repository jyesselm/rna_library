from .motif import Motif

from .enums import *

class Junction(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.JUNCTION
        self.closing_pairs_ = []
        self.gaps_ = [len(strand) - 2 for strand in self.strands_]

        self.structure_ = ["."] * len(self.sequence_)
        self.closing_pairs_.append(self.sequence_[0] + self.sequence_[-1])

        for index, nt in enumerate(self.sequence_):
            if nt == "&":
                self.structure_[index - 1] = "("
                self.structure_[index] = "&"
                self.structure_[index + 1] = ")"
                self.closing_pairs_.append(
                    self.sequence_[index - 1] + self.sequence_[index + 1]
                )

        self.structure_[0] = "("
        self.structure_[-1] = ")"

        self.structure_ = "".join(self.structure_)
        self.token_ = f"Junction{len(self.strands())}_" + "|".join(
            [str(len(strand) - 2) for strand in self.strands_]
        )
        self.num_branches_ = len(self.strands())
        self.symmetric_ = len(set([len(strand) - 2 for strand in self.strands_])) == 1

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
        secstruct_chunks = [subseq[1:-1] for subseq in self.structure_.split("&")]
        for secstruct_chunk, child in zip(secstruct_chunks, self.children()):
            result += secstruct_chunk
            result += child.recursive_structure()
        result += secstruct_chunks[-1]
        return result

    def recursive_sequence(self):
        result = str()
        seq_chunks = [subseq[1:-1] for subseq in self.sequence_.split("&")]
        for seq_chunk, child in zip(seq_chunks, self.children()):
            result += seq_chunk
            result += child.recursive_sequence()
        result += seq_chunks[-1]
        return result

    def closing_pairs(self):
        return self.closing_pairs_

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

