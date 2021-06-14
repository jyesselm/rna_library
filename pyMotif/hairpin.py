import itertools

from .util import * 
from .motif import Motif
from .enums import *

class Hairpin(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.HAIRPIN

        size = len(self.sequence_) - 2
        self.structure_ = "(" + "." * size + ")"
        self.token_ = f"Hairpin{size}"

    def buffer(self):
        return self.parent().buffer()

    def is_hairpin(self):
        return True

    def recursive_structure(self):
        return self.structure_[1:-1]

    def recursive_sequence(self):
        return self.sequence_[1:-1]

    def has_non_canonical(self):
        pair = self.sequence_[0] + self.sequence_[-1]
        return pair not in ALLOWED_PAIRS
    
    def generate_sequences( self ):
        nts = []
        for n in self.recursive_sequence():
            if n != 'N':
                nts.append( [int( NUCLEOTIDE_MAPPER[ n ] )])
            else:
                nts.append( NT_VALS )
        nt_combos = list(itertools.product( *nts ))
        sequences = list(map( nt_codes_to_sequences, nt_combos ))
        self.sequences( ['N'+seq+'N' for seq in sequences] )

