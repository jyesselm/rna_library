import itertools

from .util import * 
from .motif import Motif
from .enums import *

class Hairpin(Motif):
    """Represents a hairpin loop in an RNA structure. Inherits from :class:`Motif()`."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._Motif__type = MotifType.HAIRPIN

        size = len(self.sequence()) - 2
        self.__structure = "(" + "." * size + ")"
        self.__token = f"Hairpin{size}"

    def buffer(self) -> int:
        return self.parent().buffer()

    def is_hairpin(self) -> bool:
        return True

    def recursive_structure(self) -> str:
        return self.sequence()[1:-1]

    def recursive_sequence(self) -> str :
        return self.sequence()[1:-1]

    def has_non_canonical(self) -> bool:
        pair = self.__sequence[0] + self.__sequence[-1]
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

