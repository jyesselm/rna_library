import itertools
from plum import dispatch

from .motif import Motif
from .enums import *
from .util import *


class Helix(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._Motif__type = MotifType.HELIX
        self.__size = (len(self.sequence()) - 1) // 2
        self.__structure = f'{"("*self.__size}&{")"*self.__size}'
        self.__token = f"Helix{self.__size}"
        self.__pairs = []

        for index in range(self.__size):
            self.__pairs.append(self.sequence()[index] + self.sequence()[-index - 1])


    @dispatch
    def size(self):
        return self.__size
    
    @dispatch
    def size(self, val ):
        self.__size = val

    def buffer(self):
        return self.__size

    def pairs(self):
        return self.__pairs

    def is_helix(self):
        return True

    def recursive_structure(self):
        result = self.structure().split("&")

        result = result[0] + self.children()[0].recursive_structure() + result[1]
        if len(self.children()) == 2:
            result += self.children()[1].recursive_structure()

        return result

    def recursive_sequence(self):
        result = self.sequence().split("&")

        result = result[0] + self.children()[0].recursive_sequence() + result[1]
        if len(self.children()) == 2:
            result += self.children()[1].recursive_sequence()

        return result

    def has_non_canonical(self):
        for pair in self.pairs():
            if pair not in ALLOWED_PAIRS:
                return True
        return False

    def generate_sequences( self ):
        bp_codes = []   
        for idx, bp in enumerate( self.pairs() ):
            n_count = bp.count( 'N' )
            if n_count == 0:
                bp_codes.append( [ int(BASEPAIR_MAPPER[ bp ]) ] )
            elif n_count == 2:
                bp_codes.append( BP_VALS )
            elif n_count == 1:
                (left, right) = bp
                allowed = []
                if left != 'N':
                    allowed = list(filter( lambda pr: pr[0] == left, BPS))
                elif right != 'N':
                    allowed = list(filter( lambda pr: pr[1] == right, BPS))
                bp_codes.append( list(map( lambda code: int(BASEPAIR_MAPPER[ code ]), allowed )) )
        nt_combos = list(itertools.product( *bp_codes ))
        sequences = list(map( bp_codes_to_sequence, nt_combos ))
        self.sequences( sequences ) 

