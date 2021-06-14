from .motif import Motif

class SingleStrand(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.SINGLESTRAND

        self.structure_ = "." * len(self.sequence_)
        self.token_ = f"SingleStrand{len(self.sequence_)}"

    def buffer(self):
        return -1

    def is_singlestrand(self):
        return True

    def recursive_structure(self):
        result = self.structure_
        for child in self.children():
            result += child.recursive_structure()
        return result

    def recursive_sequence(self):
        result = self.sequence_
        for child in self.children():
            result += child.recursive_sequence()
        return result

    def has_non_canonical(self):
        return False
    
    def generate_sequences( self ):
        nts = []
        for nt in self.sequence_:
            if nt == 'N':
                nts.append( NT_VALS )
            else:
                nts.append( [ int( NUCLEOTIDE_MAPPER[ nt ] )] )
        nt_combos = list(itertools.product( *nts ))
        sequences = list(map( nt_codes_to_sequences, nt_combos ))
        self.sequences( sequences )

