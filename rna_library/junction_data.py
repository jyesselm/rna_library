"""Contains the JunctionEntry and JunctionData objects"""


class JunctionEntry:
	"""Represents a single junction entry from an RNA construct""" 
    def __init__(self, **kwargs ):
		#sequence, structure, reactivity, construct, sn, reads, score 
        # TODO add some checks here... I think it will make 
        # it better if the checks occur before the actual variable
        # assignment? not sure about that one though 
        assert len( structure ) == len( sequence ) and len( sequence ) == len( reactivity )
        assert len( sequence ) 
        self.__sequence = kwargs.get(sequence, None)
        self.__structure = kwargs.get(structure, None)
        self.__reactivity = kwargs.get(reactivity, None)
        self.__construct = kwargs.get(construct, None)
        self.__sn = kwargs.get(sn, None)
        self.__reads = kwargs.get(reads, None)
        self.__score = kwargs.get(score, None)
        self.__symmetrical = is_symmetrical( self.sequence )

    def key( self ):
        return ( self.sequence, self.structure )

    def is_symmetrical( self ) -> bool:
        return sel.symmetrical

    def __getitem__( self, idx ):
        return self.reactivity[ idx ]

	pass


class JunctionData:
	pass
