import pytest
# TypeError: Can't instantiate abstract class Motif with abstract methods buffe
#r, generate_sequences, has_non_canonical, recursive_sequence, recursive_structure


from rna_library.motif import Motif


class BaseTest( Motif ):

    def buffer( self ):
        pass
    
    def generate_sequences( self ):
        pass

    def has_non_canonical( self ):
        pass

    def recursive_sequence( self ):
        pass

    def recursive_structure( self ):
        pass



TEST_MOTIF = BaseTest()


def test_base_motif():

    assert not TEST_MOTIF.is_hairpin()
    assert not TEST_MOTIF.is_helix()
    assert not TEST_MOTIF.is_junction()
    assert not TEST_MOTIF.is_singlestrand()
