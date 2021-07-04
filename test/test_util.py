import pytest

from rna_library import util





def test_is_cicrular():
    ic = util.is_circular
    # TODO should probably make the connecvitiy list form a dot-bracket structure
    assert not util.is_circular( 0, [-1, -1] )
    assert not util.is_circular( 0, [-1, 3, -1,  1, -1] )
    assert util.is_circular( 2, [-1, 3, -1,  1, -1] )
