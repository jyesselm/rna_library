import pytest

from rna_library import util





def test_is_cicrular():
    ic = util.is_circular

    assert not ic( 0, [-1, -1] )
    assert not ic( 0, [-1, 3, -1,  1, -1] )
    assert ic( 2, [-1, 3, -1,  1, -1] )
