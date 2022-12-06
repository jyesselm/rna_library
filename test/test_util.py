import pytest

from rna_library.core import util


def test_is_cicrular():
    ic = util.is_circular
    # TODO should probably make the connecvitiy list form a dot-bracket structure
    assert not util.is_circular(0, [-1, -1])
    assert not util.is_circular(0, [-1, 3, -1, 1, -1])
    assert util.is_circular(2, [-1, 3, -1, 1, -1])


def tests_is_symmetrical():
    # good cases
    assert util.is_symmetrical("N&N")
    assert util.is_symmetrical(".&.")
    assert util.is_symmetrical("((((&))))")
    # bad cases
    assert not util.is_symmetrical("N&")
    assert not util.is_symmetrical(".&..")
    assert not util.is_symmetrical("(((&))))")
    # assert raises when no ampersand is supplied
    with pytest.raises(Exception):
        util.is_symmetrical("")
    with pytest.raises(Exception):
        util.is_symmetrical("NNN")
    with pytest.raises(Exception):
        util.is_symmetrical("&&")
