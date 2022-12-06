import sys
import pytest

from rna_library import Motif, highest_id
from rna_library import MotifType


class BaseTest(Motif):
    """Test class that inherits from :class:`Motif()` designed to test basic functionality since :class:`Motif()` has abstract methods."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def buffer(self):
        """Overrides :class:`Motif()`'s abstract method to allow instantiation. No other purpose."""
        pass

    def generate_sequences(self):
        """Overrides :class:`Motif()`'s abstract method to allow instantiation. No other purpose."""
        pass

    def has_non_canonical(self):
        """Overrides :class:`Motif()`'s abstract method to allow instantiation. No other purpose."""
        pass

    def recursive_sequence(self):
        """Overrides :class:`Motif()`'s abstract method to allow instantiation. No other purpose."""
        pass

    def recursive_structure(self):
        """Overrides :class:`Motif()`'s abstract method to allow instantiation. No other purpose."""
        pass


def test_const_data():
    """Verify that constant data is correct"""
    tm = BaseTest()

    assert not tm.is_hairpin()
    assert not tm.is_helix()
    assert not tm.is_junction()
    assert not tm.is_singlestrand()
    assert not tm.is_barcode()
    assert not tm.has_parent()
    assert not tm.has_children()

    assert tm.type() == MotifType.UNASSIGNED


def test_getters_and_setters():
    """Checking that getters and setters affect internal state properly."""
    tm = BaseTest()
    tk = "test"
    tm.token(tk)
    assert tm.token() == tk

    ss = "(((...)))"
    tm.structure(ss)

    assert tm.structure() == ss

    seq = "AAAAAA"
    tm.sequence(seq)

    assert tm.sequence() == seq

    new_id = 1
    tm.id(new_id)

    assert tm.id() == new_id

    new_depth = 1
    tm.depth(new_depth)

    assert tm.depth() == new_depth

    other = BaseTest()
    other.token("other-motif")
    other.sequence("AACCC")

    tm.parent(other)

    assert tm.parent() == other


def test_equality_op():
    """Ensure that overload for ``==`` works correctly and considers the token, sequence and type."""
    lhs, rhs = BaseTest(), BaseTest()

    assert lhs == rhs

    lhs.token("Helix2")
    assert lhs != rhs

    rhs.token("Helix2")
    assert lhs == rhs

    lhs.sequence("AAAA")
    assert lhs != rhs

    rhs.sequence("AAAA")
    assert lhs == rhs


def test_get_str():
    """Ensure that ``__str__()`` method functions correctly."""
    tm = BaseTest()
    tm.token("")
    tm.sequence("")
    tm.structure("")

    assert str(tm) == "Unassigned,,"


def test_same_pattern():
    """Checking pattern validation."""
    strands = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10]]
    tm = BaseTest(strands=strands)

    assert tm.same_pattern("NNNN&NNNN&NN")
    assert not tm.same_pattern("NNNN&NN&NN")
    assert not tm.same_pattern("NNNN&NNNNNNN")


def test_highest_id():
    """Check that the recursive ``highest_id()`` method works correctly."""
    tm0 = BaseTest()
    tm0.id(0)

    tm1 = BaseTest()
    tm1.id(1)
    tm0.add_child(tm1)

    tm2 = BaseTest()
    tm2.id(2)
    tm1.add_child(tm2)

    result = highest_id(tm0)

    assert len(tm0.children()) == 1
    assert len(tm1.children()) == 1
    assert not len(tm2.children())

    assert result == 2
