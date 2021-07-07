import pytest

from rna_library.helix import Helix
from rna_library.enums import MotifType


def test_constant_data():
    """Ensure that the constant state data is accurate."""
    th = Helix()
    assert th.is_helix()
    
    assert not th.is_singlestrand()
    assert not th.is_junction()
    assert not th.is_hairpin()
    assert th.type() == MotifType.HELIX

def test_size_buffer_methods():
    """Making sure that the buffer/size methods function properly"""
    th = Helix(sequence='GGAA&UUCC', strands=[[0, 1, 2, 3 ], [5, 6, 7, 8 ]])

    assert th.size() == 4
    assert th.buffer() == 4

def test_pair_methods():
    """Checking that pairs are iterated through properly."""
    
    EXPECTED = ['GC', 'GC', 'AU', 'AU']

    th = Helix(sequence='GGAA&UUCC', strands=[[0, 1, 2, 3 ], [5, 6, 7, 8 ]])

    assert th.pairs() == EXPECTED

    assert not th.has_non_canonical()


    for pr in ['A&U', 'U&A', 'G&C', 'C&G', 'U&G', 'G&U']:
        th = Helix(sequence=pr, strands=[[0], [2]])
        assert not th.has_non_canonical()


    for pr in ['A&A', 'A&G', 'A&C', 'C&C', 'C&A', 'C&U', 'G&G', 'G&A', 'U&U', 'U&C' ]:
        th = Helix(sequence=pr, strands=[[0], [2]])
        assert th.has_non_canonical()


