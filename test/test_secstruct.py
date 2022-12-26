from rna_library.structure.secstruct import SecStruct


def test_simple():
    """
    test simple secstruct construction
    """
    seq = "GGGGAAAACCCC"
    dot_bracket = "((((....))))"
    struct = SecStruct(dot_bracket, seq)
    assert struct.get_num_motifs() == 2

def test_replace_motif():
    """
    test replace motif
    """
    seq = "GGGGAAAACCCC"
    dot_bracket = "((((....))))"
    struct = SecStruct(dot_bracket, seq)
    struct.change_motif(1, "(....)", "CUUCGG", include_flanking=True)
    assert struct.sequence == "GGGCUUCGGCCC"
