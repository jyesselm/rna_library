import re
from typing import List, Union

from .motif import Motif
from .helix import Helix
from .hairpin import Hairpin
from .junction import Junction
from .singlestrand import SingleStrand

from rna_library.core import connectivity_list, is_circular, valid_db

# wrap these into a class or make private with __func_name
def helix_length(connections: List[int], start: int) -> int:
    """
    TODO
    """
    complement = connections[start]
    length = 0
    while (
        connections[start + length] == complement - length
        and connections[complement - length] == start + length
    ):
        length += 1
    return length


def get_junction_or_hairpin(
    sequence: str, connections: List[int], start: int
) -> Union[Junction, Hairpin]:
    """
    TODO
    """
    strands = []
    offset = 1
    it = start

    while True:
        next = [it]
        it += 1

        while connections[it] == -1:
            next.append(it)
            it += 1

        next.append(it)
        strands.append(next)
        it = connections[it]

        if it == start:
            break
    if len(strands) > 1:
        # JUNCTION
        junction = Junction(strands=strands, sequence=sequence)
        for strand in strands[:-1]:
            junction.add_child(get_motifs(sequence, connections, strand[-1]))
        return junction
    else:
        hp = Hairpin(sequence=sequence, strands=strands)
        return hp


def get_helix(
    sequence: str, connections: List[int], start: int
) -> Union[Hairpin, Helix, Junction, SingleStrand]:
    """
    TODO
    """
    # figure out how big this current helix i
    helix_len = helix_length(connections, start)
    lhs, rhs = [], []
    for index in range(start, start + helix_len):
        lhs.append(index)
        rhs.append(connections[index])
    rhs.reverse()
    helix = Helix(sequence=sequence, strands=[lhs, rhs])
    # figure out what is at the other end of the junction
    if connections[start + helix_len - 1] > start:
        helix.add_child(
            get_junction_or_hairpin(sequence, connections, start + helix_len - 1)
        )

    if not is_circular(rhs[-1], connections):
        helix.add_child(get_motifs(sequence, connections, rhs[-1] + 1))
    return helix


def get_singlestrand(sequence: str, connections: List[int], start: int) -> SingleStrand:
    """
    TODO
    """
    offset = 0
    while start + offset < len(connections) and connections[start + offset] == -1:
        offset += 1

    strand = list(range(start, start + offset))
    singlestrand = SingleStrand(sequence=sequence, strands=[strand])
    # checking here if there is an outer nt
    if start + offset < len(connections):
        singlestrand.children_.append(get_motifs(sequence, connections, start + offset))

    return singlestrand


def get_motifs(
    sequence: str, connections: List[int], start: int
) -> Union[Helix, SingleStrand]:
    """
    Helper method
    """
    if start >= len(connections):
        return []
    elif connections[start] < 0:
        return get_singlestrand(sequence, connections, start)
    return get_helix(sequence, connections, start)


def parse_to_motifs(structure: str, sequence: str) -> Motif:
    """
    Method takes a structure sequence pair and returns a root :class:`Motif()` with a complete associated graph.

    :param str structure: a valid dot-bracket structure 
    :param str sequence: the corresponding sequence composed of the alphabet [ACGUTNB] 

    :return: motif
    :rtype: Motif
    """
    # basic sanity checks
    # TODO should probably deal with pseudoknots in here
    valid_db( structure )
    assert len(structure) == len(sequence)
    assert len(re.sub("[ACGUTNB]", "", sequence)) == 0
    assert structure.count("(") == structure.count(")")
    # TODO gotta fix this so that there can be smaller hairpins
    for ii in range(3):
        invalid = "(" + "." * ii + ")"
        assert structure.find(invalid) == -1
    # generate the connections
    connections = connectivity_list(structure)
    # get the motifs
    root = get_motifs(sequence, connections, 0)
    # reverse link them to the parents
    root.link_children()
    return root
