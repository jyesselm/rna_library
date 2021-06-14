import re
from typing import List

from .motif import Motif
from .helix import Helix
from .hairpin import Hairpin
from .junction import Junction 
from .singlestrand import SingleStrand 

from .util import connectivity_list, is_circular

# wrap these into a class or make private with __func_name
def helix_length(connections, start):
    complement = connections[start]
    length = 0
    while (
        connections[start + length] == complement - length
        and connections[complement - length] == start + length
    ):
        length += 1
    return length


def get_junction_or_hairpin(sequence, connections, start):

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
            junction.children_.append(get_motifs(sequence, connections, strand[-1]))
        return junction
    else:
        hp = Hairpin(sequence=sequence, strands=strands)
        return hp


def get_helix(sequence, connections, start):
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
        helix.children_.append(
            get_junction_or_hairpin(sequence, connections, start + helix_len - 1)
        )

    if not is_circular(rhs[-1], connections):
        helix.children_.append(get_motifs(sequence, connections, rhs[-1] + 1))
    return helix


def get_singlestrand(sequence, connections, start):
    offset = 0
    while start + offset < len(connections) and connections[start + offset] == -1:
        offset += 1

    strand = list(range(start, start + offset))
    singlestrand = SingleStrand(sequence=sequence, strands=[strand])
    # checking here if there is an outer nt
    if start + offset < len(connections):
        singlestrand.children_.append(
            get_motifs(sequence, connections, start + offset)
        )

    return singlestrand

def get_motifs(sequence : str, connections: List[ int ], start : int):
    if start >= len(connections):
        return []
    elif connections[start] < 0:
        return get_singlestrand(sequence, connections, start)
    return get_helix(sequence, connections, start)


def parse_to_motifs(structure : str, sequence : str) -> Motif:
    # basic sanity checks
    assert len(structure) == len(sequence)
    assert len(re.sub("[\(\.\)]", "", structure)) == 0
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
    root.link_children_(0)
    return root
