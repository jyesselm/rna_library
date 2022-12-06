"""This file contains a class description for a SecStruct which corresponds to a non-pseudoknotted RNA
secondary structure. TODO

Author: Chris Jurich <cjurich2@huskers.unl.edu>
Date: 2022-05-17
"""
from __future__ import annotations

from typing import Tuple, List

# TODO(CJ): fix these imports or at least make them prettier
from rna_library.structure.motif import *
from rna_library.structure.hairpin import *
from rna_library.core import *
from .parser import parse_to_motifs


class SecStruct:
    """Representation of a non-pseudoknotted RNA.
    Represents a secondary structure composed of both a dot-bracket struture and a sequene.
        #TODO(CJ): finihs this

        Attributes:
                root_ :
                structure_ :
                sequence_ =
                id_mapping_ :
                total_structures_ :
                counter_ :
    """

    def __init__(self, secstruct: str, sequence: str):
        self.root_ = None
        self.structure_ = ""
        self.sequence_ = ""
        self.id_mapping_ = {}
        self.total_structures_ = 1

        assert len(secstruct) == len(sequence) and len(secstruct) > 0
        assert valid_db(secstruct)
        assert len(re.sub("[ACGUTN]", "", sequence)) == 0
        assert secstruct.count("(") == secstruct.count(")")

        self.structure_ = secstruct
        self.sequence_ = sequence
        self.root_ = parse_to_motifs(self.structure_, self.sequence_)
        # helper for recusion
        self.counter_ = 0
        self.set_ids_(self.root_)

        self.it_ = self.root_
        self.__end_id = highest_id(self.root_)

    def set_ids_(self, m: Motif) -> None:
        """Method that recursively sets ids for the Motif() objects wihtin in the structure. SHOULD NOT
        be called by users directly.

        :param: Motif m: the current motif to perform the operation on
        :rtype: None
        """
        m.id(self.counter_)
        self.id_mapping_[self.counter_] = m

        self.counter_ += 1

        for child in m.children():
            self.set_ids_(child)

    def display(self) -> None:
        """Convenience function that prints the nested Motif() graph to terminal."""
        print(self.root_.str())

    @property
    def sequence(self) -> str:
        """Getter for the underlying sequence.
        :rtype: str
        """
        return self.sequence_

    @property
    def structure(self) -> str:
        """Getter for the underlying secondary structure in dot bracket format..
        :rtype: str
        """
        return self.structure_

    def helix_replace_(self, h_id: int, secstruct: str, sequence: str) -> None:
        """Replaces the sequence and structure for a given Helix() object specified by h_id. SHOULD NOT
        be called by users directly.

        :param: int h_id: The counter_ -based index of the helix.
        :param: str secstruct: A dot bracket representation of the secondary structure replacement for the Helix().
        :param: str sequence: The sequence replacement for the Helix().
        :rtype: None
        """
        helix = self.id_mapping_[h_id]
        helix.structure(secstruct)
        helix.sequence(sequence)

    def motif_replace_(self, m_id: int, new_secstruct: str, new_sequence: str) -> None:
        """Replaces the sequence and structure for a given generic motif() object specified by m_id. SHOULD NOT
        be called by users directly.

        :param: int m_id: The counter_ -based index of the motif.
        :param: str secstruct: A dot bracket representation of the secondary structure replacement for the Motif().
        :param: str sequence: The sequence replacement for the Motif().
        :rtype: None
        """

        def is_hairpin_only(ss: str) -> bool:
            # TODO(CJ): Fix this
            total_len = len(ss)
            return (
                ss[0] == "("
                and ss[-1] == ")"
                and ss[1:-1].count(".") == (total_len - 2)
            )

        motif = self.id_mapping_[m_id]
        if motif.has_parent():
            if is_hairpin_only(new_secstruct):
                new_motif = Hairpin(
                    sequence=new_sequence, strands=[[-1] * len(new_sequence)]
                )
                new_motif.sequence_ = new_sequence
            else:
                new_motif = parse_to_motifs(new_secstruct, new_sequence)
            cleaned_children = []
            for child in motif.parent().children():
                if child.id() == m_id:
                    cleaned_children.append(new_motif)
                else:
                    cleaned_children.append(child)

            # reset the parent of the new motif
            new_motif.parent(motif.parent())
            motif.parent(None)

            new_motif.parent().set_children(cleaned_children)
            refined_seq, refined_ss = (
                self.root_.recursive_sequence(),
                self.root_.recursive_structure(),
            )
            self.root_ = parse_to_motifs(refined_ss, refined_seq)
        else:
            # this is the case where you are actually just replacing the root? seems like a dumb thing to do
            # but who knows what people will try these days
            pass

    def change_motif_sequence(
        self,
        m_id: int,
        new_secstruct: str,
        new_sequence: str,
        include_flanking: bool = False,
    ) -> None:
        # performing some basic checks
        assert len(new_sequence) == len(new_secstruct)
        assert len(re.sub("[AUCGT&]", "", new_sequence)) == 0
        assert len(re.sub("[.)(&]", "", new_secstruct)) == 0
        assert (
            new_secstruct.count("&")
            == new_sequence.count("&")
            # and new_sequence.count("&") <= 1
        )
        assert m_id in self.id_mapping_

        cpy: Motif = self.id_mapping_[m_id]

        if include_flanking and (cpy.is_junction() or cpy.is_hairpin()):
            print(new_sequence)
            exit(0)

        self.id_mapping_[m_id].sequence(new_sequence)

        full_secstruct, full_sequence = (
            self.root_.recursive_structure(),
            self.root_.recursive_sequence(),
        )

        self.structure_ = full_secstruct
        self.sequence_ = full_sequence

        self.root_ = parse_to_motifs(full_secstruct, full_sequence)
        # help for recusion
        self.counter_ = 0
        self.id_mapping_ = dict()
        self.set_ids_(self.root_)
        self.total_structures_ += 1

    def change_motif(
        self,
        m_id: int,
        new_secstruct: str,
        new_sequence: str,
        include_flanking: bool = False,
    ) -> None:
        """Client method for changing a Motif() object within a SecStruct() object. Currently supports the changing of
        a Motif()'s sequence without changing the type/size of the motif as well as changing the Motif() itself and the
        child Motif() objects that it owns. For the case of Junction()/Hairpin() objects, the flanking pairs are not updated
        by default. This behavior can be changed via the include_flanking parameter. It is set to False by default but can be
        set to True to override the outermost respective closing pairs in the adjacent Helix() object(s).


        :param: int m_id : The id of the Motif() to change.
        :param: str new_secstruct: The new secondary structure to use, in dot-bracket format.
        :param: str new_sequence: The new sequence to use.
        :param: bool include_flanking: Whether the flanking pairs of the sequence for a Junction() or Hairpin() object should be overridden on their respective Helix() parent/children.

        :rtype: None
        """
        # performing some basic checks
        # TODO(CJ): This method needs a lot of work in general. Things to fix include:
        # 1. more sophisticated checks than just assert statements
        # 2. breaking the root reset and counting into a different function
        # 3. more sophisticated breakout for changes in motif vs sequence vs tree etcs
        assert len(new_sequence) == len(new_secstruct)
        assert len(re.sub("[AUCGT&]", "", new_sequence)) == 0
        assert len(re.sub("[.)(&]", "", new_secstruct)) == 0
        assert (
            new_secstruct.count("&")
            == new_sequence.count("&")
            # and new_sequence.count("&") <= 1
        )
        assert m_id in self.id_mapping_

        helix_replace = new_secstruct.find("&") != -1

        if helix_replace:
            if include_flanking and not self.id_mapping_[m_id].is_singlestrand():
                self.id_mapping_[m_id].parent_.change_inner_flanking(
                    new_sequence[0] + new_sequence[-1]
                )
                strands = new_sequence.split("&")
                child_cps: List[str] = list()
                for s1, s2, c in zip(
                    strands[:-1], strands[1:], self.id_mapping_[m_id].children()
                ):
                    self.id_mapping_[c.id()].change_outer_flanking(s1[-1] + s2[0])
            self.helix_replace_(m_id, new_secstruct, new_sequence)
        else:
            if include_flanking:
                self.id_mapping_[m_id].parent_.change_inner_flanking(
                    new_sequence[0] + new_sequence[-1]
                )
            self.motif_replace_(m_id, new_secstruct, new_sequence)

        full_secstruct, full_sequence = (
            self.root_.recursive_structure(),
            self.root_.recursive_sequence(),
        )

        self.structure_ = full_secstruct
        self.sequence_ = full_sequence

        self.root_ = parse_to_motifs(full_secstruct, full_sequence)
        # help for recusion
        self.counter_ = 0
        self.id_mapping_ = dict()
        self.set_ids_(self.root_)
        self.total_structures_ += 1

    def get_sequence_structure(self) -> Tuple[str, str]:
        """Getter for the sequence and structure of the SecStruct() object as a Tuple().

        :rtype: Tuple[str,str]
        """
        return (self.sequence_, self.structure_)

    def _get_ids_internal(self, m, ids, mtype):
        if m.type() == mtype:
            ids.append((m.id(), m.token(), m.sequence()))

        for c in m.children():
            self._get_ids_internal(c, ids, mtype)

    def get_ids(self, motif_type):
        ids = []
        for mtype, string in TYPE_MAPPER.items():
            if string == motif_type:
                self._get_ids_internal(self.root_, ids, mtype)
                return ids
        else:
            print(
                f"{motif_type} is an invalid motif type. Allowed values are:\n{', '.join(list(TYPE_MAPPER.values()))}"
            )

    def get_motif(self, id) -> Motif:
        motif = self.id_mapping_[id]
        return motif

    def get_substructure(self, id1, id2=None):
        # get substructure between id1 and id2
        # or get all children of id1
        # make
        pass

    def get(self, m_id: int) -> Motif:
        """Allows for direct accession of the Motif() objects via their int() keys. DOES NOT CHECK FOR KEY CORRECTNESS.

        :param: int m_id: The index target that we are looking for.
        :rtype: Motif
        """
        return self.id_mapping_[m_id]

    def __iter__(self):
        for m in self.id_mapping_.values():
            yield m

    def itermotifs(self):
        for (idx, motif) in self.id_mapping_.items():
            yield (idx, motif)

    def hairpins(self, **kwargs):
        for m in self.id_mapping_.values():
            if m.is_hairpin():
                yield m

    def helix(self):
        for m in self.id_mapping_.values():
            if m.is_helix():
                yield m

    def junctions(self):
        for m in self.id_mapping_.values():
            if m.is_junction():
                yield m

    def singlestrands(self):
        for m in self.id_mapping_.values():
            if m.is_singlestrand():
                yield m

    def set_barcode(self, m_id, bc_seq):
        m: Motif
        m = self.id_mapping_[m_id]
        if not m.is_singlestrand() and not m.is_helix():
            raise TypeError(f"Barcode motif must be either a singlestrand or helix")

        if not m.same_pattern(bc_seq):
            raise TypeError(
                f'The supplied barcode "{bc_seq}" is not compatible with the selected motif'
            )

        m.sequence(bc_seq)
        self.sequence_ = self.root_.recursive_sequence()

    def __add__(self, other: SecStruct) -> SecStruct:
        """
        Operator overload to allow the concatenation of SecStruct objects using `+`
        """
        seq_total = self.sequence + other.sequence
        ss_total = self.structure + other.structure
        # TODO bring over barcode info
        return SecStruct(ss_total, seq_total)
