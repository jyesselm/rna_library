#!/usr/bin/env python3
"""
Motif classes that serve as the driver for this library's functionality.
"""
# this is so that the typehinting can reference itself
# TODO figure out which of these methods fail on after an internal change
from __future__ import annotations

import math
from rna_library.core.util import *
from abc import ABC, abstractmethod
from plum import dispatch


class Motif(ABC):
    """
    Abstract base class that :class:`Hairpin()`, :class:`Helix()`, :class:`Junction()` and :class:`SingleStrand()` all inherit from.
    
    """

    def __init__(self, **kwargs):
        """
        Constructor for ``Motif``.
        
        .. warning:: Should **NOT** be called directly. All instantiations are handled by ``rna_library.parser.parse_to_motifs()``

        """
        self.type_ = MotifType.UNASSIGNED
        self.parent_ = None
        self.sequence_ = str()
        self.children_ = []
        self.strands_ = []
        self.depth_ = int()
        self.structure_ = str()
        self.token_ = str()
        self.start_pos_ = math.inf
        self.end_pos_ = -1
        self.positions_ = set()
        self.id_ = None
        self.is_barcode_ = False
        # these are used for making barcodes... should probably be changed
        self.sequences_ = []

        if "sequence" in kwargs:
            self.sequence_ = kwargs["sequence"]

        if "strands" in kwargs:
            self.strands_ = kwargs["strands"]

        if len(self.sequence_) == 0:
            return

        strand_seqs = []
        for strand in self.strands_:
            self.start_pos_ = min(min(strand), self.start_pos_)
            self.end_pos_ = max(max(strand), self.end_pos_)
            for pos in strand:
                self.positions_.add(pos)

            strand_seqs.append("".join([self.sequence_[index] for index in strand]))

        self.sequence_ = "&".join(strand_seqs)

    def link_children(self, depth: int = 0):
        """
        Method used to link a :class:`Motif()` object to its children and vice versa. Should only be called once by the root :class:`Motif()`.

        :param int depth: depth of the current :class:`Motif()` object. defaults to 0
        :rtype: None
        """
        if depth < 0:
            raise TypeError(f"Depth supplied to Motf.link_children() must be >= 0")

        if self.type() == MotifType.SINGLESTRAND:
            depth = 0

        self.depth(depth)

        # there might be empty lists in here... remove them if that is the case
        self.children_ = [
            child for child in self.children() if not isinstance(child, list)
        ]

        for child in self.children():
            child.parent(self)

        for child in self.children():
            child.link_children(depth + 1)

    def str(self) -> str:
        """
        Creates a recursive string representation of the current :class:`Motif()` object.
        
        :rtype: :class:`str()`
        """
        if self.id_ is not None:
            identification = f"ID: {self.id_}, "
        else:
            identification = ""

        depth = "\t" * self.depth()
        if not self.has_children():
            return f"{depth}{identification}{self.token_} {self.structure_} {self.sequence_}"
        else:
            contents = [""]
            for index, child in enumerate(self.children()):
                contents.append(child.str())

                if (
                    index < len(self.children()) - 1
                    and self.children()[index + 1].is_singlestrand()
                ):
                    tks = contents[-1].splitlines()[0]
                    first_line = contents[-1].splitlines()[0]
                    length = len(first_line) - len(first_line.lstrip())

            children = "\n".join(contents)
            return f"{depth}{identification}{self.token()} {self.structure()} {self.sequence()}{children}"

    def __eq__(self, other: Motif) -> bool:
        """
        Overloaded ``==`` operator for :class:`Motif()`. Requires that type of motif, sequence and token are identical.
        
        :param Motif other: Another :class:`Motif()` to be compared against.
        
        """
        return (
            self.type() == other.type()
            and self.sequence() == other.sequence()
            and self.token() == other.token()
        )

    def __str__(self) -> str:
        """
        String representation of just the motif at hand.
        
        :return: The :class:`str()` representation of the :class:`Motif()`.
        :rtype: :class:`str()`
        """
        return f"{ TYPE_MAPPER[ self.type_ ] },{ self.sequence_ },{ self.structure_ }"

    def is_helix(self) -> bool:
        """
        If the motif is a helix or not. Overridden by child :class:`Helix()` class.

        :return: If the motif is of type :class:`Helix()`
        :rtype: :class:`bool()`
        """
        return False

    def is_singlestrand(self) -> bool:
        """
        If the motif is a singlestrand or not. Overridden by child :class:`SingleStrand()` class.

        :return: If the motif is of type :class:`SingleStrand()`
        :rtype: :class:`bool()`
        """
        return False

    def is_hairpin(self) -> bool:
        """
        If the motif is a hairpin or not. Overridden by child :class:`Hairpin()` class.

        :return: If the motif is of type :class:`Hairpin()`
        """
        return False

    def is_junction(self) -> bool:
        """
        If the motif is a junction or not. Overridden by child :class:`Junction()` class.

        :return: If the motif is of type :class:`Junction()`
        :rtype: :class:`bool()`
        """
        return False

    def type(self) -> MotifType:
        """
        Returns the :class:`MotifType()` type for the given motif.

        :return: The :class:`MotifType()` enum value for the given motif.
        :rtype: :class:`MotifType()`
        """
        return self.type_

    def children(self) -> List[Motif]:
        """
        Getter for the :class:`Motif()`'s child motifs. Returned as a list for iteration. Only returns direct children or an empty list if the motif has not children.
        
        :return: A :class:`list()` of :class:`Motif()` if the current :class:`Motif()` has any.
        :rtype: :class:`list[Motif]`
        """
        return self.children_

    def add_child(self, other: Motif) -> None:
        """
        Appends a new :class:`Motif()` to the internal list of children for the current :class:`Motif()`.

        .. warning:: Should **NOT** be called directly. Other function calls must occur to ensure that the internal graph is accurate.
        
        :param: Motif other: Another :class:`Motif()` to be appended to the internal children list.
        """
        self.children_.append(other)

    def set_children(self, other: List[Motif]) -> None:
        """
        Sets the entire list of `Motif()` to the internal list of children for the current :class:`Motif()`.

        .. warning:: Should **NOT** be called directly. Other function calls must occur to ensure that the internal graph is accurate.
        
        :param  List[Motif] other: Another :class:`Motif()` to be appended to the internal children list.
        """
        self.children_ = other

    @dispatch
    def parent(self, other):
        """ 
        Sets the :class:`Motif()`'s parent to the supplied :class:`Motif()`. 

        :param Motif other: The new parent for the current :class:`Motif()`.
        
        :return: None
        :rtype: NoneType
        """
        self.parent_ = other

    @dispatch
    def parent(self):
        """ 
        Gets the parent :class:`Motif()`'s for the current :class:`Motif()`. 
        
        :return: the parent motif 
        :rtype: :class:`Motif()`
        """
        return self.parent_

    @dispatch
    def token(self, tk):
        """
        Sets the :class:`Motif()`'s identifying token to an inputted string. Input is **NOT** validated.

        :param str tk: the new token for the :class:`Motif()`.
        
        :return: None
        :rtype: NoneType
        """
        # TODO add some kind of validation
        self.token_ = tk

    @dispatch
    def token(self):
        """
        Gets the identifying token for the :class:`Motif()`.

        :return: token
        :rtype: str
        """
        return self.token_

    @dispatch
    def structure(self, secstruct):
        """
        Sets the :class:`Motif()`'s structure to an inputted string. Input is **NOT** validated.

        :param str tk: the new structure for the :class:`Motif()`.
        
        :return: None
        :rtype: NoneType
        """
        # TODO add some kind of validation... maybe the level of validation
        # should be checked as well
        self.structure_ = secstruct

    @dispatch
    def structure(self):
        """
        Gets the secondary structure for the :class:`Motif()`.

        :return: token
        :rtype: str
        """
        return self.structure_

    def strands(self) -> List[List[int]]:
        """
        Returns a list of list of :class:`int()`'ss where each sub list contains a contiguous set of nucleotides that "belong" to the :class:`Motif()`.
        Output varies by motif type and the expected values are below:

        - :class:`Hairpin()` => 1
        - :class:`Helix()` => 2
        - :class:`SingleStrand()` => 1
        - :class:`Junction()` => number of branches in :class:`Junction()`

        :return: strands
        :rtype: List[List[int]]
        """
        return self.strands_

    @dispatch
    def sequence(self):
        """
        Gets the sequence for the :class:`Motif()`. 
        Because the nucleotides owned by the :class:`Motif()` may not be contiguous, breaks will 
        be separated by an ampersand '&'. 

        :return: sequence
        :rtype: str
        """
        return self.sequence_

    @dispatch
    def sequence(self, seq):
        """
        Sets the sequence for the :class:`Motif()` to the supplied string. Warning the input **NOT** validated.
        
        :param str seq: the new sequence for the :class:`Motif()`. 
        """
        # TODO some kind of validation
        self.sequence_ = seq

    @dispatch
    def id(self):
        """
        Gets the id :class:`int` value for the given :class:`Motif()`.
        
        :return: id
        :rtype: int
        """
        return self.id_

    @dispatch
    def id(self, new_id):
        """
        Sets the id for the :class:`Motif()`. Warning: It is **NOT** currently validated.
        
        :param int new_id: the new id for the :class:`Motif()`
        :return: none
        :rtype: NoneType
        """
        # TODO some kind of validation
        self.id_ = new_id

    @dispatch
    def depth(self):
        """
        The depth of the :class:`Motif()`, which describes how deep it is in the internal graph.

        :return: depth
        :rtype: int
        """
        return self.depth_

    @dispatch
    def depth(self, value):
        """
        Sets the depth of the current :class:`Motif()`. 

        :param int value: the new depth value for the current :class:`Motif()`.

        """
        # TODO some kind of validation
        self.depth_ = value

    @abstractmethod
    def buffer(self) -> int:
        """
        Buffer refers to the size of the closest adjacent :motif:`Helix()`.
        Varies by type of motif as seen below:

        - :class:`Helix()` => size of the helix itself
        - :class:`Hairpin()` =>  size of its parent helix
        - :class:`SingleStrand()` => -1, meaningless in this context
        - :class:`Junction()` => a :class:`list()` of the branching helices' length with the parent helix going first the in the direction of increasing nucleotide index.

        :return: buffer
        :rtype: int
        """
        pass

    def has_children(self) -> bool:
        """
        Returns whether the :class:`Motif()` has any children.

        :return:  has_children
        :rtype: bool 
        """
        return len(self.children_) > 0

    def has_parent(self) -> bool:
        """
        Returns whether the :class:`Motif()` has a parent.

        :return: has_parent
        :rtype: bool
        """
        return self.parent_ is not None

    @abstractmethod
    def recursive_sequence(self) -> str:
        """
        Builds and returns the continguous sequence of the structure viewing the current
        :class:`Motif()` as the root of the structure. The returned sequence will be part of 
        the main sequence.

        :return: sequence
        :rtype: str
        """
        pass

    @abstractmethod
    def recursive_structure(self) -> str:
        """
        Builds and returns the continguous structure of the structure viewing the current
        :class:`Motif()` as the root of the structure. The returned structure will be part of 
        the main structure.

        :return: structure
        :rtype: str
        """
        pass

    @abstractmethod
    def has_non_canonical(self) -> bool:
        """
        Checks if the :class:`Motif()` has any non-canonical (i.e. non AU/UA, UG/GU or GC/CG) pairings.

        :return: has_nc
        :rtype: bool
        """
        pass

    def same_pattern(self, sequence: str) -> bool:
        """
        Checks if a template sequence is compatible with an inputted sequence. Specifically if the length
        and placement of '&' are the same.

        :param str sequence: template string to compare against.
        
        :return: is_same
        :rtype: bool
        """
        template = "&".join(["N" * len(s) for s in self.strands()])
        if len(sequence) != len(template):
            return False

        for s, t in zip(sequence, template):
            if (s == "&" and t != "&") or (s != "&" and t == "&"):
                return False
        return True

    def start_pos(self) -> int:
        """
        Starting (or lowest) nucleotide index owned by the :class:`Motif()`.

        :return: start_pos
        :rtype: int
        """
        return self.start_pos_

    def end_pos(self) -> int:
        """
        Ending (or highest) nucleotide index owned by the :class:`Motif()`.

        :return: end_pos
        :rtype: int
        """
        return self.end_pos_

    def contains(self, pos: int) -> bool:
        """
        Indicates if a nucleotide index is contained or belongs to the current :class:`Motif()`.
        
        :param list[int] pos: the querying index

        :return: is_contained
        :rtype: bool
        """
        return pos in self.positions_

    # this is for making barcodes
    def sequences(self, seqs: List[str]) -> None:
        """
        Used to set the internal list of barcode temp sequences.

        :param List[str] seqs: the new barcode sequences to be applied to the current :class:`Motif()`.

        """
        # TODO some kind of validation
        self.sequences_ = seqs

    def number_sequences(self) -> int:
        """
        Gives the number of barcode sequences that the :class:`Motif()` currently has.

        :return: num_sequence
        :rtype: int
        """
        return len(self.sequences_)

    def set_sequence(
        self, idx: int
    ) -> None:  # TODO maybe this should be a dispatch overload?
        """
        Sets the current sequence to the sequence of the existing index from the internal barcodes list.
        Note that the `Motif.number_sequences()` method should be queried prior so that the index call will
        be known to be valid.

        :param int idx: The index to be used. 
        """
        self.sequence_ = self.sequences_[idx]

    @abstractmethod
    def generate_sequences(self):
        """
        Builds out all possible barcode sequences that fit the known constraints.
        """
        pass

    def is_barcode(self) -> bool:
        """
        Returns whether the current :class:`Motif()` serves as a barcode.

        :return: is_barcode
        :rtype: bool
        """
        return self.is_barcode_


def highest_id(m: Motif, best: int = 0) -> int:
    """
    Figures out the highest id number in a given :class:`Motif()` graph.

    :param Motif m: motif to start the query on
    :param int best: current highest or "best" motif id at that recursion level. 

    :return: highest_id
    :rtype: int
    """
    best = max(best, m.id())

    for c in m.children():
        best = highest_id(c, best)

    return best
