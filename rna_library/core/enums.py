"""Enumerated types for base pairs, nucleotides and motifs. 
Assist with ensuring type correctness and lowers overhead vs string-based
implementations."""

from enum import IntEnum

ALLOWED_PAIRS = {"AU", "UA", "GC", "CG", "UG", "GU"}
"""A :class:`set()` containing all 6 canonical and wobble basepairings."""
BPS = ("GU", "UG", "AU", "UA", "GC", "CG")
"""A :class:`tuple()` of allowed canonical and wobble basepairings.
Ordered for easy conversion by the :class:`BasePair()` class."""
LEGAL_BPS = {"A", "C", "G", "U"}
"""A :class:`set()` of all 4 allowed nucleotide types."""
NTS = ("A", "C", "G", "U")
"""A :class:`tuple()` of all 4 canonical nucleotide types. Ordered 
for each conversion by the :class:`Nucleotide()` class."""


class BasePair(IntEnum):
    """
    Enumerated type for canoncial and wobble basepairs.
    """

    GU = 0
    UG = 1
    AU = 2
    UA = 3
    GC = 4
    CG = 5

    def is_GU(self):
        """
        :return: If the instance is a UG or GU pair.
        :rtype: :class:`bool()`
        """
        return self == BasePair.GU or self == BasePair.UG

    def is_AU(self):
        """
        :return: If the instance is a UA or AU pair.
        :rtype: :class:`bool()`
        """
        return self == BasePair.AU or self == BasePair.UA

    def is_GC(self):
        """
        :return: If the instance is a CG or GC pair.
        :rtype: :class:`bool()`
        """
        return self == BasePair.CG or self == BasePair.GC

    def is_canoncial(self):
        """
        :return: If the instance is a canonical Watson-Crick basepair.
        :rtype: :class:`bool()`
        """

        return self.is_AU() or self.is_GC()

    def to_str(self):
        """
        :return: The :class:`BasePair()` instance in text form.
        :rtype: :class:`str()`
        """
        return BPS[int(self)]


BP_VALS = [int(bp) for bp in BasePair]
"""A :class:`list()` that contains the integer values for all of the 
:class:`BasePair()` enumerations."""

BASEPAIR_MAPPER = {
    "GU": BasePair.GU,
    "UG": BasePair.UG,
    "AU": BasePair.AU,
    "UA": BasePair.UA,
    "GC": BasePair.GC,
    "CG": BasePair.CG,
}

"""A :class:`dict()` object that maps a canonical basepair to its
:class:`BasePair()` value."""


class Nucleotide(IntEnum):
    """Enumerated type for all nucleotide types."""

    A = (0,)
    C = (1,)
    G = (2,)
    U = 3

    def to_str(self):
        """
        :return: The :class:`Nucleotide()` instance in text form.
        :rtype: :class:`str()`
        """
        return NTS[int(self)]


NT_VALS = [int(nt) for nt in Nucleotide]
"""A :class:`list()` that contains the integer values for all of the 
:class:`Nucleotide()` enumerations."""

NUCLEOTIDE_MAPPER = {
    "A": Nucleotide.A,
    "C": Nucleotide.C,
    "G": Nucleotide.G,
    "U": Nucleotide.U,
}
"""A :class:`dict()` object that maps a canoncial nucleotide to its
:class:`Nucleotide()` value."""


class MotifType(IntEnum):
    """Enumerated type for all motif types """

    UNASSIGNED = 0
    SINGLESTRAND = 1
    HELIX = 2
    HAIRPIN = 3
    JUNCTION = 4


TYPE_MAPPER = {  # TODO should make this a class of MotifType
    MotifType.UNASSIGNED: "Unassigned",
    MotifType.SINGLESTRAND: "SingleStrand",
    MotifType.HELIX: "Helix",
    MotifType.HAIRPIN: "Hairpin",
    MotifType.JUNCTION: "Junction",
}
"""A :class: `dict()` object that maps a :class: `MotifType` to its value
as a `str().`"""
