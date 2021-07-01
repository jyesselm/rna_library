#!/usr/bin/env python3
"""
Motif classes that serve as the driver for this library's functionality.
"""
import re
import sys
import RNA
import math
import pickle
import itertools
import pandas as pd
import editdistance
from enum import Enum, IntEnum
from .enums import * 
from .util import *
from abc import ABC, abstractmethod
from plum import dispatch


class Motif(ABC):
    """Abstract base class that :class:`Hairpin()`, 
    :class:`Helix()`,
    :class:`Junction()` and
    :class:`SingleStrand()`
    all inherit from.
    """
    def __init__(self, **kwargs):
        """
        Constructor
        """
        self.__type = MotifType.UNASSIGNED
        self.__parent = None
        self.__sequence = str()
        self.__children = []
        self.__strands = []
        self.__depth = int()
        self.__structure = str()
        self.__token = str()
        self.__start_pos = math.inf
        self.__end_pos = -1
        self.__positions = set()
        self.__id = None
        self.__is_barcode = False
        self.__sequences = []

        if "sequence" in kwargs:
            self.__sequence = kwargs["sequence"]

        if "strands" in kwargs:
            self.__strands = kwargs["strands"]

        if len(self.__sequence) == 0:
            return

        strand_seqs = []
        for strand in self.__strands:
            self.start_pos_ = min(min(strand), self.__start_pos)
            self.end_pos_ = max(max(strand), self.__end_pos) 
            for pos in strand:
                self.__positions.add(pos)

            strand_seqs.append("".join([self.__sequence[index] for index in strand]))

        self.__sequence = "&".join(strand_seqs)

    def link_children(self, depth : int = 0):
        """
        Method used to link a :class:`Motif()` object to its children and vice versa. Should only be called once by the root :class:`Motif()`.
        :param depth: depth of the current :class:`Motif()` object. defaults to 0
        :type depth: int, optional
        """            
        if depth < 0:
            raise TypeError(f"Depth supplied to Motf.link_children() must be >= 0")

        if self.type() == MotifType.SINGLESTRAND:
            depth = 0

        self.depth(depth)

        # there might be empty lists in here... remove them if that is the case
        self.__children = [
            child for child in self.children() if not isinstance(child, list)
        ]

        for child in self.children():
            child.parent(self)

        for child in self.children():
            child.link_children(depth + 1)

    def str(self) -> str:
        """
        Creates a recursive string representation of the current :class:`Motif()` object.
        :return: The :class:`Motif()` instance in text form.
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
            return f"{depth}{identification}{self.token_} {self.structure_} {self.sequence_}{children}"

    def __eq__(self, other) -> bool:
        """
        Overloaded `==` operator for :class:`Motif()`. Requires that type of motif, sequence and token are identical.
        :param: other: Another :class:`Motif()` to be compared against.
        :type: other: :class:`Motif()`
        """
        return (
            self.type_ == other.type_
            and self.sequence_ == other.sequence_
            and self.token_ == other.token_
        )

    def __str__(self) -> str:
        """
        String representation of just the motif at hand
        :return: The :class:`str()` representation of the :class:`Motif()`.
        :rtype: :class:`str()`
        """
        return f"{ TYPE_MAPPER[ self.__type ] },{ self.__sequence },{ self.__structure }"

    def is_helix(self):
        return False

    def is_singlestrand(self):
        return False

    def is_hairpin(self):
        return False

    def is_junction(self):
        return False

    def type(self):
        return self.__type

    def children(self):
        return self.__children

    def add_child(self, other):
        self.__children.append( other )

    def set_children(self, other):
        self.__children = other

    def parent(self, other=None):
        if other is not None:
            self.__parent = other
        else:
            return self.__parent

    def token(self, tk=None):
        if tk is None:
            return self.token_
        else:
            self.token_ = tk

    def structure(self, secstruct=None):
        if secstruct is None:
            return self.structure_
        else:
            self.structure_ = secstruct

    def strands(self):
        return self.__strands


    @dispatch
    def sequence(self):
        return self.__sequence
    
    @dispatch
    def sequence(self, seq ):
        self.__sequence = seq

    @dispatch
    def id(self):
        return self.id_

    @dispatch
    def id(self, new_id):
        self.id_ = new_id

    @dispatch
    def depth(self):
        return self.depth_

    @dispatch
    def depth(self, value):
        self.depth_ = value

    @abstractmethod
    def buffer(self):
        pass

    def has_children(self):
        return len(self.children_) > 0

    def has_parent(self):
        return self.parent_ is not None

    @abstractmethod
    def recursive_sequence(self):
        pass

    @abstractmethod
    def recursive_structure(self):
        pass

    @abstractmethod
    def has_non_canonical(self):
        pass

    def same_pattern(self, sequence):
        template = '&'.join(['N'*len( s ) for s in self.strands_]) 
        if len( sequence ) != len( template ):
            return False

        for s, t in zip( sequence, template ):
            if (s == '&' and t != '&') or (s != '&' and t == '&'):
                return False
        return True 

    def start_pos(self):
        return self.start_pos_

    def end_pos(self):
        return self.end_pos_

    def contains(self, pos):
        return pos in self.positions_
    

    # this is for making barcodes
    def sequences( self, seqs ):
        self.__sequences = seqs

    def number_sequences( self ):
        return len( self.__sequences )

    def set_sequence( self, idx ):
        self.sequence_ = self.__sequences[ idx ]
    
    @abstractmethod
    def generate_sequences( self ):
        pass


def traverse(motif):
    assert not motif.has_non_canonical()
    for c in motif.children():
        # print('chilren')
        # print('\t',c)
        traverse(c)

def highest_id( m : Motif, best=0 ):
    best = max( best, m.id() )

    for c in m.children():
        best = highest_id(c, best)

    return best








if __name__ == "__main__":
    ss =  "..(((...)))......"
    seq = "NNGAGAAACUCNNNNNN"
    #seq = "NNNNNUGAAACANN"
    barcodes = build_barcodes( ss, seq )
    pickle.dump( barcodes, open('barcodes-2.p', 'wb'))
    print(len(barcodes))
    exit( 0 ) 
    d1 = SecStruct( '(((...)))', 'GGGAAACCC')
    d2 = SecStruct( '....', 'AAAA')
    print( (d1 + d2).sequence )
    print( (d1 + d2).structure )

