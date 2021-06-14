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
from multipledispatch import dispatch


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
        self.__is_barcode = False
        self.__sequences = []

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

    def link_children_(self, depth : int):
            
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
            child.link_children_(depth + 1)

    def str(self):
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
                    # contents.append('<-' + '-'*len(contents[-1].splitlines()[0]) )
            # children = '\n'.join([""] + [child.str() for child in self.children()])
            children = "\n".join(contents)
            return f"{depth}{identification}{self.token_} {self.structure_} {self.sequence_}{children}"

    def __eq__(self, other) -> bool:
        return (
            self.type_ == other.type_
            and self.sequence_ == other.sequence_
            and self.token_ == other.token_
        )

    def __str__(self) -> str:
        return TYPE_MAPPER[ self.type_ ] + "," + self.sequence_ + "," + self.structure_

    def is_helix(self):
        return False

    def is_singlestrand(self):
        return False

    def is_hairpin(self):
        return False

    def is_junction(self):
        return False

    def type(self):
        return self.type_

    def children(self):
        return self.children_

    def set_children(self, other):
        self.children_ = other

    def parent(self, other=None):
        if other is not None:
            self.parent_ = other
        else:
            return self.parent_

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
        return self.strands_

    def sequence(self, seq=None):
        if seq is None:
            return self.sequence_
        else:
            self.sequence_ = seq

    @dispatch()
    def id(self):
        return self.id_

    @dispatch(int)
    def id(self, new_id):
        self.id_ = new_id

    @dispatch()
    def depth(self):
        return self.depth_

    @dispatch(int)
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

