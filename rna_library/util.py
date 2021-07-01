"""Utility functions for the rna_library module"""

import re
import sys
import RNA
import math
import pickle
import itertools
import pandas as pd
import editdistance
from .enums import * 
from abc import ABC, abstractmethod
from plum import dispatch

from typing import List, Tuple


def satisfies_constraints( sequence, template ):
    """ This is a test"""    
    for s, t in zip( sequence, template):
        if t == 'N' or t == 'B': continue
        elif t != s:
            return False
    return True


def pool_with_distance( sequences : List[ str ], min_dist : int ) -> List[ str ] :
    """
    Creates a pool of sequences where each sequence has at least a Levenshtein distance between it and all other sequences.
    :param: sequences, a list of starting RNA sequences
    :type: :sequences: list[ str ], required
    :param: min_dist, minimum edit distance between each sequence in the pool
    :type: :min_dist: int, required, must be >= 0
    :rtype: list[ str ]
    """ 
    result = []
    for seq in sorted( sequences ):
        for ii in range(len(result)-1, -1, -1):
            if editdistance.eval( seq, result[ii] ) < min_dist:
                break
        else:
            result.append( seq )

    return result


def bp_codes_to_sequence( bp_code ) -> str:
    """
    Converts a list of :class:`BasePair()`'s into a sequence string.
    :param: bp_code, a list of basepairs to be converted. Basepairs are in order of nesting.
    :type: :bp_code: list[ :class:`BasePair()` ], required
    :rtype: str
    """
    size = len( bp_code )
    left, right = ['N']*size, ['N']*size
    
    for idx, bp in enumerate( bp_code ):
        string = BasePair( bp ).to_str()
        left[ idx ] = string[0]
        right[ size - idx -1 ] = string[1]
    
    return f"{''.join(left)}&{''.join(right)}" 


def nt_codes_to_sequences( codes ) -> str:
    """
    Converts a list of :class:`Nucleotide()`'s into a sequence string.
    :param: codes, a list of nucleotides to be converted. 
    :type: :codes: list[ :class:`Nucleotide()` ], required
    :rtype: str
    """
    result = []
    for c in codes:
        result.append( Nucleotide(c).to_str())
    return ''.join( result )


def get_pair_list( secstruct : str ) -> List[ Tuple[int, int] ]:
    """
    Creates a list of pair indices from a dot-bracket secstruct string. Note 
    that the function assumes the incoming structure is valid.
    :param: secstruct, a dot-bracket structure which is assumed to be valid
    :type: :secstruct: str
    :rtype: list[ tuple( int, int ) ]
    """
    result = [] 
    lparens = [] 
    for ii, db in enumerate( secstruct ):
        if db == '(':
            lparens.append( ii )
        elif db == ')':
            result.append((lparens.pop(), ii ))
    
    if len( lparens ):
        raise TypeError("Unbalanced parentheses in structure")
    
    return result


def connectivity_list( structure : str ) -> List[ int ] :
    """
    Generates a connectivity list or pairmap from a dot-bracket secondary structure.
    :param: structure, a dot-bracket structure
    :type: :structure: str
    :rtype: list[ int ]
    """
    connections, pairs = [-1] * len(structure), []

    for index, db in enumerate(structure):
        if db == "(":
            pairs.append(index)
        elif db == ")":
            complement = pairs.pop()
            connections[complement] = index
            connections[index] = complement
    
    if len( pairs ):
        raise TypeError("Extra pairs left over in structure")

    return connections


def is_circular( start : int, connections : List[ int ] ) -> bool:
    """
    Checks if a starting point in a pairmap is in a circular portion aka a loop.
    :param: start, staring index in the pairmap
    :type: :start: int, required
    :param: connections, pairmap generated from `connectivity_list()`
    :type: :connections:, list[ int ]
    :rtype: bool
    """
    if start != -1:
        return abs( start - connections[ start] ) == 1

    it = start + 1
    while True:
        while it < len(connections) and connections[it] == -1:
            it += 1
        if it == len(connections):
            return False

        it = connections[it] + 1
        if it == start or it < start:
            return True

