"""Utility functions for the ``rna_library`` module"""

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


def satisfies_constraints( sequence : str, template : str ) -> bool:
    """
    Confirms whether a sequence and template are the same or not. A template has one of 6
    possible characters at each position: one of the normal A/C/G/U, N for "any" or "B" to indicate
    the position is meant to be part of a barcode. Will return ``False`` if sequence and template
    are not the same length.

    :param str sequence: sequence in question
    :param str template: templated sequence
    :rtype: bool
    """    
    if len( sequence ) != len( template ):
        return False

    for s, t in zip( sequence, template):
        if t == 'N' or t == 'B': 
            continue
        elif t != s:
            return False
    return True


def pool_with_distance( sequences : List[ str ], min_dist : int ) -> List[ str ] :
    """
    Creates a pool of sequences where each sequence has at least the specified Levenshtein distance between it and all other sequences.
    Method sorts sequences internally so input order is not relevant to final pool.
    
    .. warning:: This function can runs in polynomial time so large pools **WILL** take a significant amount of time to run. For reference, pools on the order of hundreds of thousands took multiple hours to run on an i7 in 2021.

    :param: list[str] sequences: A list of starting RNA sequences.
    :param: int min_dist: Minimum edit distance between each sequence in the pool. Must be >= 0.
    :rtype: list[str]
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
    
    :param list[BasePair] bp_code: a list of basepairs to be converted. Basepairs are in order of nesting.
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
    
    :param list[Nucleotide] codes: a list of nucleotides to be converted. 
    :rtype: str
    """
    resulu = []
    for c in codes:
        result.append( Nucleotide(c).to_str())
    return ''.join( result )


def get_pair_list( secstruct : str ) -> List[ Tuple[int, int] ]:
    """
    Creates a list of pairs of indices from a dot-bracket secstruct string. Note 
    that the function assumes the incoming structure is valid. 
    
    :param str secstruct: a dot-bracket structure which is assumed to be valid
    :rtype: list[tuple(int,int)]
    :raises TypeError: if the number of left parentheses exceeds the number of right parentheses
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
    The list has a value of ``-1`` for unpaired positions else has the index of a 
    positions complement.
    
    :param str structure: a dot-bracket structure 
    :rtype: list[int]
    :raises TypeError: if the number of left parentheses exceeds the number of right parentheses
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
        raise TypeError("Unbalanced parentheses in structure")

    return connections


def is_circular( start : int, connections : List[ int ] ) -> bool:
    """
    Checks if a starting point in a pairmap is in a circular portion.
    This can include the closing pairs of both hairpins and junctions.

    :param int start: staring index in the pairmap
    :param list[int] connections: pairmap generated from ``util.connectivity_list()``
    :rtype: bool
    """
    if connections[ start ] != -1:
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
