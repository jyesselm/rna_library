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
from multipledispatch import dispatch




def satisfies_constraints( sequence, template ):
    """ This is a test"""    
    for s, t in zip( sequence, template):
        if t == 'N' or t == 'B': continue
        elif t != s:
            return False
    return True


def pool_with_distance( sequences, min_dist ):
    
    result = []
    
    for seq in sorted( sequences ):
        for ii in range(len(result)-1, -1, -1):
            if editdistance.eval( seq, result[ii] ) < min_dist:
                break
        else:
            result.append( seq )

    return result


def bp_codes_to_sequence( bp_code ):
    size = len( bp_code )
    left, right = ['N']*size, ['N']*size
    
    for idx, bp in enumerate( bp_code ):
        string = BasePair( bp ).to_str()
        left[ idx ] = string[0]
        right[ size - idx -1 ] = string[1]
    
    return f"{''.join(left)}&{''.join(right)}" 


def nt_codes_to_sequences( codes ):
    result = []
    for c in codes:
        result.append( Nucleotide(c).to_str())
    return ''.join( result )


def get_pair_list( secstruct : str ): # note assumes that the incoming structure has been validated
    result = [] 
    lparens = [] 
    for ii, db in enumerate( secstruct ):
        if db == '(':
            lparens.append( ii )
        elif db == ')':
            result.append((lparens.pop(), ii ))
    assert len( lparens ) == 0
    return result


def connectivity_list(structure):
    connections, pairs = [-1] * len(structure), []

    for index, db in enumerate(structure):
        if db == "(":
            pairs.append(index)
        elif db == ")":
            complement = pairs.pop()
            connections[complement] = index
            connections[index] = complement

    assert len(pairs) == 0

    return connections


def is_circular(start, connections):
    it = start + 1
    while True:
        while it < len(connections) and connections[it] == -1:
            it += 1
        if it == len(connections):
            return False

        it = connections[it] + 1
        if it == start or it < start:
            return True

