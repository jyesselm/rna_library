import re
import itertools
from typing import Union, List

from .util import * 
from .motif import Motif
from .enums import LEGAL_BPS
from .parser import parse_to_motifs



def generate_sequences( m : Motif ):
    m.generate_sequences()
    for c in m.children():
        generate_sequences( c )

def get_idxs( m : Motif, idxs : List[ int ] ):
    idxs.append( m.number_sequences() )
    for c in m.children():
        get_idxs( c, idxs )

def make_sequences( m : Motif ):
    # 1. make the idx's
    idxs = [] 
    get_idxs( m, idxs )
    idxs = list(map(lambda n: list(range(n)), idxs))
    # 2. get the permutations
    permutations = list(itertools.product( *idxs ))
    seqs = [] 
    for perm in permutations:
        set_sequence( m, perm )
        seqs.append( m.recursive_sequence() )
    return list( set( seqs ))

def set_sequence( m : Motif, perm, ii=0):
    m.set_sequence(perm[ ii ])
    for c in m.children():
        ii += 1
        set_sequence( c, perm, ii )

def longest_repeat( nums ):
    if len( nums ) == 0: # maybe throw an error here... shouldn't be getting this
        return 0
    count, best, curr = 0, 0, math.nan
    
    for n in nums:
        if n != curr:
            best = max( best, count )
            count = 1
            curr = n
        else:
            count += 1
    
    return max( best, count )




def contains_jct( m : Motif ):
    count = m.is_junction()
    
    for c in m.children():
        count |= c.is_junction()
   
    return count


def check_pairs( sequence, indices ):
    for pair_idx in indices:
        pair = sequence[pair_idx[0]] + sequence[pair_idx[1]]
        
        if 'N' in pair or 'B' in pair:
            continue
        
        if pair not in LEGAL_BPS:
            raise TypeError(f"{pair} is not a legal basepair. Only {', '.join(BPS)} are allowed")


def validate_barcode_constraints( secstruct : str, sequence : str ):
    # must be the same length
    assert len(secstruct) == len(sequence)
    # valid characters only    
    assert len( re.findall('[^(.)]', secstruct) ) == 0
    assert len( re.findall('[^AUCGNB]', sequence) ) == 0
    # check that the basepairs are good
    pairs = get_pair_list( secstruct )
    check_pairs( sequence, pairs )


def build_barcodes( secstruct : str, start : Union[str, None] = None, distance: int =3 ):
    #TODO this code works generally well right now, but it just treats everything
    # as a graph of Motifs. It would be better to use a SecStruct object
    # and then to also make some of the above methods part of the SecStruct object
    # As a third thing, it needs to check for emptry strings and know how to deal with them

    if start is None or len(start) == 0:
        start = 'N'*len( secstruct )
    validate_barcode_constraints( secstruct, start )
    
    m : Motif
    m = parse_to_motifs( secstruct, start )
    if contains_jct( m ):
        raise TypeError(f"the build_barcodes() method cannot build barcodes containing junctions")
    generate_sequences( m )
    sequences = make_sequences( m )
    if start.count('N') != len( start ):
        filt = list(filter(lambda s: satisfies_constraints(s, start), sequences))
    else:
        filt = sequences
    # repeats  
    filt = list(filter( lambda seq: longest_repeat( seq ) < 4, filt )) 
    # get the edit distance right
    filt = pool_with_distance( filt, distance )
    folding_info = list(map( lambda seq: (seq, *RNA.fold( seq )), filt))
    # must fold correctly 
    filt = list(filter( lambda entry: entry[1] == secstruct, folding_info))
    # order by mfe
    filt = sorted( filt, key=lambda entry: entry[-1] )
    # get just the sequence
    result = list(map( lambda entry: entry[0], filt ))
    return result

