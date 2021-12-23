"""
Holds utility functions for folding RNA sequences.
"""
import os
import pickle
from collections import namedtuple
from typing import List, Tuple
from .util import safe_rm

FoldResult = namedtuple("FoldResult", "seq ss mfe ed params")
"""namedtuple that holds sequence, structure, ensembled defect and folding parameters for an RNAfold prediction."""

_CACHE_FILE = f"{os.path.expandvars('$HOME')}/.vienna_cache"
"""Value where the cached fold results will be stored. Set to $HOME/.vienna_cache in current system."""

_DEFAULT_PARAMS = ("-p", "--noLP", "-d2")
"""Default folding parameters for Vienna's RNAfold. Default values are ("-p", "--noLP", "-d2")"""

_CACHE = None
"""Cache that holds folding results."""

_CURRENT_SIZE = None
"""The current size of the cache."""

_LAST_SIZE = None
"""The size of the cache the last time it was saved."""


try:
    _CACHE = pickle.load(open(_CACHE_FILE, "rb"))
    _CURRENT_SIZE = len(_CACHE)
    _LAST_SIZE = len(_CACHE)
except FileNotFoundError:
    _CACHE = dict()
    _CURRENT_SIZE = 0
    _LAST_SIZE = 0




def folding_params() -> Tuple[str]:
    """
    See the default folding parameters being used with Vienna.
    :rtype: Tuple[str]
    """
    return _DEFAULT_PARAMS

def fold(sequence: str, params: Tuple[str] = _DEFAULT_PARAMS) -> FoldResult:
    """
    Uses RNAfold to predict the mfe for a structure. Does NOT cache the result.
    
    :param: str sequence: RNA sequence ot be folded.
    :param: Tuple[str] params: default folding params for RNAfold, defaults to ('-p','--noLP','-d2')
	:rtype: FoldResult
    """

    def get_mfe(raw):
        tk = ""
        for ch in raw:
            if ch != ")" and ch != "(":
                tk += ch
        return float(tk)


    raw_result = (
        os.popen(f"RNAfold {' '.join(list(params))} <<< {sequence}").read().splitlines()
    )
    safe_rm("dot.ps")
    safe_rm("rna.ps")
    raw_result = list(map(lambda raw: raw.strip(), raw_result))
    folded, raw_nrg = raw_result[1].split(" ", 1)
    fold_result = FoldResult(
        seq=raw_result[0],
        ss=folded,
        mfe=get_mfe(raw_nrg),
        ed=float(raw_result[4].split()[-1]),
        params=params,
    )
    return fold_result

def fold_cache(sequence: str, params: Tuple[str] = _DEFAULT_PARAMS) -> FoldResult:
    """
    Uses RNAfold to predict the mfe for a structure using a global cache of results to save time.
    
    :param: str sequence: RNA sequence ot be folded.
    :param: Tuple[str] params: default folding params for RNAfold, defaults to ('-p','--noLP','-d2')

	:rtype: FoldResult
    """
    global _LAST_SIZE
    global _CURRENT_SIZE
    global _CACHE

    def get_mfe(raw):
        tk = ""
        for ch in raw:
            if ch != ")" and ch != "(":
                tk += ch
        return float(tk)

    if sequence in _CACHE:
        cached = _CACHE[sequence]
        #if cached.params == params: TODO this seems to not work
        return cached

    raw_result = (
        os.popen(f"RNAfold {' '.join(list(params))} <<< {sequence}").read().splitlines()
    )
    safe_rm("dot.ps")
    safe_rm("rna.ps")
    raw_result = list(map(lambda raw: raw.strip(), raw_result))
    folded, raw_nrg = raw_result[1].split(" ", 1)
    fold_result = FoldResult(
        seq=raw_result[0],
        ss=folded,
        mfe=get_mfe(raw_nrg),
        ed=float(raw_result[4].split()[-1]),
        params=params,
    )
    _CACHE[sequence] = fold_result
    _CURRENT_SIZE += 1
    if abs(_CURRENT_SIZE - _LAST_SIZE) > 250:
        save_cache()
    return fold_result


def save_cache() -> None:
    """
    Saves fold results to the cache and updats the _LAST_SIZE value to _CURRENT_SIZE.
    """
    global _LAST_SIZE
    global _CURRENT_SIZE
    if _LAST_SIZE == len(_CACHE):
        return
    _LAST_SIZE = _CURRENT_SIZE
    pickle.dump(_CACHE, open(_CACHE_FILE, "wb"))
