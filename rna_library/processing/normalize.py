"""Tools for normalizing reactivity"""
import re
import numpy as np
import pandas as pd
from typing import List, Tuple
from .row_utils import row_normalize_hairpin

def normalize_hairpin( df: pd.DataFrame, seq : str, ss : str, **kwargs ) -> List[List[float]]:
    """
    Normalizes a reactivity pattern to a normalization hairpin. Creates fully normalize values
    for an entire pd.DataFrame

    :param: pd.DataFrame df: reactivity_df created from rna_library.build_react_df 
    :param: str seq: reference hairpin sequence
    :param: str ss: reference hairpin structure
    :param: float factor: factor to set the hairpin values to, is a keyword argument
    :param: str nts: string of nucleotide's to be used for calc, is a keywrod argument
    """
    #TODO add argument validation
    factor = kwargs.get('factor', 1)
    nts = list(kwargs.get('nts','A'))
	
    return df.apply( lambda row: row_normalize_hairpin( row, seq, ss, factor, nts ), axis=1 )
