"""
core module responsible for enums, nucleotide methodology, errors, utility functions and folding algorithms.

Author: Chris Jurich <cjurich2@unl.edu>
Date: 2022-04-29
"""
from .enums import (
    IntEnum,
    ALLOWED_PAIRS,
    BPS,
    LEGAL_BPS,
    NTS,
    BasePair,
    BP_VALS,
    BASEPAIR_MAPPER,
    Nucleotide,
    NT_VALS,
    NUCLEOTIDE_MAPPER,
    MotifType,
    TYPE_MAPPER,
)

from .error import InvalidDotBracket, MissingDependency, InvalidArgument

from .util import (
    valid_db,
    connectivity_list,
    is_circular,
    load_fasta,
    is_symmetrical,
    safe_mkdir,
    safe_rm,
    dsci,
)

from .folding import fold_cache, FoldResult, save_cache, folding_params
