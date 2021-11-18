"""Utility functions for transforming rows in a reactivity dataframe"""
from __future__ import annotations

import re
import numpy as np
import pandas as pd
from glob import glob
from typing import List
from rna_library.core import dsci
from collections import defaultdict
from rna_library.structure import Motif
from .junction_data import JunctionEntry


def add_reactivity(
    row: pd.Series, output_directory: str, col_name: str = "mismatches"
) -> List[float]:
    """
    Function that gets the raw reactivities for a given construct

    :param: pandas.Series row: row to get reacitivity for
    :param: str output_directory : the base directory from which to get the population files
	:param: str col_name: column name to get values from, defaults to 'mismatches'
	:rtype: List[float]
    """
    # TODO the output directory should probably be indicated somewhere as a param
    raw_data = list(glob(f'{output_directory}/{row["construct"]}_*popavg_reacts.csv'))
    assert len(raw_data) == 1, len(raw_data)
    df = pd.read_csv(raw_data[0])
    num_rows = len(df.index)
    assert num_rows == len(row["RNA"])
    return df[col_name].to_list()


def block_commons(row: pd.Series, start: str, end: str) -> str:
    """
    Function that blocks off the common start and sequence of an RNA construct with N's
	
    :param: pandas.Series row: row to block off commons for
	:param: str start: the common start sequence
	:param: str end: the common end sequence
	:rtype: str
    """
    sequence = row["RNA"]
    assert sequence.startswith(start)
    assert sequence.endswith(end)
    sequence = re.sub(f"^{start}", "N" * len(start), sequence)
    sequence = re.sub(f"{end}$", "N" * len(end), sequence)
    return sequence


def score(row: pd.Series) -> float:
    """
    Function that generates a dsci score for an RNA and its reactivity values. Value is on 
    the range [0,1] with 0.95 being a common quality cutoff.

	:param: pandas.Series row: row to get a score for  
    :rtype: float
    """
    # ok now we actually have all of the values needed to find the dsci score
    return dsci(row["blocked"], row["structure"], row["reactivity"])[0]


def signal_to_noise(row: pd.Series) -> float:
    """
    Function that calculates the signal to noise ratio for a DMS entry by using the 
    ratio of mutations for (A + C)/(G + U).

	:param: pandas.Series row: row to get sn ratio for  
    :rtype: float
    """
    reactivity = row["reactivity"]

    seq = row["blocked"]
    AC = 0
    GU = 0
    AC_count = seq.count("A") + seq.count("C")
    GU_count = seq.count("G") + seq.count("U") + seq.count("T")

    for idx, nt in enumerate(seq):
        if nt == "N":
            continue

        if nt == "A" or nt == "C":
            AC += reactivity[idx]
        else:
            GU += reactivity[idx]

    AC /= float(AC_count)
    GU /= float(GU_count)

    return round(float(AC / GU), 2)


def num_reads(row: pd.Series, histos: Dict[str, any]) -> int:
    """
    Function that finds the number of reads for a given DMS entry row.

    :param: pandas.Series row: row to get the number of reads for
    :param: Dict[str,dreem.MutationHistogram] histos: histogram dictionary that has read information
    :rytype; int
    """
    # TODO get this type hint actually right
    hist = histos[row["construct"]]
    return hist.num_reads


def collect_junction_entries(
    m: Motif,
    reactivity: List[float],
    construct: str,
    sn: float,
    reads: int,
    score: float,
    holder: Dict[str, List[JunctionEntry]],
) -> None:
    """
    Utility function that gets all JunctionEntry objects across a reactivity dataframe.

    :param: Motif m: the base motif to get data for
    :param: List[int] reactivity: list of reactivity values for the construct
    :param: str construct: name of the construct
    :param: float sn: signal to noise ratio of the construct
    :param: int reads: number of sequencer reads for the construct
    :param: float score: DSCI score for the construct
	:param: Dict[str,List[JunctionEntry]] holder: temporary holder for all of the JunctionEntry objects 
    """
    if m.is_junction():
        m: Junction
        strands = m.strands()
        assert len(strands) == 2
        reacts = (
            [reactivity[idx] for idx in strands[0]]
            + [-1]
            + [reactivity[idx] for idx in strands[1]]
        )
        je = JunctionEntry(
            sequence=m.sequence(),
            structure=m.structure(),
            reactivity=reacts,
            construct=construct,
            sn=sn,
            reads=reads,
            score=score,
        )

        holder[je.key()].append(je)

    for c in m.children():
        collect_junction_entries(c, reactivity, construct, sn, reads, score, holder)


def row_normalize_hairpin(
    row: pd.Series, norm_seq: str, norm_ss: str, factor: float, nts: List[str]
) -> List[float]:
    """
    Function that performs a hairpin normalization on a pd.Series representing a construct.
    Returns the normalized reactivity series.

    :param: pd.Series row: dataframe row describing a construct. must have 'RNA', 'structure' and 'reactivity' columns
    :param: str norm_seq: normalization hairpin sequence
    :param: str norm_ss: normalize hairpin secondary structure
    :param: float factor: factor to which the reference value will be set
    :param: List[str] nts: nucleotides to be considered in the normalization scheme. must be unpaired!
    :rtype: List[float]
	"""
    seq, ss, react = row["RNA"], row["structure"], row["reactivity"]
    assert seq.count(norm_seq)
    assert ss.count(norm_ss)

    for idx, it in enumerate(re.finditer(norm_seq, seq)):
        if ss[it.start() : it.end()] == norm_ss:
            idx = it.start()
            break

    assert idx

    vals = []

    for ii in range(idx, idx + len(norm_seq)):
        if ss[ii] == "." and seq[ii] in nts:
            vals.append(react[ii])

    assert len(vals) > 1

    vals = np.array(vals)
    react = np.array(react)
    avg = np.mean(vals)

    return factor * react / avg
