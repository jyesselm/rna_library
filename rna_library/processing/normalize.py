"""Tools for normalizing reactivity"""
import re
import numpy as np
import pandas as pd
from typing import List, Tuple
from .stats import EXCESS_MAPPER
from itertools import accumulate
from scipy.optimize import minimize
from scipy.stats import zscore
from .process_histos import build_motif_df
from .row_utils import row_normalize_hairpin
from collections import namedtuple, defaultdict
import matplotlib.pyplot as plt

# seed the rng
np.random.seed( 100 ) 

def normalize_hairpin(
    df: pd.DataFrame, seq: str, ss: str, **kwargs
) -> List[List[float]]:
    """
    Normalizes a reactivity pattern to a normalization hairpin. Creates fully normalize values
    for an entire pd.DataFrame

    :param: pd.DataFrame df: reactivity_df created from rna_library.build_react_df 
    :param: str seq: reference hairpin sequence
    :param: str ss: reference hairpin structure
    :param: float factor: factor to set the hairpin values to, is a keyword argument
    :param: str nts: string of nucleotide's to be used for calc, is a keywrod argument
    """
    # TODO add argument validation
    factor = kwargs.get("factor", 1)
    nts = list(kwargs.get("nts", "A"))

    return df.apply(
        lambda row: row_normalize_hairpin(row, seq, ss, factor, nts), axis=1
    )


class OptimizeJunction:
    def __init__(self, sequence, n_data ):
        self.sequence = sequence
        num_pts = sequence.count('A')+sequence.count('C')
        self.data = np.zeros((num_pts, n_data))
        self.n_data = n_data

    def add_data(self, row, data_idx):
        for r_idx, val in enumerate( row ):
            self.data[r_idx][data_idx] = val

    def spread( self, how ):
        func = EXCESS_MAPPER[ how ]
        total = 0
        for row in self.data:
            total += func( row ) 
        return total
        #return accumulate(func, self.data )

def remove_outliers( reactivity_df, z_cutoff=2 ):
    """
    """
    motif_df = build_motif_df( reactivity_df )
    
    con_names = reactivity_df['construct'].to_list()
    total_variance = dict(
       zip( con_names, [0]*len( con_names ))
    )
    
    for jd in motif_df.data:
        temp = jd.measure_variance() 
        for cname, val in temp.items():
            total_variance[cname] += val
    
    variances = np.array(list(total_variance.values())) 
    mask = zscore(variances)<z_cutoff

    kept_constructs = np.array(con_names)[mask]
    
    reactivity_df = reactivity_df[reactivity_df.construct.isin(kept_constructs)].reset_index(drop=True)
    motif_df = build_motif_df( reactivity_df )
    
    return (reactivity_df,motif_df)

NormJob = namedtuple("NormJob", "c_idx m_idx idxs data_idx")

def normalize_coeff_fit(reactivity_df):
    raise Exception('Don\'t use this!')
    # what do we want to do?
    def get_junctions(m, holder):
        if m.is_junction():
            holder.append(m)
        for c in m.children():
            get_junctions(c, holder)
    
    def add_jcts( row ):
        holder = []
        get_junctions(row.motif, holder)
        return holder

    def summarize_jcts( df ):
        jct_summary = defaultdict( int )
        for i, row in df.iterrows():
            for j in row.jcts:
                jct_summary[ j.sequence() ] += 1
        return jct_summary
	# print(reactivity_df)
    light_df = reactivity_df[["construct", "reactivity", "motif"]]
    light_df["jcts"] = light_df.apply( lambda row: add_jcts( row ),axis=1)
    
    jct_seqs = summarize_jcts( light_df )


    jct_mapper = dict()
    o_jcts = []
    for idx, (j_seq, j_ct) in enumerate(jct_seqs.items()):
        jct_mapper[j_seq] = idx
        o_jcts.append(OptimizeJunction(j_seq, j_ct))

    jct_idx_holder = dict(zip(list(jct_seqs.keys()),[0]*len(jct_seqs)))
    job_queue = []
    for idx, row in light_df.iterrows():
        for j in row.jcts:
            data_idx = jct_idx_holder[ j.sequence() ]
            jct_idx_holder[ j.sequence() ] += 1
            job_queue.append(
                NormJob(
                    c_idx=idx, m_idx=jct_mapper[j.sequence()], idxs=j.dms_active_idxs(), data_idx=data_idx
                )
            )


    fh = open('/Users/chrisjurich/projects/rna_library/cpp/njs.txt','w')
    
    for nj in job_queue:
        fh.write(f'{nj.c_idx},{nj.m_idx},{"|".join(map(str,nj.idxs))},{nj.data_idx}\n')
    fh.close()
    
    fh = open('/Users/chrisjurich/projects/rna_library/cpp/motifs.txt','w')
    for m in o_jcts:
        fh.write(f'{m.sequence},{m.n_data}\n')
    fh.close()

    def error_fxn( params, reacts, motifs, jq ):
        # apply the params
        for idx, rxt in enumerate( reacts ):
        	rxt *= params[idx]
        # clear out the motifs
        nj : NormJob
        # process the job queue	
        for nj in jq:
        	rxt = reacts[ nj.c_idx ]
        	motifs[ nj.m_idx ].add_data(list(map(lambda ii: rxt[ii], nj.idxs)), nj.data_idx )
        
        # undo the multiplier
        for idx, rxt in enumerate( reacts ):
        	rxt /= params[idx]
        ans = sum(list(map(lambda m: m.spread('iqr_excess'), motifs)))
        print(ans)
        return ans
	
    reactivites = [np.array(v) for v in light_df['reactivity']]
    def score_fxn( params ):
        return error_fxn( params, reactivites, o_jcts, job_queue )
    
    dParams = list(map(lambda idx: 1, range(len(light_df))))
	
    return minimize( score_fxn, dParams)
