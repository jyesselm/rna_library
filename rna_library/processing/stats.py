"""A collection of functions for measuring statistical variation in data entries"""
import numpy as np
from typing import List
#from scipy.stats import median_absolute_deviation

def iqr_excess(data: List[float]) -> float:
    """
    Finds the inter-quartiles excess, or the sum of excess range outside of the 
    first quartile or third quartile minus and plus the inter-quartile range, respectively.
    
    :param: List[float] data: data series to analyze.
    :rtype: float
    """
    if not len(data):
        return 0
    (first, third) = np.percentile(data, [25, 75])
    iqr = third - first
    lower, upper = first - iqr, third + iqr
    total = 0
    for d in data:
        if d <= lower:
            total += abs(lower - d)
        elif d >= upper:
            total += abs(upper - d)

    return total


def sigma_excess(data):
    if not len(data):
        return 0
    mean, std = np.mean(data), np.std(data)
    return np.sum(np.abs(data - mean) / std)


def sigma_excess_normed(data):
    if not len(data):
        return 0
    mean, std = np.mean(data), np.std(data)
    std /= mean
    return np.sum(np.abs(data - 1) / std)

def mad_excess(data):
    if not len(data):
        return 0
    raise NotImplemented("Does not work now")
    #return median_absolute_deviation(data)

def comparative_variance( data1, data2 ):
    assert len(data1) == len(data2)
    total = 0.
    for d1, d2 in zip(data1, data2):
        if not (d1 + d2):
            continue
        total += abs(d1-d2)/(d1+d2)
    return total


EXCESS_MAPPER = {
    "iqr_excess": iqr_excess,
    "sigma_excess": sigma_excess,
    "sigma_excess_normed": sigma_excess_normed,
    "mad_excess": mad_excess,
}
