"""
This type stub file was generated by pyright.
"""

from enum import Enum
from typing import List, Tuple
from pandas import DataFrame, Index

def chunk(nb_item: int, nb_chunks: int, start_offset=...) -> List[slice]:
    """
    Return `nb_chunks` slices of approximatively `nb_item / nb_chunks` each.

    Parameters
    ----------
    nb_item : int
        Total number of items

    nb_chunks : int
        Number of chunks to return

    start_offset : int
        Shift start of slice by this amount

    Returns
    -------
    A list of slices

    Examples
    --------
    >>> chunks = chunk(103, 4)
    >>> chunks
    [slice(0, 26, None), slice(26, 52, None), slice(52, 78, None), slice(78, 103, None)]
    """
    ...

def df_indexed_like(df: DataFrame, axes: List[Index]) -> bool:
    """
    Returns whether a data frame is indexed in the way specified by the
    provided axes.

    Used by DataFrameGroupBy to determine whether a group has been modified.

    Function adapted from pandas.core.groupby.ops._is_indexed_like

    Parameters
    ----------
    df : DataFrame
        The data frame in question

    axes : List[Index]
        The axes to which the data frame is compared

    Returns
    -------
    Whether or not the data frame is indexed in the same wa as the axes.
    """
    ...

def get_pandas_version() -> Tuple[int, int]:
    ...

class WorkerStatus(int, Enum):
    Running = ...
    Success = ...
    Error = ...

