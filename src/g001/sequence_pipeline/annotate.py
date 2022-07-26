from __future__ import annotations
from ntpath import join
from typing import Any, List

import numpy as np
from g001.data import Data
import pandas as pd
import logging
from sadie.airr import Airr

logger = logging.getLogger("Annotate")


def slice_phred_row(row: pd.Series[Any]) -> float | List[int]:
    """
    Slice the phreds per row
    """
    phred_list: List[int] = row[0]  # type: ignore
    sequence: str = row[1]  # type: ignore
    sequence_alignment: str = row[2]  # type: ignore
    rev_comp: bool = row[3]  # type: ignore
    if rev_comp:
        # If we have rev comp reverse the phred
        phred_list = phred_list[::-1]

    if pd.isna(sequence_alignment):  # type: ignore
        return np.nan

    start = sequence.index(sequence_alignment.replace("-", ""))
    end = start + len(sequence_alignment.replace("-", ""))
    # of the phred list, return the slice that is just the sequence alignment
    return phred_list[start:end]


def annotate(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Correct the data, unused..precalculated
    """
    logger.info(f"attempting to annotate {len(working_dataframe):,}")
    airr_api = Airr("human", adaptable=True)
    airr_results = airr_api.run_dataframe(working_dataframe, "fastq_sequence_id", "fastq_sequence")

    # rename fastq_sequence_id to sequence_id to make it joinable
    joined_dataframe = working_dataframe.rename({"fastq_sequence_id": "sequence_id"}, axis=1).merge(
        airr_results, on="sequence_id"
    )

    if len(joined_dataframe) != len(working_dataframe):
        raise ValueError(f"dataframes are not equal size, {joined_dataframe.shape} != {working_dataframe.shape}")

    logger.info("Slicing phreds")
    targeted_phreds = joined_dataframe[["fastq_phred", "sequence", "sequence_alignment", "rev_comp"]].apply(
        lambda x: slice_phred_row(x), axis=1
    )
    # Now insert those sliced_abi_phred after abi_phred
    joined_dataframe.insert(joined_dataframe.columns.get_loc("fastq_phred") + 1, "sliced_phred", targeted_phreds)

    # Now grab the means of just the sliced portion and put it abi mean
    logger.info("Getting means of sliced phreds")
    joined_dataframe.insert(
        joined_dataframe.columns.get_loc("fastq_phred") + 1,
        "mean_sliced_phred",
        joined_dataframe["sliced_phred"].apply(np.mean),
    )
    return joined_dataframe
