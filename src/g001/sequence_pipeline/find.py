# -*- coding: utf-8 -*-
from __future__ import annotations

import gzip
import logging
import multiprocessing
from pathlib import Path
from typing import Any, List

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from g001.data import Data

logger = logging.getLogger("Find")


def parse_single_fastq_entry(seq_record: SeqRecord) -> pd.Series[Any]:
    """parse an individual sequnce into a dataseries. These build up dataframes

    Parameters
    ----------
    seq_record : SeqRecord
        biopython sequence record

    Returns
    -------
    pd.Series[Any]
        dataframe row of sequence information
    """
    quality_records: List[int] = seq_record.letter_annotations["phred_quality"]  # type: ignore
    _series: pd.Series[Any] = pd.Series(
        {
            "fastq_sequence_id": seq_record.id,  # type: ignore
            "fastq_sequence": str(seq_record.seq),  # type: ignore
            "fastq_phred": quality_records,
            "fastq_mean_phred": np.mean(quality_records),  # type: ignore
        }
    )
    return _series


def parse_single_fastq_file(file: Path) -> pd.DataFrame:
    """Parse a single fastq file into a dataframe

    Parameters
    ----------
    file : Path
        Path to fastq file

    Returns
    -------
    pd.DataFrame
        single dataframe fastq
    """
    logger.debug(f"parse_single_fastq_file: {file}")
    single_fastq: List[SeqRecord] = list(SeqIO.parse(gzip.open(file, "rt"), "fastq"))  # type: ignore
    _dataframe = pd.DataFrame(list(map(parse_single_fastq_entry, single_fastq)))
    _dataframe["fastq_file"] = file.name
    return _dataframe


def find(data: Data) -> pd.DataFrame:
    """Parse all fastq files in the dataframe

    Parameters
    ----------
    data : Data
        Data object

    Returns
    -------
    pd.DataFrame
        dataframe with all fastq files parsed
    """
    all_sequence_files = data.get_fastq_files()
    logger.info(f"Found {len(all_sequence_files):,} fastq files...parsing")

    # using multiprocessing to speed up parsing
    pool = multiprocessing.Pool()

    # concat many dataframes together
    results = pd.concat(pool.map(parse_single_fastq_file, all_sequence_files)).reset_index(drop=True)  # type: ignore
    logger.info(f"Parsed {len(results):,} sequences into dataframe")
    return results
