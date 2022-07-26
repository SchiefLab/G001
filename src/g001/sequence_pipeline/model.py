# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
from typing import Any, List

import pandas as pd
from pydantic import BaseModel, validator
from g001.data import Data

logger = logging.getLogger("Model")


valid_pub_ids = [
    "PubID001",
    "PubID005",
    "PubID009",
    "PubID014",
    "PubID016",
    "PubID023",
    "PubID028",
    "PubID030",
    "PubID032",
    "PubID036",
    "PubID046",
    "PubID047",
    "PubID051",
    "PubID056",
    "PubID059",
    "PubID060",
    "PubID062",
    "PubID064",
    "PubID068",
    "PubID070",
    "PubID077",
    "PubID079",
    "PubID080",
    "PubID088",
    "PubID092",
    "PubID100",
    "PubID110",
    "PubID112",
    "PubID113",
    "PubID114",
    "PubID116",
    "PubID117",
    "PubID121",
    "PubID151",
    "PubID152",
    "PubID153",
    "PubID154",
    "PubID163",
    "PubID164",
    "PubID165",
    "PubID172",
    "PubID177",
    "PubID180",
    "PubID185",
    "PubID187",
    "PubID191",
    "PubID193",
    "PubID198",
]
valid_timepoints = ["V02", "V05", "V06", "V07", "V07A", "V08", "V09", "V10"]
valid_row_wells = ["A", "B", "C", "D", "E", "F", "G", "H"]


class SequenceIDModel(BaseModel):
    """
    Validate the model
    """

    pub_id: str
    timepoint: str
    plate: str
    chain: str
    well: str
    replicate: str

    @validator("pub_id")
    @classmethod
    def validate_pub_id(cls, v: str) -> str:
        if v not in valid_pub_ids:
            raise ValueError(f"{v} is not a valid pub_id must be one of {valid_pub_ids}")
        return v

    @validator("timepoint")
    @classmethod
    def validate_timepoint(cls, v: str) -> str:
        if v not in valid_timepoints:
            raise ValueError(f"{v} is not a valid timpoint must be one of {valid_timepoints}")
        return v

    @validator("plate")
    @classmethod
    def validate_plate(cls, v: str) -> str:
        p = v[0]
        if p != "P":
            raise ValueError(f"{v} is not a valid plate must start with P")
        integer_part = int(v[1:])
        if integer_part > 18:
            raise ValueError(f"{v} is not a valid plate integer must be less than 18")
        return p + str(integer_part).zfill(2)  # pad with a 0

    @validator("chain")
    @classmethod
    def validate_chain(cls, v: str) -> str:
        if v not in ["HEAVY", "KAPPA", "LAMBDA"]:
            raise ValueError(f"{v} is not a valid chain must be one of HEAVY, KAPPA, LAMBDA")
        return v

    @validator("well")
    @classmethod
    def validate_well(cls, v: str) -> str:
        row_part = v[0]
        if row_part not in valid_row_wells:
            raise ValueError(f"{v} is not a valid well row must be one of {valid_row_wells}")
        column_part = int(v[1:])
        if column_part > 13:
            raise ValueError(f"{v} is not a valid row column integer must be less than 13")
        return row_part + str(column_part).zfill(2)  # pad with a 0

    @validator("replicate")
    @classmethod
    def validate_replicate(cls, v: str) -> str:
        integer_part = int(v)
        if integer_part < 0 or integer_part > 2:
            raise ValueError(f"{v} is not a valid replicate integer must be between 0 and 2")
        return v


def apply_model(sequence_id: str) -> pd.Series[Any]:
    """
    Apply the model to a single row
    """
    tokens: List[str] = sequence_id.split("_")
    model = SequenceIDModel(
        **{
            "pub_id": tokens[0],
            "timepoint": tokens[1],
            "plate": tokens[2],
            "chain": tokens[3],
            "well": tokens[4],
            "replicate": tokens[5],
        }
    )
    return pd.Series(model.dict())


def model(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Model the data
    """
    logger.info(f"Validating {len(working_dataframe):,} sequence_id data into model")
    return working_dataframe["fastq_sequence_id"].apply(apply_model).join(working_dataframe)  # type: ignore
