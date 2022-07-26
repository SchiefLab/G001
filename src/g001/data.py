# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from pydantic import BaseModel, validator

logger = logging.getLogger("DataManage")


class SequenceDataPaths(BaseModel):
    base_path: Path
    base_sequence_path: Optional[Path] = None
    fastq_path: Optional[Path] = None
    manifest_path: Optional[Path] = None
    ngs_path: Optional[Path] = None
    controls_path: Optional[Path] = None
    dose_path: Optional[Path] = None
    haplotype_path: Optional[Path] = None

    @validator("base_path")
    def validate_base_path(cls, v: Path, values: Dict[str, Path]) -> Path:
        if not v.exists():
            raise FileNotFoundError(f"{v} does not exist")
        if "base_sequence_path" in values or "fastq_path" in values or "manifest_path" in values:
            raise ValueError("base_sequence_path,fastq_path/manifest_path cannot be set if base_path is set")
        values["base_sequence_path"] = v / Path("sequence")
        values["fastq_path"] = values["base_sequence_path"] / Path("fastq/")
        values["manifest_path"] = values["base_sequence_path"] / Path("manifests/")
        values["ngs_path"] = values["base_sequence_path"] / Path("ngs/")
        values["controls_path"] = values["base_sequence_path"] / Path("controls/")
        values["dose_path"] = values["base_sequence_path"] / Path("dose_groups/")
        values["haplotype_path"] = values["base_sequence_path"] / Path("haplotype/")
        return v

    @validator("base_sequence_path", always=True)
    def validate_base_sequence_path(cls, v: Path | None, values: Dict[str, Path]) -> Path | None:
        if v is None:
            return values["base_sequence_path"]

    @validator("fastq_path", always=True)
    def validate_fastq_path(cls, v: Path | None, values: Dict[str, Path]) -> Path | None:
        if v is None:
            return Path(values["fastq_path"])

    @validator("manifest_path", always=True)
    def validate_manifest_path(cls, v: Path | None, values: Dict[str, Path]) -> Path | None:
        if v is None:
            return Path(values["manifest_path"])

    @validator("ngs_path", always=True)
    def validate_ngs_path(cls, v: Path | None, values: Dict[str, Path]) -> Path | None:
        if v is None:
            return Path(values["ngs_path"])

    @validator("controls_path", always=True)
    def validate_control_path(cls, v: Path | None, values: Dict[str, Path]) -> Path | None:
        if v is None:
            return Path(values["controls_path"])

    @validator("dose_path", always=True)
    def validate_dose_path(cls, v: Path | None, values: Dict[str, Path]) -> Path | None:
        if v is None:
            return Path(values["dose_path"])

    @validator("haplotype_path", always=True)
    def validate_haplotype_path(cls, v: Path | None, values: Dict[str, Path]) -> Path | None:
        if v is None:
            return Path(values["haplotype_path"])


class Data:
    def __init__(self, base_path: Path):
        self.sequence_data_paths = SequenceDataPaths(base_path=base_path)
        self.intra_plate_cluster_distance: int = 2

    def get_fastq_files(self) -> List[Path]:
        """
        Parse all fastq files in the fastq directory to a list of SeqRecord objects
        """
        return list(self.sequence_data_paths.fastq_path.glob("*.fastq.gz"))  # type: ignore

    def get_vrc_flow_manifest(self) -> pd.DataFrame:
        """Get the VRC flow manifest file. This matches sequences to flow properties

        Returns
        -------
        pd.DataFrame
            The VRC flow manifest file as a dataframe
        """
        return pd.read_csv(self.sequence_data_paths.manifest_path / Path("VRC_flow_manifest.csv.gz"), index_col=0)  # type: ignore

    def get_fhcrc_flow_manifest(self) -> pd.DataFrame:
        """Get the Fred Hutch Cancer Research Center flow manifest file. This matches sequences to flow properties

        Returns
        -------
        pd.DataFrame
            The FHCRC flow manifest file as a dataframe
        """
        return pd.read_csv(self.sequence_data_paths.manifest_path / Path("FHCRC_flow_manifest.csv.gz"), index_col=0)  # type: ignore

    def get_ngs_resequencing_file(self) -> pd.DataFrame:
        """Get the NGS resequencing file. This is sanger sequecnes that were verified by NGS

        Returns
        -------
        pd.DataFrame
            The NGS resequencing file as a dataframe
        """
        return pd.read_csv(self.sequence_data_paths.ngs_path / Path("ngs_dataset.csv.gz"), index_col=0)  # type: ignore

    def get_controls_datafrmae(self) -> pd.DataFrame:
        """Get the controls dataframe. This is a dataframe of controls that were used in the resequencing process

        Returns
        -------
        pd.DataFrame
            The controls dataframe as a dataframe
        """
        return pd.read_csv(self.sequence_data_paths.controls_path / Path("controls_airr.csv.gz"), index_col=0)  # type: ignore

    def get_intra_plate_cluster_distance(self) -> int:
        """Get the intra plate cluster distance to be used in the intra plate clustering algorithm

        Returns
        -------
        int
            The intra plate cluster distance
        """
        return self.intra_plate_cluster_distance

    def get_dose_group(self) -> pd.DataFrame:
        """Get the dose group csv that maps pubids to dose and placebo broups

        Returns
        -------
        pd.DataFrame
            a pandas dataframe with the dose groups, placebo and pub id
        """
        return pd.read_csv(self.sequence_data_paths.dose_path / Path("dose_groups.csv.gz"), index_col=0)  # type: ignore

    def get_haplotype_data(self) -> pd.DataFrame:
        """Get the dose group csv that maps pubids to dose and placebo broups

        Returns
        -------
        pd.DataFrame
            a pandas dataframe with the dose groups, placebo and pub id
        """
        return pd.read_csv(self.sequence_data_paths.haplotype_path / Path("haplotype.csv.gz"), index_col=0)  # type: ignore
