# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from pydantic import BaseModel, validator

logger = logging.getLogger("DataManage")

class PlotParameters:
    def __init__(self):
        # define visit lookups for
        self.memory_visit_lookup = {
            "V02": -4,
            "V06": 4,
            "V07": 8,
            "V08": 10,
            "V10": 16,
        }
        # germinal center timepoints
        self.gc_visit_lookup = {"V05": 3, "V09": 11}

        # plasmablast timepoints
        self.pb_visit_lookup = {"V07A": 9}

        # and combine all timepoitns
        self.combined_timepoints = {**self.memory_visit_lookup, **self.gc_visit_lookup, **self.pb_visit_lookup}

        # any annotations we make will ahve these arguments
        self.small_annotate_args = dict(
            xycoords="axes fraction",
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=9,
            fontweight="bold",
        )

        self.placebo_color = "#A6A6A6"
        self.low_dose_vrc01_color = "#6C65FF"
        self.high_dose_vrc01_color = "#B885EE"
        self.low_dose_nonvrc01_color = "#F7C3B1"
        self.high_dose_nonvrc01_color = "#FFF059"

        self.dose_palette = dict(
            zip(["Low Dose", "High Dose"], [self.low_dose_vrc01_color, self.high_dose_vrc01_color])
        )
        self.non_vrv01_dose_palette = dict(
            zip(["Low Dose", "High Dose"], [self.low_dose_nonvrc01_color, self.high_dose_nonvrc01_color])
        )

        self.treatment_palette = dict(
            zip(
                ["DPBS sucrose", "20 µg eOD-GT8 60mer + AS01B", "100 µg eOD-GT8 60mer + AS01B"],
                [self.placebo_color, self.low_dose_vrc01_color, self.high_dose_vrc01_color],
            )
        )
        self.treatment_palette_nonvrc01 = dict(
            zip(
                ["DPBS sucrose", "20 µg eOD-GT8 60mer + AS01B", "100 µg eOD-GT8 60mer + AS01B"],
                [self.placebo_color, self.low_dose_nonvrc01_color, self.high_dose_nonvrc01_color],
            )
        )
        self.full_palette = {
            "Low Dose": {
                "vrc01_class": self.low_dose_vrc01_color,
                "vrc01_class_kappa": self.low_dose_vrc01_color,
                "vrc01_class_lambda": "#65ABFF",
                "nonvrc01_class": self.low_dose_nonvrc01_color,
                "nonvrc01_class_kappa": self.low_dose_nonvrc01_color,
                "nonvrc01_class_lambda": "#B1E5F7",
            },
            "High Dose": {
                "vrc01_class": self.high_dose_vrc01_color,
                "vrc01_class_kappa": self.high_dose_vrc01_color,
                "vrc01_class_lambda": "#ED85EE",
                "nonvrc01_class": self.high_dose_nonvrc01_color,
                "nonvrc01_class_kappa": self.high_dose_nonvrc01_color,
                "nonvrc01_class_lambda": "#F059FF",
            },
        }

    def get_memory_visit_lookup(self) -> Dict[str, int]:
        """The timepoint and week visit dictionry of memory b cell associated samples

        Returns
        -------
        Dict[str, int]
        """
        return self.memory_visit_lookup

    def get_gc_visit_lookup(self) -> Dict[str, int]:
        """The timepoint and week visit dictionry of germinal center derived (FNA) associated samples

        Returns
        -------
        Dict[str, int]
        """
        return self.gc_visit_lookup

    def get_pb_visit_lookup(self) -> Dict[str, int]:
        """The timepoint and week visit dictionry of plasmablasts associated samples

        Returns
        -------
        Dict[str, int]
        """
        return self.pb_visit_lookup

    def get_combined_timepoints(self) -> Dict[str, int]:
        """The timepoint and week visit dictionry of all associated samples

        Returns
        -------
        Dict[str, int]
        """
        return self.combined_timepoints

    def get_pallete(self) -> Dict[str, Dict[str, str]]:
        """The color palette for the plots

        Low Dose:
            vrc01_class: #6C65FF
            vrc01_class_kappa: #6C65FF
            vrc01_class_lambda: #65ABFF
            nonvrc01_class: #F7C3B1
            nonvrc01_class_kappa: #F7C3B1
            nonvrc01_class_lambda: #B1E5F7
        High Dose:
            vrc01_class: #B885EE
            vrc01_class_kappa: #B885EE
            vrc01_class_lambda: #ED85EE
            nonvrc01_class: #FFF059
            nonvrc01_class_kappa: #FFF059
            nonvrc01_class_lambda: #D8FF65

        Returns
        -------
        Dict[str, Dict[str, str]]
        """
        return self.full_palette

    def get_treatment_pallete(self) -> Dict[str, str]:
        """Color palette but lookup by Treatment, e.g. DPBS sucrose, 20 µg eOD-GT8 60mer + AS01B

        DPBS sucrose: #A6A6A6
        20 µg eOD-GT8 60mer + AS01B: #6C65FF
        100 µg eOD-GT8 60mer + AS01B: #B885EE

        Returns
        -------
        Dict[str, str]
        """
        return self.treatment_palette

    def get_nonvrc01_treatment_pallete(self) -> Dict[str, str]:
        """Color palette but lookup by Treatment, e.g. DPBS sucrose, 20 µg eOD-GT8 60mer + AS01B

        DPBS sucrose: #A6A6A6
        20 µg eOD-GT8 60mer + AS01B: #6C65FF
        100 µg eOD-GT8 60mer + AS01B: #B885EE

        Returns
        -------
        Dict[str, str]
        """
        return self.treatment_palette_nonvrc01

    def get_small_annotate_args(self) -> Dict[str, str | int]:
        """The arguments for the small annotations on the plots

        Returns
        -------
        Dict[str, str]
        """
        return self.small_annotate_args

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
class FigureDataPaths(BaseModel):
    base_path: Path = Path('data')
    flow_and_frequency_path : Path = base_path / Path('figures/flow_summary/flow_and_sequences.feather')

    @validator("*",always=True)
    def validate_paths(cls, v:Path):
        if not v.exists():
            raise FileNotFoundError(f"{v} is not found")
        return v

class Data:
    def __init__(self, base_path: Path):
        self.sequence_data_paths = SequenceDataPaths(base_path=base_path)
        self.figure_data_paths = FigureDataPaths(base_path=base_path)
        self.intra_plate_cluster_distance: int = 2
        self.plot_parameters = PlotParameters()

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

    def get_flow_and_frequency_data(self) -> pd.DataFrame:
        """Returns the combined flow and frequency dataframe

        Returns
        -------
        pd.DataFrame
            a pandas dataframe with all computed values from the combining flow and sequencing module
        """
        return pd.read_feather(self.figure_data_paths.flow_and_frequency_path)