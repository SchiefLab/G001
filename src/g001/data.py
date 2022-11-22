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
    base_path: Path = Path("data")
    flow_and_frequency_path: Path = base_path / Path("figures/flow_summary/flow_and_sequences.feather")
    sequence_path: Path = base_path / Path("figures/sequences/unblinded_sequences.feather")
    distance_path_vrc01_path: Path = base_path / Path("figures/cluster/distance_df_vrc01.feather")
    distance_path_nonvrc01_path: Path = base_path / Path("figures/cluster/distance_df_nonvrc01.feather")
    dekosky_vh12_path: Path = base_path / Path("controls/dekosky_vh12.feather")
    kappa_vrc01_aa_path: Path = base_path / Path("controls/kappa_vrc01_aa.feather")
    lambda_vrc01_aa_path: Path = base_path / Path("controls/lambda_vrc01_aa.feather")
    oas_five_len_path: Path = base_path / Path("controls/oas_5_lc_5_len.feather")
    oas_vh12_path: Path = base_path / Path("controls/oas_vh12.feather")
    vrc01_reference_airr_path: Path = base_path / Path("controls/vrc01_class_bnabs_ref_airr.feather")
    human_naive_path: Path = base_path / Path("controls/human_naive_vrc01_class.feather")
    cottrell_path: Path = base_path / Path("controls/cottrell.json")
    igl_spr_path: Path = base_path / Path("figures/binding/igl_spr.feather")
    clk_spr_path: Path = base_path / Path("figures/binding/clk_spr.feather")
    mature_spr_path: Path = base_path / Path("figures/binding/mature_spr.feather")
    boost_v_gt8_trend_path: Path = base_path / Path("figures/trends/boost_v_gt8.csv")

    @validator("*", always=True)
    def validate_paths(cls, v: Path):
        if not v.exists():
            raise FileNotFoundError(f"{v} is not found")
        return v


class FlowDataPaths(BaseModel):
    base_path: Path = Path("data/flow")
    processed_flow_path: Path = base_path / Path("flow_processed_out")

    @validator("*", always=True)
    def validate_paths(cls, v: Path):
        if not v.exists():
            raise FileNotFoundError(f"{v} is not found")
        return v


class TableDataPaths(BaseModel):
    base_path: Path = Path("data/tables")
    table_src: Path = base_path / Path("input")
    table_out: Path = base_path / Path("output")

    @validator("*", always=True)
    def validate_paths(cls, v: Path):
        if not v.exists():
            raise FileNotFoundError(f"{v} is not found")
        return v


class Data:
    def __init__(self, base_path: Path):
        self.sequence_data_paths = SequenceDataPaths(base_path=base_path)
        self.flow_data_paths = FlowDataPaths(base_path=base_path / Path("flow"))
        self.figure_data_paths = FigureDataPaths(base_path=base_path)
        self.intra_plate_cluster_distance: int = 2
        self.plot_parameters = PlotParameters()
        self.table_data_paths = TableDataPaths(base_path=base_path)

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

    def get_unblinded_sequences(self) -> pd.DataFrame:
        """Returns the unblinded sequences dataframe

        Returns
        -------
        pd.DataFrame
            a pandas dataframe with all computed values from the combining flow and sequencing module
        """
        return pd.read_feather(self.figure_data_paths.sequence_path)

    def get_distance_df(self, is_vrc01_class: bool) -> pd.DataFrame:
        """Get the pre-computed distance dataframe for all vrc01 and non-vrc01 sequences

        Returns
        -------
        is_vrc01_class: bool
            if True, return the vrc01 distance dataframe, otherwise return the non-vrc01 distance dataframe
        """
        if is_vrc01_class:
            return pd.read_feather(self.figure_data_paths.distance_path_vrc01_path)
        return pd.read_feather(self.figure_data_paths.distance_path_nonvrc01_path)

    def get_dekosky(self) -> pd.DataFrame:
        """Get VH1-2 pairs from Dekosky et al

        Returns
        -------
        pd.DataFrame
            the dekosky vh12 dataframe
        """
        return pd.read_feather(self.figure_data_paths.dekosky_vh12_path)

    def get_kappa_vrc01_aa(self) -> pd.DataFrame:
        """Get the kappa vrc01 aa dataframe

        Returns
        -------
        pd.DataFrame
            the kappa vrc01 aa dataframe from an ANARCI alignment
        """
        return pd.read_feather(self.figure_data_paths.kappa_vrc01_aa_path)

    def get_lambda_vrc01_aa(self) -> pd.DataFrame:
        """Get the lambda vrc01 aa dataframe

        Returns
        -------
        pd.DataFrame
            the lambda vrc01 aa dataframe from an ANARCI alignment
        """
        return pd.read_feather(self.figure_data_paths.lambda_vrc01_aa_path)

    def get_oas_five_len(self) -> pd.DataFrame:
        """Get the oas five len dataframe

        A collection of HD light chains that have a LCDR3 len of 5

        Returns
        -------
        pd.DataFrame
            the oas five len dataframe
        """
        return pd.read_feather(self.figure_data_paths.oas_five_len_path)

    def get_oas_vh12(self) -> pd.DataFrame:
        """Get the oas vh12 dataframe

        A collection of HD VH1-2 heavy chains

        Returns
        -------
        pd.DataFrame
            the oas vh12 dataframe
        """
        return pd.read_feather(self.figure_data_paths.oas_vh12_path)

    def get_vr01_ref_airr(self) -> pd.DataFrame:
        """Get the vrc01 reference airr dataframe

        Returns
        -------
        pd.DataFrame
            the vrc01 reference airr dataframe
        """
        return pd.read_feather(self.figure_data_paths.vrc01_reference_airr_path)

    def get_human_naive_airr(self) -> pd.DataFrame:
        """Get the human naive airr dataframe

        Returns
        -------
        pd.DataFrame
            the human naive airr dataframe
        """
        return pd.read_feather(self.figure_data_paths.human_naive_path)

    def get_cottrell_path(self) -> Path:
        """Get the cottrell path

        Returns
        -------
        Path
            the cottrell path
        """
        return self.figure_data_paths.cottrell_path

    def get_igl_spr_df(self) -> pd.DataFrame:
        """Get the igl spr dataframe

        Returns
        -------
        pd.DataFrame
            the igl spr dataframe
        """
        return pd.read_feather(self.figure_data_paths.igl_spr_path)

    def get_clk_spr_df(self) -> pd.DataFrame:
        """Get the clk/hugl (human naives) spr dataframe

        Returns
        -------
        pd.DataFrame
            the clk/hugl spr dataframe
        """
        return pd.read_feather(self.figure_data_paths.clk_spr_path)

    def get_mature_spr_df(self) -> pd.DataFrame:
        """Get the mature spr dataframe

        Returns
        -------
        pd.DataFrame
            the mature spr dataframe
        """
        return pd.read_feather(self.figure_data_paths.mature_spr_path)

    def get_boost_v_gt8_trend(self) -> pd.DataFrame:
        """Get the boost v gt8 trend dataframe

        Returns
        -------
        pd.DataFrame
            the boost v gt8 trend dataframe
        """
        return pd.read_csv(self.figure_data_paths.boost_v_gt8_trend_path, index_col=0)

    def get_processed_flow_paths(self) -> Path:
        return self.flow_data_paths.processed_flow_path
