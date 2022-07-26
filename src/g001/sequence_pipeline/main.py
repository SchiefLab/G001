from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import shutil
from typing import Callable, Optional, Tuple
from g001.sequence_pipeline import (
    find,
    model,
    ngs,
    split,
    correct,
    annotate,
    join,
    tag,
    swap,
    pair,
    unblind,
    personalize,
    mutation,
)
from g001.data import Data
import logging
import pandas as pd


logger = logging.getLogger("G001")


@dataclass
class ModuleRunner:
    cached_output_name: str
    cached_output_path: Path
    data: Data
    callable_function: Callable[[Data], pd.DataFrame] | Callable[[Data, pd.DataFrame], pd.DataFrame]
    resume: bool
    working_dataframe: Optional[pd.DataFrame] = None
    secondary_working_dataframe: Optional[pd.DataFrame] = None
    secondary_cache_output_name: Optional[str] = None

    def run(self) -> pd.DataFrame | Tuple[pd.DataFrame, pd.DataFrame]:
        """Run the sequence module and return a working dataframe

        Returns
        -------
        pd.DataFrame
            The working dataframe
        """
        self.outpath: Path = self.cached_output_path / Path(self.cached_output_name)
        if self.outpath.exists() and self.resume:
            logger.info(f"Loading cached output from {self.outpath}")
            if self.secondary_cache_output_name:
                self.secondary_outpath: Path = self.cached_output_path / Path(self.secondary_cache_output_name)
                return pd.read_feather(self.outpath), pd.read_feather(self.secondary_outpath)
            else:
                return pd.read_feather(self.outpath)

        # two types of run, accepting current working dataframe and not.
        if isinstance(self.working_dataframe, pd.DataFrame):
            if self.secondary_cache_output_name:
                self.secondary_outpath: Path = self.cached_output_path / Path(self.secondary_cache_output_name)
                self.working_dataframe, self.secondary_working_dataframe = self.callable_function(self.data, self.working_dataframe)  # type: ignore
            else:
                self.working_dataframe = self.callable_function(self.data, self.working_dataframe)  # type: ignore
        else:
            self.working_dataframe = self.callable_function(self.data)  # type: ignore
        logger.info(f"Writing output to {self.outpath}")

        # always will rewrite the dataframe
        if isinstance(self.secondary_working_dataframe, pd.DataFrame):
            logger.info(f"Writing output to {self.secondary_outpath}")
            self.secondary_working_dataframe.to_feather(self.secondary_outpath)  # type: ignore
            self.working_dataframe.to_feather(self.outpath)  # type: ignore
            return self.working_dataframe, self.secondary_working_dataframe
        else:
            self.working_dataframe.to_feather(self.outpath)
            return self.working_dataframe  # type: ignore


def run_sequence_analysis(data: Data, outpath: Path, resume: bool, skip_csv: bool) -> None:
    """
    Run the sequence analysis pipeline
    """

    # Find module - Break fastq into dataframes
    find_module = ModuleRunner(
        cached_output_name="find.feather", cached_output_path=outpath, data=data, callable_function=find, resume=resume
    )
    working_dataframe: pd.DataFrame = find_module.run()
    # Model module - Model the data
    model_module = ModuleRunner(
        cached_output_name="model.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=model,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = model_module.run()
    logger.info(f"Finished validating {len(working_dataframe)} sequence_id data into model")

    # Split module - split the data
    split_module = ModuleRunner(
        cached_output_name="split.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=split,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = split_module.run()
    logger.info(f"Finished split {len(working_dataframe)} dataframe")

    # Correct module - correct the data
    correct_module = ModuleRunner(
        cached_output_name="correct.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=correct,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = correct_module.run()
    logger.info(f"Finished correcting {len(working_dataframe)} dataframe")

    # annotation module - annotate the data
    annotate_module = ModuleRunner(
        cached_output_name="annotate.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=annotate,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = annotate_module.run()
    logger.info(f"Finished annotate {len(working_dataframe)} dataframe")

    # join module - join the data
    join_module = ModuleRunner(
        cached_output_name="join.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=join,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = join_module.run()

    # ngs module - correct with ngs where applicable
    ngs_module = ModuleRunner(
        cached_output_name="ngs.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=ngs,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = ngs_module.run()

    # tag module - tag the data with possible control contamination
    tag_module = ModuleRunner(
        cached_output_name="tag.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=tag,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = tag_module.run()

    # swap module - correct the sample swap
    swap_module = ModuleRunner(
        cached_output_name="swap.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=swap,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe: pd.DataFrame = swap_module.run()

    # pair module - make heavy and light pairs
    pair_module = ModuleRunner(
        cached_output_name="pair.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=pair,
        working_dataframe=working_dataframe,
        secondary_cache_output_name="pairing_candidates.feather",
        resume=resume,
    )

    # we don't do anything with unpaired dataframe so we don't need it...we write it out in the module class
    working_dataframe, _ = pair_module.run()

    # unblind module - add dose and vaccine group
    unblind_module = ModuleRunner(
        cached_output_name="unblind.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=unblind,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe = unblind_module.run()

    # personalize module - rerun annotation with personalized alleles
    personalize_module = ModuleRunner(
        cached_output_name="personalize.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=personalize,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe = personalize_module.run()

    # mutation module - annotate the data with mutations
    mutation_module = ModuleRunner(
        cached_output_name="mutation.feather",
        cached_output_path=outpath,
        data=data,
        callable_function=mutation,
        working_dataframe=working_dataframe,
        resume=resume,
    )
    working_dataframe = mutation_module.run()

    # Complete pipeline
    logger.info(f"Pipeline complete, moving files to {outpath.parent}")
    for file in outpath.glob("*.feather"):
        full_path = outpath.parent / file.name
        logger.info(f"Moving {file} to {full_path}")
        # uncomment after testing
        shutil.move(file, full_path)
        if not skip_csv:
            pd.read_feather(full_path).to_csv(full_path.with_suffix(".csv.gz"))
            logger.info(f"Making {full_path} to {full_path.with_suffix('.csv.gz')}")

    # uncomment after testing
    # shutil.rmtree(outpath)
    logger.info(f"Results in {outpath.parent}")
