# -*- coding: utf-8 -*-
from __future__ import annotations

# stdlib
import logging
from pathlib import Path
import subprocess

# third party
import click
import seaborn as sns
import numpy as np

# G001 package
from g001.data import Data
from g001.sequence_pipeline.main import run_sequence_analysis
from g001.figures import plot_flow_frequencies


@click.group()
@click.pass_context
def main(ctx: click.Context) -> None:
    """
    Generate G001
    """
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%m-%d-%Y %I:%M%p",
        format="%(name)s:%(levelname)s:%(asctime)s--> %(message)s",
    )
    data = Data(Path("data"))

    # make empty context object
    ctx.obj = {"data": data}
    pass


@main.command("sequence_analysis")
@click.pass_context
@click.option(
    "--data-path",
    "-d",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=True,
    help="Path to the sequence data directory",
)
@click.option(
    "--outpath",
    "-o",
    type=click.Path(dir_okay=True, readable=True, resolve_path=True),
    required=True,
    help="Outpath to sequence_analysis",
    default="sequence_analysis_output",
)
@click.option("--resume", "-r", is_flag=True, help="""Resume pipeline from .cache path""")
@click.option("--skip-csv", is_flag=True, help="""Skip copying feathers to .csv in final path""")
def click_run_sa(ctx: click.Context, data_path: str | Path, outpath: Path, resume: bool, skip_csv: bool) -> None:
    """
    Run sequence analysis
    """
    logger = logging.getLogger("G001")
    data = Data(Path(data_path))
    logger.info(
        f"Running complete sequence analysis pipeline using data: {data.sequence_data_paths.base_sequence_path}"
    )

    # create output directory and cache path
    outpath = Path(outpath)
    cache_path = outpath / Path(".cache")

    # if they dont exist, create them
    if not outpath.exists():
        logger.info(f"Creating output directory: {outpath}")
        outpath.mkdir(parents=True)
    if not cache_path.exists():
        logger.info(f"Creating cache directory: {cache_path}")
        cache_path.mkdir(parents=True)

    # run sequene analysis pipeline
    run_sequence_analysis(data, cache_path, resume, skip_csv)


@main.command("fhcrc")
@click.pass_context
@click.option(
    "--flow-output",
    "-o",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=False,
    default="flow_output",
    help="Path to create the flow_output directory",
)
def fhcrc(ctx: click.Context, flow_output: Path | str) -> None:
    """
    Run for FHCRC
    """
    cmd = [
        "Rscript",
        "src/g001/R/Flow_Processing.R",
        "FHCRC",
        "flow_input/fhcrc/FHCRC_Flow_Manifest.csv",
        "flow_input/fhcrc/",
        "yes",
        flow_output,
    ]
    stdout = subprocess.run(cmd, capture_output=True)
    print(stdout.stdout.decode("utf-8"))


@main.command("vrc")
@click.pass_context
def vcr(ctx: click.Context) -> None:
    """
    Run for VRC
    """
    cmd = [
        "Rscript",
        "src/g001/R/Flow_Processing.R",
        "VRC",
        "flow_input/vrc/VRC_Flow_Manifest.csv",
        "flow_input/vrc/",
        "yes",
        "flow_output",
    ]
    stdout = subprocess.run(cmd, capture_output=True)
    print(stdout.stdout.decode("utf-8"))


@main.command("collate")
@click.pass_context
def collate(ctx: click.Context) -> None:
    """
    Collation of Flow Data
    """
    cmd = [
        "Rscript",
        "src/g001/R/Collate_Flow_Data.R",
        "data/flow/processed_flow/",
    ]
    stdout = subprocess.run(cmd, capture_output=True)
    # TODO: collate doesnt seem to have an stdout
    print(stdout.stdout.decode("utf-8"))


if __name__ == "__main__":
    main()
