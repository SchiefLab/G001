# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
from __future__ import annotations
from itertools import product

# stdlib
import logging
from pathlib import Path
import subprocess

# third party
import click
import seaborn as sns
import numpy as np
from g001.figures.binding import plot_boosting_binding, plot_gt8_binding
from g001.figures.features import plot_sequence_features
from g001.figures.polyclonality import plot_class_clonality

# G001 package
from g001.data import Data
from g001.figures.frequency import (
    plot_count_frequency_response,
    plot_flow_frequencies,
    plot_vrc01_among_epitope_specific,
)
from g001.figures.overview import plot_overview
from g001.figures.somatic import plot_somatic_mutation_frequencies, plot_somatic_mutation_frequencies_violin
from g001.sequence_pipeline.main import run_sequence_analysis
from g001.utils import RScript


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


@main.command("process-flow")
@click.pass_context
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    required=False,
    default=False,
    help="Prints logs to console",
)
@click.argument(
    "gate",
    type=str,
    required=True,
)
@click.option(
    "--manifest",
    "-m",
    type=click.Path(dir_okay=False, readable=True, exists=None, resolve_path=True),
    required=False,
    default=Path(__file__).parent.parent.parent / "flow_input/{gate_lower}/{gate_upper}_Flow_Manifest.csv",
    show_default=True,
    help="Path to read the flow manifest for specified gate",
)
@click.option(
    "--flow-input-dir",
    "-i",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=False,
    default=Path(__file__).parent.parent.parent / f"flow_input",
    show_default=True,
    help="Path to read the flow_input directory",
)
@click.option(
    "--flow-output-dir",
    "-o",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=False,
    default=Path(__file__).parent.parent.parent / f"flow_output",
    show_default=True,
    help="Path to create the flow_output directory",
)
@click.option("-f", "--force_overwrite_output_dir", default=True, is_flag=True, help="Force overwrite output directory")
def fhcrc(
    ctx: click.Context,
    verbose: bool,
    gate: str,
    manifest: Path | str,
    flow_input_dir: Path | str,
    flow_output_dir: Path | str,
    force_overwrite_output_dir: str,
) -> None:
    """
    GATE: FHCRC or VRC
    """
    manifest = Path(str(manifest).format(gate_lower=gate.lower(), gate_upper=gate.upper()))
    RScript(verbose=verbose).flow_processing(
        gate=gate,
        manifest=manifest,
        flow_input_dir=flow_input_dir,
        flow_output_dir=flow_output_dir,
        force_overwrite_output_dir=force_overwrite_output_dir,
    )


@main.command("collate")
@click.pass_context
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    required=False,
    default=False,
    help="Prints final command that was run to console",
)
@click.option(
    "--collate-output-dir",
    "-o",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=False),
    required=False,
    default=Path(__file__).parent.parent.parent / "data/flow/processed_flow",
    show_default=True,
    help="Path to create the collated output directory",
)
def collate(
    ctx: click.Context,
    verbose: bool,
    collate_output_dir: Path | str,
) -> None:
    """ """
    RScript(verbose=verbose).collate_flow(
        collate_output_dir=collate_output_dir,
    )


if __name__ == "__main__":
    main()
