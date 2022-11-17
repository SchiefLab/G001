# -*- coding: utf-8 -*-
from __future__ import annotations

# stdlib
import logging
from pathlib import Path

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


@main.group()
def figures():
    """
    Generate for figures.
    """
    sns.set_context("paper", font_scale=1)
    sns.set_style("ticks")
    sns.set_style({"font.family": "Arial"})
    np.random.seed(1000)
    pass


@figures.command("fig1")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure1",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_1(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_flow_frequencies(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 1 generated to {outpath}.png")


if __name__ == "__main__":
    main()
