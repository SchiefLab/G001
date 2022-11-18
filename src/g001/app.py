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
@click.option(
    "--data-path", "-d", type=click.Path(exists=True, dir_okay=True), default=Path("./data"), show_default=True
)
def main(ctx: click.Context, data_path: Path) -> None:
    """
    Generate G001
    """
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%m-%d-%Y %I:%M%p",
        format="%(name)s:%(levelname)s:%(asctime)s--> %(message)s",
    )
    data = Data(data_path)

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
@click.option(
    "--site",
    "-s",
    required=True,
    type=click.Choice(["FHCRC", "VRC"]),
)
@click.option(
    "--manifest",
    "-m",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=True,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--flow-input-dir",
    "-i",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=False,
    default="flow_input",
    show_default=True,
    help="Path to read the flow_input directory",
)
@click.option(
    "--flow-output-dir",
    "-o",
    type=click.Path(dir_okay=True, readable=True, resolve_path=True),
    required=False,
    default=Path("flow_output"),
    show_default=True,
    help="Path to create the flow_output directory",
)
@click.option("-f", "--force_overwrite_output_dir", default=True, is_flag=True, help="Force overwrite output directory")
def process_flow(
    ctx: click.Context,
    verbose: bool,
    site: str,
    manifest: Path | str,
    flow_input_dir: Path | str,
    flow_output_dir: Path | str,
    force_overwrite_output_dir: str,
) -> None:
    """
    Process flow data
    """
    if not Path(flow_output_dir).exists():
        Path(flow_output_dir).mkdir()
    RScript(verbose=verbose).flow_processing(
        gate=site,
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
    "--flow-output-dir",
    "-f",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=False),
    required=False,
    default=None,
    help="Path to the processed flow directory",
)
def collate(
    ctx: click.Context,
    verbose: bool,
    flow_output_dir: Path | str,
) -> None:
    """ """
    if not flow_output_dir:
        data = ctx.obj["data"]
        flow_output_dir = data.get_processed_flow_paths()
    RScript(verbose=verbose).collate_flow(
        collate_output_dir=flow_output_dir,
    )


@main.group()
def figures():
    """
    Generate for figures.
    """
    sns.set_context("paper", font_scale=1)
    sns.set_style("ticks")
    sns.set_style({"font.family": "Arial"})
    import warnings

    warnings.simplefilter("ignore")
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


@figures.command("fig2")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure2",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_2(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_count_frequency_response(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 2A-C generated to {outpath}.png")

    figure_d = plot_overview(data)
    figure_d.savefig(outpath + "_d.png", dpi=300)
    click.echo(f"Figure 2D generated to {outpath}_d.png")


@figures.command("fig3")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure3",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_3(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_class_clonality(data, True)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 3 generated to {outpath}.png")


@figures.command("fig4")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure4",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_4(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_vrc01_among_epitope_specific(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 4 generated to {outpath}.png")


@figures.command("fig5")
@click.pass_context
@click.option(
    "--aa-outpath",
    "-o",
    default="figures/figure5",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
@click.option(
    "--nt-outpath",
    "-n",
    default="figures/figureS20",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure S5",
    show_default=True,
)
def plot_figure_5_and_s26(ctx: click.Context, aa_outpath: str, nt_outpath: str) -> None:
    data: Data = ctx.obj["data"]
    plot_params = data.plot_parameters
    unblind_df = data.get_unblinded_sequences()
    pallete = plot_params.get_pallete()
    vrc01_dose_palette = {
        "Low Dose": pallete["Low Dose"]["vrc01_class"],
        "High Dose": pallete["High Dose"]["vrc01_class"],
    }
    nonvrc01_dose_palette = {
        "Low Dose": pallete["Low Dose"]["nonvrc01_class"],
        "High Dose": pallete["High Dose"]["nonvrc01_class"],
    }

    for ab_type, molecule in list(product(["vrc01", "nonvrc01"], ["nt", "aa"])):
        if ab_type == "vrc01":
            sup_title = "VRC01-class BCRs"
            palette = vrc01_dose_palette
        else:
            sup_title = "non-VRC01-class BCRs"
            palette = nonvrc01_dose_palette
        if molecule == "nt":
            outpath_ = nt_outpath
        else:
            outpath_ = aa_outpath
        fig = plot_somatic_mutation_frequencies(unblind_df, palette, ab_type, molecule, sup_title)
        new_outpath = outpath_ + f"_{ab_type}_{molecule}.png"
        fig.savefig(new_outpath, dpi=300)
        if molecule == "nt":
            _figure = "S20"
        else:
            _figure = "2"
        click.echo(f"Figure {_figure} generated to {new_outpath}")
        fig = plot_somatic_mutation_frequencies_violin(unblind_df, palette, ab_type, molecule, sup_title)
        new_outpath = outpath_ + f"_{ab_type}_{molecule}_violin.png"
        click.echo(f"Figure {_figure} generated to {new_outpath}")
        fig.savefig(new_outpath, dpi=300)


@figures.command("fig6")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure6",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_6(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_sequence_features(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 6 generated to {outpath}.png")


@figures.command("fig7")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure7",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_6(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_gt8_binding(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 7 generated to {outpath}.png")


@figures.command("fig8")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure8",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_8(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_boosting_binding(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 8 generated to {outpath}.png")


if __name__ == "__main__":
    main()
