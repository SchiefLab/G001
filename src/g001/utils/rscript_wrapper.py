from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import subprocess
from subprocess import CalledProcessError

from rich.console import Console


@dataclass
class RScript:
    """RScript Wrapper"""

    rpath = Path(__file__).parent.parent / "R"
    sites = ["FHCRC", "VRC"]
    verbose: bool = False
    console = Console(record=False)

    def __run_cmd(self, cmd: list[str]) -> None:
        if self.verbose:
            self.console.print(f"Running command:\n{' '.join(cmd)}")
        try:
            proc = subprocess.run(cmd, capture_output=True, encoding="utf-8")
        except CalledProcessError:
            self.console.print_exception(show_locals=True)
        else:
            # ignore ugly R warnings from cytoqc :: even they are still trying to figure out how to remove them.
            stderr = "\n".join([e for e in proc.stderr.replace(":\n", ":").split("\n") if not e.startswith("Warning")])
            # ignore ugly R column renames
            stderr = "\n".join([e for e in stderr.split("\n") if not e.startswith("New names:")])
            if stderr:
                self.console.print(stderr)
            self.console.print(proc.stdout)

    def flow_processing(
        self,
        site: str,
        manifest: Path | str,
        flow_input_dir: Path | str,
        flow_output_dir: Path | str,
        force_overwrite_output_dir: bool | str = True,
    ) -> None:
        """Flow Processing for FHCRC or VRC

        Parameters
        ----------
        site : str
            FHCRC or VRC
        flow_input_dir : Path | str
            Path to the directory containing the flow data
        flow_output_dir : Path | str
            Path to the directory where the flow data will be written
        manifest : Path | str | None, optional
            Manifest linked to site, by default flow_input/{site.lower()}/{site}_Flow_Manifest.csv
        force_overwrite_output_dir : bool | str, optional
            overwrite output folder, by default True

        Raises
        ------
        ValueError
            site must be FHCRC or VRC
        ValueError
            flow_input_dir must be a Path or str and must exist
        ValueError
            flow_output_dir must be a Path or str and must exist
        ValueError
            manifest must be a Path or str
        """
        site = site.upper()
        if site not in self.sites:
            raise ValueError(f"site {site} not found, must be one of {self.sites}")

        flow_input_dir = Path(flow_input_dir)
        if not flow_input_dir.exists():
            raise ValueError(f"flow_input_dir {flow_input_dir} does not exist")

        flow_output_dir = Path(flow_output_dir)
        if not flow_output_dir.parent.exists():
            raise ValueError(f"flow_output_dir parent {flow_output_dir.parent} does not exist")

        manifest = (
            Path(manifest)
            if manifest
            else Path(__file__).parent.parent / f"flow_input/{site.lower()}/{site}_Flow_Manifest.csv"
        )

        if not manifest.exists():
            raise ValueError(f"manifest {manifest} does not exist")

        force_overwrite_output_dir = "yes" if force_overwrite_output_dir in ["yes", True] else "no"

        cmd = [
            "Rscript",
            str(self.rpath / "Flow_Processing.R"),
            site,
            str(manifest),
            str(flow_input_dir),
            str(force_overwrite_output_dir),
            str(flow_output_dir),
        ]
        self.__run_cmd(cmd)

    def collate_flow(self, fhcrc_manifest: Path | str, vrc_manifest: Path | str, flow_output_dir: Path | str) -> None:
        """Collate Flow

        Parameters
        ----------
        collate_output_dir : Path | str
            Path to the directory where the collated flow data will be written
        """
        flow_output_dir = Path(flow_output_dir)
        if not flow_output_dir.exists():
            raise ValueError(f"collate_output_dir {flow_output_dir} does not exist")
        cmd = [
            "Rscript",
            str(self.rpath / "Collate_Flow_Data.R"),
            str(fhcrc_manifest),
            str(vrc_manifest),
            str(flow_output_dir),
        ]
        self.__run_cmd(cmd)
