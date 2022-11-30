from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import subprocess
from subprocess import CalledProcessError

from rich.console import Console


class RProcessError(CalledProcessError):
    pass


@dataclass
class RScript:
    """RScript Wrapper"""

    rpath = Path(__file__).parent.parent / "R"
    sites = ["fhcrc", "vrc"]
    verbose: bool = False
    console = Console(record=False)

    def __run_cmd(self, cmd: list[str], cap_output: bool = True) -> None:
        self.console.print(f"Running command:\n{' '.join(cmd)}")
        try:
            if not cap_output:
                proc = subprocess.run(cmd, stderr=subprocess.STDOUT)
            else:
                proc = subprocess.run(cmd, capture_output=True, encoding="utf-8")
        except CalledProcessError:
            self.console.print_exception(show_locals=True)
        else:
            # ignore ugly R warnings from cytoqc :: even they are still trying to figure out how to remove them.
            if proc.stderr:
                # stderr: str = "\n".join(
                #     [e for e in str(proc.stderr).replace(":\n", ":").split("\n") if not e.startswith("Warning")]
                # )
                # # ignore ugly R column renames
                if str(proc.stderr).startswith("Warning message:"):
                    self.console.print(proc.stderr)
                else:
                    raise RProcessError(proc.returncode, " ".join(cmd), proc.stderr)
            if proc.stdout:
                self.console.print(f"Output from R:\n{proc.stdout}")

    def flow_processing(
        self,
        site: str,
        manifest: Path | str,
        flow_input_dir: Path | str,
        flow_output_dir: Path | str,
        force_overwrite_output_dir: bool,
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
        site = site.lower()
        if site not in self.sites:
            raise ValueError(f"site {site} not found, must be one of {self.sites}")

        flow_input_dir = Path(flow_input_dir)
        if not flow_input_dir.exists():
            raise ValueError(f"flow_input_dir {flow_input_dir} does not exist")

        flow_output_dir = Path(flow_output_dir)
        if not flow_output_dir.parent.exists():
            raise ValueError(f"flow_output_dir parent {flow_output_dir.parent} does not exist")

        if not Path(manifest).exists():
            raise ValueError(f"manifest {manifest} does not exist")

        force_overwrite_output_dir: str = "yes" if force_overwrite_output_dir else "no"

        cmd = [
            "Rscript",
            str(self.rpath / "Flow_Processing.R"),
            site,
            str(manifest),
            str(flow_input_dir),
            str(force_overwrite_output_dir),
            str(flow_output_dir),
        ]
        print("Running...\n", " ".join(cmd))

        # I don't want to cap the output here since people will think the script is not working
        self.__run_cmd(cmd, cap_output=False)

    def collate_flow(
        self,
        fhcrc_manifest: Path,
        vrc_manifest: Path,
        flow_processed_dir: Path,
        collated_output_dir: Path,
        swap_file: Path,
    ) -> None:
        """Collate Flow

        Parameters
        ----------
        fhcrc_manifest : Path | str
            Path to the FHCRC manifest
        vrc_manifest : Path | str
            Path to the VRC manifest
        flow_processed_dir : Path | str
            Path to the directory containing the flow data for FHCRC and VRC
        collate_output_dir : Path | str
            Path to the directory where the collated flow data will be written
        """
        flow_processed_dir = Path(flow_processed_dir)

        # these generally shouldn't be called since click will handle theem
        if not flow_processed_dir.exists():
            raise ValueError(f"flow_processed_dir {flow_processed_dir} does not exist")
        if not fhcrc_manifest.exists():
            raise ValueError(f"fhcrc_manifest {fhcrc_manifest} does not exist")
        if not vrc_manifest.exists():
            raise ValueError(f"vrc_manifest {vrc_manifest} does not exist")
        if not swap_file.exists():
            raise ValueError(f"swap_file {swap_file} does not exist")
        cmd = [
            "Rscript",
            str(self.rpath / "Collate_Flow_Data.R"),
            str(fhcrc_manifest),
            str(vrc_manifest),
            str(flow_processed_dir),
            str(swap_file),
            str(collated_output_dir),
        ]
        print("Running...\n", " ".join(cmd))
        self.__run_cmd(cmd)

    def combine_flow_and_sequence(
        self,
        fhcrc_manifest: Path,
        vrc_manifest: Path,
        sequence_dir: Path,
        collated_dir: Path,
        combined_output_dir: Path,
        treatment_path: Path,
    ) -> None:
        """Combine Flow Data with Sequence Data

        Parameters
        ----------
        fhcrc_manifest : Path | str
            Path to the FHCRC manifest
        vrc_manifest : Path | str
            Path to the VRC manifest
        sequence_dir : Path | str
            Path to the directory containing the sequence data
        collate_dir : Path | str
            Path to the directory where the collated flow data was be written
        combined_output_dir : Path | str
            Path to the directory where the combined data will be written
        """
        if not fhcrc_manifest.exists():
            raise ValueError(f"fhcrc_manifest {fhcrc_manifest} does not exist")
        if not vrc_manifest.exists():
            raise ValueError(f"vrc_manifest {vrc_manifest} does not exist")
        if not collated_dir.exists():
            raise ValueError(f"collated folder {collated_dir} does not exist")
        if not sequence_dir.exists():
            raise ValueError(f"sequence folder {sequence_dir} does not exist")
        cmd = [
            "Rscript",
            str(self.rpath / "Combine_Flow_and_Seq_Results.R"),
            str(fhcrc_manifest),
            str(vrc_manifest),
            str(sequence_dir),
            str(collated_dir),
            str(combined_output_dir),
            str(treatment_path),
        ]
        print("Running...\n", " ".join(cmd))
        self.__run_cmd(cmd)
