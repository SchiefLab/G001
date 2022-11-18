from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import subprocess
from subprocess import CalledProcessError

from rich.console import Console


@dataclass
class RScript:
    """RScript Wrapper"""

    rpath = Path(__file__).parent / "R"
    gates = ["FHCRC", "VRC"]
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
        gate: str,
        manifest: Path | str,
        flow_input_dir: Path | str,
        flow_output_dir: Path | str,
        force_overwrite_output_dir: bool | str = True,
    ) -> None:
        """Flow Processing for FHCRC or VRC

        Parameters
        ----------
        gate : str
            FHCRC or VRC
        flow_input_dir : Path | str
            Path to the directory containing the flow data
        flow_output_dir : Path | str
            Path to the directory where the flow data will be written
        manifest : Path | str | None, optional
            Manifest linked to gate, by default flow_input/{gate.lower()}/{gate}_Flow_Manifest.csv
        force_overwrite_output_dir : bool | str, optional
            overwrite output folder, by default True

        Raises
        ------
        ValueError
            gate must be FHCRC or VRC
        ValueError
            flow_input_dir must be a Path or str and must exist
        ValueError
            flow_output_dir must be a Path or str and must exist
        ValueError
            manifest must be a Path or str
        """
        gate = gate.upper()
        if gate not in self.gates:
            raise ValueError(f"gate {gate} not found, must be one of {self.gates}")

        flow_input_dir = Path(flow_input_dir)
        if not flow_input_dir.exists():
            raise ValueError(f"flow_input_dir {flow_input_dir} does not exist")

        flow_output_dir = Path(flow_output_dir)
        if not flow_output_dir.parent.exists():
            raise ValueError(f"flow_output_dir parent {flow_output_dir.parent} does not exist")

        manifest = (
            Path(manifest)
            if manifest
            else Path(__file__).parent.parent / f"flow_input/{gate.lower()}/{gate}_Flow_Manifest.csv"
        )

        if not manifest.exists():
            raise ValueError(f"manifest {manifest} does not exist")

        force_overwrite_output_dir = "yes" if force_overwrite_output_dir in ["yes", True] else "no"

        cmd = [
            "Rscript",
            str(self.rpath / "Flow_Processing.R"),
            gate,
            str(manifest),
            str(flow_input_dir),
            str(force_overwrite_output_dir),
            str(flow_output_dir),
        ]
        self.__run_cmd(cmd)

    def collate_flow(self, collate_output_dir: Path | str) -> None:
        """Collate Flow

        Parameters
        ----------
        collate_output_dir : Path | str
            Path to the directory where the collated flow data will be written
        """
        collate_output_dir = Path(collate_output_dir)
        if not collate_output_dir.exists():
            raise ValueError(f"collate_output_dir {collate_output_dir} does not exist")
        cmd = [
            "Rscript",
            str(self.rpath / "Collate_Flow_Data.R"),
            "data/flow/processed_flow",  # BUG: hard coded since
            # str(collate_output_dir),
        ]
        self.__run_cmd(cmd)
