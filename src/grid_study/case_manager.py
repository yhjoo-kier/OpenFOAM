"""OpenFOAM case management for grid study."""

import shutil
import subprocess
from pathlib import Path
from typing import Optional

from .config import GridStudyConfig, MeshLevel


class CaseManager:
    """Manages OpenFOAM case creation and configuration for grid studies."""

    def __init__(self, config: GridStudyConfig):
        self.config = config
        self.cases_dir = config.output_dir / "cases"

    def setup_study_directory(self) -> None:
        """Create study output directory structure."""
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        self.cases_dir.mkdir(exist_ok=True)
        (self.config.output_dir / "results").mkdir(exist_ok=True)
        (self.config.output_dir / "logs").mkdir(exist_ok=True)

    def create_case(self, mesh_level: MeshLevel) -> Path:
        """
        Create a new case for a specific mesh level.

        Args:
            mesh_level: Mesh refinement configuration

        Returns:
            Path to the created case directory
        """
        case_name = f"{self.config.study_name}_{mesh_level.name}"
        case_dir = self.cases_dir / case_name

        # Remove existing case if present
        if case_dir.exists():
            shutil.rmtree(case_dir)

        # Copy base case structure (0, constant, system)
        self._copy_case_template(case_dir)

        return case_dir

    def _copy_case_template(self, target_dir: Path) -> None:
        """Copy essential OpenFOAM case files from base case."""
        base = self.config.base_case_path

        # Create target directory
        target_dir.mkdir(parents=True, exist_ok=True)

        # Copy 0/ directory (initial conditions)
        src_0 = base / "0"
        if src_0.exists():
            shutil.copytree(src_0, target_dir / "0")

        # Copy constant/ directory (but not polyMesh - will be regenerated)
        src_constant = base / "constant"
        if src_constant.exists():
            dst_constant = target_dir / "constant"
            dst_constant.mkdir(exist_ok=True)

            # Copy files and subdirectories except polyMesh
            for item in src_constant.iterdir():
                if item.name == "polyMesh":
                    continue
                if item.is_dir():
                    # For region directories (fluid, solid), copy but skip polyMesh
                    dst_region = dst_constant / item.name
                    dst_region.mkdir(exist_ok=True)
                    for subitem in item.iterdir():
                        if subitem.name == "polyMesh":
                            continue
                        if subitem.is_file():
                            shutil.copy2(subitem, dst_region)
                        elif subitem.is_dir():
                            shutil.copytree(subitem, dst_region / subitem.name)
                else:
                    shutil.copy2(item, dst_constant)

        # Copy system/ directory
        src_system = base / "system"
        if src_system.exists():
            shutil.copytree(src_system, target_dir / "system")

    def get_case_path(self, mesh_level: MeshLevel) -> Path:
        """Get the path to a case for a specific mesh level."""
        case_name = f"{self.config.study_name}_{mesh_level.name}"
        return self.cases_dir / case_name

    def clean_case(self, case_dir: Path) -> None:
        """Remove time directories and logs from a case."""
        # Remove time directories (except 0)
        for item in case_dir.iterdir():
            if item.is_dir():
                try:
                    time_val = float(item.name)
                    if time_val > 0:
                        shutil.rmtree(item)
                except ValueError:
                    pass

        # Remove processor directories
        for proc_dir in case_dir.glob("processor*"):
            shutil.rmtree(proc_dir)

        # Remove log files
        for log_file in case_dir.glob("log.*"):
            log_file.unlink()

    def update_control_dict(self, case_dir: Path, end_time: Optional[int] = None) -> None:
        """
        Update controlDict for the simulation.

        Args:
            case_dir: Path to case directory
            end_time: Number of iterations (for steady-state)
        """
        if end_time is None:
            end_time = self.config.num_iterations

        control_dict_path = case_dir / "system" / "controlDict"

        if not control_dict_path.exists():
            return

        content = control_dict_path.read_text()

        # Update endTime
        import re
        content = re.sub(
            r"endTime\s+\d+;",
            f"endTime         {end_time};",
            content
        )

        control_dict_path.write_text(content)
