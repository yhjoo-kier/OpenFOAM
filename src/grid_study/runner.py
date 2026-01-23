"""Simulation runner for grid study."""

import subprocess
import time
from pathlib import Path
from typing import Optional
from dataclasses import dataclass


@dataclass
class RunResult:
    """Result of a simulation run."""
    success: bool
    elapsed_time: float  # seconds
    log_file: Path
    error_message: Optional[str] = None
    final_iteration: Optional[int] = None


class SimulationRunner:
    """Runs OpenFOAM simulations."""

    def __init__(self, solver_app: str = "chtMultiRegionSimpleFoam"):
        self.solver_app = solver_app

    def run(
        self,
        case_dir: Path,
        log_dir: Optional[Path] = None,
        case_name: Optional[str] = None,
    ) -> RunResult:
        """
        Run simulation for a case.

        Args:
            case_dir: Path to OpenFOAM case directory
            log_dir: Directory to store log files (default: case_dir)
            case_name: Name for log file (default: case directory name)

        Returns:
            RunResult with execution details
        """
        if log_dir is None:
            log_dir = case_dir

        if case_name is None:
            case_name = case_dir.name

        log_file = log_dir / f"log.{self.solver_app}_{case_name}"

        start_time = time.time()

        try:
            # Run solver - use absolute path to avoid path doubling issues
            case_dir_abs = case_dir.resolve()
            with open(log_file, "w") as f:
                result = subprocess.run(
                    [self.solver_app, "-case", str(case_dir_abs)],
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    timeout=7200,  # 2 hour timeout
                )

            elapsed_time = time.time() - start_time

            if result.returncode != 0:
                return RunResult(
                    success=False,
                    elapsed_time=elapsed_time,
                    log_file=log_file,
                    error_message=f"Solver exited with code {result.returncode}",
                )

            # Check for convergence/completion
            final_iter = self._get_final_iteration(log_file)

            return RunResult(
                success=True,
                elapsed_time=elapsed_time,
                log_file=log_file,
                final_iteration=final_iter,
            )

        except subprocess.TimeoutExpired:
            elapsed_time = time.time() - start_time
            return RunResult(
                success=False,
                elapsed_time=elapsed_time,
                log_file=log_file,
                error_message="Simulation timed out",
            )
        except Exception as e:
            elapsed_time = time.time() - start_time
            return RunResult(
                success=False,
                elapsed_time=elapsed_time,
                log_file=log_file,
                error_message=str(e),
            )

    def _get_final_iteration(self, log_file: Path) -> Optional[int]:
        """Extract final iteration number from log file."""
        try:
            content = log_file.read_text()
            # Look for "Time = <number>" patterns
            import re
            matches = re.findall(r"Time = (\d+)", content)
            if matches:
                return int(matches[-1])
        except Exception:
            pass
        return None


def run_mesh_conversion(case_dir: Path, msh_file: Path, log_dir: Path) -> bool:
    """
    Run gmshToFoam and splitMeshRegions.

    Args:
        case_dir: OpenFOAM case directory
        msh_file: Path to Gmsh .msh file
        log_dir: Directory for log files

    Returns:
        True if successful
    """
    # gmshToFoam
    log_gmsh = log_dir / f"log.gmshToFoam_{case_dir.name}"
    # Use absolute path for msh_file since cwd will be case_dir
    msh_file_abs = msh_file.resolve()
    with open(log_gmsh, "w") as f:
        result = subprocess.run(
            ["gmshToFoam", str(msh_file_abs)],
            stdout=f,
            stderr=subprocess.STDOUT,
            cwd=case_dir,
        )

    if result.returncode != 0:
        print(f"gmshToFoam failed. See {log_gmsh}")
        return False

    # splitMeshRegions
    log_split = log_dir / f"log.splitMeshRegions_{case_dir.name}"
    with open(log_split, "w") as f:
        result = subprocess.run(
            ["splitMeshRegions", "-cellZones", "-overwrite"],
            stdout=f,
            stderr=subprocess.STDOUT,
            cwd=case_dir,
        )

    if result.returncode != 0:
        print(f"splitMeshRegions failed. See {log_split}")
        return False

    return True
