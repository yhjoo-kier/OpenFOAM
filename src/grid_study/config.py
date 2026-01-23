"""Configuration schema for grid independence study."""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from pathlib import Path
import json


@dataclass
class MeshLevel:
    """Single mesh refinement level configuration."""
    name: str
    mesh_factor: float  # Global mesh size multiplier (smaller = finer)
    bl_first_height: float  # Boundary layer first cell height [m]
    bl_growth_ratio: float  # Boundary layer growth ratio
    bl_num_layers: int  # Number of boundary layer cells


@dataclass
class GridStudyConfig:
    """Configuration for grid independence study."""

    # Case identification
    base_case_path: Path
    study_name: str = "grid_study"
    output_dir: Path = field(default_factory=lambda: Path("grid_study_output"))

    # Mesh refinement levels
    mesh_levels: List[MeshLevel] = field(default_factory=list)

    # Monitoring metric
    metric_name: str = "T_base_avg"  # Name for reporting
    metric_patch: str = "heat_source"  # OpenFOAM patch name
    metric_field: str = "T"  # Field to extract
    metric_region: str = "solid"  # Region (for multi-region)

    # Convergence criterion
    convergence_threshold: float = 1.0  # Percentage change threshold

    # Solver settings
    solver_app: str = "chtMultiRegionSimpleFoam"
    num_iterations: int = 2000

    # Geometry script
    geometry_script: Optional[Path] = None

    # Adaptive refinement mode
    adaptive_mode: bool = False  # Continue refining until convergence
    refinement_ratio: float = 0.7  # mesh_factor multiplier for each new level
    bl_refinement_ratio: float = 0.7  # boundary layer refinement ratio

    # Stopping criteria for adaptive mode
    max_cells: int = 2_000_000  # Maximum total cells
    max_levels: int = 10  # Maximum number of refinement levels
    max_runtime_per_level: Optional[int] = None  # Max seconds per level (None=unlimited)

    def __post_init__(self):
        """Initialize default mesh levels if not provided."""
        if not self.mesh_levels:
            self.mesh_levels = self.get_default_mesh_levels()

        # Convert paths
        if isinstance(self.base_case_path, str):
            self.base_case_path = Path(self.base_case_path)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)
        if isinstance(self.geometry_script, str):
            self.geometry_script = Path(self.geometry_script)

    @staticmethod
    def get_default_mesh_levels() -> List[MeshLevel]:
        """Return default 4-level mesh refinement configuration."""
        return [
            MeshLevel(
                name="L1_coarse",
                mesh_factor=2.0,
                bl_first_height=0.0005,  # 0.5mm
                bl_growth_ratio=1.2,
                bl_num_layers=3,
            ),
            MeshLevel(
                name="L2_medium",
                mesh_factor=1.0,
                bl_first_height=0.00025,  # 0.25mm
                bl_growth_ratio=1.2,
                bl_num_layers=4,
            ),
            MeshLevel(
                name="L3_fine",
                mesh_factor=0.7,
                bl_first_height=0.00015,  # 0.15mm
                bl_growth_ratio=1.2,
                bl_num_layers=5,
            ),
            MeshLevel(
                name="L4_very_fine",
                mesh_factor=0.5,
                bl_first_height=0.0001,  # 0.1mm
                bl_growth_ratio=1.2,
                bl_num_layers=6,
            ),
        ]

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary for serialization."""
        return {
            "base_case_path": str(self.base_case_path),
            "study_name": self.study_name,
            "output_dir": str(self.output_dir),
            "mesh_levels": [
                {
                    "name": ml.name,
                    "mesh_factor": ml.mesh_factor,
                    "bl_first_height": ml.bl_first_height,
                    "bl_growth_ratio": ml.bl_growth_ratio,
                    "bl_num_layers": ml.bl_num_layers,
                }
                for ml in self.mesh_levels
            ],
            "metric_name": self.metric_name,
            "metric_patch": self.metric_patch,
            "metric_field": self.metric_field,
            "metric_region": self.metric_region,
            "convergence_threshold": self.convergence_threshold,
            "solver_app": self.solver_app,
            "num_iterations": self.num_iterations,
            "geometry_script": str(self.geometry_script) if self.geometry_script else None,
            "adaptive_mode": self.adaptive_mode,
            "refinement_ratio": self.refinement_ratio,
            "bl_refinement_ratio": self.bl_refinement_ratio,
            "max_cells": self.max_cells,
            "max_levels": self.max_levels,
            "max_runtime_per_level": self.max_runtime_per_level,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GridStudyConfig":
        """Create config from dictionary."""
        mesh_levels = [
            MeshLevel(**ml) for ml in data.get("mesh_levels", [])
        ]
        return cls(
            base_case_path=Path(data["base_case_path"]),
            study_name=data.get("study_name", "grid_study"),
            output_dir=Path(data.get("output_dir", "grid_study_output")),
            mesh_levels=mesh_levels,
            metric_name=data.get("metric_name", "T_base_avg"),
            metric_patch=data.get("metric_patch", "heat_source"),
            metric_field=data.get("metric_field", "T"),
            metric_region=data.get("metric_region", "solid"),
            convergence_threshold=data.get("convergence_threshold", 1.0),
            solver_app=data.get("solver_app", "chtMultiRegionSimpleFoam"),
            num_iterations=data.get("num_iterations", 2000),
            geometry_script=Path(data["geometry_script"]) if data.get("geometry_script") else None,
            adaptive_mode=data.get("adaptive_mode", False),
            refinement_ratio=data.get("refinement_ratio", 0.7),
            bl_refinement_ratio=data.get("bl_refinement_ratio", 0.7),
            max_cells=data.get("max_cells", 2_000_000),
            max_levels=data.get("max_levels", 10),
            max_runtime_per_level=data.get("max_runtime_per_level"),
        )

    def save(self, filepath: Path) -> None:
        """Save configuration to JSON file."""
        with open(filepath, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, filepath: Path) -> "GridStudyConfig":
        """Load configuration from JSON file."""
        with open(filepath, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)

    def generate_next_level(self, previous_level: MeshLevel, level_number: int) -> MeshLevel:
        """
        Generate the next finer mesh level based on the previous level.
        Used in adaptive refinement mode.

        Args:
            previous_level: The previous mesh level configuration
            level_number: The level number (for naming)

        Returns:
            New MeshLevel with refined parameters
        """
        return MeshLevel(
            name=f"L{level_number}_adaptive",
            mesh_factor=previous_level.mesh_factor * self.refinement_ratio,
            bl_first_height=previous_level.bl_first_height * self.bl_refinement_ratio,
            bl_growth_ratio=previous_level.bl_growth_ratio,
            bl_num_layers=previous_level.bl_num_layers,
        )
