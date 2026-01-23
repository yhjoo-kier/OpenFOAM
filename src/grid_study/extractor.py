"""Result extractor for grid study metrics."""

import subprocess
import re
from pathlib import Path
from typing import Optional, Dict, Any
from dataclasses import dataclass


@dataclass
class MetricResult:
    """Extracted metric from simulation."""
    value: float
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    patch_name: str = ""
    field_name: str = ""
    region: str = ""
    time_step: Optional[int] = None


class ResultExtractor:
    """Extracts metrics from OpenFOAM simulation results."""

    def __init__(self, case_dir: Path):
        self.case_dir = case_dir

    def get_patch_average(
        self,
        field: str,
        patch: str,
        region: Optional[str] = None,
        time_step: Optional[int] = None,
    ) -> MetricResult:
        """
        Get area-weighted average of a field on a patch.

        Args:
            field: Field name (e.g., "T", "p", "U")
            patch: Patch name (e.g., "heat_source", "inlet")
            region: Region name for multi-region cases (e.g., "fluid", "solid")
            time_step: Time step to extract (default: latest)

        Returns:
            MetricResult with extracted values
        """
        if time_step is None:
            time_step = self._get_latest_time(region)

        if time_step is None:
            raise ValueError("No time directories found")

        # Build path to field file
        if region:
            field_path = self.case_dir / str(time_step) / region / field
        else:
            field_path = self.case_dir / str(time_step) / field

        if not field_path.exists():
            raise FileNotFoundError(f"Field file not found: {field_path}")

        # Parse the field file
        content = field_path.read_text()

        # Extract boundary field values for the specified patch
        avg_value, min_val, max_val = self._parse_boundary_field(content, patch)

        return MetricResult(
            value=avg_value,
            min_value=min_val,
            max_value=max_val,
            patch_name=patch,
            field_name=field,
            region=region or "",
            time_step=time_step,
        )

    def _get_latest_time(self, region: Optional[str] = None) -> Optional[int]:
        """Find the latest time directory."""
        time_dirs = []

        for item in self.case_dir.iterdir():
            if item.is_dir():
                try:
                    time_val = float(item.name)
                    # For multi-region, check if region subdirectory exists
                    if region:
                        if (item / region).exists():
                            time_dirs.append((time_val, item.name))
                    else:
                        time_dirs.append((time_val, item.name))
                except ValueError:
                    pass

        if not time_dirs:
            return None

        # Sort by time value and return the latest
        time_dirs.sort(key=lambda x: x[0], reverse=True)
        return int(float(time_dirs[0][1]))

    def _parse_boundary_field(
        self, content: str, patch: str
    ) -> tuple[float, Optional[float], Optional[float]]:
        """
        Parse boundary field values from OpenFOAM field file.

        Returns: (average, min, max)
        """
        # Find the boundary field section
        boundary_match = re.search(r"boundaryField\s*\{(.*?)\n\}", content, re.DOTALL)
        if not boundary_match:
            raise ValueError("Could not find boundaryField section")

        boundary_content = boundary_match.group(1)

        # Find the specific patch
        # Pattern: patch_name { ... value uniform <number>; ... } or value nonuniform ...
        patch_pattern = rf"{re.escape(patch)}\s*\{{([^}}]+)\}}"
        patch_match = re.search(patch_pattern, boundary_content, re.DOTALL)

        if not patch_match:
            raise ValueError(f"Could not find patch '{patch}' in boundary field")

        patch_content = patch_match.group(1)

        # Try to extract value
        # Case 1: uniform value
        uniform_match = re.search(r"value\s+uniform\s+([\d.eE+-]+)", patch_content)
        if uniform_match:
            val = float(uniform_match.group(1))
            return val, val, val

        # Case 2: nonuniform List<scalar>
        nonuniform_match = re.search(
            r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(([\d\s.eE+-]+)\)",
            patch_content,
        )
        if nonuniform_match:
            values_str = nonuniform_match.group(2)
            values = [float(v) for v in values_str.split()]
            if values:
                avg = sum(values) / len(values)
                return avg, min(values), max(values)

        raise ValueError(f"Could not parse value for patch '{patch}'")

    def get_cell_count(self, region: Optional[str] = None) -> int:
        """Get the number of cells in the mesh."""
        if region:
            mesh_path = self.case_dir / "constant" / region / "polyMesh"
        else:
            mesh_path = self.case_dir / "constant" / "polyMesh"

        owner_file = mesh_path / "owner"
        if not owner_file.exists():
            raise FileNotFoundError(f"Mesh owner file not found: {owner_file}")

        content = owner_file.read_text()

        # Find the number of cells from the header or count unique cell IDs
        # OpenFOAM format has nCells in the FoamFile header or we count entries
        n_cells_match = re.search(r"nCells:\s*(\d+)", content)
        if n_cells_match:
            return int(n_cells_match.group(1))

        # Alternative: count from the list size
        size_match = re.search(r"^\s*(\d+)\s*\(", content, re.MULTILINE)
        if size_match:
            # This is number of faces, need to find max cell index + 1
            values_match = re.search(r"\(\s*([\d\s]+)\s*\)", content, re.DOTALL)
            if values_match:
                values = [int(v) for v in values_match.group(1).split()]
                return max(values) + 1 if values else 0

        return 0

    def get_mesh_info(self) -> Dict[str, Any]:
        """Get mesh statistics for both regions."""
        info = {}

        # Try to get cell counts for fluid and solid regions
        for region in ["fluid", "solid"]:
            try:
                info[f"{region}_cells"] = self.get_cell_count(region)
            except Exception:
                info[f"{region}_cells"] = None

        # Total cells
        total = sum(v for v in [info.get("fluid_cells"), info.get("solid_cells")] if v)
        info["total_cells"] = total

        return info
