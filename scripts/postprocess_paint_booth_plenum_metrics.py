#!/usr/bin/env python3
"""Post-process paint-booth plenum/filter OpenFOAM VTK output.

The script expects `foamToVTK -latestTime` output inside an OpenFOAM case and
computes engineering metrics useful for filter/plenum development:

- supply/floor patch flow rates from boundary VTP files
- region-averaged velocity and pressure metrics
- filter pressure drop proxy from slabs above/below the filter
- downdraft uniformity metrics below the filter and in the work zone

Run inside an environment with PyVista/VTK available, e.g. the project Docker
image `openfoam-pipeline-local:latest`.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv


PATCH_NORMALS = {
    "supplyInlet": np.array([0.0, 0.0, 1.0]),   # top patch outward normal
    "floorExhaust": np.array([0.0, 0.0, -1.0]), # bottom patch outward normal
}

DEFAULT_REGIONS = {
    "near_car": {"x": [-0.25, 4.75], "y": [-1.25, 1.25], "z": [0.10, 2.00]},
    "roof_region": {"x": [0.0, 4.5], "y": [-1.0, 1.0], "z": [1.25, 1.90]},
    "work_zone": {"x": [-1.0, 6.0], "y": [-1.6, 1.6], "z": [0.20, 2.70]},
    "filter_below": {"x": [-1.0, 6.0], "y": [-1.6, 1.6], "z": [2.78, 2.94]},
    "filter_layer_panel": {"x": [-1.0, 6.0], "y": [-1.6, 1.6], "z": [2.95, 3.05]},
    "plenum": {"x": [-1.0, 6.0], "y": [-1.6, 1.6], "z": [3.08, 3.70]},
    "pressure_above_filter": {"x": [-1.0, 6.0], "y": [-1.6, 1.6], "z": [3.08, 3.18]},
    "pressure_below_filter": {"x": [-1.0, 6.0], "y": [-1.6, 1.6], "z": [2.82, 2.92]},
}


def latest_vtk_dir(case_dir: Path) -> Path:
    vtk_root = case_dir / "VTK"
    if not vtk_root.exists():
        raise FileNotFoundError(f"VTK directory not found: {vtk_root}. Run foamToVTK -latestTime first.")
    candidates = [p for p in vtk_root.iterdir() if p.is_dir() and (p / "internal.vtu").exists()]
    if not candidates:
        raise FileNotFoundError(f"No foamToVTK time directory with internal.vtu found under {vtk_root}")

    def time_key(path: Path) -> tuple[float, str]:
        # OpenFOAM foamToVTK dirs usually end with _<time>.
        token = path.name.rsplit("_", 1)[-1]
        try:
            return (float(token), path.name)
        except ValueError:
            return (-math.inf, path.name)

    return sorted(candidates, key=time_key)[-1]


def cell_centers(mesh: pv.DataSet) -> np.ndarray:
    return np.asarray(mesh.cell_centers().points)


def select_box_indices(centers: np.ndarray, box: dict[str, list[float]]) -> np.ndarray:
    mask = np.ones(len(centers), dtype=bool)
    for axis, idx in (("x", 0), ("y", 1), ("z", 2)):
        lo, hi = box[axis]
        mask &= centers[:, idx] >= lo
        mask &= centers[:, idx] <= hi
    return np.where(mask)[0]


def summarize_values(values: np.ndarray) -> dict[str, Any]:
    if values.size == 0:
        return {"count": 0}
    return {
        "count": int(values.size),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "min": float(np.min(values)),
        "p05": float(np.percentile(values, 5)),
        "p50": float(np.percentile(values, 50)),
        "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def region_metrics(mesh: pv.UnstructuredGrid, regions: dict[str, dict[str, list[float]]]) -> dict[str, Any]:
    centers = cell_centers(mesh)
    U = np.asarray(mesh.cell_data["U"])
    p = np.asarray(mesh.cell_data["p"])
    U_mag = np.linalg.norm(U, axis=1)
    out: dict[str, Any] = {}
    for name, box in regions.items():
        idx = select_box_indices(centers, box)
        if idx.size == 0:
            out[name] = {"box": box, "n_cells": 0}
            continue
        uz = U[idx, 2]
        mean_uz = float(np.mean(uz))
        out[name] = {
            "box": box,
            "n_cells": int(idx.size),
            "U_mag": summarize_values(U_mag[idx]),
            "Ux": summarize_values(U[idx, 0]),
            "Uy": summarize_values(U[idx, 1]),
            "Uz": summarize_values(uz),
            "p": summarize_values(p[idx]),
            "downdraft_uniformity": {
                "Uz_std_over_abs_mean": float(np.std(uz) / max(abs(mean_uz), 1e-12)),
                "reverse_flow_fraction_Uz_positive": float(np.mean(uz > 0.0)),
                "low_speed_fraction_Umag_lt_0p05": float(np.mean(U_mag[idx] < 0.05)),
            },
        }
    return out


def poly_cell_areas(poly: pv.PolyData) -> np.ndarray:
    sized = poly.compute_cell_sizes(length=False, area=True, volume=False)
    return np.asarray(sized.cell_data["Area"])


def patch_flow_metrics(boundary_dir: Path, patch_name: str) -> dict[str, Any]:
    path = boundary_dir / f"{patch_name}.vtp"
    if not path.exists():
        return {"patch": patch_name, "exists": False, "path": str(path)}
    mesh = pv.read(path)
    if "U" not in mesh.cell_data:
        raise KeyError(f"Patch {patch_name} has no cell_data['U'] in {path}")
    if "p" not in mesh.cell_data:
        raise KeyError(f"Patch {patch_name} has no cell_data['p'] in {path}")
    areas = poly_cell_areas(mesh)
    U = np.asarray(mesh.cell_data["U"])
    p = np.asarray(mesh.cell_data["p"])
    normal = PATCH_NORMALS[patch_name]
    un = U @ normal
    signed_flow = float(np.sum(un * areas))
    area = float(np.sum(areas))
    return {
        "patch": patch_name,
        "exists": True,
        "path": str(path),
        "n_cells": int(mesh.n_cells),
        "area_m2": area,
        "normal": normal.tolist(),
        "signed_outward_flow_m3s": signed_flow,
        "inflow_m3s": float(max(-signed_flow, 0.0)),
        "outflow_m3s": float(max(signed_flow, 0.0)),
        "area_weighted_normal_velocity_mps": float(signed_flow / max(area, 1e-12)),
        "area_weighted_abs_normal_velocity_mps": float(np.sum(np.abs(un) * areas) / max(area, 1e-12)),
        "area_weighted_p": float(np.sum(p * areas) / max(area, 1e-12)),
        "normal_velocity": summarize_values(un),
    }




def virtual_horizontal_plane_flux(
    mesh: pv.UnstructuredGrid,
    name: str,
    z: float,
    x_range: list[float],
    y_range: list[float],
    normal: Iterable[float] = (0.0, 0.0, -1.0),
) -> dict[str, Any]:
    """Integrate velocity through a virtual horizontal plane slice.

    The default normal points downward, so positive `normal_flow_m3s` means
    downward flow through the plane.
    """
    sl = mesh.slice(normal="z", origin=(0.0, 0.0, z))
    if sl.n_cells == 0:
        return {"name": name, "z": z, "n_cells": 0}
    centers = np.asarray(sl.cell_centers().points)
    idx = np.where(
        (centers[:, 0] >= x_range[0])
        & (centers[:, 0] <= x_range[1])
        & (centers[:, 1] >= y_range[0])
        & (centers[:, 1] <= y_range[1])
    )[0]
    if idx.size == 0:
        return {"name": name, "z": z, "n_cells": 0, "x": x_range, "y": y_range}
    if "U" in sl.cell_data:
        U = np.asarray(sl.cell_data["U"])[idx]
    elif "U" in sl.point_data:
        sl = sl.point_data_to_cell_data()
        U = np.asarray(sl.cell_data["U"])[idx]
    else:
        raise KeyError("Slice has no U data")
    areas = poly_cell_areas(sl)[idx]
    n = np.asarray(list(normal), dtype=float)
    un = U @ n
    area = float(np.sum(areas))
    normal_flow = float(np.sum(un * areas))
    return {
        "name": name,
        "z": z,
        "x": x_range,
        "y": y_range,
        "normal": n.tolist(),
        "n_cells": int(idx.size),
        "area_m2": area,
        "normal_flow_m3s": normal_flow,
        "area_weighted_normal_velocity_mps": float(normal_flow / max(area, 1e-12)),
        "normal_velocity": summarize_values(un),
    }


def load_case_metadata(case_dir: Path) -> dict[str, Any] | None:
    path = case_dir / "constant" / "geometry_metadata.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path, help="OpenFOAM case directory")
    parser.add_argument("--output", type=Path, default=None, help="Output JSON path. Default: <case>/post_plenum_metrics.json")
    parser.add_argument("--vtk-dir", type=Path, default=None, help="Specific foamToVTK time directory")
    args = parser.parse_args()

    case_dir = args.case_dir.resolve()
    vtk_dir = args.vtk_dir.resolve() if args.vtk_dir else latest_vtk_dir(case_dir)
    internal_path = vtk_dir / "internal.vtu"
    boundary_dir = vtk_dir / "boundary"
    if not internal_path.exists():
        raise FileNotFoundError(internal_path)
    if not boundary_dir.exists():
        raise FileNotFoundError(boundary_dir)

    mesh = pv.read(internal_path)
    regions = region_metrics(mesh, DEFAULT_REGIONS)
    patches = {
        name: patch_flow_metrics(boundary_dir, name)
        for name in ["supplyInlet", "floorExhaust"]
    }
    virtual_planes = {
        "filter_panel_below_z2p94": virtual_horizontal_plane_flux(mesh, "filter_panel_below_z2p94", 2.94, [-1.0, 6.0], [-1.6, 1.6]),
        "work_zone_mid_z1p50": virtual_horizontal_plane_flux(mesh, "work_zone_mid_z1p50", 1.50, [-1.0, 6.0], [-1.6, 1.6]),
    }

    supply_in = patches["supplyInlet"].get("inflow_m3s", 0.0)
    floor_out = patches["floorExhaust"].get("outflow_m3s", 0.0)
    above_p = regions["pressure_above_filter"].get("p", {}).get("mean")
    below_p = regions["pressure_below_filter"].get("p", {}).get("mean")
    pressure_drop = None
    if above_p is not None and below_p is not None:
        pressure_drop = float(above_p - below_p)

    metrics: dict[str, Any] = {
        "case_dir": str(case_dir),
        "vtk_dir": str(vtk_dir),
        "metadata": load_case_metadata(case_dir),
        "patches": patches,
        "virtual_planes": virtual_planes,
        "mass_balance": {
            "supply_inflow_m3s": float(supply_in),
            "floor_outflow_m3s": float(floor_out),
            "net_out_minus_in_m3s": float(floor_out - supply_in),
            "relative_imbalance_out_minus_in": float((floor_out - supply_in) / max(abs(supply_in), 1e-12)),
        },
        "regions": regions,
        "filter_pressure_drop_proxy_m2_per_s2": pressure_drop,
        "notes": [
            "OpenFOAM incompressible kinematic pressure p has units m2/s2; multiply by rho for Pa.",
            "Filter pressure drop is a slab-average proxy, not a patch-integrated porous-source pressure loss.",
            "Patch fluxes use expected outward normals for supplyInlet (+z) and floorExhaust (-z).",
            "Virtual plane fluxes use downward normal (0,0,-1), so positive values mean downdraft.",
        ],
    }

    output = args.output or (case_dir / "post_plenum_metrics.json")
    output.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    print(f"Wrote {output}")
    sys.stdout.flush()
    sys.stderr.flush()
    # Some VTK/PyVista builds can leave non-daemon cleanup threads alive at
    # interpreter shutdown inside Docker. The script has already written all
    # outputs, so exit immediately to keep large sweeps from hanging.
    os._exit(0)


if __name__ == "__main__":
    main()
