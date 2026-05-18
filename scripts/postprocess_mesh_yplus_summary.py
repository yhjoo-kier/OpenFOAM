#!/usr/bin/env python3
"""Summarize snappyHexMesh resolution and wall y+ for a paint-booth OpenFOAM case.

Expected workflow before running:
  blockMesh
  snappyHexMesh -overwrite
  checkMesh
  simpleFoam
  simpleFoam -postProcess -func yPlus -latestTime
  foamToVTK -latestTime

The script is intentionally lightweight and relies on the Docker image's PyVista/VTK stack.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pyvista as pv


def percentile_stats(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {"count": 0}
    return {
        "count": int(values.size),
        "min": float(np.min(values)),
        "p05": float(np.percentile(values, 5)),
        "p50": float(np.percentile(values, 50)),
        "mean": float(np.mean(values)),
        "p95": float(np.percentile(values, 95)),
        "p99": float(np.percentile(values, 99)),
        "max": float(np.max(values)),
    }


def latest_vtk_dir(case_dir: Path) -> Path:
    vtk_root = case_dir / "VTK"
    if not vtk_root.exists():
        raise FileNotFoundError(f"VTK directory not found: {vtk_root}")
    dirs = [p for p in vtk_root.iterdir() if p.is_dir()]
    if not dirs:
        raise FileNotFoundError(f"No VTK time directories found under {vtk_root}")

    def key(p: Path) -> tuple[int, float | str]:
        try:
            return (1, float(p.name))
        except ValueError:
            return (0, p.name)

    return sorted(dirs, key=key)[-1]


def parse_check_mesh(case_dir: Path) -> dict[str, float | int | str | bool]:
    out: dict[str, float | int | str | bool] = {}
    log = case_dir / "log.checkMesh_snappy"
    if not log.exists():
        return out
    text = log.read_text(errors="replace")
    patterns = {
        "points": r"points:\s+(\d+)",
        "faces": r"faces:\s+(\d+)",
        "internal_faces": r"internal faces:\s+(\d+)",
        "cells": r"cells:\s+(\d+)",
        "max_aspect_ratio": r"Max aspect ratio =\s+([0-9.eE+-]+)",
        "max_non_orthogonality": r"Max non-orthogonality =\s+([0-9.eE+-]+)",
        "average_non_orthogonality": r"average:\s+([0-9.eE+-]+)",
        "max_skewness": r"Max skewness =\s+([0-9.eE+-]+)",
    }
    for name, pat in patterns.items():
        m = re.search(pat, text)
        if not m:
            continue
        value = m.group(1)
        out[name] = int(value) if value.isdigit() else float(value)
    out["mesh_ok"] = "Mesh OK" in text
    return out


def parse_layer_summary(case_dir: Path) -> dict[str, float | int]:
    log = case_dir / "log.snappyHexMesh"
    if not log.exists():
        return {}
    text = log.read_text(errors="replace")
    # Use the final layer report, e.g.:
    # carBody 9935     3.97     0.0155    99.6
    matches = re.findall(r"^carBody\s+(\d+)\s+([0-9.]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)", text, flags=re.M)
    if not matches:
        return {}
    faces, layers, thickness, percent = matches[-1]
    return {
        "carBody_faces_at_layering": int(faces),
        "carBody_average_layers_added": float(layers),
        "carBody_overall_layer_thickness_m": float(thickness),
        "carBody_layer_coverage_percent": float(percent),
    }


def summarize_case(case_dir: Path) -> dict:
    vtk_dir = latest_vtk_dir(case_dir)
    internal_files = list(vtk_dir.glob("*internal*.vtu")) + list(vtk_dir.glob("internal.vtu"))
    if not internal_files:
        raise FileNotFoundError(f"No internal .vtu found under {vtk_dir}")
    mesh = pv.read(internal_files[0])
    sized = mesh.compute_cell_sizes(length=False, area=False, volume=True)
    volumes = np.asarray(sized.cell_data["Volume"])
    leq = np.cbrt(volumes)
    centers = np.asarray(mesh.cell_centers().points)

    regions = {
        "whole_domain": ((-1.5, 6.5), (-2.0, 2.0), (0.0, 3.8)),
        "work_zone": ((-1.0, 6.0), (-1.6, 1.6), (0.2, 2.7)),
        "near_car": ((-0.25, 4.75), (-1.25, 1.25), (0.1, 2.0)),
        "car_close_box": ((0.0, 4.5), (-1.0, 1.0), (0.18, 1.8)),
        "filter_panel": ((-1.0, 6.0), (-1.6, 1.6), (2.95, 3.05)),
        "plenum": ((-1.0, 6.0), (-1.6, 1.6), (3.08, 3.7)),
    }
    region_stats = {}
    for name, (xr, yr, zr) in regions.items():
        mask = (
            (centers[:, 0] >= xr[0]) & (centers[:, 0] <= xr[1])
            & (centers[:, 1] >= yr[0]) & (centers[:, 1] <= yr[1])
            & (centers[:, 2] >= zr[0]) & (centers[:, 2] <= zr[1])
        )
        region_stats[name] = percentile_stats(leq[mask])

    boundary_dir = vtk_dir / "boundary"
    car_files = list(boundary_dir.glob("carBody.vtp")) + list(boundary_dir.glob("carBody*.vtp"))
    if not car_files:
        raise FileNotFoundError(f"No carBody boundary VTP found under {boundary_dir}")
    car = pv.read(car_files[0])
    car_area = np.asarray(car.compute_cell_sizes(length=False, area=True, volume=False).cell_data["Area"])
    face_l = np.sqrt(car_area)

    yplus = None
    if "yPlus" in car.cell_data:
        yplus = np.asarray(car.cell_data["yPlus"])
    elif "yPlus" in car.point_data:
        yplus = np.asarray(car.point_data["yPlus"])

    yplus_stats = {}
    if yplus is not None:
        yplus_stats = percentile_stats(yplus)
        yplus_stats.update(
            {
                "fraction_lt_1": float(np.mean(yplus < 1.0)),
                "fraction_1_to_5": float(np.mean((yplus >= 1.0) & (yplus <= 5.0))),
                "fraction_5_to_30": float(np.mean((yplus > 5.0) & (yplus <= 30.0))),
                "fraction_gt_30": float(np.mean(yplus > 30.0)),
            }
        )

    return {
        "case": str(case_dir),
        "vtk_time_dir": str(vtk_dir.relative_to(case_dir)),
        "checkMesh": parse_check_mesh(case_dir),
        "snappy_layers": parse_layer_summary(case_dir),
        "cell_equivalent_size_m": region_stats,
        "carBody_face_size_m": {
            **percentile_stats(face_l),
            "surface_area_m2": float(np.sum(car_area)),
        },
        "carBody_yPlus": yplus_stats,
        "notes": {
            "yplus_target_interpretation": "For kOmegaSST with integrated/low-y+ treatment, carBody y+ around 1 is desirable. This summary reports whether the car body is mostly in y+<5 and flags y+>30 outliers.",
            "wall_scope": "Prism layers are currently applied to carBody only; booth walls still use wall functions and show high y+ in OpenFOAM postProcess logs.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    result = summarize_case(args.case_dir)
    out = args.output or args.case_dir / "mesh_yplus_summary.json"
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
