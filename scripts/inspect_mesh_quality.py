#!/usr/bin/env python3
"""Inspect mesh quality for CFD preprocessing.

Usage inside the existing Docker image, for example:

    docker run --rm -v "$PWD:$PWD" -w "$PWD" openfoam-pipeline-local:latest \
      python scripts/inspect_mesh_quality.py data/geometry/raw/body_shell.stl \
      --out docs/mesh_quality_report.json

The script uses PyVista/VTK when available and reports basic mesh statistics,
bounds, triangle/non-triangle cell counts, connectivity, open edges, and a few
cell-quality indicators relevant to OpenFOAM preprocessing.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pyvista as pv


def _safe_float(x):
    try:
        y = float(x)
        if math.isnan(y) or math.isinf(y):
            return None
        return y
    except Exception:
        return None


def inspect_mesh(path: Path) -> dict:
    mesh = pv.read(path)
    surf = mesh.extract_surface().triangulate()

    bounds = surf.bounds
    lengths = [bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]]

    # Boundary/open edges. For a watertight closed surface, this should be zero.
    try:
        feature_edges = surf.extract_feature_edges(
            boundary_edges=True,
            non_manifold_edges=True,
            feature_edges=False,
            manifold_edges=False,
        )
        open_or_nonmanifold_edges = int(feature_edges.n_cells)
    except Exception as exc:
        open_or_nonmanifold_edges = None
        feature_edges_error = repr(exc)
    else:
        feature_edges_error = None

    # Connected components / bodies.
    try:
        conn = surf.connectivity()
        region_ids = np.asarray(conn.cell_data.get("RegionId", []))
        n_regions = int(region_ids.max() + 1) if region_ids.size else 0
        region_cell_counts = np.bincount(region_ids.astype(int)).tolist() if region_ids.size else []
    except Exception as exc:
        n_regions = None
        region_cell_counts = []
        connectivity_error = repr(exc)
    else:
        connectivity_error = None

    # Cell area quality.
    try:
        sizes = surf.compute_cell_sizes(length=False, area=True, volume=False)
        areas = np.asarray(sizes.cell_data["Area"])
        positive = areas[areas > 0]
        area_stats = {
            "min": _safe_float(np.min(positive)) if positive.size else None,
            "p01": _safe_float(np.percentile(positive, 1)) if positive.size else None,
            "median": _safe_float(np.median(positive)) if positive.size else None,
            "p99": _safe_float(np.percentile(positive, 99)) if positive.size else None,
            "max": _safe_float(np.max(positive)) if positive.size else None,
            "zero_or_negative_count": int(np.count_nonzero(areas <= 0)),
        }
    except Exception as exc:
        area_stats = {"error": repr(exc)}

    # Volume only meaningful if closed/watertight and consistently oriented.
    try:
        volume = _safe_float(surf.volume)
    except Exception:
        volume = None

    return {
        "path": str(path),
        "file_size_bytes": path.stat().st_size,
        "n_points": int(surf.n_points),
        "n_cells": int(surf.n_cells),
        "bounds": {
            "xmin": _safe_float(bounds[0]),
            "xmax": _safe_float(bounds[1]),
            "ymin": _safe_float(bounds[2]),
            "ymax": _safe_float(bounds[3]),
            "zmin": _safe_float(bounds[4]),
            "zmax": _safe_float(bounds[5]),
        },
        "lengths": {
            "x": _safe_float(lengths[0]),
            "y": _safe_float(lengths[1]),
            "z": _safe_float(lengths[2]),
        },
        "open_or_nonmanifold_edges": open_or_nonmanifold_edges,
        "feature_edges_error": feature_edges_error,
        "connected_regions": n_regions,
        "connected_region_cell_counts_top10": sorted(region_cell_counts, reverse=True)[:10],
        "connectivity_error": connectivity_error,
        "area_stats": area_stats,
        "signed_volume_if_closed": volume,
        "heuristic_cfd_readiness": {
            "watertight_likely": open_or_nonmanifold_edges == 0 if open_or_nonmanifold_edges is not None else None,
            "single_component_likely": n_regions == 1 if n_regions is not None else None,
            "very_large_for_snappyhexmesh": int(surf.n_cells) > 1_000_000,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    result = inspect_mesh(args.mesh)
    text = json.dumps(result, indent=2, ensure_ascii=False)
    print(text)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
