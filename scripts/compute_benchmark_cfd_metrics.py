#!/usr/bin/env python3
"""Compute normalized-grid CFD comparison metrics between two indoor cases.

This script compares reference vs predicted OpenFOAM outputs by sampling both
VTK volumes on the same normalized room-coordinate lattice. It is intended for
benchmark evaluation, where geometry can differ across runs but we still want a
coarse flow-field similarity score in a room-relative frame.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True


DEFAULT_GRID = (18, 18, 10)
DEFAULT_EPSILON = (0.05, 0.05, 0.08)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def round_or_none(value: float | None, digits: int = 6) -> float | None:
    return round(float(value), digits) if value is not None else None


def room_blocks(scene: dict[str, Any]) -> list[dict[str, Any]]:
    room = scene.get("room", {})
    if "blocks" in room:
        return room["blocks"]
    size = room.get("size", {})
    return [{
        "origin": {"x": 0.0, "y": 0.0, "z": 0.0},
        "size": {"dx": size.get("Lx", 0.0), "dy": size.get("Ly", 0.0), "dz": size.get("Lz", 0.0)},
    }]


def room_bounds(scene: dict[str, Any]) -> dict[str, float]:
    blocks = room_blocks(scene)
    return {
        "Lx": max(block["origin"]["x"] + block["size"]["dx"] for block in blocks),
        "Ly": max(block["origin"]["y"] + block["size"]["dy"] for block in blocks),
        "Lz": max(block["origin"]["z"] + block["size"]["dz"] for block in blocks),
    }


def find_latest_vtk(case_dir: Path) -> Path:
    vtk_dir = case_dir / "VTK"
    candidates: list[tuple[int, Path]] = []
    for child in vtk_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            timestep = int(child.name.split("_")[-1])
        except ValueError:
            continue
        vtk_path = child / "internal.vtu"
        if vtk_path.exists():
            candidates.append((timestep, vtk_path))
    if not candidates:
        raise FileNotFoundError(f"No internal.vtu found under {vtk_dir}")
    candidates.sort(key=lambda item: item[0])
    return candidates[-1][1]


def load_flow_mesh(case_dir: Path) -> tuple[pv.DataSet, Path]:
    vtk_path = find_latest_vtk(case_dir)
    mesh = pv.read(vtk_path)
    if "U" not in mesh.point_data and "U" in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()
    if "p" not in mesh.point_data and "p" in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()
    if "U" not in mesh.point_data:
        raise KeyError(f"Velocity field U not found in {vtk_path}")
    u = np.asarray(mesh.point_data["U"], dtype=float)
    mesh.point_data["Umag"] = np.linalg.norm(u, axis=1)
    return mesh, vtk_path


def make_normalized_lattice(nx: int, ny: int, nz: int, epsx: float, epsy: float, epsz: float) -> np.ndarray:
    xs = np.linspace(epsx, 1.0 - epsx, nx)
    ys = np.linspace(epsy, 1.0 - epsy, ny)
    zs = np.linspace(epsz, 1.0 - epsz, nz)
    grid = np.array([(x, y, z) for z in zs for y in ys for x in xs], dtype=float)
    return grid


def scale_points(norm_pts: np.ndarray, bounds: dict[str, float]) -> np.ndarray:
    pts = np.empty_like(norm_pts)
    pts[:, 0] = norm_pts[:, 0] * bounds["Lx"]
    pts[:, 1] = norm_pts[:, 1] * bounds["Ly"]
    pts[:, 2] = norm_pts[:, 2] * bounds["Lz"]
    return pts


def sample_case(mesh: pv.DataSet, bounds: dict[str, float], norm_pts: np.ndarray) -> dict[str, Any]:
    physical_pts = scale_points(norm_pts, bounds)
    sampled = pv.PolyData(physical_pts).sample(mesh)
    valid = np.asarray(sampled["vtkValidPointMask"], dtype=bool) if "vtkValidPointMask" in sampled.array_names else np.ones(len(norm_pts), dtype=bool)

    payload: dict[str, Any] = {
        "valid_mask": valid,
        "valid_count": int(valid.sum()),
        "valid_ratio": float(valid.mean()) if len(valid) else 0.0,
    }

    if "U" in sampled.array_names:
        payload["U"] = np.asarray(sampled["U"], dtype=float)
        payload["Umag"] = np.linalg.norm(payload["U"], axis=1)
    if "p" in sampled.array_names:
        payload["p"] = np.asarray(sampled["p"], dtype=float).reshape(-1)
    return payload


def vector_metrics(ref: np.ndarray, pred: np.ndarray) -> dict[str, float | None]:
    diff = pred - ref
    l2 = np.linalg.norm(diff, axis=1)
    rmse = math.sqrt(float(np.mean(np.sum(diff ** 2, axis=1)))) if len(diff) else None
    mae = float(np.mean(l2)) if len(l2) else None

    ref_mag = np.linalg.norm(ref, axis=1)
    pred_mag = np.linalg.norm(pred, axis=1)
    denom = ref_mag * pred_mag
    direction_mask = denom > 1e-12
    cosines = np.sum(ref[direction_mask] * pred[direction_mask], axis=1) / denom[direction_mask] if np.any(direction_mask) else np.array([])
    ref_rms = math.sqrt(float(np.mean(ref_mag ** 2))) if len(ref_mag) else None

    return {
        "mae_l2": round_or_none(mae),
        "rmse_l2": round_or_none(rmse),
        "reference_rms": round_or_none(ref_rms),
        "relative_rmse": round_or_none((rmse / ref_rms) if (rmse is not None and ref_rms and ref_rms > 0) else None),
        "mean_direction_cosine": round_or_none(float(np.mean(cosines)) if len(cosines) else None),
        "direction_comparable_count": int(direction_mask.sum()),
    }


def scalar_metrics(ref: np.ndarray, pred: np.ndarray) -> dict[str, float | None]:
    diff = pred - ref
    mae = float(np.mean(np.abs(diff))) if len(diff) else None
    rmse = math.sqrt(float(np.mean(diff ** 2))) if len(diff) else None
    ref_rms = math.sqrt(float(np.mean(ref ** 2))) if len(ref) else None
    return {
        "mae": round_or_none(mae),
        "rmse": round_or_none(rmse),
        "reference_rms": round_or_none(ref_rms),
        "relative_rmse": round_or_none((rmse / ref_rms) if (rmse is not None and ref_rms and ref_rms > 0) else None),
        "reference_mean": round_or_none(float(np.mean(ref)) if len(ref) else None),
        "predicted_mean": round_or_none(float(np.mean(pred)) if len(pred) else None),
    }


def compare_cases(
    reference_scene: Path,
    reference_case: Path,
    predicted_scene: Path,
    predicted_case: Path,
    *,
    nx: int,
    ny: int,
    nz: int,
    epsx: float,
    epsy: float,
    epsz: float,
) -> dict[str, Any]:
    ref_scene = load_json(reference_scene)
    pred_scene = load_json(predicted_scene)
    ref_bounds = room_bounds(ref_scene)
    pred_bounds = room_bounds(pred_scene)
    norm_pts = make_normalized_lattice(nx, ny, nz, epsx, epsy, epsz)

    ref_mesh, ref_vtk = load_flow_mesh(reference_case)
    pred_mesh, pred_vtk = load_flow_mesh(predicted_case)
    ref_sample = sample_case(ref_mesh, ref_bounds, norm_pts)
    pred_sample = sample_case(pred_mesh, pred_bounds, norm_pts)

    ref_valid = ref_sample["valid_mask"]
    pred_valid = pred_sample["valid_mask"]
    overlap = ref_valid & pred_valid
    union = ref_valid | pred_valid

    overlap_count = int(overlap.sum())
    union_count = int(union.sum())
    total_count = int(len(norm_pts))

    summary: dict[str, Any] = {
        "ok": True,
        "grid": {
            "nx": nx,
            "ny": ny,
            "nz": nz,
            "point_count": total_count,
            "epsilon": {"x": epsx, "y": epsy, "z": epsz},
        },
        "reference": {
            "scene": str(reference_scene),
            "case": str(reference_case),
            "vtk": str(ref_vtk),
            "room_bounds": ref_bounds,
            "valid_count": ref_sample["valid_count"],
            "valid_ratio": round_or_none(ref_sample["valid_ratio"]),
        },
        "predicted": {
            "scene": str(predicted_scene),
            "case": str(predicted_case),
            "vtk": str(pred_vtk),
            "room_bounds": pred_bounds,
            "valid_count": pred_sample["valid_count"],
            "valid_ratio": round_or_none(pred_sample["valid_ratio"]),
        },
        "overlap": {
            "count": overlap_count,
            "ratio_vs_total": round_or_none(overlap_count / total_count if total_count else 0.0),
            "ratio_vs_union": round_or_none(overlap_count / union_count if union_count else 1.0),
            "union_count": union_count,
        },
    }

    if overlap_count == 0:
        summary["field_metrics"] = {
            "velocity_vector": None,
            "velocity_magnitude": None,
            "pressure": None,
        }
        summary["aggregate_score"] = {
            "cfd_score": 0.0,
            "notes": ["No overlapping valid sampled points between reference and prediction."],
        }
        return summary

    field_metrics: dict[str, Any] = {}

    if "U" in ref_sample and "U" in pred_sample:
        ref_u = ref_sample["U"][overlap]
        pred_u = pred_sample["U"][overlap]
        ref_umag = ref_sample["Umag"][overlap]
        pred_umag = pred_sample["Umag"][overlap]
        field_metrics["velocity_vector"] = vector_metrics(ref_u, pred_u)
        field_metrics["velocity_magnitude"] = scalar_metrics(ref_umag, pred_umag)
    else:
        field_metrics["velocity_vector"] = None
        field_metrics["velocity_magnitude"] = None

    if "p" in ref_sample and "p" in pred_sample:
        field_metrics["pressure"] = scalar_metrics(ref_sample["p"][overlap], pred_sample["p"][overlap])
    else:
        field_metrics["pressure"] = None

    overlap_ratio = overlap_count / union_count if union_count else 1.0
    components = [overlap_ratio]

    umag_rel_rmse = None
    if field_metrics["velocity_magnitude"]:
        umag_rel_rmse = field_metrics["velocity_magnitude"].get("relative_rmse")
        if umag_rel_rmse is not None:
            components.append(max(0.0, 1.0 - min(1.0, umag_rel_rmse)))

    dir_cos = None
    if field_metrics["velocity_vector"]:
        dir_cos = field_metrics["velocity_vector"].get("mean_direction_cosine")
        if dir_cos is not None:
            components.append(max(0.0, min(1.0, 0.5 * (dir_cos + 1.0))))

    p_rel_rmse = None
    if field_metrics["pressure"]:
        p_rel_rmse = field_metrics["pressure"].get("relative_rmse")
        if p_rel_rmse is not None:
            components.append(max(0.0, 1.0 - min(1.0, p_rel_rmse)))

    summary["field_metrics"] = field_metrics
    summary["aggregate_score"] = {
        "cfd_score": round_or_none(sum(components) / len(components) if components else None),
        "components": {
            "overlap_ratio_vs_union": round_or_none(overlap_ratio),
            "velocity_magnitude_similarity": round_or_none(max(0.0, 1.0 - min(1.0, umag_rel_rmse)) if umag_rel_rmse is not None else None),
            "velocity_direction_similarity": round_or_none(max(0.0, min(1.0, 0.5 * (dir_cos + 1.0))) if dir_cos is not None else None),
            "pressure_similarity": round_or_none(max(0.0, 1.0 - min(1.0, p_rel_rmse)) if p_rel_rmse is not None else None),
        },
    }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute normalized-grid CFD comparison metrics between two indoor cases")
    parser.add_argument("--reference-scene", type=Path, required=True)
    parser.add_argument("--reference-case", type=Path, required=True)
    parser.add_argument("--predicted-scene", type=Path, required=True)
    parser.add_argument("--predicted-case", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--nx", type=int, default=DEFAULT_GRID[0])
    parser.add_argument("--ny", type=int, default=DEFAULT_GRID[1])
    parser.add_argument("--nz", type=int, default=DEFAULT_GRID[2])
    parser.add_argument("--epsx", type=float, default=DEFAULT_EPSILON[0])
    parser.add_argument("--epsy", type=float, default=DEFAULT_EPSILON[1])
    parser.add_argument("--epsz", type=float, default=DEFAULT_EPSILON[2])
    args = parser.parse_args()

    payload = compare_cases(
        args.reference_scene.expanduser().resolve(),
        args.reference_case.expanduser().resolve(),
        args.predicted_scene.expanduser().resolve(),
        args.predicted_case.expanduser().resolve(),
        nx=args.nx,
        ny=args.ny,
        nz=args.nz,
        epsx=args.epsx,
        epsy=args.epsy,
        epsz=args.epsz,
    )
    write_json(args.output.expanduser().resolve(), payload)
    print(json.dumps(payload, indent=2))
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":
    raise SystemExit(main())
