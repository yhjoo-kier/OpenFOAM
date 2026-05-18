#!/usr/bin/env python3
"""Export completed paint-booth OpenFOAM cases to ML-ready NPZ samples.

The first target is a steady cell-centered point-cloud/graph representation for
neural-operator/GINO experiments.  The script intentionally does not run CFD; it
only consumes completed cases with foamToVTK output and postprocessing metadata.

Expected runtime: project Docker image `openfoam-pipeline-local:latest` where
PyVista/VTK and NumPy are available.
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv


REQUIRED_FIELDS = ("U", "p")
OPTIONAL_FIELDS = ("k", "omega", "nut", "yPlus")


def latest_vtk_dir(case_dir: Path) -> Path:
    """Return latest foamToVTK time directory containing internal.vtu."""
    vtk_root = case_dir / "VTK"
    if not vtk_root.exists():
        raise FileNotFoundError(f"VTK directory not found: {vtk_root}")
    candidates = [p for p in vtk_root.iterdir() if p.is_dir() and (p / "internal.vtu").exists()]
    if not candidates:
        raise FileNotFoundError(f"No internal.vtu found below {vtk_root}")

    def time_key(path: Path) -> tuple[float, str]:
        token = path.name.rsplit("_", 1)[-1]
        try:
            return (float(token), path.name)
        except ValueError:
            return (-math.inf, path.name)

    return sorted(candidates, key=time_key)[-1]


def load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def get_nested(data: dict[str, Any] | None, path: str, default: Any = None) -> Any:
    cur: Any = data
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def supply_velocity_from_case(case_dir: Path, metadata: dict[str, Any] | None) -> float | None:
    """Infer supply velocity from metadata or case name such as *_u4p36."""
    for key_path in (
        "supply_velocity_mps",
        "flow.supply_velocity_mps",
        "boundary_conditions.supply_velocity_mps",
        "inputs.supply_velocity_mps",
    ):
        value = get_nested(metadata, key_path)
        if value is not None:
            try:
                return float(value)
            except (TypeError, ValueError):
                pass
    match = re.search(r"_u([0-9]+p[0-9]+|[0-9]+(?:\.[0-9]+)?)$", case_dir.name)
    if match:
        return float(match.group(1).replace("p", "."))
    return None


def read_flow_summary(case_root: Path) -> dict[str, dict[str, Any]]:
    """Load existing QA rows keyed by case name if flow_qa_summary.json exists."""
    summary = load_json(case_root / "flow_qa_summary.json")
    if not summary:
        return {}
    rows = summary.get("flow_conditions", [])
    if not isinstance(rows, list):
        return {}
    return {str(row.get("case")): row for row in rows if isinstance(row, dict) and row.get("case")}


def discover_cases(case_root: Path, include: str, exclude: str | None) -> list[Path]:
    if not case_root.exists():
        raise FileNotFoundError(case_root)
    cases: list[Path] = []
    for path in sorted(p for p in case_root.iterdir() if p.is_dir()):
        if not fnmatch.fnmatch(path.name, include):
            continue
        if exclude and fnmatch.fnmatch(path.name, exclude):
            continue
        if (path / "VTK").exists():
            cases.append(path)
    return cases


def extract_arrays(mesh: pv.UnstructuredGrid, fields: tuple[str, ...]) -> dict[str, np.ndarray]:
    missing = [name for name in REQUIRED_FIELDS if name not in mesh.cell_data]
    if missing:
        raise KeyError(f"Required cell fields missing from internal.vtu: {missing}")
    arrays: dict[str, np.ndarray] = {
        "cell_centers": np.asarray(mesh.cell_centers().points, dtype=np.float32),
    }
    for name in fields:
        if name in mesh.cell_data:
            data = np.asarray(mesh.cell_data[name])
            if data.dtype.kind == "f":
                data = data.astype(np.float32)
            arrays[name] = data
    return arrays


def build_metadata(
    *,
    case_dir: Path,
    case_root: Path,
    vtk_dir: Path,
    mesh: pv.UnstructuredGrid,
    arrays: dict[str, np.ndarray],
    qa_row: dict[str, Any] | None,
) -> dict[str, Any]:
    geometry_metadata = load_json(case_dir / "constant" / "geometry_metadata.json")
    plenum_metrics = load_json(case_dir / "post_plenum_metrics.json")
    yplus_summary = load_json(case_dir / "mesh_yplus_summary.json")
    supply_velocity = supply_velocity_from_case(case_dir, geometry_metadata)
    if qa_row and qa_row.get("supply_velocity_mps") is not None:
        supply_velocity = float(qa_row["supply_velocity_mps"])

    metadata: dict[str, Any] = {
        "case": case_dir.name,
        "case_path": str(case_dir.relative_to(case_root.parent)),
        "vtk_dir": str(vtk_dir.relative_to(case_root.parent)),
        "sample_id": case_dir.name.replace("l2_flow_qa_", "l2_"),
        "supply_velocity_mps": supply_velocity,
        "n_cells": int(mesh.n_cells),
        "n_points": int(mesh.n_points),
        "bounds": [float(v) for v in mesh.bounds],
        "arrays": {name: list(value.shape) for name, value in arrays.items()},
        "qa": qa_row or {},
        "qa_pass_basic": None if not qa_row else bool(qa_row.get("qa_pass_basic")),
        "mesh_yplus_summary_path": str((case_dir / "mesh_yplus_summary.json").relative_to(case_root.parent)) if (case_dir / "mesh_yplus_summary.json").exists() else None,
        "post_plenum_metrics_path": str((case_dir / "post_plenum_metrics.json").relative_to(case_root.parent)) if (case_dir / "post_plenum_metrics.json").exists() else None,
        "key_metrics": {
            "relative_imbalance": get_nested(plenum_metrics, "mass_balance.relative_imbalance_out_minus_in"),
            "supply_inflow_m3s": get_nested(plenum_metrics, "mass_balance.supply_inflow_m3s"),
            "work_zone_Uz_mean": get_nested(plenum_metrics, "regions.work_zone.Uz.mean"),
            "near_car_Uz_mean": get_nested(plenum_metrics, "regions.near_car.Uz.mean"),
            "work_zone_reverse_fraction": get_nested(plenum_metrics, "regions.work_zone.downdraft_uniformity.reverse_flow_fraction_Uz_positive"),
            "near_car_reverse_fraction": get_nested(plenum_metrics, "regions.near_car.downdraft_uniformity.reverse_flow_fraction_Uz_positive"),
            "car_yplus_p95": get_nested(yplus_summary, "carBody_yPlus.p95"),
            "car_yplus_max": get_nested(yplus_summary, "carBody_yPlus.max"),
            "mesh_ok": get_nested(yplus_summary, "checkMesh.mesh_ok"),
        },
    }
    return metadata


def write_manifest(output_root: Path, rows: list[dict[str, Any]]) -> None:
    manifest = {
        "dataset_name": output_root.name,
        "representation": "cell-centered OpenFOAM foamToVTK point cloud / graph sample",
        "sample_count": len(rows),
        "required_arrays": list(REQUIRED_FIELDS),
        "optional_arrays": list(OPTIONAL_FIELDS),
        "samples": rows,
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    csv_fields = [
        "sample_id",
        "case",
        "supply_velocity_mps",
        "n_cells",
        "qa_pass_basic",
        "sample_path",
        "metadata_path",
        "relative_imbalance",
        "work_zone_Uz_mean",
        "near_car_Uz_mean",
        "car_yplus_p95",
        "car_yplus_max",
    ]
    with (output_root / "manifest.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in csv_fields})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, default=Path("cases/paint_booth_l2_flow_qa"))
    parser.add_argument("--output-root", type=Path, default=Path("data/paint_booth_ml_dataset/l2_flow_qa_seed_v0"))
    parser.add_argument("--include", default="l2_flow_qa_*", help="fnmatch pattern for case directories")
    parser.add_argument("--exclude", default=None, help="optional fnmatch exclude pattern")
    parser.add_argument("--fields", nargs="+", default=list(REQUIRED_FIELDS + OPTIONAL_FIELDS), help="cell fields to export when present")
    parser.add_argument("--dry-run", action="store_true", help="discover cases and print manifest preview without writing NPZ files")
    parser.add_argument("--force", action="store_true", help="replace existing output-root")
    args = parser.parse_args()

    repo = Path.cwd().resolve()
    case_root = args.case_root.resolve() if args.case_root.is_absolute() else (repo / args.case_root).resolve()
    output_root = args.output_root.resolve() if args.output_root.is_absolute() else (repo / args.output_root).resolve()

    cases = discover_cases(case_root, args.include, args.exclude)
    if not cases:
        raise SystemExit(f"No cases matching {args.include!r} under {case_root}")

    if output_root.exists() and args.force and not args.dry_run:
        shutil.rmtree(output_root)
    if not args.dry_run:
        (output_root / "samples").mkdir(parents=True, exist_ok=True)
        (output_root / "metadata").mkdir(parents=True, exist_ok=True)

    qa_rows = read_flow_summary(case_root)
    manifest_rows: list[dict[str, Any]] = []

    for case_dir in cases:
        vtk_dir = latest_vtk_dir(case_dir)
        mesh = pv.read(vtk_dir / "internal.vtu")
        arrays = extract_arrays(mesh, tuple(args.fields))
        metadata = build_metadata(
            case_dir=case_dir,
            case_root=case_root,
            vtk_dir=vtk_dir,
            mesh=mesh,
            arrays=arrays,
            qa_row=qa_rows.get(case_dir.name),
        )
        sample_id = metadata["sample_id"]
        sample_rel = Path("samples") / f"{sample_id}.npz"
        meta_rel = Path("metadata") / f"{sample_id}.json"
        metadata["sample_path"] = str(sample_rel)
        metadata["metadata_path"] = str(meta_rel)

        if not args.dry_run:
            np.savez_compressed(output_root / sample_rel, **arrays)
            (output_root / meta_rel).write_text(json.dumps(metadata, indent=2), encoding="utf-8")

        key_metrics = metadata["key_metrics"]
        manifest_rows.append({
            "sample_id": sample_id,
            "case": case_dir.name,
            "supply_velocity_mps": metadata["supply_velocity_mps"],
            "n_cells": metadata["n_cells"],
            "qa_pass_basic": metadata["qa_pass_basic"],
            "sample_path": str(sample_rel),
            "metadata_path": str(meta_rel),
            "relative_imbalance": key_metrics.get("relative_imbalance"),
            "work_zone_Uz_mean": key_metrics.get("work_zone_Uz_mean"),
            "near_car_Uz_mean": key_metrics.get("near_car_Uz_mean"),
            "car_yplus_p95": key_metrics.get("car_yplus_p95"),
            "car_yplus_max": key_metrics.get("car_yplus_max"),
        })

    manifest_rows = sorted(manifest_rows, key=lambda row: float(row["supply_velocity_mps"] or 0.0))

    if args.dry_run:
        print(json.dumps({"sample_count": len(manifest_rows), "samples": manifest_rows}, indent=2), flush=True)
    else:
        write_manifest(output_root, manifest_rows)
        print(f"Exported {len(manifest_rows)} samples to {output_root}", flush=True)
        print(f"Manifest: {output_root / 'manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
