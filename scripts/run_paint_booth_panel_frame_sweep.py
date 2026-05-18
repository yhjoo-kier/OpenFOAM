#!/usr/bin/env python3
"""Run a small calibration sweep for the paint-booth panel-frame case.

This host-side script creates panel-frame cases with different supply velocities
and filter Forchheimer coefficients, runs the standard OpenFOAM smoke workflow in
Docker, post-processes metrics, and writes an aggregate JSON/CSV table.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any


def run(cmd: list[str], *, cwd: Path, log: Path | None = None, check: bool = True) -> subprocess.CompletedProcess[str]:
    print("+", " ".join(map(str, cmd)), flush=True)
    cp = subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(cp.stdout, encoding="utf-8")
    if check and cp.returncode != 0:
        tail = "\n".join(cp.stdout.splitlines()[-80:])
        raise RuntimeError(f"Command failed with exit {cp.returncode}: {' '.join(cmd)}\n--- tail ---\n{tail}")
    return cp


def docker_case_cmd(repo: Path, image: str, case_dir: Path, shell_script: str) -> list[str]:
    return [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{repo}:{repo}",
        "-w",
        str(case_dir),
        image,
        "bash",
        "-lc",
        shell_script,
    ]


def flatten_metrics(case_name: str, supply_velocity: float, filter_f: float, metrics: dict[str, Any]) -> dict[str, Any]:
    def get(path: str, default: Any = None) -> Any:
        cur: Any = metrics
        for part in path.split("."):
            if not isinstance(cur, dict) or part not in cur:
                return default
            cur = cur[part]
        return cur

    dp_kin = get("filter_pressure_drop_proxy_m2_per_s2")
    return {
        "case": case_name,
        "supply_velocity_mps": supply_velocity,
        "filter_forchheimer": filter_f,
        "supply_inflow_m3s": get("mass_balance.supply_inflow_m3s"),
        "floor_outflow_m3s": get("mass_balance.floor_outflow_m3s"),
        "relative_imbalance": get("mass_balance.relative_imbalance_out_minus_in"),
        "filter_dp_kinematic_m2s2": dp_kin,
        "filter_dp_pa_rho1p2": None if dp_kin is None else dp_kin * 1.2,
        "filter_below_Uz_mean": get("regions.filter_below.Uz.mean"),
        "filter_below_Uz_cv": get("regions.filter_below.downdraft_uniformity.Uz_std_over_abs_mean"),
        "work_zone_Uz_mean": get("regions.work_zone.Uz.mean"),
        "work_zone_Uz_cv": get("regions.work_zone.downdraft_uniformity.Uz_std_over_abs_mean"),
        "work_zone_reverse_fraction": get("regions.work_zone.downdraft_uniformity.reverse_flow_fraction_Uz_positive"),
        "near_car_Uz_mean": get("regions.near_car.Uz.mean"),
        "near_car_Uz_cv": get("regions.near_car.downdraft_uniformity.Uz_std_over_abs_mean"),
        "near_car_reverse_fraction": get("regions.near_car.downdraft_uniformity.reverse_flow_fraction_Uz_positive"),
        "filter_plane_downflow_m3s": get("virtual_planes.filter_panel_below_z2p94.normal_flow_m3s"),
        "filter_plane_Udown_mps": get("virtual_planes.filter_panel_below_z2p94.area_weighted_normal_velocity_mps"),
        "work_midplane_downflow_m3s": get("virtual_planes.work_zone_mid_z1p50.normal_flow_m3s"),
        "work_midplane_Udown_mps": get("virtual_planes.work_zone_mid_z1p50.area_weighted_normal_velocity_mps"),
    }


def score_row(row: dict[str, Any], target_near_car: float, target_dp_pa: float) -> float:
    near = abs(float(row["near_car_Uz_mean"]))
    work = abs(float(row["work_zone_Uz_mean"]))
    dp = float(row["filter_dp_pa_rho1p2"])
    rev = float(row["near_car_reverse_fraction"])
    # Weighted heuristic: near-car target first, then work-zone reasonableness,
    # pressure-drop target, and reverse flow penalty.
    return (
        3.0 * abs(near - target_near_car) / max(target_near_car, 1e-12)
        + 1.0 * abs(work - target_near_car) / max(target_near_car, 1e-12)
        + 1.2 * abs(dp - target_dp_pa) / max(target_dp_pa, 1e-12)
        + 2.0 * rev
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--image", default="openfoam-pipeline-local:latest")
    ap.add_argument("--supply", type=float, nargs="+", default=[2.3, 3.0, 3.5, 4.36])
    ap.add_argument("--filter-f", type=float, nargs="+", default=[6800.0, 10000.0, 20000.0, 40000.0])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--target-near-car", type=float, default=0.30)
    ap.add_argument("--target-dp-pa", type=float, default=50.0)
    args = ap.parse_args()

    repo = args.repo.resolve()
    sweep_root = repo / "cases" / "paint_booth_panel_frame_sweep"
    sweep_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []

    for supply_velocity in args.supply:
        for filter_f in args.filter_f:
            tag = f"u{supply_velocity:.2f}_f{filter_f:.0f}".replace(".", "p")
            case_name = f"panel_frame_{tag}"
            case_dir = sweep_root / case_name
            if args.force and case_dir.exists():
                shutil.rmtree(case_dir)
            metrics_path = case_dir / "post_plenum_metrics.json"
            if metrics_path.exists() and not args.force:
                print(f"Skipping existing completed case: {case_dir}", flush=True)
            else:
                run(
                    [
                        "python3",
                        "scripts/create_paint_booth_plenum_filter_case.py",
                        "--central-filter-panel-frame",
                        "--case-dir",
                        str(case_dir.relative_to(repo)),
                        "--supply-velocity",
                        str(supply_velocity),
                        "--filter-forchheimer",
                        str(filter_f),
                        "--force",
                    ],
                    cwd=repo,
                    log=case_dir / "log.generate",
                )
                foam_script = """
set -e
surfaceCheck constant/triSurface/simplified_car_shell.stl > log.surfaceCheck 2>&1
blockMesh > log.blockMesh 2>&1
snappyHexMesh -overwrite > log.snappyHexMesh 2>&1
checkMesh > log.checkMesh_snappy 2>&1
topoSet > log.topoSet 2>&1
timeout 240 simpleFoam > log.simpleFoam 2>&1
foamToVTK -latestTime > log.foamToVTK 2>&1
"""
                run(docker_case_cmd(repo, args.image, case_dir, foam_script), cwd=repo, log=case_dir / "log.docker_pipeline")
                run(
                    [
                        "docker",
                        "run",
                        "--rm",
                        "-v",
                        f"{repo}:{repo}",
                        "-w",
                        str(repo),
                        args.image,
                        "bash",
                        "-lc",
                        f"python3 scripts/postprocess_paint_booth_plenum_metrics.py {case_dir.relative_to(repo)} > {case_dir.relative_to(repo)}/log.postprocess 2>&1",
                    ],
                    cwd=repo,
                    log=case_dir / "log.postprocess_docker",
                )
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            row = flatten_metrics(case_name, supply_velocity, filter_f, metrics)
            row["score"] = score_row(row, args.target_near_car, args.target_dp_pa)
            rows.append(row)
            print(json.dumps(row, indent=2), flush=True)

    rows.sort(key=lambda r: r["score"])
    out_json = sweep_root / "sweep_summary.json"
    out_csv = sweep_root / "sweep_summary.csv"
    out_json.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print("\nBest candidates:")
    print(json.dumps(rows[:5], indent=2))
    print(f"Wrote {out_json}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
