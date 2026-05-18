#!/usr/bin/env python3
"""Run a small L2-quality flow-rate QA sweep for the paint-booth case.

This is not a large training-data generator.  It is a gatekeeping sweep meant to
answer whether the current L2 CFD baseline behaves consistently across a few
representative supply-flow conditions before producing GINO/neural-operator
training snapshots.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
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
        tail = "\n".join(cp.stdout.splitlines()[-120:])
        raise RuntimeError(f"Command failed with exit {cp.returncode}: {' '.join(cmd)}\n--- tail ---\n{tail}")
    return cp


def docker_cmd(repo: Path, image: str, workdir: Path, shell_script: str) -> list[str]:
    return [
        "docker", "run", "--rm",
        "-v", f"{repo}:{repo}",
        "-w", str(workdir),
        image,
        "bash", "-lc", shell_script,
    ]


def get_nested(data: dict[str, Any], path: str, default: Any = None) -> Any:
    cur: Any = data
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def parse_residual_tail(case_dir: Path) -> dict[str, Any]:
    log = case_dir / "log.simpleFoam"
    if not log.exists():
        return {"exists": False}
    text = log.read_text(errors="replace")
    out: dict[str, Any] = {
        "exists": True,
        "ended": "End" in text,
        "final_time": None,
    }
    times = re.findall(r"^Time =\s*([0-9.eE+-]+)", text, flags=re.M)
    if times:
        try:
            out["final_time"] = float(times[-1])
        except ValueError:
            out["final_time"] = times[-1]
    for field in ["Ux", "Uy", "Uz", "p", "k", "omega"]:
        matches = re.findall(
            rf"Solving for {field}, Initial residual = ([0-9.eE+-]+), Final residual = ([0-9.eE+-]+), No Iterations ([0-9]+)",
            text,
        )
        if matches:
            init, final, niter = matches[-1]
            out[f"{field}_initial_residual_last"] = float(init)
            out[f"{field}_final_residual_last"] = float(final)
            out[f"{field}_iterations_last"] = int(niter)
    cont = re.findall(r"time step continuity errors : sum local = ([0-9.eE+-]+), global = ([0-9.eE+-]+), cumulative = ([0-9.eE+-]+)", text)
    if cont:
        local, global_, cumulative = cont[-1]
        out["continuity_sum_local_last"] = float(local)
        out["continuity_global_last"] = float(global_)
        out["continuity_cumulative_last"] = float(cumulative)
    return out


def flatten_row(case_name: str, supply_velocity: float, metrics: dict[str, Any], yplus: dict[str, Any], residuals: dict[str, Any]) -> dict[str, Any]:
    dp_kin = get_nested(metrics, "filter_pressure_drop_proxy_m2_per_s2")
    mesh = yplus.get("checkMesh", {}) if isinstance(yplus, dict) else {}
    layers = yplus.get("snappy_layers", {}) if isinstance(yplus, dict) else {}
    car_y = yplus.get("carBody_yPlus", {}) if isinstance(yplus, dict) else {}
    row = {
        "case": case_name,
        "supply_velocity_mps": supply_velocity,
        "supply_inflow_m3s": get_nested(metrics, "mass_balance.supply_inflow_m3s"),
        "floor_outflow_m3s": get_nested(metrics, "mass_balance.floor_outflow_m3s"),
        "relative_imbalance": get_nested(metrics, "mass_balance.relative_imbalance_out_minus_in"),
        "filter_dp_pa_rho1p2": None if dp_kin is None else float(dp_kin) * 1.2,
        "filter_below_Uz_mean": get_nested(metrics, "regions.filter_below.Uz.mean"),
        "filter_below_Uz_cv": get_nested(metrics, "regions.filter_below.downdraft_uniformity.Uz_std_over_abs_mean"),
        "work_zone_Uz_mean": get_nested(metrics, "regions.work_zone.Uz.mean"),
        "work_zone_Uz_cv": get_nested(metrics, "regions.work_zone.downdraft_uniformity.Uz_std_over_abs_mean"),
        "work_zone_reverse_fraction": get_nested(metrics, "regions.work_zone.downdraft_uniformity.reverse_flow_fraction_Uz_positive"),
        "near_car_Uz_mean": get_nested(metrics, "regions.near_car.Uz.mean"),
        "near_car_Uz_cv": get_nested(metrics, "regions.near_car.downdraft_uniformity.Uz_std_over_abs_mean"),
        "near_car_reverse_fraction": get_nested(metrics, "regions.near_car.downdraft_uniformity.reverse_flow_fraction_Uz_positive"),
        "filter_plane_downflow_m3s": get_nested(metrics, "virtual_planes.filter_panel_below_z2p94.normal_flow_m3s"),
        "work_midplane_downflow_m3s": get_nested(metrics, "virtual_planes.work_zone_mid_z1p50.normal_flow_m3s"),
        "mesh_ok": mesh.get("mesh_ok"),
        "cells": mesh.get("cells"),
        "max_aspect_ratio": mesh.get("max_aspect_ratio"),
        "max_skewness": mesh.get("max_skewness"),
        "car_layer_coverage_percent": layers.get("carBody_layer_coverage_percent"),
        "car_yplus_median": car_y.get("p50"),
        "car_yplus_p95": car_y.get("p95"),
        "car_yplus_max": car_y.get("max"),
        "simpleFoam_ended": residuals.get("ended"),
        "simpleFoam_final_time": residuals.get("final_time"),
        "continuity_global_last": residuals.get("continuity_global_last"),
        "p_final_residual_last": residuals.get("p_final_residual_last"),
        "Ux_final_residual_last": residuals.get("Ux_final_residual_last"),
        "Uy_final_residual_last": residuals.get("Uy_final_residual_last"),
        "Uz_final_residual_last": residuals.get("Uz_final_residual_last"),
    }
    row["qa_pass_basic"] = bool(
        row["mesh_ok"] is True
        and row["simpleFoam_ended"] is True
        and row["relative_imbalance"] is not None
        and abs(float(row["relative_imbalance"])) < 0.01
        and row["car_yplus_p95"] is not None
        and float(row["car_yplus_p95"]) < 5.0
        and row["car_layer_coverage_percent"] is not None
        and float(row["car_layer_coverage_percent"]) > 95.0
    )
    return row


def monotonic_checks(rows: list[dict[str, Any]]) -> dict[str, Any]:
    ordered = sorted(rows, key=lambda r: float(r["supply_velocity_mps"]))

    def values(key: str) -> list[float]:
        return [float(r[key]) for r in ordered if r.get(key) is not None]

    def strictly_increasing(vals: list[float]) -> bool:
        return all(b > a for a, b in zip(vals, vals[1:]))

    def strictly_decreasing(vals: list[float]) -> bool:
        return all(b < a for a, b in zip(vals, vals[1:]))

    checks = {
        "supply_inflow_increases_with_supply_velocity": strictly_increasing(values("supply_inflow_m3s")),
        "filter_dp_increases_with_supply_velocity": strictly_increasing(values("filter_dp_pa_rho1p2")),
        "abs_work_zone_downdraft_increases_with_supply_velocity": strictly_increasing([abs(v) for v in values("work_zone_Uz_mean")]),
        "abs_near_car_downdraft_increases_with_supply_velocity": strictly_increasing([abs(v) for v in values("near_car_Uz_mean")]),
        "work_zone_Uz_more_negative_with_supply_velocity": strictly_decreasing(values("work_zone_Uz_mean")),
        "near_car_Uz_more_negative_with_supply_velocity": strictly_decreasing(values("near_car_Uz_mean")),
    }
    return checks


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--image", default="openfoam-pipeline-local:latest")
    ap.add_argument("--root", type=Path, default=Path("cases/paint_booth_l2_flow_qa"))
    ap.add_argument("--supply", type=float, nargs="+", default=[2.18, 4.36, 6.54], help="Representative supply velocities [m/s]. Default: 0.5x/1.0x/1.5x nominal.")
    ap.add_argument("--filter-forchheimer", type=float, default=10000.0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    repo = args.repo.resolve()
    root = (repo / args.root).resolve() if not args.root.is_absolute() else args.root
    root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []

    for supply in args.supply:
        tag = f"u{supply:.2f}".replace(".", "p")
        case_name = f"l2_flow_qa_{tag}"
        case_dir = root / case_name
        if args.force and case_dir.exists():
            shutil.rmtree(case_dir)
        metrics_path = case_dir / "post_plenum_metrics.json"
        yplus_path = case_dir / "mesh_yplus_summary.json"
        if not (metrics_path.exists() and yplus_path.exists()) or args.force:
            run([
                "python3", "scripts/create_paint_booth_plenum_filter_case.py",
                "--central-filter-panel-frame",
                "--case-dir", str(case_dir.relative_to(repo)),
                "--supply-velocity", str(supply),
                "--filter-forchheimer", str(args.filter_forchheimer),
                "--cell-size", "0.125",
                "--filter-z-cells", "6",
                "--car-refinement-min", "3",
                "--car-refinement-max", "4",
                "--add-layers",
                "--n-surface-layers", "5",
                "--absolute-layer-sizes",
                "--final-layer-thickness", "0.004",
                "--min-layer-thickness", "0.0003",
                "--expansion-ratio", "1.2",
                "--force",
            ], cwd=repo, log=case_dir / "log.generate")
            foam_script = """
set -e
surfaceCheck constant/triSurface/simplified_car_shell.stl > log.surfaceCheck 2>&1
blockMesh > log.blockMesh 2>&1
snappyHexMesh -overwrite > log.snappyHexMesh 2>&1
checkMesh > log.checkMesh_snappy 2>&1
topoSet > log.topoSet 2>&1
simpleFoam > log.simpleFoam 2>&1
simpleFoam -postProcess -func yPlus -latestTime > log.yPlus 2>&1
foamToVTK -latestTime > log.foamToVTK 2>&1
"""
            run(docker_cmd(repo, args.image, case_dir, foam_script), cwd=repo, log=case_dir / "log.docker_pipeline")
            run(docker_cmd(repo, args.image, repo, f"python3 scripts/postprocess_paint_booth_plenum_metrics.py {case_dir.relative_to(repo)} > {case_dir.relative_to(repo)}/log.postprocess 2>&1"), cwd=repo, log=case_dir / "log.postprocess_docker")
            run(docker_cmd(repo, args.image, repo, f"python3 scripts/postprocess_mesh_yplus_summary.py {case_dir.relative_to(repo)} > {case_dir.relative_to(repo)}/log.mesh_yplus_summary 2>&1"), cwd=repo, log=case_dir / "log.mesh_yplus_summary_docker")
        else:
            print(f"Skipping completed case: {case_dir}", flush=True)

        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        yplus = json.loads(yplus_path.read_text(encoding="utf-8"))
        residuals = parse_residual_tail(case_dir)
        row = flatten_row(case_name, supply, metrics, yplus, residuals)
        rows.append(row)
        print(json.dumps(row, indent=2), flush=True)

    rows = sorted(rows, key=lambda r: float(r["supply_velocity_mps"]))
    summary = {
        "purpose": "Three-point L2 flow-rate QA before large CFD data generation for neural-operator training.",
        "flow_conditions": rows,
        "monotonic_checks": monotonic_checks(rows),
        "all_basic_qa_pass": all(bool(r.get("qa_pass_basic")) for r in rows),
    }
    (root / "flow_qa_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (root / "flow_qa_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2), flush=True)
    print(f"Wrote {root / 'flow_qa_summary.json'}")
    print(f"Wrote {root / 'flow_qa_summary.csv'}")


if __name__ == "__main__":
    main()
