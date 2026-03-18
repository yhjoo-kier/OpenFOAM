#!/usr/bin/env python3
"""Run or collect the benchmark's known stress subset using the existing reference batch flow.

Purpose:
- keep the canonical frozen-bundle stress cases easy to rerun after meshing/solver changes
- avoid hand-maintaining scene lists for regression checks
- optionally refresh aggregate manifests / bundle verification after the subset pass

By default this script reads `benchmark/manifests/dataset_integrity_summary.json`, extracts
cases tagged with `stress_case_*` signals, maps them back to benchmark scene JSONs, and then
invokes `scripts/run_benchmark_reference_batch.py` on just those scenes.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK = PROJECT_ROOT / "benchmark"
MANIFESTS = BENCHMARK / "manifests"
DEFAULT_INTEGRITY = MANIFESTS / "dataset_integrity_summary.json"
DEFAULT_SCENE_MANIFEST = MANIFESTS / "scene_manifest.json"
SCRIPTS = PROJECT_ROOT / "scripts"
DEFAULT_VIEWS = ["perspective", "birdseye", "floorplan", "wireframe", "section"]


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True)


def tail(text: str, lines: int = 40) -> str:
    return "\n".join((text or "").splitlines()[-lines:])


def resolve_scene_selection(integrity_path: Path, scene_manifest_path: Path) -> tuple[list[Path], list[dict[str, Any]]]:
    integrity = load_json(integrity_path)
    scene_manifest = load_json(scene_manifest_path)
    scene_by_case = {row["case_name"]: Path(row["scene_file"]) for row in scene_manifest}

    stress_rows = integrity.get("stress_cases", [])
    scene_paths: list[Path] = []
    selected_rows: list[dict[str, Any]] = []
    for row in stress_rows:
        case_name = row["case_name"]
        scene_path = scene_by_case.get(case_name)
        if scene_path is None:
            raise FileNotFoundError(f"Stress case {case_name} missing from {scene_manifest_path}")
        scene_paths.append(scene_path)
        selected_rows.append(row)
    return scene_paths, selected_rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the benchmark stress subset through the reference CFD batch flow")
    parser.add_argument("--integrity-summary", type=Path, default=DEFAULT_INTEGRITY)
    parser.add_argument("--scene-manifest", type=Path, default=DEFAULT_SCENE_MANIFEST)
    parser.add_argument("--mesh-size", type=float, default=0.35)
    parser.add_argument("--end-time", type=int, default=250)
    parser.add_argument("--solver-timeout", type=int, default=600)
    parser.add_argument("--name-prefix", default="bench")
    parser.add_argument("--views", nargs="+", default=DEFAULT_VIEWS, choices=DEFAULT_VIEWS)
    parser.add_argument("--collect-only", action="store_true", help="Refresh subset bookkeeping without rerunning CFD")
    parser.add_argument("--skip-existing-success", action="store_true", help="Pass through to batch runner")
    parser.add_argument("--refresh-full-manifests", action="store_true", help="After the subset pass, rebuild full frozen-set aggregate manifests via collect-only + evaluation scaffold refresh")
    parser.add_argument("--verify-after", action="store_true", help="Run verify_benchmark_bundle.py after the subset/full refresh sequence")
    parser.add_argument("--output", type=Path, default=MANIFESTS / "stress_subset_run_summary.json")
    args = parser.parse_args()

    scene_paths, stress_rows = resolve_scene_selection(args.integrity_summary, args.scene_manifest)
    if not scene_paths:
        raise SystemExit("No stress cases found in the integrity summary")

    batch_cmd = [
        "python3", str(SCRIPTS / "run_benchmark_reference_batch.py"),
        "--scenes", *[str(path) for path in scene_paths],
        "--mesh-size", str(args.mesh_size),
        "--end-time", str(args.end_time),
        "--solver-timeout", str(args.solver_timeout),
        "--name-prefix", args.name_prefix,
        "--views", *args.views,
    ]
    if args.collect_only:
        batch_cmd.append("--collect-only")
    if args.skip_existing_success:
        batch_cmd.append("--skip-existing-success")

    subset_proc = run(batch_cmd)

    refresh_proc = None
    scaffold_proc = None
    verify_proc = None
    if args.refresh_full_manifests:
        refresh_cmd = ["python3", str(SCRIPTS / "run_benchmark_reference_batch.py"), "--collect-only"]
        refresh_proc = run(refresh_cmd)
        scaffold_cmd = ["python3", str(SCRIPTS / "scaffold_benchmark_evaluations.py")]
        scaffold_proc = run(scaffold_cmd)

    if args.verify_after:
        verify_cmd = ["python3", str(SCRIPTS / "verify_benchmark_bundle.py")]
        verify_proc = run(verify_cmd)

    payload = {
        "ok": (
            subset_proc.returncode == 0
            and (refresh_proc is None or refresh_proc.returncode == 0)
            and (scaffold_proc is None or scaffold_proc.returncode == 0)
            and (verify_proc is None or verify_proc.returncode == 0)
        ),
        "selected_case_count": len(stress_rows),
        "selected_cases": stress_rows,
        "selected_scene_paths": [str(path) for path in scene_paths],
        "subset_command": batch_cmd,
        "subset_returncode": subset_proc.returncode,
        "subset_stdout_tail": tail(subset_proc.stdout),
        "subset_stderr_tail": tail(subset_proc.stderr),
        "refreshed_full_manifests": bool(args.refresh_full_manifests),
        "refresh_returncode": None if refresh_proc is None else refresh_proc.returncode,
        "refresh_stdout_tail": None if refresh_proc is None else tail(refresh_proc.stdout),
        "refresh_stderr_tail": None if refresh_proc is None else tail(refresh_proc.stderr),
        "scaffold_returncode": None if scaffold_proc is None else scaffold_proc.returncode,
        "scaffold_stdout_tail": None if scaffold_proc is None else tail(scaffold_proc.stdout),
        "scaffold_stderr_tail": None if scaffold_proc is None else tail(scaffold_proc.stderr),
        "verified_after": bool(args.verify_after),
        "verify_returncode": None if verify_proc is None else verify_proc.returncode,
        "verify_stdout_tail": None if verify_proc is None else tail(verify_proc.stdout),
        "verify_stderr_tail": None if verify_proc is None else tail(verify_proc.stderr),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
