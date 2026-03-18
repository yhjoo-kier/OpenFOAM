#!/usr/bin/env python3
"""Batch runner for benchmark scene -> OpenFOAM reference CFD -> visualization.

This wraps the project-local indoor CFD pipeline for pre-generated benchmark scene JSONs.
It runs each scene through:
- stabilization pipeline (`run_indoor_stabilized.py`)
- benchmark input-view rendering (`render_benchmark_views.py`)
- summary collection
- artifact indexing

Outputs:
- `benchmark/reference_cfd/<case_name>/...` via symlink/copy-friendly canonical paths
- `benchmark/visualizations/<case_name>/...`
- `benchmark/renderings/<case_name>/...`
- `benchmark/manifests/reference_batch_summary.json`
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
BENCHMARK = PROJECT_ROOT / "benchmark"
SCENES_DIR = BENCHMARK / "scenes"
REFERENCE_DIR = BENCHMARK / "reference_cfd"
VIS_DIR = BENCHMARK / "visualizations"
RENDER_DIR = BENCHMARK / "renderings"
MANIFEST_DIR = BENCHMARK / "manifests"


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True)


def ensure_link(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    dst.symlink_to(src)


def render_input_views(scene_path: Path, views: list[str]) -> subprocess.CompletedProcess[str]:
    cmd = [
        "python3", str(SCRIPTS / "render_benchmark_views.py"),
        "--scenes", str(scene_path),
        "--render-root", str(RENDER_DIR),
        "--views", *views,
    ]
    return run(cmd)


def collect_case(case_name: str) -> dict[str, Any]:
    results_dir = PROJECT_ROOT / "results" / case_name
    case_dir = PROJECT_ROOT / "cases" / case_name
    summary_path = results_dir / "stabilization_summary.json"
    if not summary_path.exists():
        return {
            "case_name": case_name,
            "success": False,
            "reason": "missing_summary",
            "case_dir": str(case_dir),
            "results_dir": str(results_dir),
        }

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    ref_link = REFERENCE_DIR / case_name
    vis_link = VIS_DIR / case_name
    if case_dir.exists():
        ensure_link(case_dir, ref_link)
    if results_dir.exists():
        ensure_link(results_dir, vis_link)

    comparison = results_dir / "comparison_1x2.png"
    render3d = results_dir / "indoor_pipeline_3d_comparison.png"
    panel_geo = results_dir / "panel_geometry_3d.png"
    panel_flow = results_dir / "panel_flow_3d.png"
    render_manifest = RENDER_DIR / case_name / "manifest.json"

    return {
        "case_name": case_name,
        "success": bool(summary.get("success")),
        "scene_json": summary.get("scene_json"),
        "used_scene_json": summary.get("used_scene_json"),
        "successful_preset": summary.get("successful_preset"),
        "successful_mode": summary.get("successful_mode"),
        "successful_mesh_size": summary.get("successful_mesh_size"),
        "successful_scene_variant": summary.get("successful_scene_variant"),
        "case_dir": str(case_dir),
        "results_dir": str(results_dir),
        "reference_link": str(ref_link),
        "visualization_link": str(vis_link),
        "comparison_1x2": str(comparison) if comparison.exists() else None,
        "render3d": str(render3d) if render3d.exists() else None,
        "panel_geometry_3d": str(panel_geo) if panel_geo.exists() else None,
        "panel_flow_3d": str(panel_flow) if panel_flow.exists() else None,
        "benchmark_render_manifest": str(render_manifest) if render_manifest.exists() else None,
        "attempt_count": len(summary.get("attempts", [])),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run benchmark scenes through the indoor CFD reference pipeline")
    parser.add_argument("--scenes", nargs="*", default=None, help="Scene JSON paths; default is benchmark/scenes/*.json")
    parser.add_argument("--mesh-size", type=float, default=0.35)
    parser.add_argument("--end-time", type=int, default=250)
    parser.add_argument("--solver-timeout", type=int, default=600)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--name-prefix", default="bench")
    parser.add_argument("--collect-only", action="store_true", help="Refresh manifests/symlinks from existing case/results artifacts without rerunning CFD")
    parser.add_argument("--skip-existing-success", action="store_true", help="Skip CFD reruns for scenes whose stabilization summary already reports success")
    parser.add_argument(
        "--views",
        nargs="+",
        default=["perspective", "birdseye", "floorplan", "wireframe"],
        choices=["perspective", "birdseye", "floorplan", "wireframe"],
        help="Benchmark input-view families to export alongside the CFD reference run",
    )
    args = parser.parse_args()

    scene_paths = [Path(p) for p in args.scenes] if args.scenes else sorted(SCENES_DIR.glob("*.json"))
    if args.limit is not None:
        scene_paths = scene_paths[: args.limit]

    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    VIS_DIR.mkdir(parents=True, exist_ok=True)
    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)

    batch_results: list[dict[str, Any]] = []
    for scene_path in scene_paths:
        case_name = f"{args.name_prefix}_{scene_path.stem}"
        existing = collect_case(case_name)
        already_success = bool(existing.get("success"))
        should_run = not args.collect_only and not (args.skip_existing_success and already_success)

        if should_run:
            cmd = [
                "python3", str(SCRIPTS / "run_indoor_stabilized.py"),
                "--scenario", str(scene_path),
                "--name", case_name,
                "--mesh-size", str(args.mesh_size),
                "--end-time", str(args.end_time),
                "--solver-timeout", str(args.solver_timeout),
            ]
            proc = run(cmd)
            render_proc = render_input_views(scene_path, args.views)
        else:
            cmd = None
            proc = subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr="")
            render_proc = subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr="")

        entry = collect_case(case_name)
        entry.update({
            "command": cmd,
            "returncode": proc.returncode,
            "stdout_tail": "\n".join((proc.stdout or "").splitlines()[-20:]),
            "stderr_tail": "\n".join((proc.stderr or "").splitlines()[-20:]),
            "source_scene": str(scene_path),
            "render_views": args.views,
            "render_returncode": render_proc.returncode,
            "render_stdout_tail": "\n".join((render_proc.stdout or "").splitlines()[-20:]),
            "render_stderr_tail": "\n".join((render_proc.stderr or "").splitlines()[-20:]),
            "collect_only": args.collect_only,
            "skipped_existing_success": bool(args.skip_existing_success and already_success and not args.collect_only),
        })
        batch_results.append(entry)

    payload = {
        "ok": True,
        "count": len(batch_results),
        "success_count": sum(1 for x in batch_results if x.get("success")),
        "views": args.views,
        "results": batch_results,
    }
    out = MANIFEST_DIR / "reference_batch_summary.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
