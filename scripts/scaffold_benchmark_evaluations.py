#!/usr/bin/env python3
"""Scaffold benchmark evaluation tasks from the frozen reference bundle.

This prepares canonical per-case/per-view evaluation folders so image-conditioned
pipeline runs can be tracked against the frozen benchmark reference set.

Inputs:
- benchmark/manifests/scene_manifest.json
- benchmark/manifests/reference_batch_summary.json (or compatible aggregate reference-status manifest)
- benchmark/renderings/renderings_manifest.json

Outputs:
- benchmark/evaluations/manifest.json
- benchmark/evaluations/summary.json
- benchmark/evaluations/<case>/<view>/task.json
- benchmark/evaluations/<case>/<view>/{input.png,reference_scene.json,reference_case,reference_results}
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK = PROJECT_ROOT / "benchmark"
MANIFESTS = BENCHMARK / "manifests"
RENDERINGS = BENCHMARK / "renderings"
EVALUATIONS = BENCHMARK / "evaluations"

SCENE_MANIFEST = MANIFESTS / "scene_manifest.json"
DEFAULT_REFERENCE_STATUS = MANIFESTS / "reference_batch_summary.json"
RENDERINGS_MANIFEST = RENDERINGS / "renderings_manifest.json"

DEFAULT_VIEWS = ["perspective", "birdseye", "floorplan", "wireframe", "section"]


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def ensure_link(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    dst.symlink_to(src)


def main() -> int:
    parser = argparse.ArgumentParser(description="Scaffold benchmark evaluation tasks from the frozen reference bundle")
    parser.add_argument("--scene-manifest", type=Path, default=SCENE_MANIFEST)
    parser.add_argument("--reference-status", type=Path, default=DEFAULT_REFERENCE_STATUS)
    parser.add_argument("--renderings-manifest", type=Path, default=RENDERINGS_MANIFEST)
    args = parser.parse_args()

    scene_rows = load_json(args.scene_manifest)
    ref_payload = load_json(args.reference_status)
    render_payload = load_json(args.renderings_manifest)

    ref_by_case = {row["case_name"]: row for row in ref_payload["results"]}
    render_by_case = {row["case_name"]: row for row in render_payload["cases"]}

    EVALUATIONS.mkdir(parents=True, exist_ok=True)

    tasks: list[dict[str, Any]] = []
    case_summaries: list[dict[str, Any]] = []

    for scene_row in scene_rows:
        case_name = scene_row["case_name"]
        ref_row = ref_by_case.get(case_name)
        render_row = render_by_case.get(case_name)
        if ref_row is None or render_row is None:
            raise RuntimeError(f"Missing reference/render manifest data for {case_name}")
        if not ref_row.get("success"):
            raise RuntimeError(f"Reference CFD is not successful for {case_name}; refusing to scaffold evaluation tasks")

        case_dir = EVALUATIONS / case_name
        case_dir.mkdir(parents=True, exist_ok=True)

        scene_src = Path(scene_row["scene_file"])
        if not scene_src.is_absolute():
            scene_src = PROJECT_ROOT / scene_src
        ref_case = PROJECT_ROOT / "cases" / case_name
        ref_results = PROJECT_ROOT / "results" / case_name

        views_present = []
        for view in DEFAULT_VIEWS:
            input_image = render_row["views"].get(view)
            if not input_image:
                continue
            views_present.append(view)
            task_dir = case_dir / view
            task_dir.mkdir(parents=True, exist_ok=True)

            input_path = Path(input_image)
            ensure_link(input_path, task_dir / "input.png")
            ensure_link(scene_src, task_dir / "reference_scene.json")
            ensure_link(ref_case, task_dir / "reference_case")
            ensure_link(ref_results, task_dir / "reference_results")

            task_path = task_dir / "task.json"
            existing_task = load_json(task_path) if task_path.exists() else {}
            task = {
                "case_name": case_name,
                "category": scene_row["category"],
                "variant_index": scene_row["variant_index"],
                "room_kind": scene_row["room_kind"],
                "obstacle_count": scene_row["obstacle_count"],
                "view": view,
                "status": existing_task.get("status", "pending"),
                "input_image": str(task_dir / "input.png"),
                "reference_scene": str(task_dir / "reference_scene.json"),
                "reference_case": str(task_dir / "reference_case"),
                "reference_results": str(task_dir / "reference_results"),
                "expected_outputs": {
                    "predicted_scene_json": str(task_dir / "predicted_scene.json"),
                    "predicted_case_dir": str(task_dir / "predicted_case"),
                    "predicted_results_dir": str(task_dir / "predicted_results"),
                    "evaluation_summary_json": str(task_dir / "evaluation_summary.json"),
                },
                "reference_summary": {
                    "success": ref_row["success"],
                    "successful_preset": ref_row.get("preset"),
                    "successful_mode": ref_row.get("mode"),
                    "successful_mesh_size": ref_row.get("mesh_size"),
                    "attempt_count": ref_row.get("attempt_count"),
                    "risk": ref_row.get("risk"),
                    "risk_score": ref_row.get("risk_score"),
                    "risk_reasons": ref_row.get("risk_reasons"),
                },
            }
            for key in ("last_started_at", "last_finished_at", "last_run_name", "evaluation_summary", "actual_outputs", "run_request"):
                if key in existing_task:
                    task[key] = existing_task[key]
            task_path.write_text(json.dumps(task, indent=2), encoding="utf-8")
            tasks.append(task)

        case_summary = {
            "case_name": case_name,
            "category": scene_row["category"],
            "room_kind": scene_row["room_kind"],
            "obstacle_count": scene_row["obstacle_count"],
            "views": views_present,
            "reference_scene": str(scene_src),
            "reference_case": str(ref_case),
            "reference_results": str(ref_results),
            "task_count": len(views_present),
            "case_root": str(case_dir),
        }
        (case_dir / "manifest.json").write_text(json.dumps(case_summary, indent=2), encoding="utf-8")
        case_summaries.append(case_summary)

    status_counts: dict[str, int] = {}
    for task in tasks:
        status = task.get("status", "unknown")
        status_counts[status] = status_counts.get(status, 0) + 1

    summary = {
        "ok": True,
        "benchmark_root": str(BENCHMARK),
        "evaluation_root": str(EVALUATIONS),
        "views": DEFAULT_VIEWS,
        "case_count": len(case_summaries),
        "task_count": len(tasks),
        "status_counts": status_counts,
        "cases": case_summaries,
    }
    manifest = {
        "ok": True,
        "cases": case_summaries,
        "tasks": tasks,
    }

    (EVALUATIONS / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (EVALUATIONS / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
