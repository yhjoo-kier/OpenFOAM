#!/usr/bin/env python3
"""Scaffold benchmark evaluation tasks from the frozen reference bundle.

This prepares canonical per-case/per-view evaluation folders so image-conditioned
pipeline runs can be tracked against the frozen benchmark reference set.

Inputs:
- benchmark/manifests/scene_manifest.json
- benchmark/manifests/reference_batch_summary.json (or compatible aggregate reference-status manifest)
- benchmark/renderings/renderings_manifest.json

Outputs:
- <evaluation-root>/manifest.json
- <evaluation-root>/summary.json
- <evaluation-root>/<case>/<view>/task.json
- <evaluation-root>/<case>/<view>/{input.png,reference_scene.json,reference_case,reference_results}
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
DEFAULT_SETTING = "no_scale_hint_baseline"
SCALE_HINTED_SETTING = "scale_hinted_longest_horizontal_span_v1"
SCALE_HINTED_DUAL_SETTING = "scale_hinted_longest_span_plus_height_v1"
SCALE_HINTED_LAYOUT_PROTECTED_SETTING = "scale_hinted_longest_span_layout_protected_v1"
SCALE_HINTED_VIEW_GUARDED_SETTING = "scale_hinted_longest_span_view_guarded_v1"
SCALE_HINTED_GUARD_WEIGHTED_SETTING = "scale_hinted_longest_span_guard_weighted_v1"

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


def room_blocks(scene: dict[str, Any]) -> list[dict[str, Any]]:
    room = scene.get("room", {})
    if "blocks" in room:
        return room["blocks"]
    size = room.get("size", {})
    return [{
        "origin": {"x": 0.0, "y": 0.0, "z": 0.0},
        "size": {"dx": size.get("Lx", 0.0), "dy": size.get("Ly", 0.0), "dz": size.get("Lz", 0.0)},
    }]


def box_bounds(box: dict[str, Any]) -> tuple[float, float, float, float, float, float]:
    origin = box.get("origin") or box.get("min") or {"x": 0.0, "y": 0.0, "z": 0.0}
    size = box.get("size", {})
    dx = size.get("dx", size.get("Lx", 0.0))
    dy = size.get("dy", size.get("Ly", 0.0))
    dz = size.get("dz", size.get("Lz", 0.0))
    x0 = float(origin.get("x", 0.0))
    y0 = float(origin.get("y", 0.0))
    z0 = float(origin.get("z", 0.0))
    return x0, y0, z0, x0 + float(dx), y0 + float(dy), z0 + float(dz)


def overall_bbox(boxes: list[dict[str, Any]]) -> dict[str, float]:
    bounds = [box_bounds(box) for box in boxes]
    if not bounds:
        return {"Lx": 0.0, "Ly": 0.0, "Lz": 0.0}
    return {
        "Lx": max(b[3] for b in bounds) - min(b[0] for b in bounds),
        "Ly": max(b[4] for b in bounds) - min(b[1] for b in bounds),
        "Lz": max(b[5] for b in bounds) - min(b[2] for b in bounds),
    }


def compute_scale_hint(scene_path: Path, room_kind: str, setting: str) -> dict[str, Any]:
    scene = load_json(scene_path)
    bbox = overall_bbox(room_blocks(scene))
    axis = "x" if bbox["Lx"] >= bbox["Ly"] else "y"
    span = max(bbox["Lx"], bbox["Ly"])
    height = bbox["Lz"]
    if setting == SCALE_HINTED_DUAL_SETTING:
        kind = "longest_horizontal_span_plus_room_height"
        prompt_text = (
            f"Scale hint: the longest horizontal span of the room is approximately {span:.2f} m, "
            f"and the ceiling height is approximately {height:.2f} m. "
            "Use the span as the primary global metric anchor, and keep the ceiling height close to the hinted value when choosing room dimensions, opening sizes, and obstacle sizes. "
            "Preserve the qualitative layout from the image instead of treating these hints as an exact full bounding box."
        )
    elif setting == SCALE_HINTED_LAYOUT_PROTECTED_SETTING:
        kind = "longest_horizontal_span_layout_protected"
        prompt_text = (
            f"Scale hint: the longest horizontal span of the room is approximately {span:.2f} m. "
            "First infer the room topology, opening-wall placement, and obstacle layout from the image. "
            "Then use this number only as a soft global metric anchor for the overall room span and proportionally scale the rest of the geometry. "
            "Do not move openings to different walls or collapse a clearly composite room just to satisfy the hinted span. "
            "If exact scale agreement conflicts with the image, preserve the image-consistent layout/topology first and treat the hint as approximate."
        )
    elif setting == SCALE_HINTED_VIEW_GUARDED_SETTING:
        kind = "longest_horizontal_span_view_guarded"
        prompt_text = (
            f"Scale hint: the longest horizontal span of the room is approximately {span:.2f} m. "
            "First preserve room topology, opening-wall identity, and the main flow path supported by the image. "
            "Then use this number only as a soft global anchor for the overall room span and proportionally scale the rest of the geometry. "
            "Do not move openings to different walls, invent unsupported hidden depth, or collapse a clearly composite room just to satisfy the hinted span. "
            "If exact scale agreement conflicts with the image evidence, preserve the image-supported layout/topology first and treat the hint as approximate."
        )
    elif setting == SCALE_HINTED_GUARD_WEIGHTED_SETTING:
        kind = "longest_horizontal_span_guard_weighted"
        prompt_text = (
            f"Scale hint: the longest horizontal span of the room is approximately {span:.2f} m. "
            "Keep the image-supported room topology, opening-wall identity, and dominant flow path first. "
            "Use this number only as a soft metric anchor for the dominant horizontal span, not as an exact full bounding box and not as a reason to regularize every view in the same way. "
            "Apply the hint conservatively on views that already expose layout well, and prefer conservative unsupported geometry over aggressive hidden-depth or unseen-height completion. "
            "If the scene is dense or composite, preserve connected-room topology and opening placement before refining obstacle detail. "
            "If exact scale agreement conflicts with the visible evidence, preserve the image-supported layout/topology first and treat the hint as approximate."
        )
    else:
        kind = "longest_horizontal_span"
        prompt_text = (
            f"Scale hint: the longest horizontal span of the room is approximately {span:.2f} m. "
            "Use this as a global metric anchor when choosing room dimensions, opening sizes, and obstacle sizes. "
            "Preserve the qualitative layout from the image instead of treating this hint as an exact full bounding box."
        )
    return {
        "kind": kind,
        "axis": axis,
        "span_m": round(span, 3),
        "room_height_m": round(height, 3),
        "room_kind": room_kind,
        "setting": setting,
        "source_scene": str(scene_path),
        "prompt_text": prompt_text,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Scaffold benchmark evaluation tasks from the frozen reference bundle")
    parser.add_argument("--scene-manifest", type=Path, default=SCENE_MANIFEST)
    parser.add_argument("--reference-status", type=Path, default=DEFAULT_REFERENCE_STATUS)
    parser.add_argument("--renderings-manifest", type=Path, default=RENDERINGS_MANIFEST)
    parser.add_argument("--evaluation-root", type=Path, default=EVALUATIONS)
    parser.add_argument("--setting", choices=[DEFAULT_SETTING, SCALE_HINTED_SETTING, SCALE_HINTED_DUAL_SETTING, SCALE_HINTED_LAYOUT_PROTECTED_SETTING, SCALE_HINTED_VIEW_GUARDED_SETTING, SCALE_HINTED_GUARD_WEIGHTED_SETTING], default=DEFAULT_SETTING)
    args = parser.parse_args()

    scene_rows = load_json(args.scene_manifest)
    ref_payload = load_json(args.reference_status)
    render_payload = load_json(args.renderings_manifest)

    ref_by_case = {row["case_name"]: row for row in ref_payload["results"]}
    render_by_case = {row["case_name"]: row for row in render_payload["cases"]}

    evaluation_root = args.evaluation_root.expanduser().resolve()
    evaluation_root.mkdir(parents=True, exist_ok=True)

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

        case_dir = evaluation_root / case_name
        case_dir.mkdir(parents=True, exist_ok=True)

        scene_src = Path(scene_row["scene_file"])
        if not scene_src.is_absolute():
            scene_src = PROJECT_ROOT / scene_src
        ref_case = PROJECT_ROOT / "cases" / case_name
        ref_results = PROJECT_ROOT / "results" / case_name
        scale_hint = None if args.setting == DEFAULT_SETTING else compute_scale_hint(scene_src, scene_row["room_kind"], args.setting)

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
                "setting": args.setting,
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
                "scale_hint": scale_hint,
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
            "setting": args.setting,
            "scale_hint": scale_hint,
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
        "evaluation_root": str(evaluation_root),
        "setting": args.setting,
        "views": DEFAULT_VIEWS,
        "case_count": len(case_summaries),
        "task_count": len(tasks),
        "status_counts": status_counts,
        "cases": case_summaries,
    }
    manifest = {
        "ok": True,
        "setting": args.setting,
        "cases": case_summaries,
        "tasks": tasks,
    }

    (evaluation_root / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (evaluation_root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
