#!/usr/bin/env python3
"""Verify frozen benchmark bundle integrity and write a canonical summary manifest.

Checks the benchmark dataset end-to-end across:
- scene manifest rows
- reference CFD aggregate manifest + linked case/results artifacts
- benchmark 2D render bundle (5-view input images)
- evaluation scaffold completeness

Outputs:
- benchmark/manifests/dataset_integrity_summary.json
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK = PROJECT_ROOT / "benchmark"
MANIFESTS = BENCHMARK / "manifests"
RENDERINGS = BENCHMARK / "renderings"
EVALUATIONS = BENCHMARK / "evaluations"
DEFAULT_VIEWS = ["perspective", "birdseye", "floorplan", "wireframe", "section"]


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_project_path(value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path.resolve()


def path_exists(value: str | None) -> bool:
    path = resolve_project_path(value)
    return bool(path and path.exists())


def verify_case(
    scene_row: dict[str, Any],
    ref_row: dict[str, Any] | None,
    render_row: dict[str, Any] | None,
    eval_row: dict[str, Any] | None,
    views: list[str],
) -> dict[str, Any]:
    case_name = scene_row["case_name"]
    issues: list[str] = []

    scene_path = resolve_project_path(scene_row.get("scene_file"))
    if not scene_path or not scene_path.exists():
        issues.append("missing_scene_file")

    if ref_row is None:
        issues.append("missing_reference_manifest_row")
    else:
        if not ref_row.get("success"):
            issues.append("reference_not_successful")
        for key in ("case_dir", "results_dir", "comparison_1x2", "render3d", "benchmark_render_manifest"):
            if ref_row.get(key) and not path_exists(ref_row.get(key)):
                issues.append(f"missing_reference_artifact:{key}")
        if ref_row.get("attempt_count", 0) > 1:
            issues.append("stress_case_multi_attempt")
        if ref_row.get("successful_mesh_size") not in (None, 0.35):
            issues.append("stress_case_mesh_fallback")
        if ref_row.get("successful_mode") not in (None, "RAS"):
            issues.append("stress_case_mode_fallback")
        if ref_row.get("successful_preset") not in (None, "robust"):
            issues.append("stress_case_preset_fallback")

    section_axis = None
    section_plane = None
    if render_row is None:
        issues.append("missing_render_manifest_row")
    else:
        source_scene = resolve_project_path(render_row.get("source_scene"))
        if not source_scene or not source_scene.exists():
            issues.append("render_manifest_source_scene_missing")
        for view in views:
            if view not in render_row.get("views", {}):
                issues.append(f"missing_render_view:{view}")
            elif not path_exists(render_row["views"].get(view)):
                issues.append(f"missing_render_image:{view}")
        section_meta = render_row.get("section_view") or {}
        section_axis = section_meta.get("axis")
        section_plane = section_meta.get("plane_coordinate")
        if "section" in views:
            if section_axis not in {"x", "y"}:
                issues.append("invalid_section_axis")
            if section_plane is None:
                issues.append("missing_section_plane")

    if eval_row is None:
        issues.append("missing_evaluation_case_manifest")
    else:
        expected_task_count = len(views)
        if eval_row.get("task_count") != expected_task_count:
            issues.append("evaluation_task_count_mismatch")
        present_views = set(eval_row.get("views", []))
        missing_eval_views = [view for view in views if view not in present_views]
        for view in missing_eval_views:
            issues.append(f"missing_evaluation_view:{view}")
        case_root = resolve_project_path(eval_row.get("case_root"))
        if not case_root or not case_root.exists():
            issues.append("evaluation_case_root_missing")
        else:
            for view in views:
                task_path = case_root / view / "task.json"
                if not task_path.exists():
                    issues.append(f"missing_task_json:{view}")

    return {
        "case_name": case_name,
        "category": scene_row.get("category"),
        "room_kind": scene_row.get("room_kind"),
        "obstacle_count": scene_row.get("obstacle_count"),
        "scene_file": str(scene_path) if scene_path else scene_row.get("scene_file"),
        "reference_success": bool(ref_row and ref_row.get("success")),
        "attempt_count": ref_row.get("attempt_count") if ref_row else None,
        "successful_preset": ref_row.get("successful_preset") if ref_row else None,
        "successful_mode": ref_row.get("successful_mode") if ref_row else None,
        "successful_mesh_size": ref_row.get("successful_mesh_size") if ref_row else None,
        "render_views_present": sorted(render_row.get("views", {}).keys()) if render_row else [],
        "evaluation_views_present": sorted(eval_row.get("views", [])) if eval_row else [],
        "section_axis": section_axis,
        "section_plane_coordinate": section_plane,
        "ok": not any(not issue.startswith("stress_case_") for issue in issues),
        "issues": issues,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify the frozen benchmark bundle and write an integrity summary")
    parser.add_argument("--scene-manifest", type=Path, default=MANIFESTS / "scene_manifest.json")
    parser.add_argument("--reference-manifest", type=Path, default=MANIFESTS / "reference_batch_summary.json")
    parser.add_argument("--renderings-manifest", type=Path, default=RENDERINGS / "renderings_manifest.json")
    parser.add_argument("--evaluations-summary", type=Path, default=EVALUATIONS / "summary.json")
    parser.add_argument("--output", type=Path, default=MANIFESTS / "dataset_integrity_summary.json")
    args = parser.parse_args()

    scene_rows = load_json(args.scene_manifest)
    ref_payload = load_json(args.reference_manifest)
    render_payload = load_json(args.renderings_manifest)
    eval_payload = load_json(args.evaluations_summary)

    ref_by_case = {row["case_name"]: row for row in ref_payload.get("results", [])}
    render_by_case = {row["case_name"]: row for row in render_payload.get("cases", [])}
    eval_by_case = {row["case_name"]: row for row in eval_payload.get("cases", [])}

    case_rows = [
        verify_case(scene_row, ref_by_case.get(scene_row["case_name"]), render_by_case.get(scene_row["case_name"]), eval_by_case.get(scene_row["case_name"]), DEFAULT_VIEWS)
        for scene_row in scene_rows
    ]

    hard_issue_counter: Counter[str] = Counter()
    soft_issue_counter: Counter[str] = Counter()
    status_by_category: dict[str, dict[str, int]] = defaultdict(lambda: {"total": 0, "ok": 0})
    section_axis_counts: Counter[str] = Counter()
    stress_cases: list[dict[str, Any]] = []

    for row in case_rows:
        category = row.get("category") or "unknown"
        status_by_category[category]["total"] += 1
        if row.get("ok"):
            status_by_category[category]["ok"] += 1
        if row.get("section_axis"):
            section_axis_counts[row["section_axis"]] += 1
        hard_issues = [issue for issue in row["issues"] if not issue.startswith("stress_case_")]
        soft_issues = [issue for issue in row["issues"] if issue.startswith("stress_case_")]
        hard_issue_counter.update(hard_issues)
        soft_issue_counter.update(soft_issues)
        if soft_issues:
            stress_cases.append(
                {
                    "case_name": row["case_name"],
                    "category": row.get("category"),
                    "attempt_count": row.get("attempt_count"),
                    "successful_preset": row.get("successful_preset"),
                    "successful_mode": row.get("successful_mode"),
                    "successful_mesh_size": row.get("successful_mesh_size"),
                    "signals": soft_issues,
                }
            )

    expected_task_count = len(scene_rows) * len(DEFAULT_VIEWS)
    payload = {
        "ok": not hard_issue_counter,
        "project_root": str(PROJECT_ROOT),
        "benchmark_root": str(BENCHMARK),
        "views": DEFAULT_VIEWS,
        "summary": {
            "scene_count": len(scene_rows),
            "reference_case_count": ref_payload.get("count"),
            "reference_success_count": ref_payload.get("success_count"),
            "render_bundle_count": render_payload.get("count"),
            "evaluation_case_count": eval_payload.get("case_count"),
            "evaluation_task_count": eval_payload.get("task_count"),
            "expected_evaluation_task_count": expected_task_count,
            "hard_issue_count": sum(hard_issue_counter.values()),
            "soft_stress_signal_count": sum(soft_issue_counter.values()),
        },
        "status_by_category": status_by_category,
        "section_axis_counts": dict(section_axis_counts),
        "hard_issue_counts": dict(hard_issue_counter),
        "soft_stress_signal_counts": dict(soft_issue_counter),
        "stress_cases": stress_cases,
        "cases": case_rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
