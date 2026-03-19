#!/usr/bin/env python3
"""Scaffold post-hoc scaled evaluation tasks from an existing baseline evaluation root."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK = PROJECT_ROOT / "benchmark"
SOURCE_EVALUATION_ROOT = BENCHMARK / "evaluations"
TARGET_EVALUATION_ROOT = BENCHMARK / "evaluations_posthoc_scaled_longest_span"
SETTING = "posthoc_uniform_longest_span_v1"
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
    parser = argparse.ArgumentParser(description="Scaffold post-hoc scaled benchmark evaluation tasks")
    parser.add_argument("--source-evaluation-root", type=Path, default=SOURCE_EVALUATION_ROOT)
    parser.add_argument("--target-evaluation-root", type=Path, default=TARGET_EVALUATION_ROOT)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--views", nargs="*", choices=DEFAULT_VIEWS, default=None)
    args = parser.parse_args()

    source_root = args.source_evaluation_root.expanduser().resolve()
    target_root = args.target_evaluation_root.expanduser().resolve()
    target_root.mkdir(parents=True, exist_ok=True)

    cases_payload: list[dict[str, Any]] = []
    tasks_payload: list[dict[str, Any]] = []

    case_dirs = sorted(p for p in source_root.iterdir() if p.is_dir() and p.name.startswith("bench_"))
    for case_dir in case_dirs:
        if args.cases and case_dir.name not in args.cases:
            continue
        view_dirs = sorted(p for p in case_dir.iterdir() if p.is_dir() and p.name in (args.views or DEFAULT_VIEWS))
        if not view_dirs:
            continue

        case_tasks: list[dict[str, Any]] = []
        target_case_dir = target_root / case_dir.name
        target_case_dir.mkdir(parents=True, exist_ok=True)

        for view_dir in view_dirs:
            baseline_task_path = view_dir / "task.json"
            baseline_eval_summary_path = view_dir / "evaluation_summary.json"
            baseline_predicted_scene = view_dir / "predicted_scene.json"
            if not baseline_task_path.exists() or not baseline_eval_summary_path.exists() or not baseline_predicted_scene.exists():
                continue
            baseline_task = load_json(baseline_task_path)
            if baseline_task.get("status") != "success":
                continue

            target_view_dir = target_case_dir / view_dir.name
            target_view_dir.mkdir(parents=True, exist_ok=True)
            for name in ["input.png", "reference_scene.json", "reference_case", "reference_results"]:
                ensure_link(view_dir / name, target_view_dir / name)
            ensure_link(baseline_predicted_scene, target_view_dir / "baseline_predicted_scene.json")
            ensure_link(baseline_eval_summary_path, target_view_dir / "baseline_evaluation_summary.json")
            ensure_link(baseline_task_path, target_view_dir / "baseline_task.json")

            task = {
                "case_name": baseline_task["case_name"],
                "category": baseline_task["category"],
                "variant_index": baseline_task["variant_index"],
                "room_kind": baseline_task["room_kind"],
                "obstacle_count": baseline_task["obstacle_count"],
                "view": baseline_task["view"],
                "setting": SETTING,
                "status": "pending",
                "input_image": str(target_view_dir / "input.png"),
                "reference_scene": str(target_view_dir / "reference_scene.json"),
                "reference_case": str(target_view_dir / "reference_case"),
                "reference_results": str(target_view_dir / "reference_results"),
                "baseline_predicted_scene": str(target_view_dir / "baseline_predicted_scene.json"),
                "baseline_evaluation_summary": str(target_view_dir / "baseline_evaluation_summary.json"),
                "baseline_task_json": str(target_view_dir / "baseline_task.json"),
                "posthoc_scaling": {
                    "kind": "uniform_scaling",
                    "characteristic": "longest_horizontal_span",
                    "source_evaluation_root": str(source_root),
                },
                "expected_outputs": {
                    "scaled_scene_json": str(target_view_dir / "scaled_scene.json"),
                    "predicted_scene_json": str(target_view_dir / "predicted_scene.json"),
                    "predicted_case_dir": str(target_view_dir / "predicted_case"),
                    "predicted_results_dir": str(target_view_dir / "predicted_results"),
                    "evaluation_summary_json": str(target_view_dir / "evaluation_summary.json"),
                },
                "reference_summary": baseline_task.get("reference_summary"),
            }
            task_path = target_view_dir / "task.json"
            task_path.write_text(json.dumps(task, indent=2), encoding="utf-8")
            case_tasks.append(task)
            tasks_payload.append(task)

        if case_tasks:
            case_manifest = {
                "case_name": case_tasks[0]["case_name"],
                "category": case_tasks[0]["category"],
                "room_kind": case_tasks[0]["room_kind"],
                "obstacle_count": case_tasks[0]["obstacle_count"],
                "setting": SETTING,
                "views": [task["view"] for task in case_tasks],
                "task_count": len(case_tasks),
                "case_root": str(target_case_dir),
                "source_case_root": str(case_dir),
            }
            (target_case_dir / "manifest.json").write_text(json.dumps(case_manifest, indent=2), encoding="utf-8")
            cases_payload.append(case_manifest)

    summary = {
        "ok": True,
        "source_evaluation_root": str(source_root),
        "evaluation_root": str(target_root),
        "setting": SETTING,
        "case_count": len(cases_payload),
        "task_count": len(tasks_payload),
        "status_counts": {"pending": len(tasks_payload)},
        "cases": cases_payload,
    }
    manifest = {"ok": True, "setting": SETTING, "cases": cases_payload, "tasks": tasks_payload}
    (target_root / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (target_root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
