#!/usr/bin/env python3
"""Build a compact dataset card for the frozen OpenFOAM benchmark bundle.

This script consolidates the benchmark's existing manifests into two publication-friendly
artifacts:
- benchmark/manifests/dataset_card.json
- docs/<date>_benchmark_dataset_card.md

The goal is to keep the frozen bundle easy to inspect/release without hand-copying numbers
from multiple manifests.
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
DOCS = PROJECT_ROOT / "docs"

DEFAULT_SCENE_MANIFEST = MANIFESTS / "scene_manifest.json"
DEFAULT_REFERENCE_SUMMARY = MANIFESTS / "reference_batch_summary.json"
DEFAULT_RENDERINGS_MANIFEST = BENCHMARK / "renderings" / "renderings_manifest.json"
DEFAULT_EVALUATION_SUMMARY = BENCHMARK / "evaluations" / "summary.json"
DEFAULT_INTEGRITY_SUMMARY = MANIFESTS / "dataset_integrity_summary.json"
DEFAULT_JSON_OUTPUT = MANIFESTS / "dataset_card.json"
DEFAULT_MD_OUTPUT = DOCS / "26-03-18_benchmark_dataset_card.md"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def format_counter(counter: Counter[str]) -> dict[str, int]:
    return {key: counter[key] for key in sorted(counter)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Build benchmark dataset card artifacts")
    parser.add_argument("--scene-manifest", type=Path, default=DEFAULT_SCENE_MANIFEST)
    parser.add_argument("--reference-summary", type=Path, default=DEFAULT_REFERENCE_SUMMARY)
    parser.add_argument("--renderings-manifest", type=Path, default=DEFAULT_RENDERINGS_MANIFEST)
    parser.add_argument("--evaluation-summary", type=Path, default=DEFAULT_EVALUATION_SUMMARY)
    parser.add_argument("--integrity-summary", type=Path, default=DEFAULT_INTEGRITY_SUMMARY)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON_OUTPUT)
    parser.add_argument("--markdown-output", type=Path, default=DEFAULT_MD_OUTPUT)
    args = parser.parse_args()

    scene_manifest = load_json(args.scene_manifest)
    reference_summary = load_json(args.reference_summary)
    renderings_manifest = load_json(args.renderings_manifest)
    evaluation_summary = load_json(args.evaluation_summary)
    integrity_summary = load_json(args.integrity_summary)

    scene_rows = list(scene_manifest)
    reference_rows = list(reference_summary.get("results", []))
    rendering_rows = list(renderings_manifest.get("cases", []))
    evaluation_rows = list(evaluation_summary.get("cases", []))
    stress_rows = list(integrity_summary.get("stress_cases", []))

    by_case_scene = {row["case_name"]: row for row in scene_rows}
    by_case_ref = {row["case_name"]: row for row in reference_rows}
    by_case_render = {row["case_name"]: row for row in rendering_rows}
    by_case_eval = {row["case_name"]: row for row in evaluation_rows}

    category_counts: Counter[str] = Counter()
    room_kind_counts: Counter[str] = Counter()
    obstacle_hist: Counter[str] = Counter()
    section_axis_counts: Counter[str] = Counter()
    successful_preset_counts: Counter[str] = Counter()
    successful_mode_counts: Counter[str] = Counter()
    successful_mesh_size_counts: Counter[str] = Counter()
    evaluation_status_counts: Counter[str] = Counter(evaluation_summary.get("status_counts", {}))
    per_category_view_coverage: dict[str, Counter[str]] = defaultdict(Counter)

    missing_reference_cases: list[str] = []
    missing_render_cases: list[str] = []
    missing_evaluation_cases: list[str] = []
    case_cards: list[dict[str, Any]] = []

    for case_name in sorted(by_case_scene):
        scene = by_case_scene[case_name]
        ref = by_case_ref.get(case_name)
        render = by_case_render.get(case_name)
        eval_case = by_case_eval.get(case_name)

        category = str(scene.get("category"))
        room_kind = str(scene.get("room_kind"))
        obstacle_count = int(scene.get("obstacle_count", 0))
        category_counts[category] += 1
        room_kind_counts[room_kind] += 1
        obstacle_hist[str(obstacle_count)] += 1

        if ref is None:
            missing_reference_cases.append(case_name)
        else:
            if ref.get("successful_preset"):
                successful_preset_counts[str(ref["successful_preset"])] += 1
            if ref.get("successful_mode"):
                successful_mode_counts[str(ref["successful_mode"])] += 1
            if ref.get("successful_mesh_size") is not None:
                successful_mesh_size_counts[str(ref["successful_mesh_size"])] += 1

        render_views = sorted((render or {}).get("views", {}).keys())
        if render is None:
            missing_render_cases.append(case_name)
        else:
            section_view = render.get("section_view") or {}
            axis = section_view.get("axis")
            if axis:
                section_axis_counts[str(axis)] += 1
            for view in render_views:
                per_category_view_coverage[category][view] += 1

        if eval_case is None:
            missing_evaluation_cases.append(case_name)

        case_cards.append(
            {
                "case_name": case_name,
                "category": category,
                "room_kind": room_kind,
                "obstacle_count": obstacle_count,
                "seed": scene.get("seed"),
                "scene_file": scene.get("scene_file"),
                "reference_success": None if ref is None else bool(ref.get("success")),
                "attempt_count": None if ref is None else ref.get("attempt_count"),
                "successful_preset": None if ref is None else ref.get("successful_preset"),
                "successful_mode": None if ref is None else ref.get("successful_mode"),
                "successful_mesh_size": None if ref is None else ref.get("successful_mesh_size"),
                "render_views": render_views,
                "section_axis": None if render is None else (render.get("section_view") or {}).get("axis"),
                "evaluation_task_count": None if eval_case is None else eval_case.get("task_count"),
                "evaluation_case_root": None if eval_case is None else eval_case.get("case_root"),
            }
        )

    stress_case_names = [row["case_name"] for row in stress_rows]
    payload = {
        "ok": True,
        "project_root": str(PROJECT_ROOT),
        "benchmark_root": str(BENCHMARK),
        "summary": {
            "scene_count": len(scene_rows),
            "reference_case_count": len(reference_rows),
            "render_bundle_count": len(rendering_rows),
            "evaluation_case_count": len(evaluation_rows),
            "evaluation_task_count": int(evaluation_summary.get("task_count", 0)),
            "reference_success_count": int(reference_summary.get("success_count", 0)),
            "stress_case_count": len(stress_rows),
            "hard_issue_count": int(integrity_summary.get("summary", {}).get("hard_issue_count", 0)),
        },
        "coverage": {
            "category_counts": format_counter(category_counts),
            "room_kind_counts": format_counter(room_kind_counts),
            "obstacle_histogram": format_counter(obstacle_hist),
            "section_axis_counts": format_counter(section_axis_counts),
            "views": list(renderings_manifest.get("views", [])),
            "per_category_view_coverage": {
                category: format_counter(counter)
                for category, counter in sorted(per_category_view_coverage.items())
            },
        },
        "stabilization_profile": {
            "successful_preset_counts": format_counter(successful_preset_counts),
            "successful_mode_counts": format_counter(successful_mode_counts),
            "successful_mesh_size_counts": format_counter(successful_mesh_size_counts),
            "stress_case_names": stress_case_names,
            "stress_cases": stress_rows,
        },
        "evaluation_profile": {
            "status_counts": format_counter(evaluation_status_counts),
            "blocked_reason_note": "Blocked tasks in the current cron shell are backend-environment blocks, not benchmark reference failures.",
        },
        "integrity": {
            "summary": integrity_summary.get("summary", {}),
            "soft_stress_signal_counts": integrity_summary.get("soft_stress_signal_counts", {}),
            "missing_reference_cases": missing_reference_cases,
            "missing_render_cases": missing_render_cases,
            "missing_evaluation_cases": missing_evaluation_cases,
        },
        "paths": {
            "scene_manifest": str(args.scene_manifest),
            "reference_summary": str(args.reference_summary),
            "renderings_manifest": str(args.renderings_manifest),
            "evaluation_summary": str(args.evaluation_summary),
            "integrity_summary": str(args.integrity_summary),
        },
        "cases": case_cards,
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines: list[str] = []
    lines.append("# Benchmark Dataset Card")
    lines.append("")
    lines.append("> Date: 2026-03-18")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Frozen benchmark bundle: **{payload['summary']['scene_count']} scenes**")
    lines.append(f"- Reference CFD success: **{payload['summary']['reference_success_count']}/{payload['summary']['reference_case_count']}**")
    lines.append(f"- Render bundles: **{payload['summary']['render_bundle_count']}**")
    lines.append(f"- Evaluation scaffold: **{payload['summary']['evaluation_task_count']} tasks** across **{payload['summary']['evaluation_case_count']} cases**")
    lines.append(f"- Stress subset: **{payload['summary']['stress_case_count']} cases**")
    lines.append("")
    lines.append("## Bundle composition")
    lines.append("")
    lines.append("### Categories")
    lines.append("")
    for category, count in sorted(payload["coverage"]["category_counts"].items()):
        lines.append(f"- {category}: {count}")
    lines.append("")
    lines.append("### Room kinds")
    lines.append("")
    for room_kind, count in sorted(payload["coverage"]["room_kind_counts"].items()):
        lines.append(f"- {room_kind}: {count}")
    lines.append("")
    lines.append("### Obstacle histogram")
    lines.append("")
    for obstacle_count, count in sorted(payload["coverage"]["obstacle_histogram"].items(), key=lambda item: int(item[0])):
        lines.append(f"- {obstacle_count} obstacles: {count}")
    lines.append("")
    lines.append("### View coverage")
    lines.append("")
    lines.append(f"- Views: {', '.join(payload['coverage']['views'])}")
    lines.append(f"- Section axis coverage: {payload['coverage']['section_axis_counts']}")
    lines.append("")
    lines.append("Per-category view counts:")
    lines.append("")
    for category, counts in sorted(payload["coverage"]["per_category_view_coverage"].items()):
        lines.append(f"- {category}: {counts}")
    lines.append("")
    lines.append("## Stabilization profile")
    lines.append("")
    lines.append(f"- Successful preset counts: {payload['stabilization_profile']['successful_preset_counts']}")
    lines.append(f"- Successful mode counts: {payload['stabilization_profile']['successful_mode_counts']}")
    lines.append(f"- Successful mesh size counts: {payload['stabilization_profile']['successful_mesh_size_counts']}")
    lines.append(f"- Stress cases: {', '.join(stress_case_names)}")
    lines.append("")
    lines.append("## Evaluation profile")
    lines.append("")
    lines.append(f"- Status counts: {payload['evaluation_profile']['status_counts']}")
    lines.append(f"- Note: {payload['evaluation_profile']['blocked_reason_note']}")
    lines.append("")
    lines.append("## Integrity / release notes")
    lines.append("")
    lines.append(f"- Integrity summary: {payload['integrity']['summary']}")
    lines.append(f"- Soft stress signal counts: {payload['integrity']['soft_stress_signal_counts']}")
    lines.append(f"- Missing reference cases: {payload['integrity']['missing_reference_cases']}")
    lines.append(f"- Missing render cases: {payload['integrity']['missing_render_cases']}")
    lines.append(f"- Missing evaluation cases: {payload['integrity']['missing_evaluation_cases']}")
    lines.append("")
    lines.append("## Artifact paths")
    lines.append("")
    for key, value in payload["paths"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- This dataset card is generated from the benchmark manifests rather than maintained manually.")
    lines.append("- It is intended to be a stable paper/release snapshot for the current frozen benchmark bundle.")
    lines.append("- Re-run this script after future tranche expansions or partial-rerun manifest refreshes.")
    lines.append("")

    args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown_output.write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({
        "ok": True,
        "json_output": str(args.json_output),
        "markdown_output": str(args.markdown_output),
        "scene_count": payload["summary"]["scene_count"],
        "stress_case_count": payload["summary"]["stress_case_count"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
