#!/usr/bin/env python3
"""Run one scaffolded benchmark evaluation task end-to-end.

This consumes a task scaffold under benchmark/evaluations/<case>/<view>/task.json,
reuses the existing indoor pipeline entrypoint, writes predicted outputs back into
that scaffold, and refreshes aggregate benchmark evaluation manifests.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
BENCHMARK = PROJECT_ROOT / "benchmark"
EVALUATIONS = BENCHMARK / "evaluations"
DEFAULT_TASK = EVALUATIONS / "bench_a1_01" / "perspective" / "task.json"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def ensure_link(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    dst.symlink_to(src)


def room_kind(scene: dict) -> str:
    return "composite" if "blocks" in scene.get("room", {}) else "rectangular"


def build_scenario_prompt(view: str) -> str:
    view_hint = {
        "perspective": "The image is a perspective benchmark rendering with visible depth cues.",
        "birdseye": "The image is a bird's-eye benchmark rendering emphasizing overall layout.",
        "floorplan": "The image is a floor-plan style benchmark rendering emphasizing the horizontal layout.",
        "wireframe": "The image is a wireframe benchmark rendering emphasizing structural edges.",
        "section": "The image is a section-view benchmark rendering emphasizing one vertical cut through the space.",
    }.get(view, "The image is a benchmark rendering of an indoor ventilated space.")
    return (
        "Generate a simulation-ready indoor ventilation scene from the provided image. "
        "Use the image as the primary cue and infer a solver-friendly abstraction of room shape, openings, and only the largest flow-relevant obstacles. "
        "Do not assume hidden geometry beyond what is reasonably supported by the image. "
        "Favor simple box obstacles and use two joined room blocks only if the visible layout is clearly non-rectangular. "
        f"{view_hint}"
    )


def summarize_prediction(reference_scene_path: Path, predicted_scene_path: Path) -> dict[str, Any]:
    reference_scene = load_json(reference_scene_path)
    predicted_scene = load_json(predicted_scene_path)
    ref_openings = sorted(op.get("wall") for op in reference_scene.get("openings", []))
    pred_openings = sorted(op.get("wall") for op in predicted_scene.get("openings", []))
    return {
        "reference_room_kind": room_kind(reference_scene),
        "predicted_room_kind": room_kind(predicted_scene),
        "room_kind_match": room_kind(reference_scene) == room_kind(predicted_scene),
        "reference_obstacle_count": len(reference_scene.get("obstacles", [])),
        "predicted_obstacle_count": len(predicted_scene.get("obstacles", [])),
        "obstacle_count_delta": len(predicted_scene.get("obstacles", [])) - len(reference_scene.get("obstacles", [])),
        "reference_opening_walls": ref_openings,
        "predicted_opening_walls": pred_openings,
        "opening_wall_match": ref_openings == pred_openings,
    }


def refresh_evaluation_index() -> None:
    task_paths = sorted(EVALUATIONS.glob("*/ */task.json".replace(" ", "")))
    # Fallback if the above glob behaves unexpectedly across platforms.
    if not task_paths:
        task_paths = sorted(EVALUATIONS.glob("*/*/task.json"))

    tasks: list[dict[str, Any]] = []
    case_map: dict[str, dict[str, Any]] = {}
    status_counts: Counter[str] = Counter()

    for task_path in task_paths:
        task = load_json(task_path)
        tasks.append(task)
        status_counts[task.get("status", "unknown")] += 1
        case_name = task["case_name"]
        case_entry = case_map.setdefault(
            case_name,
            {
                "case_name": case_name,
                "category": task.get("category"),
                "room_kind": task.get("room_kind"),
                "obstacle_count": task.get("obstacle_count"),
                "views": [],
                "reference_scene": task.get("reference_scene"),
                "reference_case": task.get("reference_case"),
                "reference_results": task.get("reference_results"),
                "task_count": 0,
                "status_counts": Counter(),
                "case_root": str(task_path.parent.parent),
            },
        )
        case_entry["views"].append(task["view"])
        case_entry["task_count"] += 1
        case_entry["status_counts"][task.get("status", "unknown")] += 1

    case_summaries: list[dict[str, Any]] = []
    for case_name in sorted(case_map):
        entry = case_map[case_name]
        entry["views"] = sorted(entry["views"])
        entry["status_counts"] = dict(entry["status_counts"])
        case_summaries.append(entry)
        write_json(Path(entry["case_root"]) / "manifest.json", entry)

    summary = {
        "ok": True,
        "benchmark_root": str(BENCHMARK),
        "evaluation_root": str(EVALUATIONS),
        "views": ["perspective", "birdseye", "floorplan", "wireframe", "section"],
        "case_count": len(case_summaries),
        "task_count": len(tasks),
        "status_counts": dict(status_counts),
        "cases": case_summaries,
    }
    manifest = {"ok": True, "cases": case_summaries, "tasks": tasks}
    write_json(EVALUATIONS / "summary.json", summary)
    write_json(EVALUATIONS / "manifest.json", manifest)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one benchmark evaluation task end-to-end")
    parser.add_argument("--task", type=Path, default=DEFAULT_TASK, help="Path to task.json")
    parser.add_argument("--backend", choices=["cli", "api"], default="api")
    parser.add_argument("--model", default="gemini-3.1-pro-preview")
    parser.add_argument("--mesh-size", type=float, default=0.35)
    parser.add_argument("--end-time", type=int, default=250)
    parser.add_argument("--solver-timeout", type=int, default=600)
    parser.add_argument("--skip-existing-success", action="store_true")
    args = parser.parse_args()

    task_path = args.task.expanduser().resolve()
    if not task_path.exists():
        raise FileNotFoundError(f"Task file does not exist: {task_path}")

    task = load_json(task_path)
    task_dir = task_path.parent
    output_paths = {k: Path(v) for k, v in task["expected_outputs"].items()}
    evaluation_summary_path = output_paths["evaluation_summary_json"]

    if args.skip_existing_success and task.get("status") == "success" and evaluation_summary_path.exists():
        refresh_evaluation_index()
        print(json.dumps({"ok": True, "skipped": True, "reason": "existing_success", "task": str(task_path)}, indent=2))
        return 0

    task["status"] = "running"
    task["last_started_at"] = utc_now()
    task["run_request"] = {
        "backend": args.backend,
        "model": args.model,
        "mesh_size": args.mesh_size,
        "end_time": args.end_time,
        "solver_timeout": args.solver_timeout,
    }
    write_json(task_path, task)
    refresh_evaluation_index()

    input_image = Path(task["input_image"])
    reference_scene_path = Path(task["reference_scene"])
    run_name = f"eval_{task['case_name']}_{task['view']}"
    run_cmd = [
        "python3", str(SCRIPTS / "run_indoor_stabilized.py"),
        "--backend", args.backend,
        "--model", args.model,
        "--image", str(input_image),
        "--scenario", build_scenario_prompt(task["view"]),
        "--name", run_name,
        "--mesh-size", str(args.mesh_size),
        "--end-time", str(args.end_time),
        "--solver-timeout", str(args.solver_timeout),
    ]

    proc = subprocess.run(run_cmd, cwd=PROJECT_ROOT, text=True, capture_output=True)
    results_dir = PROJECT_ROOT / "results" / run_name
    case_dir = PROJECT_ROOT / "cases" / run_name
    stabilization_summary_path = results_dir / "stabilization_summary.json"
    stabilization_summary = load_json(stabilization_summary_path) if stabilization_summary_path.exists() else None

    predicted_scene_src = None
    if stabilization_summary and stabilization_summary.get("used_scene_json"):
        predicted_scene_src = Path(stabilization_summary["used_scene_json"])
        if not predicted_scene_src.is_absolute():
            predicted_scene_src = PROJECT_ROOT / predicted_scene_src

    predicted_scene_out = output_paths["predicted_scene_json"]
    predicted_case_out = output_paths["predicted_case_dir"]
    predicted_results_out = output_paths["predicted_results_dir"]

    success = bool(stabilization_summary and stabilization_summary.get("success"))
    if success and predicted_scene_src and predicted_scene_src.exists():
        predicted_scene_out.write_text(predicted_scene_src.read_text(encoding="utf-8"), encoding="utf-8")
        if case_dir.exists():
            ensure_link(case_dir, predicted_case_out)
        if results_dir.exists():
            ensure_link(results_dir, predicted_results_out)

    comparison = summarize_prediction(reference_scene_path, predicted_scene_out) if success and predicted_scene_out.exists() else None
    evaluation_summary = {
        "ok": success,
        "task": {
            "case_name": task["case_name"],
            "view": task["view"],
            "task_json": str(task_path),
            "input_image": str(input_image),
            "reference_scene": str(reference_scene_path),
        },
        "run": {
            "name": run_name,
            "command": run_cmd,
            "returncode": proc.returncode,
            "started_at": task.get("last_started_at"),
            "finished_at": utc_now(),
        },
        "outputs": {
            "predicted_scene_json": str(predicted_scene_out),
            "predicted_case_dir": str(predicted_case_out),
            "predicted_results_dir": str(predicted_results_out),
            "stabilization_summary": str(stabilization_summary_path) if stabilization_summary_path.exists() else None,
        },
        "reference_summary": task.get("reference_summary"),
        "prediction_summary": comparison,
        "pipeline_summary": stabilization_summary,
        "stdout_tail": "\n".join((proc.stdout or "").splitlines()[-40:]),
        "stderr_tail": "\n".join((proc.stderr or "").splitlines()[-40:]),
    }
    write_json(evaluation_summary_path, evaluation_summary)

    task["status"] = "success" if success else "failed"
    task["last_finished_at"] = evaluation_summary["run"]["finished_at"]
    task["last_run_name"] = run_name
    task["evaluation_summary"] = str(evaluation_summary_path)
    task["actual_outputs"] = {
        "predicted_scene_json": str(predicted_scene_out) if predicted_scene_out.exists() else None,
        "predicted_case_dir": str(predicted_case_out) if predicted_case_out.exists() or predicted_case_out.is_symlink() else None,
        "predicted_results_dir": str(predicted_results_out) if predicted_results_out.exists() or predicted_results_out.is_symlink() else None,
    }
    write_json(task_path, task)
    refresh_evaluation_index()

    print(json.dumps(evaluation_summary, indent=2))
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
