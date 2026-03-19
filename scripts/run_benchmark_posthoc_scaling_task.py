#!/usr/bin/env python3
"""Run one scaffolded post-hoc scaling benchmark evaluation task.

This reuses an existing baseline predicted scene, applies a minimal uniform
post-hoc scaling layer, then reruns only the downstream meshing/solver/eval path.
No Gemini regeneration is performed.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import run_benchmark_evaluation_task as base_eval

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
DEFAULT_TASK = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span" / "bench_a1_01" / "perspective" / "task.json"
SETTING = "posthoc_uniform_longest_span_v1"


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


def cleanup_outputs(output_paths: dict[str, Path]) -> None:
    for key in ("scaled_scene_json", "predicted_scene_json", "predicted_case_dir", "predicted_results_dir", "evaluation_summary_json"):
        path = output_paths.get(key)
        if not path:
            continue
        if path.exists() or path.is_symlink():
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink()


def infer_evaluation_root(task_path: Path) -> Path:
    return task_path.parent.parent.parent


def compute_cfd_metrics(reference_scene_path: Path, reference_case_path: Path, predicted_scene_path: Path, predicted_case_path: Path, cfd_summary_path: Path) -> dict[str, Any] | None:
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{PROJECT_ROOT}:{PROJECT_ROOT}",
        "-w", str(PROJECT_ROOT),
        "openfoam-pipeline-local:latest",
        "python3", str(SCRIPTS / "compute_benchmark_cfd_metrics.py"),
        "--reference-scene", str(reference_scene_path),
        "--reference-case", str(reference_case_path),
        "--predicted-scene", str(predicted_scene_path),
        "--predicted-case", str(predicted_case_path),
        "--output", str(cfd_summary_path),
    ]
    proc = subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True)
    if proc.returncode != 0 or not cfd_summary_path.exists():
        return {
            "ok": False,
            "command": cmd,
            "returncode": proc.returncode,
            "stdout_tail": "\n".join((proc.stdout or "").splitlines()[-30:]),
            "stderr_tail": "\n".join((proc.stderr or "").splitlines()[-30:]),
        }
    payload = load_json(cfd_summary_path)
    payload.setdefault("runner", {})
    payload["runner"].update({
        "command": cmd,
        "returncode": proc.returncode,
        "stdout_tail": "\n".join((proc.stdout or "").splitlines()[-30:]),
        "stderr_tail": "\n".join((proc.stderr or "").splitlines()[-30:]),
    })
    write_json(cfd_summary_path, payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one post-hoc scaling benchmark evaluation task")
    parser.add_argument("--task", type=Path, default=DEFAULT_TASK)
    parser.add_argument("--mesh-size", type=float, default=0.35)
    parser.add_argument("--end-time", type=int, default=250)
    parser.add_argument("--solver-timeout", type=int, default=600)
    parser.add_argument("--skip-existing-success", action="store_true")
    args = parser.parse_args()

    task_path = args.task.expanduser().resolve()
    if not task_path.exists():
        raise FileNotFoundError(f"Task file does not exist: {task_path}")

    evaluation_root = infer_evaluation_root(task_path)
    task = load_json(task_path)
    output_paths = {k: Path(v) for k, v in task["expected_outputs"].items()}
    evaluation_summary_path = output_paths["evaluation_summary_json"]

    if args.skip_existing_success and task.get("status") == "success" and evaluation_summary_path.exists():
        print(json.dumps({"ok": True, "skipped": True, "reason": "existing_success", "task": str(task_path)}, indent=2))
        return 0

    cleanup_outputs(output_paths)
    cfd_summary_path = task_path.parent / "cfd_metrics.json"
    if cfd_summary_path.exists():
        cfd_summary_path.unlink()

    task["status"] = "running"
    task["setting"] = SETTING
    task["last_started_at"] = utc_now()
    task["run_request"] = {
        "mesh_size": args.mesh_size,
        "end_time": args.end_time,
        "solver_timeout": args.solver_timeout,
        "setting": SETTING,
        "posthoc_scaling": task.get("posthoc_scaling"),
    }
    write_json(task_path, task)

    reference_scene_path = Path(task["reference_scene"])
    reference_case_path = Path(task["reference_case"])
    baseline_predicted_scene = Path(task["baseline_predicted_scene"])
    scaled_scene_path = output_paths["scaled_scene_json"]

    scale_cmd = [
        "python3", str(SCRIPTS / "posthoc_scale_indoor_scene.py"),
        "--reference-scene", str(reference_scene_path),
        "--predicted-scene", str(baseline_predicted_scene),
        "--output", str(scaled_scene_path),
        "--characteristic", str((task.get("posthoc_scaling") or {}).get("characteristic", "longest_horizontal_span")),
    ]
    scale_proc = subprocess.run(scale_cmd, cwd=PROJECT_ROOT, text=True, capture_output=True)
    if scale_proc.returncode != 0 or not scaled_scene_path.exists():
        evaluation_summary = {
            "ok": False,
            "setting": SETTING,
            "task": {
                "case_name": task["case_name"],
                "view": task["view"],
                "task_json": str(task_path),
                "reference_scene": str(reference_scene_path),
                "evaluation_root": str(evaluation_root),
            },
            "run": {
                "name": None,
                "command": scale_cmd,
                "returncode": scale_proc.returncode,
                "started_at": task.get("last_started_at"),
                "finished_at": utc_now(),
            },
            "baseline": {
                "predicted_scene_json": str(baseline_predicted_scene),
                "evaluation_summary_json": str(task.get("baseline_evaluation_summary")),
            },
            "outputs": {
                "scaled_scene_json": str(scaled_scene_path),
                "predicted_scene_json": str(output_paths["predicted_scene_json"]),
                "predicted_case_dir": str(output_paths["predicted_case_dir"]),
                "predicted_results_dir": str(output_paths["predicted_results_dir"]),
                "cfd_metrics_json": None,
            },
            "scale_summary": None,
            "prediction_summary": None,
            "cfd_summary": None,
            "pipeline_summary": None,
            "stdout_tail": "\n".join((scale_proc.stdout or "").splitlines()[-40:]),
            "stderr_tail": "\n".join((scale_proc.stderr or "").splitlines()[-40:]),
        }
        write_json(evaluation_summary_path, evaluation_summary)
        task["status"] = "failed"
        task["last_finished_at"] = evaluation_summary["run"]["finished_at"]
        task["evaluation_summary"] = str(evaluation_summary_path)
        write_json(task_path, task)
        print(json.dumps(evaluation_summary, indent=2))
        return 1

    scale_summary = json.loads(scale_proc.stdout)

    run_name = f"eval_{SETTING}_{task['case_name']}_{task['view']}"
    run_cmd = [
        "python3", str(SCRIPTS / "run_indoor_stabilized.py"),
        "--scenario", str(scaled_scene_path),
        "--name", run_name,
        "--mesh-size", str(args.mesh_size),
        "--end-time", str(args.end_time),
        "--solver-timeout", str(args.solver_timeout),
        "--disable-repair",
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

    baseline_evaluation_summary = load_json(Path(task["baseline_evaluation_summary"])) if task.get("baseline_evaluation_summary") else None
    baseline_prediction_summary = (baseline_evaluation_summary or {}).get("prediction_summary") or {}
    comparison = None
    if success and predicted_scene_out.exists():
        posthoc_prediction_summary = base_eval.summarize_prediction(reference_scene_path, predicted_scene_out)
        comparison = dict(posthoc_prediction_summary)
        comparison["baseline_comparison"] = {
            "baseline_prediction_summary": baseline_prediction_summary,
            "room_bbox_relative_error_delta": {
                axis: round(((posthoc_prediction_summary.get("room_bbox_relative_error") or {}).get(axis) or 0.0) - (((baseline_prediction_summary.get("room_bbox_relative_error") or {}).get(axis) or 0.0)), 4)
                for axis in ["Lx", "Ly", "Lz"]
            },
            "opening_wall_match_changed": posthoc_prediction_summary.get("opening_wall_match") != baseline_prediction_summary.get("opening_wall_match"),
            "room_kind_match_changed": posthoc_prediction_summary.get("room_kind_match") != baseline_prediction_summary.get("room_kind_match"),
            "obstacle_count_delta_changed": posthoc_prediction_summary.get("obstacle_count_delta") != baseline_prediction_summary.get("obstacle_count_delta"),
        }

    cfd_summary = None
    if success and predicted_scene_out.exists() and (predicted_case_out.exists() or predicted_case_out.is_symlink()):
        cfd_summary = compute_cfd_metrics(reference_scene_path, reference_case_path, predicted_scene_out, predicted_case_out, cfd_summary_path)

    evaluation_summary = {
        "ok": success,
        "setting": SETTING,
        "task": {
            "case_name": task["case_name"],
            "view": task["view"],
            "task_json": str(task_path),
            "reference_scene": str(reference_scene_path),
            "evaluation_root": str(evaluation_root),
        },
        "baseline": {
            "predicted_scene_json": str(baseline_predicted_scene),
            "evaluation_summary_json": str(task.get("baseline_evaluation_summary")),
        },
        "run": {
            "name": run_name,
            "command": run_cmd,
            "returncode": proc.returncode,
            "started_at": task.get("last_started_at"),
            "finished_at": utc_now(),
        },
        "outputs": {
            "scaled_scene_json": str(scaled_scene_path),
            "predicted_scene_json": str(predicted_scene_out),
            "predicted_case_dir": str(predicted_case_out),
            "predicted_results_dir": str(predicted_results_out),
            "stabilization_summary": str(stabilization_summary_path) if stabilization_summary_path.exists() else None,
            "cfd_metrics_json": str(cfd_summary_path) if cfd_summary_path.exists() else None,
        },
        "scale_summary": scale_summary,
        "prediction_summary": comparison,
        "cfd_summary": cfd_summary,
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
        "scaled_scene_json": str(scaled_scene_path) if scaled_scene_path.exists() else None,
        "predicted_scene_json": str(predicted_scene_out) if predicted_scene_out.exists() else None,
        "predicted_case_dir": str(predicted_case_out) if predicted_case_out.exists() or predicted_case_out.is_symlink() else None,
        "predicted_results_dir": str(predicted_results_out) if predicted_results_out.exists() or predicted_results_out.is_symlink() else None,
        "cfd_metrics_json": str(cfd_summary_path) if cfd_summary_path.exists() else None,
    }
    write_json(task_path, task)
    print(json.dumps(evaluation_summary, indent=2))
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
