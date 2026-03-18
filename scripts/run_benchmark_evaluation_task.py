#!/usr/bin/env python3
"""Run one scaffolded benchmark evaluation task end-to-end.

This consumes a task scaffold under benchmark/evaluations/<case>/<view>/task.json,
reuses the existing indoor pipeline entrypoint, writes predicted outputs back into
that scaffold, and refreshes aggregate benchmark evaluation manifests.
"""

from __future__ import annotations

import argparse
import json
import math
import os
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


def room_blocks(scene: dict) -> list[dict[str, Any]]:
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


def union_volume(boxes: list[dict[str, Any]]) -> float:
    if not boxes:
        return 0.0
    bounds = [box_bounds(box) for box in boxes]
    xs = sorted({v for b in bounds for v in (b[0], b[3])})
    ys = sorted({v for b in bounds for v in (b[1], b[4])})
    zs = sorted({v for b in bounds for v in (b[2], b[5])})
    total = 0.0
    for ix in range(len(xs) - 1):
        for iy in range(len(ys) - 1):
            for iz in range(len(zs) - 1):
                x0, x1 = xs[ix], xs[ix + 1]
                y0, y1 = ys[iy], ys[iy + 1]
                z0, z1 = zs[iz], zs[iz + 1]
                cx = 0.5 * (x0 + x1)
                cy = 0.5 * (y0 + y1)
                cz = 0.5 * (z0 + z1)
                if any(bx0 <= cx <= bx1 and by0 <= cy <= by1 and bz0 <= cz <= bz1 for bx0, by0, bz0, bx1, by1, bz1 in bounds):
                    total += (x1 - x0) * (y1 - y0) * (z1 - z0)
    return total


def overall_bbox(boxes: list[dict[str, Any]]) -> dict[str, float]:
    bounds = [box_bounds(box) for box in boxes]
    if not bounds:
        return {"Lx": 0.0, "Ly": 0.0, "Lz": 0.0}
    return {
        "Lx": max(b[3] for b in bounds) - min(b[0] for b in bounds),
        "Ly": max(b[4] for b in bounds) - min(b[1] for b in bounds),
        "Lz": max(b[5] for b in bounds) - min(b[2] for b in bounds),
    }


def round_or_none(value: float | None, digits: int = 4) -> float | None:
    return round(value, digits) if value is not None else None


def box_iou(a: dict[str, Any], b: dict[str, Any]) -> float:
    ax0, ay0, az0, ax1, ay1, az1 = box_bounds(a)
    bx0, by0, bz0, bx1, by1, bz1 = box_bounds(b)
    ix = max(0.0, min(ax1, bx1) - max(ax0, bx0))
    iy = max(0.0, min(ay1, by1) - max(ay0, by0))
    iz = max(0.0, min(az1, bz1) - max(az0, bz0))
    inter = ix * iy * iz
    if inter <= 0.0:
        return 0.0
    av = max(0.0, (ax1 - ax0) * (ay1 - ay0) * (az1 - az0))
    bv = max(0.0, (bx1 - bx0) * (by1 - by0) * (bz1 - bz0))
    union = av + bv - inter
    return inter / union if union > 0 else 0.0


def box_center(box: dict[str, Any]) -> tuple[float, float, float]:
    x0, y0, z0, x1, y1, z1 = box_bounds(box)
    return (0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1))


def box_size(box: dict[str, Any]) -> tuple[float, float, float]:
    x0, y0, z0, x1, y1, z1 = box_bounds(box)
    return (max(0.0, x1 - x0), max(0.0, y1 - y0), max(0.0, z1 - z0))


def l2_distance(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def greedy_box_match_summary(
    reference: list[dict[str, Any]],
    predicted: list[dict[str, Any]],
    *,
    label: str,
    iou_match_threshold: float = 0.1,
) -> dict[str, Any]:
    pairs: list[dict[str, Any]] = []
    remaining_pred = list(range(len(predicted)))
    for ref_idx, ref_box in enumerate(reference):
        best_idx = None
        best_iou = -1.0
        for pred_idx in remaining_pred:
            iou = box_iou(ref_box, predicted[pred_idx])
            if iou > best_iou:
                best_iou = iou
                best_idx = pred_idx
        if best_idx is not None:
            remaining_pred.remove(best_idx)
            pred_box = predicted[best_idx]
            center_error = l2_distance(box_center(ref_box), box_center(pred_box))
            size_error = l2_distance(box_size(ref_box), box_size(pred_box))
            pairs.append(
                {
                    "reference_index": ref_idx,
                    "predicted_index": best_idx,
                    "iou": round(best_iou, 4),
                    "center_error_l2": round(center_error, 4),
                    "size_error_l2": round(size_error, 4),
                    "matched_above_threshold": best_iou >= iou_match_threshold,
                }
            )
        else:
            pairs.append(
                {
                    "reference_index": ref_idx,
                    "predicted_index": None,
                    "iou": 0.0,
                    "center_error_l2": None,
                    "size_error_l2": None,
                    "matched_above_threshold": False,
                }
            )

    matched_pairs = [pair for pair in pairs if pair["predicted_index"] is not None]
    threshold_matches = [pair for pair in matched_pairs if pair["matched_above_threshold"]]
    mean_iou = sum(pair["iou"] for pair in matched_pairs) / len(matched_pairs) if matched_pairs else 0.0
    mean_center_error = (
        sum(pair["center_error_l2"] for pair in threshold_matches if pair["center_error_l2"] is not None) / len(threshold_matches)
        if threshold_matches
        else None
    )
    mean_size_error = (
        sum(pair["size_error_l2"] for pair in threshold_matches if pair["size_error_l2"] is not None) / len(threshold_matches)
        if threshold_matches
        else None
    )
    matched_count = len(threshold_matches)
    precision = matched_count / len(predicted) if predicted else 1.0
    recall = matched_count / len(reference) if reference else 1.0
    f1 = (2.0 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
    return {
        "label": label,
        "iou_match_threshold": iou_match_threshold,
        "matched_pairs": pairs,
        "mean_iou": round(mean_iou, 4),
        "mean_center_error_l2": round_or_none(mean_center_error),
        "mean_size_error_l2": round_or_none(mean_size_error),
        "matched_above_threshold_count": matched_count,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "unmatched_reference_count": max(0, len(reference) - matched_count),
        "unmatched_predicted_count": max(0, len(predicted) - matched_count),
    }


def opening_signature(scene: dict) -> dict[str, dict[str, Any]]:
    signature: dict[str, dict[str, Any]] = {}
    for idx, op in enumerate(scene.get("openings", [])):
        key = str(op.get("type") or f"opening_{idx}")
        center = op.get("center", {})
        size = op.get("size", {})
        signature[key] = {
            "wall": op.get("wall"),
            "center_u": float(center.get("u", 0.0)),
            "center_v": float(center.get("v", 0.0)),
            "size_du": float(size.get("du", 0.0)),
            "size_dv": float(size.get("dv", 0.0)),
        }
    return signature


def opening_metrics(reference_scene: dict, predicted_scene: dict) -> dict[str, Any]:
    ref = opening_signature(reference_scene)
    pred = opening_signature(predicted_scene)
    keys = sorted(set(ref) | set(pred))
    per_opening: dict[str, Any] = {}
    wall_match_count = 0
    comparable_count = 0
    center_errors = []
    size_errors = []
    type_precision = comparable_count / len(pred) if pred else 1.0
    type_recall = comparable_count / len(ref) if ref else 1.0
    for key in keys:
        r = ref.get(key)
        p = pred.get(key)
        if r and p:
            comparable_count += 1
            wall_match = r["wall"] == p["wall"]
            wall_match_count += int(wall_match)
            center_error = math.hypot(r["center_u"] - p["center_u"], r["center_v"] - p["center_v"])
            size_error = math.hypot(r["size_du"] - p["size_du"], r["size_dv"] - p["size_dv"])
            center_errors.append(center_error)
            size_errors.append(size_error)
            per_opening[key] = {
                "present_in_reference": True,
                "present_in_prediction": True,
                "wall_match": wall_match,
                "center_error_l2": round(center_error, 4),
                "size_error_l2": round(size_error, 4),
            }
        else:
            per_opening[key] = {
                "present_in_reference": bool(r),
                "present_in_prediction": bool(p),
                "wall_match": False,
                "center_error_l2": None,
                "size_error_l2": None,
            }
    type_precision = comparable_count / len(pred) if pred else 1.0
    type_recall = comparable_count / len(ref) if ref else 1.0
    type_f1 = (2.0 * type_precision * type_recall / (type_precision + type_recall)) if (type_precision + type_recall) > 0 else 0.0
    return {
        "reference_count": len(ref),
        "predicted_count": len(pred),
        "matched_type_count": comparable_count,
        "type_precision": round(type_precision, 4),
        "type_recall": round(type_recall, 4),
        "type_f1": round(type_f1, 4),
        "wall_match_ratio": round(wall_match_count / comparable_count, 4) if comparable_count else 0.0,
        "mean_center_error_l2": round_or_none(sum(center_errors) / len(center_errors) if center_errors else None),
        "mean_size_error_l2": round_or_none(sum(size_errors) / len(size_errors) if size_errors else None),
        "per_opening": per_opening,
    }


def detect_backend_blocker(backend: str) -> dict[str, str] | None:
    if backend == "api" and not os.environ.get("GEMINI_API_KEY"):
        return {
            "reason": "missing_gemini_api_key",
            "message": "GEMINI_API_KEY is not set for Gemini API backend",
        }
    if backend == "cli" and shutil.which("gemini") is None:
        return {
            "reason": "missing_gemini_cli",
            "message": "gemini CLI is not available in PATH for CLI backend",
        }
    return None


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

    ref_room_blocks = room_blocks(reference_scene)
    pred_room_blocks = room_blocks(predicted_scene)
    ref_obstacles = reference_scene.get("obstacles", [])
    pred_obstacles = predicted_scene.get("obstacles", [])
    ref_openings = sorted(op.get("wall") for op in reference_scene.get("openings", []))
    pred_openings = sorted(op.get("wall") for op in predicted_scene.get("openings", []))

    ref_room_bbox = overall_bbox(ref_room_blocks)
    pred_room_bbox = overall_bbox(pred_room_blocks)
    ref_room_volume = union_volume(ref_room_blocks)
    pred_room_volume = union_volume(pred_room_blocks)
    room_volume_rel_error = abs(pred_room_volume - ref_room_volume) / ref_room_volume if ref_room_volume > 0 else None
    room_bbox_abs_error = {axis: abs(pred_room_bbox[axis] - ref_room_bbox[axis]) for axis in ("Lx", "Ly", "Lz")}
    room_bbox_rel_error = {
        axis: (room_bbox_abs_error[axis] / ref_room_bbox[axis]) if ref_room_bbox[axis] > 0 else None
        for axis in ("Lx", "Ly", "Lz")
    }
    room_block_match = greedy_box_match_summary(ref_room_blocks, pred_room_blocks, label="room_blocks", iou_match_threshold=0.2)
    obstacle_match = greedy_box_match_summary(ref_obstacles, pred_obstacles, label="obstacles", iou_match_threshold=0.1)
    openings = opening_metrics(reference_scene, predicted_scene)

    components = [
        room_block_match["f1"],
        obstacle_match["f1"],
        openings["type_f1"],
        openings["wall_match_ratio"],
    ]
    structural_score = sum(components) / len(components) if components else None

    return {
        "reference_room_kind": room_kind(reference_scene),
        "predicted_room_kind": room_kind(predicted_scene),
        "room_kind_match": room_kind(reference_scene) == room_kind(predicted_scene),
        "reference_room_block_count": len(ref_room_blocks),
        "predicted_room_block_count": len(pred_room_blocks),
        "reference_room_bbox": ref_room_bbox,
        "predicted_room_bbox": pred_room_bbox,
        "room_bbox_absolute_error": {k: round(v, 4) for k, v in room_bbox_abs_error.items()},
        "room_bbox_relative_error": {k: round_or_none(v) for k, v in room_bbox_rel_error.items()},
        "reference_room_volume": round(ref_room_volume, 4),
        "predicted_room_volume": round(pred_room_volume, 4),
        "room_volume_relative_error": round_or_none(room_volume_rel_error),
        "room_block_match": room_block_match,
        "reference_obstacle_count": len(ref_obstacles),
        "predicted_obstacle_count": len(pred_obstacles),
        "obstacle_count_delta": len(pred_obstacles) - len(ref_obstacles),
        "obstacle_match": obstacle_match,
        "reference_opening_walls": ref_openings,
        "predicted_opening_walls": pred_openings,
        "opening_wall_match": ref_openings == pred_openings,
        "opening_metrics": openings,
        "structural_score": round_or_none(structural_score),
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


def cleanup_outputs(output_paths: dict[str, Path]) -> None:
    for path in output_paths.values():
        if path.exists() or path.is_symlink():
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink()


def compute_cfd_metrics(
    reference_scene_path: Path,
    reference_case_path: Path,
    predicted_scene_path: Path,
    predicted_case_path: Path,
    cfd_summary_path: Path,
) -> dict[str, Any] | None:
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
    payload["runner"].update(
        {
            "command": cmd,
            "returncode": proc.returncode,
            "stdout_tail": "\n".join((proc.stdout or "").splitlines()[-30:]),
            "stderr_tail": "\n".join((proc.stderr or "").splitlines()[-30:]),
        }
    )
    write_json(cfd_summary_path, payload)
    return payload



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
    output_paths = {k: Path(v) for k, v in task["expected_outputs"].items()}
    evaluation_summary_path = output_paths["evaluation_summary_json"]

    if args.skip_existing_success and task.get("status") == "success" and evaluation_summary_path.exists():
        refresh_evaluation_index()
        print(json.dumps({"ok": True, "skipped": True, "reason": "existing_success", "task": str(task_path)}, indent=2))
        return 0

    blocker = detect_backend_blocker(args.backend)
    if blocker:
        task["status"] = "blocked"
        task["last_finished_at"] = utc_now()
        task["run_request"] = {
            "backend": args.backend,
            "model": args.model,
            "mesh_size": args.mesh_size,
            "end_time": args.end_time,
            "solver_timeout": args.solver_timeout,
        }
        task["blocked_reason"] = blocker["reason"]
        task["evaluation_summary"] = str(evaluation_summary_path)
        task["actual_outputs"] = {
            "predicted_scene_json": None,
            "predicted_case_dir": None,
            "predicted_results_dir": None,
        }
        write_json(task_path, task)
        evaluation_summary = {
            "ok": False,
            "blocked": True,
            "blocker": blocker,
            "task": {
                "case_name": task["case_name"],
                "view": task["view"],
                "task_json": str(task_path),
                "input_image": str(task["input_image"]),
                "reference_scene": str(task["reference_scene"]),
            },
            "run": {
                "name": None,
                "command": None,
                "returncode": None,
                "started_at": None,
                "finished_at": task["last_finished_at"],
            },
            "outputs": {
                "predicted_scene_json": str(output_paths["predicted_scene_json"]),
                "predicted_case_dir": str(output_paths["predicted_case_dir"]),
                "predicted_results_dir": str(output_paths["predicted_results_dir"]),
                "stabilization_summary": None,
            },
            "reference_summary": task.get("reference_summary"),
            "prediction_summary": None,
            "pipeline_summary": None,
            "stdout_tail": "",
            "stderr_tail": blocker["message"],
        }
        write_json(evaluation_summary_path, evaluation_summary)
        refresh_evaluation_index()
        print(json.dumps(evaluation_summary, indent=2))
        return 2

    cleanup_outputs(output_paths)
    cfd_summary_path = task_path.parent / "cfd_metrics.json"
    if cfd_summary_path.exists():
        cfd_summary_path.unlink()

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
    cfd_summary_path = task_path.parent / "cfd_metrics.json"

    success = bool(stabilization_summary and stabilization_summary.get("success"))
    if success and predicted_scene_src and predicted_scene_src.exists():
        predicted_scene_out.write_text(predicted_scene_src.read_text(encoding="utf-8"), encoding="utf-8")
        if case_dir.exists():
            ensure_link(case_dir, predicted_case_out)
        if results_dir.exists():
            ensure_link(results_dir, predicted_results_out)

    comparison = summarize_prediction(reference_scene_path, predicted_scene_out) if success and predicted_scene_out.exists() else None
    cfd_summary = None
    if success and predicted_scene_out.exists() and (predicted_case_out.exists() or predicted_case_out.is_symlink()):
        cfd_summary = compute_cfd_metrics(
            reference_scene_path,
            Path(task["reference_case"]),
            predicted_scene_out,
            predicted_case_out,
            cfd_summary_path,
        )
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
            "cfd_metrics_json": str(cfd_summary_path) if cfd_summary_path.exists() else None,
        },
        "reference_summary": task.get("reference_summary"),
        "prediction_summary": comparison,
        "cfd_summary": cfd_summary,
        "pipeline_summary": stabilization_summary,
        "stdout_tail": "\n".join((proc.stdout or "").splitlines()[-40:]),
        "stderr_tail": "\n".join((proc.stderr or "").splitlines()[-40:]),
    }
    write_json(evaluation_summary_path, evaluation_summary)

    task["status"] = "success" if success else "failed"
    task.pop("blocked_reason", None)
    task["last_finished_at"] = evaluation_summary["run"]["finished_at"]
    task["last_run_name"] = run_name
    task["evaluation_summary"] = str(evaluation_summary_path)
    task["actual_outputs"] = {
        "predicted_scene_json": str(predicted_scene_out) if predicted_scene_out.exists() else None,
        "predicted_case_dir": str(predicted_case_out) if predicted_case_out.exists() or predicted_case_out.is_symlink() else None,
        "predicted_results_dir": str(predicted_results_out) if predicted_results_out.exists() or predicted_results_out.is_symlink() else None,
        "cfd_metrics_json": str(cfd_summary_path) if cfd_summary_path.exists() else None,
    }
    write_json(task_path, task)
    refresh_evaluation_index()

    print(json.dumps(evaluation_summary, indent=2))
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
