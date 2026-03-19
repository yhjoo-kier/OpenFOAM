#!/usr/bin/env python3
"""Compare baseline vs post-hoc-scaled prediction against the reference scene."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def room_kind(scene: dict) -> str:
    return "composite" if "blocks" in scene.get("room", {}) else "rectangular"


def room_blocks(scene: dict) -> list[dict[str, Any]]:
    room = scene.get("room", {})
    if "blocks" in room:
        return room["blocks"]
    size = room.get("size", {})
    return [{"origin": {"x": 0.0, "y": 0.0, "z": 0.0}, "size": {"dx": size.get("Lx", 0.0), "dy": size.get("Ly", 0.0), "dz": size.get("Lz", 0.0)}}]


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


def rel_error(pred: float, ref: float) -> float | None:
    if ref == 0:
        return None
    return abs(pred - ref) / abs(ref)


def opening_signature(scene: dict) -> list[tuple[str, str]]:
    return sorted((str(op.get("type")), str(op.get("wall"))) for op in scene.get("openings", []))


def opening_metrics(reference_scene: dict, predicted_scene: dict) -> dict[str, Any]:
    ref = opening_signature(reference_scene)
    pred = opening_signature(predicted_scene)
    ref_set = set(ref)
    pred_set = set(pred)
    intersection = ref_set & pred_set
    precision = len(intersection) / len(pred_set) if pred_set else 0.0
    recall = len(intersection) / len(ref_set) if ref_set else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    wall_ref = [wall for _, wall in ref]
    wall_pred = [wall for _, wall in pred]
    matched = sum(1 for a, b in zip(wall_ref, wall_pred) if a == b)
    ratio = matched / len(wall_ref) if wall_ref else None
    return {"signature_match": ref == pred, "wall_match_ratio": round(ratio, 4) if ratio is not None else None, "type_f1": round(f1, 4)}


def scene_summary(reference_scene: dict, predicted_scene: dict) -> dict[str, Any]:
    ref_bbox = overall_bbox(room_blocks(reference_scene))
    pred_bbox = overall_bbox(room_blocks(predicted_scene))
    return {
        "room_kind": room_kind(predicted_scene),
        "room_kind_match": room_kind(reference_scene) == room_kind(predicted_scene),
        "room_bbox": pred_bbox,
        "room_bbox_relative_error": {axis: round(rel_error(pred_bbox[axis], ref_bbox[axis]), 4) if rel_error(pred_bbox[axis], ref_bbox[axis]) is not None else None for axis in ["Lx", "Ly", "Lz"]},
        "opening_metrics": opening_metrics(reference_scene, predicted_scene),
        "opening_wall_match": opening_signature(reference_scene) == opening_signature(predicted_scene),
        "reference_obstacle_count": len(reference_scene.get("obstacles", []) or []),
        "predicted_obstacle_count": len(predicted_scene.get("obstacles", []) or []),
        "obstacle_count_delta": len(predicted_scene.get("obstacles", []) or []) - len(reference_scene.get("obstacles", []) or []),
    }


def diff(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    out = {
        "room_kind_changed": before["room_kind"] != after["room_kind"],
        "opening_wall_match_changed": before["opening_wall_match"] != after["opening_wall_match"],
        "obstacle_count_delta_changed": before["obstacle_count_delta"] != after["obstacle_count_delta"],
        "room_bbox_relative_error_delta": {},
    }
    for axis in ["Lx", "Ly", "Lz"]:
        b = before["room_bbox_relative_error"].get(axis)
        a = after["room_bbox_relative_error"].get(axis)
        out["room_bbox_relative_error_delta"][axis] = round(a - b, 4) if b is not None and a is not None else None
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare baseline vs post-hoc scaled scene summaries")
    parser.add_argument("--reference-scene", type=Path, required=True)
    parser.add_argument("--baseline-scene", type=Path, required=True)
    parser.add_argument("--posthoc-scene", type=Path, required=True)
    args = parser.parse_args()

    reference_scene = load_json(args.reference_scene.expanduser().resolve())
    baseline_scene = load_json(args.baseline_scene.expanduser().resolve())
    posthoc_scene = load_json(args.posthoc_scene.expanduser().resolve())

    baseline = scene_summary(reference_scene, baseline_scene)
    posthoc = scene_summary(reference_scene, posthoc_scene)
    print(json.dumps({
        "reference_room_kind": room_kind(reference_scene),
        "baseline": baseline,
        "posthoc": posthoc,
        "delta": diff(baseline, posthoc),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
