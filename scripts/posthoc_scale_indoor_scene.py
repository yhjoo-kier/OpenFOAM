#!/usr/bin/env python3
"""Apply a minimal post-hoc uniform scaling layer to an indoor scene JSON.

This is intentionally downstream-only: it reads an already-generated predicted
scene, matches its room longest horizontal span to the reference scene, and
uniformly rescales room blocks, obstacles, and openings without changing wall
identity, topology, or orientation.
"""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


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


def longest_horizontal_span(scene: dict[str, Any]) -> float:
    bbox = overall_bbox(room_blocks(scene))
    return max(float(bbox.get("Lx", 0.0)), float(bbox.get("Ly", 0.0)))


def scale_mapping(obj: dict[str, Any], keys: tuple[str, ...], factor: float) -> None:
    for key in keys:
        if key in obj and obj[key] is not None:
            obj[key] = round(float(obj[key]) * factor, 6)


def scale_scene(scene: dict[str, Any], factor: float) -> dict[str, Any]:
    scaled = deepcopy(scene)
    room = scaled.setdefault("room", {})

    if "blocks" in room:
        for block in room.get("blocks", []):
            scale_mapping(block.setdefault("origin", {}), ("x", "y", "z"), factor)
            scale_mapping(block.setdefault("size", {}), ("dx", "dy", "dz", "Lx", "Ly", "Lz"), factor)
    elif "size" in room:
        scale_mapping(room.setdefault("size", {}), ("Lx", "Ly", "Lz", "dx", "dy", "dz"), factor)

    for obstacle in scaled.get("obstacles", []) or []:
        scale_mapping(obstacle.setdefault("min", {}), ("x", "y", "z"), factor)
        scale_mapping(obstacle.setdefault("origin", {}), ("x", "y", "z"), factor)
        scale_mapping(obstacle.setdefault("size", {}), ("dx", "dy", "dz", "Lx", "Ly", "Lz"), factor)

    for opening in scaled.get("openings", []) or []:
        scale_mapping(opening.setdefault("center", {}), ("u", "v"), factor)
        scale_mapping(opening.setdefault("size", {}), ("du", "dv"), factor)
        if "min" in opening and isinstance(opening["min"], dict):
            scale_mapping(opening["min"], ("u", "v"), factor)

    return scaled


def main() -> int:
    parser = argparse.ArgumentParser(description="Uniformly scale a predicted scene to match reference longest span")
    parser.add_argument("--reference-scene", type=Path, required=True)
    parser.add_argument("--predicted-scene", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--characteristic", choices=["longest_horizontal_span"], default="longest_horizontal_span")
    args = parser.parse_args()

    reference_scene = load_json(args.reference_scene.expanduser().resolve())
    predicted_scene = load_json(args.predicted_scene.expanduser().resolve())

    ref_span = longest_horizontal_span(reference_scene)
    pred_span = longest_horizontal_span(predicted_scene)
    if pred_span <= 0:
        raise ValueError(f"Predicted characteristic length must be positive, got {pred_span}")
    factor = ref_span / pred_span

    scaled_scene = scale_scene(predicted_scene, factor)
    scaled_scene.setdefault("meta", {})["posthoc_scaling"] = {
        "applied": True,
        "kind": "uniform_scaling",
        "characteristic": args.characteristic,
        "reference_value_m": round(ref_span, 6),
        "predicted_value_before_m": round(pred_span, 6),
        "factor": round(factor, 6),
        "reference_scene": str(args.reference_scene.expanduser().resolve()),
        "source_predicted_scene": str(args.predicted_scene.expanduser().resolve()),
        "applied_at": utc_now(),
    }
    write_json(args.output.expanduser().resolve(), scaled_scene)
    print(json.dumps({
        "ok": True,
        "reference_scene": str(args.reference_scene.expanduser().resolve()),
        "predicted_scene": str(args.predicted_scene.expanduser().resolve()),
        "output": str(args.output.expanduser().resolve()),
        "characteristic": args.characteristic,
        "reference_value_m": round(ref_span, 6),
        "predicted_value_before_m": round(pred_span, 6),
        "factor": round(factor, 6),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
