#!/usr/bin/env python3
"""Compute a naive geometric baseline for the image-to-CFD benchmark.

The naive baseline predicts the same default room for every case:
  - Room: single rectangular block, 5.0 x 4.0 x 3.0 m
  - Obstacles: none
  - Openings: one inlet on west wall, one outlet on east wall
    (each centered at u=2.0, with du=0.5, dv=0.5)

It evaluates each prediction against the reference scene using the same
structural_score formula as the real evaluation pipeline.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Reuse evaluation helpers from the benchmark evaluation script
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_benchmark_evaluation_task import (
    greedy_box_match_summary,
    load_json,
    opening_metrics,
    room_blocks,
    room_kind,
    round_or_none,
    write_json,
)

EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUTPUT_PATH = PROJECT_ROOT / "benchmark" / "manifests" / "naive_baseline_results.json"

VIEWS = ["birdseye", "floorplan", "perspective", "section", "wireframe"]


def naive_predicted_scene() -> dict[str, Any]:
    """Return the fixed naive-baseline predicted scene."""
    return {
        "schema_version": "indoor_cfd_scene_v1",
        "units": "m",
        "coordinate_system": {
            "origin": "room_min_corner",
            "x": "west_to_east",
            "y": "south_to_north",
            "z": "floor_to_ceiling",
        },
        "room": {
            "size": {
                "Lx": 5.0,
                "Ly": 4.0,
                "Lz": 3.0,
            }
        },
        "obstacles": [],
        "openings": [
            {
                "id": "inlet",
                "type": "inlet",
                "wall": "west",
                "center": {"u": 2.0, "v": 1.5},
                "size": {"du": 0.5, "dv": 0.5},
            },
            {
                "id": "outlet",
                "type": "outlet",
                "wall": "east",
                "center": {"u": 2.0, "v": 1.5},
                "size": {"du": 0.5, "dv": 0.5},
            },
        ],
        "meta": {
            "description": "Naive geometric baseline: fixed 5x4x3 room, no obstacles, west inlet / east outlet",
        },
    }


def evaluate_one(reference_scene: dict, predicted_scene: dict) -> dict[str, Any]:
    """Compute structural score for a single (reference, predicted) pair."""
    ref_room = room_blocks(reference_scene)
    pred_room = room_blocks(predicted_scene)
    ref_obstacles = reference_scene.get("obstacles", [])
    pred_obstacles = predicted_scene.get("obstacles", [])

    room_block_match = greedy_box_match_summary(
        ref_room, pred_room, label="room_blocks", iou_match_threshold=0.2
    )
    obstacle_match = greedy_box_match_summary(
        ref_obstacles, pred_obstacles, label="obstacles", iou_match_threshold=0.1
    )
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
        "room_block_match_f1": room_block_match["f1"],
        "obstacle_match_f1": obstacle_match["f1"],
        "openings_type_f1": openings["type_f1"],
        "openings_wall_match_ratio": openings["wall_match_ratio"],
        "structural_score": round_or_none(structural_score),
    }


def main() -> None:
    predicted = naive_predicted_scene()

    case_dirs = sorted(
        d for d in EVAL_ROOT.iterdir()
        if d.is_dir() and d.name.startswith("bench_")
    )

    per_case_view: list[dict[str, Any]] = []
    scores_all: list[float] = []
    scores_by_view: defaultdict[str, list[float]] = defaultdict(list)
    scores_by_category: defaultdict[str, list[float]] = defaultdict(list)
    room_kind_matches = 0
    opening_wall_matches = 0
    total_count = 0

    for case_dir in case_dirs:
        case_name = case_dir.name
        manifest_path = case_dir / "manifest.json"
        manifest = load_json(manifest_path) if manifest_path.exists() else {}
        category = manifest.get("category", case_name.split("_")[1].upper())

        for view in VIEWS:
            ref_path = case_dir / view / "reference_scene.json"
            if not ref_path.exists():
                continue

            reference_scene = load_json(ref_path)
            result = evaluate_one(reference_scene, predicted)

            total_count += 1
            score = result["structural_score"]
            if score is not None:
                scores_all.append(score)
                scores_by_view[view].append(score)
                scores_by_category[category].append(score)

            if result["room_kind_match"]:
                room_kind_matches += 1
            if result["openings_wall_match_ratio"] > 0:
                opening_wall_matches += 1

            per_case_view.append({
                "case_name": case_name,
                "category": category,
                "view": view,
                **result,
            })

    mean_all = sum(scores_all) / len(scores_all) if scores_all else 0.0
    by_view = {
        v: round(sum(s) / len(s), 4) if s else 0.0
        for v, s in sorted(scores_by_view.items())
    }
    by_category = {
        c: round(sum(s) / len(s), 4) if s else 0.0
        for c, s in sorted(scores_by_category.items())
    }

    payload = {
        "baseline": "naive_geometric",
        "description": "Fixed 5x4x3 room, no obstacles, west inlet / east outlet",
        "timestamp": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "evaluation_root": str(EVAL_ROOT),
        "total_cases": total_count,
        "overall_mean_structural_score": round(mean_all, 4),
        "room_kind_match_rate": round(room_kind_matches / total_count, 4) if total_count else 0.0,
        "opening_wall_match_rate": round(opening_wall_matches / total_count, 4) if total_count else 0.0,
        "by_view": by_view,
        "by_category": by_category,
        "per_case_view": per_case_view,
    }

    write_json(OUTPUT_PATH, payload)

    # Print summary
    print(f"Naive baseline evaluated on {total_count} case-view pairs")
    print(f"Overall mean structural score: {mean_all:.4f}")
    print()
    print("By view:")
    for v, s in by_view.items():
        print(f"  {v:12s}  {s:.4f}")
    print()
    print("By category:")
    for c, s in by_category.items():
        print(f"  {c:4s}  {s:.4f}")
    print()
    print(f"Room-kind match rate:    {room_kind_matches}/{total_count} = {room_kind_matches/total_count:.4f}")
    print(f"Opening wall match rate: {opening_wall_matches}/{total_count} = {opening_wall_matches/total_count:.4f}")
    print(f"\nResults written to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
