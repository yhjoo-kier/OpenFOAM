#!/usr/bin/env python3
"""Compute how the structural score changes when IoU matching threshold is varied.

Sweeps a single IoU threshold (used for both room_block_match and obstacle_match)
across [0.05, 0.1, 0.15, 0.2, 0.3, 0.5] and reports mean structural score at each.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_benchmark_evaluation_task import (
    greedy_box_match_summary,
    load_json,
    opening_metrics,
    room_blocks,
)

EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUTPUT_PATH = PROJECT_ROOT / "benchmark" / "manifests" / "iou_sensitivity.json"

IOU_THRESHOLDS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5]


def structural_score_at_threshold(
    ref_scene: dict, pred_scene: dict, iou_threshold: float
) -> float:
    ref_room = room_blocks(ref_scene)
    pred_room = room_blocks(pred_scene)
    ref_obs = ref_scene.get("obstacles", [])
    pred_obs = pred_scene.get("obstacles", [])

    room_match = greedy_box_match_summary(
        ref_room, pred_room, label="room_blocks", iou_match_threshold=iou_threshold
    )
    obs_match = greedy_box_match_summary(
        ref_obs, pred_obs, label="obstacles", iou_match_threshold=iou_threshold
    )
    openings = opening_metrics(ref_scene, pred_scene)

    components = [
        room_match["f1"],
        obs_match["f1"],
        openings["type_f1"],
        openings["wall_match_ratio"],
    ]
    return sum(components) / len(components) if components else 0.0


def main() -> None:
    # Collect all valid (case, view) pairs
    cases: list[tuple[str, str, dict, dict]] = []
    skipped = 0

    for case_dir in sorted(EVAL_ROOT.iterdir()):
        if not case_dir.is_dir():
            continue
        for view_dir in sorted(case_dir.iterdir()):
            if not view_dir.is_dir():
                continue
            ref_path = view_dir / "reference_scene.json"
            pred_path = view_dir / "predicted_scene.json"
            eval_summary_path = view_dir / "evaluation_summary.json"

            if not pred_path.exists():
                skipped += 1
                continue

            if eval_summary_path.exists():
                summary = load_json(eval_summary_path)
                if summary.get("ok") is False:
                    skipped += 1
                    continue

            ref_scene = load_json(ref_path)
            pred_scene = load_json(pred_path)
            cases.append((case_dir.name, view_dir.name, ref_scene, pred_scene))

    print(f"Loaded {len(cases)} valid evaluation pairs (skipped {skipped})")

    # Sweep thresholds
    results: dict[str, dict] = {}
    for thr in IOU_THRESHOLDS:
        scores = [
            structural_score_at_threshold(ref, pred, thr)
            for _, _, ref, pred in cases
        ]
        mean_score = sum(scores) / len(scores) if scores else 0.0
        results[str(thr)] = {
            "iou_threshold": thr,
            "mean_structural_score": round(mean_score, 4),
            "n_cases": len(scores),
        }

    # Print table
    print()
    print(f"{'IoU threshold':>15s}  {'Mean structural score':>22s}  {'N':>4s}")
    print("-" * 46)
    for thr in IOU_THRESHOLDS:
        r = results[str(thr)]
        print(f"{thr:>15.2f}  {r['mean_structural_score']:>22.4f}  {r['n_cases']:>4d}")

    # Save
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nSaved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
