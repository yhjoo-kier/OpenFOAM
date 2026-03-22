#!/usr/bin/env python3
"""VLM repeatability test: call Gemini multiple times on the same images and measure structural score variation."""

from __future__ import annotations

import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from generate_indoor_scene_with_gemini import (
    build_prompt,
    ensure_gemini_available,
    generate_with_fallback,
    parse_scene,
    apply_uniform_scale,
    room_horizontal_bbox,
)
from run_benchmark_evaluation_task import (
    summarize_prediction,
    load_json,
    overall_bbox,
    room_blocks,
    write_json,
    build_scenario_prompt,
)

EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
MODEL = "gemini-3.1-pro-preview"
BACKEND = "api"
N_REPEATS = 3

CASES = [
    ("bench_a1_01", "floorplan"),
    ("bench_a3_03", "floorplan"),
    ("bench_a4_03", "floorplan"),
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def compute_scale_factor(reference_scene: dict, predicted_scene: dict) -> float:
    """Compute the post-hoc uniform scale factor to match the reference longest horizontal span."""
    ref_bbox = overall_bbox(room_blocks(reference_scene))
    target_span = max(ref_bbox["Lx"], ref_bbox["Ly"])
    lx, ly = room_horizontal_bbox(predicted_scene)
    predicted_span = max(lx, ly)
    if predicted_span > 0 and target_span > 0:
        return target_span / predicted_span
    return 1.0


def run_single_generation(
    case_name: str,
    view: str,
    image_path: Path,
    reference_scene_path: Path,
    run_index: int,
) -> dict:
    """Run one VLM generation, apply post-hoc scaling, compute structural score."""
    task_json = load_json(EVAL_ROOT / case_name / view / "task.json")
    room_kind = task_json.get("room_kind")
    obstacle_count = task_json.get("obstacle_count")

    scenario = build_scenario_prompt(
        view=view,
        setting="no_scale_hint_baseline",
        room_kind=room_kind,
        obstacle_count=obstacle_count,
    )
    prompt = build_prompt(scenario, scale_hint=None)

    print(f"  [{case_name}/{view}] Run {run_index+1}/{N_REPEATS}: calling Gemini...", flush=True)
    t0 = time.time()
    raw_output, used_model, history = generate_with_fallback(
        prompt,
        MODEL,
        backend=BACKEND,
        image_paths=[image_path],
        allow_fallback=True,
    )
    elapsed = time.time() - t0
    print(f"    Model: {used_model}, elapsed: {elapsed:.1f}s", flush=True)

    scene = parse_scene(raw_output)

    # Apply post-hoc uniform scaling
    reference_scene = load_json(reference_scene_path)
    factor = compute_scale_factor(reference_scene, scene)
    scaled_scene = apply_uniform_scale(scene, factor)

    # Save to temp file for summarize_prediction
    out_dir = EVAL_ROOT / case_name / view
    temp_path = out_dir / f"repeatability_run_{run_index}.json"
    temp_path.write_text(json.dumps(scaled_scene, indent=2), encoding="utf-8")

    # Compute structural score
    summary = summarize_prediction(reference_scene_path, temp_path)

    return {
        "run_index": run_index,
        "used_model": used_model,
        "elapsed_seconds": round(elapsed, 2),
        "scale_factor": round(factor, 6),
        "structural_score": summary["structural_score"],
        "room_kind_match": summary["room_kind_match"],
        "opening_wall_match": summary["opening_wall_match"],
        "obstacle_f1": summary["obstacle_match"]["f1"],
        "room_block_f1": summary["room_block_match"]["f1"],
        "opening_type_f1": summary["opening_metrics"]["type_f1"],
        "opening_wall_match_ratio": summary["opening_metrics"]["wall_match_ratio"],
        "predicted_room_bbox": summary["predicted_room_bbox"],
        "predicted_obstacle_count": summary["predicted_obstacle_count"],
        "temp_scene_path": str(temp_path),
    }


def main():
    ensure_gemini_available(BACKEND)

    results = {
        "experiment": "vlm_repeatability_test",
        "model": MODEL,
        "backend": BACKEND,
        "temperature": 0.2,
        "n_repeats": N_REPEATS,
        "started_at": utc_now(),
        "cases": {},
    }

    for case_name, view in CASES:
        case_dir = EVAL_ROOT / case_name / view
        image_path = case_dir / "input.png"
        reference_scene_path = case_dir / "reference_scene.json"

        if not image_path.exists():
            print(f"SKIP {case_name}/{view}: input.png not found", flush=True)
            continue
        if not reference_scene_path.exists():
            print(f"SKIP {case_name}/{view}: reference_scene.json not found", flush=True)
            continue

        print(f"\n{'='*60}", flush=True)
        print(f"Case: {case_name}/{view}", flush=True)
        print(f"{'='*60}", flush=True)

        runs = []
        for i in range(N_REPEATS):
            try:
                run_result = run_single_generation(
                    case_name, view, image_path, reference_scene_path, i
                )
                runs.append(run_result)
                print(f"    structural_score = {run_result['structural_score']}", flush=True)
            except Exception as e:
                print(f"    ERROR on run {i}: {e}", flush=True, file=sys.stderr)
                runs.append({"run_index": i, "error": str(e)})

            # Brief pause between calls to avoid rate limiting
            if i < N_REPEATS - 1:
                time.sleep(2)

        scores = [r["structural_score"] for r in runs if "structural_score" in r and r["structural_score"] is not None]
        if scores:
            mean_score = sum(scores) / len(scores)
            std_score = math.sqrt(sum((s - mean_score) ** 2 for s in scores) / len(scores))
            min_score = min(scores)
            max_score = max(scores)
        else:
            mean_score = std_score = min_score = max_score = None

        case_key = f"{case_name}/{view}"
        results["cases"][case_key] = {
            "case_name": case_name,
            "view": view,
            "runs": runs,
            "scores": scores,
            "mean_structural_score": round(mean_score, 4) if mean_score is not None else None,
            "std_structural_score": round(std_score, 4) if std_score is not None else None,
            "min_structural_score": round(min_score, 4) if min_score is not None else None,
            "max_structural_score": round(max_score, 4) if max_score is not None else None,
            "range": round(max_score - min_score, 4) if (min_score is not None and max_score is not None) else None,
            "is_deterministic": std_score == 0.0 if std_score is not None else None,
        }

        print(f"\n  Summary: mean={mean_score:.4f} +/- {std_score:.4f}  range=[{min_score:.4f}, {max_score:.4f}]", flush=True)

    results["finished_at"] = utc_now()

    # Overall summary
    all_stds = [c["std_structural_score"] for c in results["cases"].values() if c["std_structural_score"] is not None]
    results["overall"] = {
        "mean_std_across_cases": round(sum(all_stds) / len(all_stds), 4) if all_stds else None,
        "max_std_across_cases": round(max(all_stds), 4) if all_stds else None,
        "any_deterministic": any(c.get("is_deterministic") for c in results["cases"].values()),
        "all_deterministic": all(c.get("is_deterministic") for c in results["cases"].values()),
        "note": "Temperature is set to 0.2 in the Gemini API generationConfig. Non-zero temperature means outputs are stochastic.",
    }

    out_path = PROJECT_ROOT / "benchmark" / "manifests" / "vlm_repeatability_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nResults saved to {out_path}", flush=True)

    # Print final summary table
    print(f"\n{'='*70}")
    print(f"{'Case':<30} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8} {'Range':>8}")
    print(f"{'-'*70}")
    for case_key, case_data in results["cases"].items():
        m = case_data["mean_structural_score"]
        s = case_data["std_structural_score"]
        mn = case_data["min_structural_score"]
        mx = case_data["max_structural_score"]
        rng = case_data["range"]
        print(f"{case_key:<30} {m:>8.4f} {s:>8.4f} {mn:>8.4f} {mx:>8.4f} {rng:>8.4f}")
    print(f"{'='*70}")
    print(f"VLM deterministic: {results['overall']['all_deterministic']}")
    print(f"Temperature: {results['temperature']}")


if __name__ == "__main__":
    main()
