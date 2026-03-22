#!/usr/bin/env python3
"""Compute the "no scale calibration" ablation baseline.

Reads baseline_evaluation_summary.json from each of the 100 evaluation
directories (20 cases x 5 views) under evaluations_posthoc_scaled_longest_span,
extracts structural_score / room_kind_match / opening_wall_match, and produces
aggregate statistics (overall, by-view, by-category).

Output: benchmark/manifests/noscaling_ablation_results.json
"""

import json
import statistics
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUTPUT_PATH = PROJECT_ROOT / "benchmark" / "manifests" / "noscaling_ablation_results.json"

CASES = [
    f"bench_a{a}_{c:02d}"
    for a in range(1, 5)
    for c in range(1, 6)
]
VIEWS = ["birdseye", "floorplan", "perspective", "section", "wireframe"]


def category_from_case(case_name: str) -> str:
    """bench_a1_01 -> A1"""
    return case_name.split("_")[1].upper()


def aggregate(records: list[dict]) -> dict:
    """Compute aggregate stats from a list of per-case records."""
    n = len(records)
    if n == 0:
        return {"n": 0}
    scores = [r["structural_score"] for r in records]
    rk = [r["room_kind_match"] for r in records]
    ow = [r["opening_wall_match"] for r in records]
    return {
        "n": n,
        "mean_structural_score": round(statistics.mean(scores), 4),
        "median_structural_score": round(statistics.median(scores), 4),
        "std_structural_score": round(statistics.stdev(scores), 4) if n > 1 else 0.0,
        "min_structural_score": round(min(scores), 4),
        "max_structural_score": round(max(scores), 4),
        "room_kind_match_rate": round(sum(rk) / n, 4),
        "opening_wall_match_rate": round(sum(ow) / n, 4),
    }


def main():
    records = []
    missing = []

    for case in CASES:
        cat = category_from_case(case)
        for view in VIEWS:
            path = EVAL_ROOT / case / view / "baseline_evaluation_summary.json"
            if not path.exists():
                missing.append(str(path))
                continue
            with open(path) as f:
                data = json.load(f)

            ps = data.get("prediction_summary", {})
            records.append({
                "case": case,
                "view": view,
                "category": cat,
                "structural_score": ps.get("structural_score", 0.0),
                "room_kind_match": bool(ps.get("room_kind_match", False)),
                "opening_wall_match": bool(ps.get("opening_wall_match", False)),
            })

    # Overall
    overall = aggregate(records)

    # By view
    by_view = {}
    for v in VIEWS:
        by_view[v] = aggregate([r for r in records if r["view"] == v])

    # By category
    by_category = {}
    for cat in sorted({r["category"] for r in records}):
        by_category[cat] = aggregate([r for r in records if r["category"] == cat])

    result = {
        "ok": True,
        "description": "No-scaling ablation baseline: structural scores from VLM raw output without post-hoc scale calibration",
        "evaluation_root": str(EVAL_ROOT),
        "source_file": "baseline_evaluation_summary.json",
        "total_expected": len(CASES) * len(VIEWS),
        "total_found": len(records),
        "missing_count": len(missing),
        "overall": overall,
        "by_view": by_view,
        "by_category": by_category,
        "missing_files": missing,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(result, f, indent=2)

    # Print summary
    print(f"No-scaling ablation baseline computed.")
    print(f"  Records: {len(records)} / {len(CASES) * len(VIEWS)}")
    if missing:
        print(f"  Missing: {len(missing)}")
    print()
    print(f"  Overall mean structural score (no scaling): {overall['mean_structural_score']}")
    print(f"  Overall room_kind_match rate:               {overall['room_kind_match_rate']}")
    print(f"  Overall opening_wall_match rate:             {overall['opening_wall_match_rate']}")
    print()
    print("  By view:")
    for v, s in by_view.items():
        print(f"    {v:12s}  structural={s['mean_structural_score']:.4f}  n={s['n']}")
    print()
    print("  By category:")
    for c, s in by_category.items():
        print(f"    {c:4s}  structural={s['mean_structural_score']:.4f}  n={s['n']}")
    print()
    print(f"  Comparison: scaled = 0.7813, no-scaling = {overall['mean_structural_score']}")
    delta = overall["mean_structural_score"] - 0.7813
    print(f"  Delta (no-scaling - scaled): {delta:+.4f}")
    print()
    print(f"  Output: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
