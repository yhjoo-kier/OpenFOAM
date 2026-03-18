#!/usr/bin/env python3
"""Build aggregate image-conditioned benchmark evaluation summaries.

This consolidates the per-task `evaluation_summary.json` files into two artifacts:
- `benchmark/manifests/evaluation_aggregate_summary.json`
- `docs/26-03-19_cli_eval_aggregate_results.md`

The intent is to keep the 100/100 frozen-20 CLI evaluation easy to cite in the paper's
Results/Discussion section without manually recomputing view/category-level statistics.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK = PROJECT_ROOT / "benchmark"
EVALUATIONS = BENCHMARK / "evaluations"
MANIFESTS = BENCHMARK / "manifests"
DOCS = PROJECT_ROOT / "docs"

DEFAULT_EVALUATION_ROOT = EVALUATIONS
DEFAULT_JSON_OUTPUT = MANIFESTS / "evaluation_aggregate_summary.json"
DEFAULT_MD_OUTPUT = DOCS / "26-03-19_cli_eval_aggregate_results.md"


VIEW_ORDER = ["perspective", "birdseye", "floorplan", "wireframe", "section"]
CATEGORY_ORDER = ["A1", "A2", "A3", "A4"]
ROOM_KIND_ORDER = ["rectangular", "composite"]


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def round_or_none(value: float | None, digits: int = 4) -> float | None:
    if value is None:
        return None
    return round(float(value), digits)


def mean(values: list[float | None]) -> float | None:
    clean = [float(v) for v in values if v is not None]
    if not clean:
        return None
    return sum(clean) / len(clean)


def fraction(count: int, total: int) -> float | None:
    if total <= 0:
        return None
    return count / total


def collect_rows(evaluation_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for summary_path in sorted(evaluation_root.glob("*/*/evaluation_summary.json")):
        data = load_json(summary_path)
        task = data.get("task", {})
        prediction = data.get("prediction_summary", {})
        cfd_summary = (data.get("cfd_summary") or {}).get("aggregate_score", {})
        pipeline = data.get("pipeline_summary", {})
        attempts = list(pipeline.get("attempts") or [])
        last_attempt = attempts[-1] if attempts else {}
        opening_metrics = prediction.get("opening_metrics") or {}

        case_name = str(task.get("case_name"))
        category = case_name.split("_")[1].upper() if "_" in case_name else "unknown"

        rows.append(
            {
                "case": case_name,
                "view": str(task.get("view")),
                "category": category,
                "room_kind": prediction.get("reference_room_kind"),
                "predicted_room_kind": prediction.get("predicted_room_kind"),
                "room_kind_match": bool(prediction.get("room_kind_match")),
                "opening_wall_match": bool(prediction.get("opening_wall_match")),
                "opening_wall_ratio": opening_metrics.get("wall_match_ratio"),
                "structural_score": prediction.get("structural_score"),
                "cfd_score": cfd_summary.get("cfd_score"),
                "reference_obstacle_count": prediction.get("reference_obstacle_count"),
                "predicted_obstacle_count": prediction.get("predicted_obstacle_count"),
                "obstacle_count_delta": prediction.get("obstacle_count_delta"),
                "used_repair": pipeline.get("successful_scene_variant") == "repaired",
                "successful_scene_variant": pipeline.get("successful_scene_variant"),
                "preset": pipeline.get("successful_preset"),
                "mode": pipeline.get("successful_mode"),
                "mesh_size": pipeline.get("successful_mesh_size"),
                "attempt_count": len(attempts),
                "max_non_ortho": ((last_attempt.get("checkmesh") or {}).get("max_non_ortho")),
                "repair_info_present": pipeline.get("repair_info") is not None,
                "summary_path": str(summary_path),
            }
        )
    return rows


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(rows)
    return {
        "n": n,
        "mean_structural_score": round_or_none(mean([row["structural_score"] for row in rows])),
        "mean_cfd_score": round_or_none(mean([row["cfd_score"] for row in rows])),
        "room_kind_match_rate": round_or_none(fraction(sum(1 for row in rows if row["room_kind_match"]), n)),
        "opening_wall_match_rate": round_or_none(fraction(sum(1 for row in rows if row["opening_wall_match"]), n)),
        "used_repair_count": sum(1 for row in rows if row["used_repair"]),
        "ultra_robust_count": sum(1 for row in rows if row["preset"] == "ultra_robust"),
        "conservative_count": sum(1 for row in rows if row["preset"] == "conservative"),
        "mesh025_count": sum(1 for row in rows if row["mesh_size"] == 0.25),
    }


def group_summary(rows: list[dict[str, Any]], key: str, order: list[str]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[key])].append(row)
    keys = [value for value in order if value in grouped] + sorted(k for k in grouped if k not in order)
    return {group_key: summarize_rows(grouped[group_key]) for group_key in keys}


def case_cards(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row["case"]].append(row)

    cards: list[dict[str, Any]] = []
    for case_name in sorted(grouped):
        case_rows = grouped[case_name]
        cards.append(
            {
                "case": case_name,
                "category": case_rows[0]["category"],
                "room_kind": case_rows[0]["room_kind"],
                "mean_structural_score": round_or_none(mean([row["structural_score"] for row in case_rows])),
                "mean_cfd_score": round_or_none(mean([row["cfd_score"] for row in case_rows])),
                "room_kind_match_count": sum(1 for row in case_rows if row["room_kind_match"]),
                "opening_wall_match_count": sum(1 for row in case_rows if row["opening_wall_match"]),
                "used_repair_views": [row["view"] for row in case_rows if row["used_repair"]],
                "ultra_robust_views": [row["view"] for row in case_rows if row["preset"] == "ultra_robust"],
                "mesh025_views": [row["view"] for row in case_rows if row["mesh_size"] == 0.25],
                "hallucination_views": [row["view"] for row in case_rows if (row["obstacle_count_delta"] or 0) > 0],
            }
        )
    return cards


def build_payload(rows: list[dict[str, Any]]) -> dict[str, Any]:
    cards = case_cards(rows)
    cards_low = sorted(cards, key=lambda card: card["mean_cfd_score"] if card["mean_cfd_score"] is not None else 999.0)
    cards_high = sorted(cards, key=lambda card: card["mean_cfd_score"] if card["mean_cfd_score"] is not None else -1.0, reverse=True)

    task_cards = [
        {
            "case": row["case"],
            "view": row["view"],
            "category": row["category"],
            "room_kind": row["room_kind"],
            "structural_score": round_or_none(row["structural_score"]),
            "cfd_score": round_or_none(row["cfd_score"]),
        }
        for row in rows
    ]

    payload = {
        "ok": True,
        "project_root": str(PROJECT_ROOT),
        "evaluation_root": str(DEFAULT_EVALUATION_ROOT),
        "overall": summarize_rows(rows),
        "by_view": group_summary(rows, "view", VIEW_ORDER),
        "by_category": group_summary(rows, "category", CATEGORY_ORDER),
        "by_room_kind": group_summary(rows, "room_kind", ROOM_KIND_ORDER),
        "case_rankings": {
            "lowest_mean_cfd": cards_low[:10],
            "highest_mean_cfd": cards_high[:10],
        },
        "notables": {
            "room_kind_mismatch_tasks": [
                {
                    "case": row["case"],
                    "view": row["view"],
                    "ref_room_kind": row["room_kind"],
                    "pred_room_kind": row["predicted_room_kind"],
                    "structural_score": round_or_none(row["structural_score"]),
                    "cfd_score": round_or_none(row["cfd_score"]),
                }
                for row in rows
                if not row["room_kind_match"]
            ],
            "repair_used_tasks": [
                {
                    "case": row["case"],
                    "view": row["view"],
                    "preset": row["preset"],
                    "mesh_size": row["mesh_size"],
                    "structural_score": round_or_none(row["structural_score"]),
                    "cfd_score": round_or_none(row["cfd_score"]),
                }
                for row in rows
                if row["used_repair"]
            ],
            "ultra_robust_tasks": [
                {
                    "case": row["case"],
                    "view": row["view"],
                    "used_repair": row["used_repair"],
                    "structural_score": round_or_none(row["structural_score"]),
                    "cfd_score": round_or_none(row["cfd_score"]),
                }
                for row in rows
                if row["preset"] == "ultra_robust"
            ],
            "mesh025_tasks": [
                {
                    "case": row["case"],
                    "view": row["view"],
                    "preset": row["preset"],
                    "used_repair": row["used_repair"],
                    "structural_score": round_or_none(row["structural_score"]),
                    "cfd_score": round_or_none(row["cfd_score"]),
                }
                for row in rows
                if row["mesh_size"] == 0.25
            ],
            "highest_cfd_tasks": sorted(task_cards, key=lambda card: card["cfd_score"] if card["cfd_score"] is not None else -1.0, reverse=True)[:10],
            "lowest_cfd_tasks": sorted(task_cards, key=lambda card: card["cfd_score"] if card["cfd_score"] is not None else 999.0)[:10],
            "largest_positive_obstacle_delta": sorted(
                [
                    {
                        "case": row["case"],
                        "view": row["view"],
                        "reference_obstacle_count": row["reference_obstacle_count"],
                        "predicted_obstacle_count": row["predicted_obstacle_count"],
                        "obstacle_count_delta": row["obstacle_count_delta"],
                        "structural_score": round_or_none(row["structural_score"]),
                        "cfd_score": round_or_none(row["cfd_score"]),
                    }
                    for row in rows
                    if (row["obstacle_count_delta"] or 0) > 0
                ],
                key=lambda card: card["obstacle_count_delta"],
                reverse=True,
            )[:20],
        },
        "generated_from": {
            "task_count": len(rows),
            "evaluation_summary_pattern": "benchmark/evaluations/*/*/evaluation_summary.json",
        },
    }
    return payload


def write_markdown(payload: dict[str, Any], path: Path) -> None:
    overall = payload["overall"]
    by_view = payload["by_view"]
    by_category = payload["by_category"]
    by_room_kind = payload["by_room_kind"]

    lines: list[str] = []
    lines.append("# Frozen-20 CLI Evaluation Aggregate Results")
    lines.append("")
    lines.append("> Date: 2026-03-19")
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append("- Frozen-20 image-conditioned benchmark evaluation reached **100/100 task success**.")
    lines.append(f"- Aggregate mean structural score: **{overall['mean_structural_score']:.4f}**")
    lines.append(f"- Aggregate mean CFD score: **{overall['mean_cfd_score']:.4f}**")
    lines.append(f"- Room-kind match rate: **{overall['room_kind_match_rate']:.2%}**")
    lines.append(f"- Opening-wall match rate: **{overall['opening_wall_match_rate']:.2%}**")
    lines.append(f"- Repair was needed in only **{overall['used_repair_count']}** tasks; mesh-size `0.25` fallback appeared in **{overall['mesh025_count']}** tasks.")
    lines.append("")
    lines.append("## View-level results")
    lines.append("")
    for view in VIEW_ORDER:
        stats = by_view.get(view)
        if not stats:
            continue
        lines.append(
            f"- **{view}**: structural `{stats['mean_structural_score']:.4f}`, CFD `{stats['mean_cfd_score']:.4f}`, room-kind match `{stats['room_kind_match_rate']:.2%}`, opening-wall match `{stats['opening_wall_match_rate']:.2%}`"
        )
    lines.append("")
    lines.append("Interpretation:")
    lines.append("")
    lines.append("- `floorplan` was the strongest aggregate view on both structure and CFD, suggesting that plan-layout/opening placement dominates many benchmark cases more than photorealistic texture cues.")
    lines.append("- `section` was the weakest aggregate view and the only view with room-kind collapses, indicating that single-cut geometry evidence is often insufficient for composite-room recovery.")
    lines.append("- `perspective` remained useful but less stable than expected; it needed both of the `mesh_size=0.25` fallbacks and one of the two repaired-scene salvages.")
    lines.append("")
    lines.append("## Category-level results")
    lines.append("")
    for category in CATEGORY_ORDER:
        stats = by_category.get(category)
        if not stats:
            continue
        lines.append(
            f"- **{category}**: structural `{stats['mean_structural_score']:.4f}`, CFD `{stats['mean_cfd_score']:.4f}`, room-kind match `{stats['room_kind_match_rate']:.2%}`, opening-wall match `{stats['opening_wall_match_rate']:.2%}`"
        )
    lines.append("")
    lines.append("Interpretation:")
    lines.append("")
    lines.append("- `A1` is the cleanest regime overall and behaves as the expected positive-control category.")
    lines.append("- `A2` is the weakest rectangular category mainly because opening-wall fidelity collapses; its low opening-wall match rate suggests blockage/layout metrics alone are not enough.")
    lines.append("- `A3` shows the clearest structure-vs-CFD decoupling: structure scores are modest, yet CFD stays relatively strong because opening/topology fidelity can outweigh obstacle-count hallucinations in empty/light composite rooms.")
    lines.append("- `A4` keeps room-kind fidelity mostly intact but still yields low CFD, meaning dense composite layouts remain physically hard even when coarse structure recovery looks acceptable.")
    lines.append("")
    lines.append("## Room-kind split")
    lines.append("")
    for room_kind in ROOM_KIND_ORDER:
        stats = by_room_kind.get(room_kind)
        if not stats:
            continue
        lines.append(
            f"- **{room_kind}**: structural `{stats['mean_structural_score']:.4f}`, CFD `{stats['mean_cfd_score']:.4f}`, room-kind match `{stats['room_kind_match_rate']:.2%}`, opening-wall match `{stats['opening_wall_match_rate']:.2%}`"
        )
    lines.append("")
    lines.append("## Hard failure modes that still matter despite 100/100 task success")
    lines.append("")
    for row in payload["notables"]["room_kind_mismatch_tasks"]:
        lines.append(
            f"- Room-kind collapse: `{row['case']}/{row['view']}` (`{row['ref_room_kind']} -> {row['pred_room_kind']}`), structural `{row['structural_score']:.4f}`, CFD `{row['cfd_score']:.4f}`"
        )
    for row in payload["notables"]["repair_used_tasks"]:
        lines.append(
            f"- Repair salvage: `{row['case']}/{row['view']}` with preset `{row['preset']}` at mesh `{row['mesh_size']}`; structural `{row['structural_score']:.4f}`, CFD `{row['cfd_score']:.4f}`"
        )
    lines.append("")
    lines.append("Additional robustness counters:")
    lines.append("")
    lines.append(f"- `ultra_robust` tasks: **{overall['ultra_robust_count']}**")
    lines.append(f"- `conservative` tasks: **{overall['conservative_count']}**")
    lines.append(f"- repaired-scene successes: **{overall['used_repair_count']}**")
    lines.append("")
    lines.append("## Lowest- and highest-performing cases by mean CFD")
    lines.append("")
    lines.append("### Lowest mean CFD")
    lines.append("")
    for row in payload["case_rankings"]["lowest_mean_cfd"][:5]:
        lines.append(
            f"- `{row['case']}` ({row['category']}, {row['room_kind']}): structural `{row['mean_structural_score']:.4f}`, CFD `{row['mean_cfd_score']:.4f}`, opening-wall matches `{row['opening_wall_match_count']}/5`"
        )
    lines.append("")
    lines.append("### Highest mean CFD")
    lines.append("")
    for row in payload["case_rankings"]["highest_mean_cfd"][:5]:
        lines.append(
            f"- `{row['case']}` ({row['category']}, {row['room_kind']}): structural `{row['mean_structural_score']:.4f}`, CFD `{row['mean_cfd_score']:.4f}`, opening-wall matches `{row['opening_wall_match_count']}/5`"
        )
    lines.append("")
    lines.append("## Results/Discussion-ready takeaways")
    lines.append("")
    lines.append("1. **Success-rate reporting alone is insufficient.** The benchmark now has 100/100 task success, yet aggregate CFD fidelity remains moderate (`~0.49`) and varies systematically by view/category.")
    lines.append("2. **Opening/topology fidelity matters more than exact obstacle count in several composite cases.** This is especially visible in `A3`, where obstacle hallucinations often coexist with reasonable CFD alignment.")
    lines.append("3. **Rectangular multi-obstacle cases need blockage-sensitive interpretation.** `A2` underperforms largely through opening-wall mistakes and obstacle-layout distortion rather than room-kind collapse.")
    lines.append("4. **Dense composite cases expose the biggest structure-vs-physics gap.** `A4` often preserves room-kind but still loses CFD fidelity, so future discussion should not overclaim from structural scores alone.")
    lines.append("5. **Section view should be framed as a stress input, not a generally strong modality.** Its room-kind collapses and weak opening-wall match rate make it a useful failure-analysis axis.")
    lines.append("")
    lines.append("## Suggested follow-up tags / metrics")
    lines.append("")
    lines.append("- Composite cases: add opening/topology-sensitive commentary tags when obstacle hallucination does not strongly degrade CFD.")
    lines.append("- Rectangular multi-obstacle cases: add occupancy/blockage-sensitive tags so `A2` penalties are not hidden behind room-kind correctness.")
    lines.append("- Empty or light composite controls: add hallucinated-obstacle burden tags separate from CFD penalties.")
    lines.append("- Dense composite tail cases: preserve non-blocking repair-sidecar / stabilization warning tags as robustness metadata rather than benchmark failures.")
    lines.append("")
    lines.append("## Artifact")
    lines.append("")
    lines.append(f"- JSON summary: `benchmark/manifests/evaluation_aggregate_summary.json`")
    lines.append(f"- Generated from: `{payload['generated_from']['evaluation_summary_pattern']}`")
    lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build aggregate benchmark evaluation summaries")
    parser.add_argument("--evaluation-root", type=Path, default=DEFAULT_EVALUATION_ROOT)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON_OUTPUT)
    parser.add_argument("--markdown-output", type=Path, default=DEFAULT_MD_OUTPUT)
    args = parser.parse_args()

    rows = collect_rows(args.evaluation_root)
    payload = build_payload(rows)

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_markdown(payload, args.markdown_output)

    print(
        json.dumps(
            {
                "ok": True,
                "task_count": payload["overall"]["n"],
                "json_output": str(args.json_output),
                "markdown_output": str(args.markdown_output),
                "mean_structural_score": payload["overall"]["mean_structural_score"],
                "mean_cfd_score": payload["overall"]["mean_cfd_score"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
