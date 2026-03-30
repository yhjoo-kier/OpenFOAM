#!/usr/bin/env python3
"""Build aggregate image-conditioned benchmark evaluation summaries.

This consolidates the per-task `evaluation_summary.json` files into two artifacts.
By default it targets the original no-scale-hint benchmark evaluation root, but the
CLI can point it at scale-hinted roots as well.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
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
TAG_ORDER = [
    "opening_topology_preserved_with_obstacle_hallucination",
    "rectangular_blockage_layout_failure",
    "dense_composite_structure_physics_gap",
    "section_room_kind_collapse",
    "repair_salvaged_success",
    "mesh025_fallback",
    "ultra_robust_escalation",
    "laminar_fallback_success",
    "nonblocking_repair_sidecar_warning",
]
TAG_DESCRIPTIONS = {
    "opening_topology_preserved_with_obstacle_hallucination": "Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.",
    "rectangular_blockage_layout_failure": "Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.",
    "dense_composite_structure_physics_gap": "Dense composite task looked structurally acceptable but still lost CFD fidelity.",
    "section_room_kind_collapse": "Section-view task collapsed a composite room into a rectangular prediction.",
    "repair_salvaged_success": "Task required repaired-scene salvage to finish successfully.",
    "mesh025_fallback": "Task needed the smaller 0.25 mesh-size fallback to succeed.",
    "ultra_robust_escalation": "Task required ultra_robust solver escalation.",
    "laminar_fallback_success": "Task finished only after falling back to a laminar run.",
    "nonblocking_repair_sidecar_warning": "Original scene succeeded, but a repair-sidecar attempt still left warning traces worth preserving as robustness metadata.",
}


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
        prediction = data.get("prediction_summary") or {}
        cfd_agg = (data.get("cfd_summary") or {}).get("aggregate_score") or {}
        cfd_summary = cfd_agg  # alias for backward compat
        pipeline = data.get("pipeline_summary") or {}
        attempts = list(pipeline.get("attempts") or [])
        opening_metrics = prediction.get("opening_metrics") or {}
        obstacle_match = prediction.get("obstacle_match") or {}
        room_block_match = prediction.get("room_block_match") or {}
        repair_info = pipeline.get("repair_info") or {}

        case_name = str(task.get("case_name"))
        category = case_name.split("_")[1].upper() if "_" in case_name else "unknown"
        used_repair = pipeline.get("successful_scene_variant") == "repaired"
        repair_attempted = bool(repair_info.get("repair_attempted"))
        repair_available = bool(repair_info.get("repair_available"))
        nonblocking_repair_sidecar_warning = (
            repair_attempted and not repair_available and not used_repair and bool(repair_info.get("repair_error"))
        )

        rows.append(
            {
                "case": case_name,
                "view": str(task.get("view")),
                "category": category,
                "setting": data.get("setting"),
                "scale_hint": data.get("scale_hint"),
                "room_kind": prediction.get("reference_room_kind"),
                "predicted_room_kind": prediction.get("predicted_room_kind"),
                "room_kind_match": bool(prediction.get("room_kind_match")),
                "opening_wall_match": bool(prediction.get("opening_wall_match")),
                "opening_wall_ratio": opening_metrics.get("wall_match_ratio"),
                "opening_type_f1": opening_metrics.get("type_f1"),
                "opening_center_error": opening_metrics.get("mean_center_error_l2"),
                "structural_score": prediction.get("structural_score"),
                "cfd_score": cfd_summary.get("cfd_score"),
                "cfd_agreement_score": cfd_summary.get("cfd_agreement_score"),
                "cfd_components": cfd_summary.get("components"),
                "room_bbox_relative_error": prediction.get("room_bbox_relative_error") or {},
                "room_volume_relative_error": prediction.get("room_volume_relative_error"),
                "reference_obstacle_count": prediction.get("reference_obstacle_count"),
                "predicted_obstacle_count": prediction.get("predicted_obstacle_count"),
                "obstacle_count_delta": prediction.get("obstacle_count_delta"),
                "obstacle_f1": obstacle_match.get("f1"),
                "obstacle_mean_iou": obstacle_match.get("mean_iou"),
                "room_block_mean_iou": room_block_match.get("mean_iou"),
                "used_repair": used_repair,
                "repair_attempted": repair_attempted,
                "repair_available": repair_available,
                "repair_error": repair_info.get("repair_error"),
                "nonblocking_repair_sidecar_warning": nonblocking_repair_sidecar_warning,
                "successful_scene_variant": pipeline.get("successful_scene_variant"),
                "preset": pipeline.get("successful_preset"),
                "mode": pipeline.get("successful_mode"),
                "mesh_size": pipeline.get("successful_mesh_size"),
                "attempt_count": len(attempts),
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
        "mean_cfd_agreement_score": round_or_none(mean([row.get("cfd_agreement_score") for row in rows])),
        "mean_cfd_overlap": round_or_none(mean([(row.get("cfd_components") or {}).get("overlap_ratio_vs_union") for row in rows])),
        "mean_cfd_vel_mag": round_or_none(mean([(row.get("cfd_components") or {}).get("velocity_magnitude_similarity") for row in rows])),
        "mean_cfd_vel_dir": round_or_none(mean([(row.get("cfd_components") or {}).get("velocity_direction_similarity") for row in rows])),
        "mean_cfd_pressure": round_or_none(mean([(row.get("cfd_components") or {}).get("pressure_similarity") for row in rows])),
        "mean_room_bbox_rel_error_Lx": round_or_none(mean([(row["room_bbox_relative_error"] or {}).get("Lx") for row in rows])),
        "mean_room_bbox_rel_error_Ly": round_or_none(mean([(row["room_bbox_relative_error"] or {}).get("Ly") for row in rows])),
        "mean_room_bbox_rel_error_Lz": round_or_none(mean([(row["room_bbox_relative_error"] or {}).get("Lz") for row in rows])),
        "mean_room_volume_relative_error": round_or_none(mean([row["room_volume_relative_error"] for row in rows])),
        "room_kind_match_rate": round_or_none(fraction(sum(1 for row in rows if row["room_kind_match"]), n)),
        "opening_wall_match_rate": round_or_none(fraction(sum(1 for row in rows if row["opening_wall_match"]), n)),
        "used_repair_count": sum(1 for row in rows if row["used_repair"]),
        "ultra_robust_count": sum(1 for row in rows if row["preset"] == "ultra_robust"),
        "laminar_fallback_count": sum(1 for row in rows if row["preset"] == "laminar_fallback"),
        "conservative_count": sum(1 for row in rows if row["preset"] == "conservative"),
        "mesh025_count": sum(1 for row in rows if row["mesh_size"] == 0.25),
        "nonblocking_repair_sidecar_warning_count": sum(1 for row in rows if row["nonblocking_repair_sidecar_warning"]),
    }


def group_summary(rows: list[dict[str, Any]], key: str, order: list[str]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[key])].append(row)
    keys = [value for value in order if value in grouped] + sorted(k for k in grouped if k not in order)
    return {group_key: summarize_rows(grouped[group_key]) for group_key in keys}


def derive_task_tags(row: dict[str, Any]) -> list[str]:
    tags: list[str] = []
    structural = row.get("structural_score")
    cfd = row.get("cfd_score")
    obstacle_delta = row.get("obstacle_count_delta") or 0
    obstacle_f1 = row.get("obstacle_f1")
    opening_ratio = row.get("opening_wall_ratio")

    if (
        row.get("room_kind") == "composite"
        and row.get("room_kind_match")
        and obstacle_delta > 0
        and (opening_ratio is not None and opening_ratio >= 0.5)
        and (cfd is not None and cfd >= 0.55)
    ):
        tags.append("opening_topology_preserved_with_obstacle_hallucination")

    if (
        row.get("room_kind") == "rectangular"
        and (row.get("reference_obstacle_count") or 0) >= 2
        and ((obstacle_f1 is not None and obstacle_f1 <= 0.5) or not row.get("opening_wall_match"))
        and (cfd is not None and cfd <= 0.5)
    ):
        tags.append("rectangular_blockage_layout_failure")

    if (
        row.get("category") == "A4"
        and row.get("room_kind") == "composite"
        and (structural is not None and structural >= 0.75)
        and (cfd is not None and cfd <= 0.5)
    ):
        tags.append("dense_composite_structure_physics_gap")

    if row.get("view") == "section" and not row.get("room_kind_match"):
        tags.append("section_room_kind_collapse")
    if row.get("used_repair"):
        tags.append("repair_salvaged_success")
    if row.get("mesh_size") == 0.25:
        tags.append("mesh025_fallback")
    if row.get("preset") == "ultra_robust":
        tags.append("ultra_robust_escalation")
    if row.get("preset") == "laminar_fallback":
        tags.append("laminar_fallback_success")
    if row.get("nonblocking_repair_sidecar_warning"):
        tags.append("nonblocking_repair_sidecar_warning")
    return tags


def case_cards(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row["case"]].append(row)

    cards: list[dict[str, Any]] = []
    for case_name in sorted(grouped):
        case_rows = grouped[case_name]
        tag_counts = Counter(tag for row in case_rows for tag in row.get("derived_tags", []))
        cards.append(
            {
                "case": case_name,
                "category": case_rows[0]["category"],
                "room_kind": case_rows[0]["room_kind"],
                "mean_structural_score": round_or_none(mean([row["structural_score"] for row in case_rows])),
                "mean_cfd_score": round_or_none(mean([row["cfd_score"] for row in case_rows])),
                "mean_room_bbox_rel_error_Lx": round_or_none(mean([(row["room_bbox_relative_error"] or {}).get("Lx") for row in case_rows])),
                "mean_room_bbox_rel_error_Ly": round_or_none(mean([(row["room_bbox_relative_error"] or {}).get("Ly") for row in case_rows])),
                "mean_room_bbox_rel_error_Lz": round_or_none(mean([(row["room_bbox_relative_error"] or {}).get("Lz") for row in case_rows])),
                "room_kind_match_count": sum(1 for row in case_rows if row["room_kind_match"]),
                "opening_wall_match_count": sum(1 for row in case_rows if row["opening_wall_match"]),
                "used_repair_views": [row["view"] for row in case_rows if row["used_repair"]],
                "ultra_robust_views": [row["view"] for row in case_rows if row["preset"] == "ultra_robust"],
                "laminar_fallback_views": [row["view"] for row in case_rows if row["preset"] == "laminar_fallback"],
                "mesh025_views": [row["view"] for row in case_rows if row["mesh_size"] == 0.25],
                "hallucination_views": [row["view"] for row in case_rows if (row["obstacle_count_delta"] or 0) > 0],
                "nonblocking_repair_sidecar_warning_views": [
                    row["view"] for row in case_rows if row["nonblocking_repair_sidecar_warning"]
                ],
                "derived_tag_counts": dict(sorted(tag_counts.items())),
                "dominant_tags": [tag for tag, count in tag_counts.items() if count >= 2],
            }
        )
    return cards


def build_tag_summary(rows: list[dict[str, Any]], case_cards_payload: list[dict[str, Any]]) -> dict[str, Any]:
    task_index: dict[str, list[dict[str, Any]]] = {}
    for tag in TAG_ORDER:
        tagged_rows = [row for row in rows if tag in row.get("derived_tags", [])]
        task_index[tag] = [
            {
                "case": row["case"],
                "view": row["view"],
                "category": row["category"],
                "room_kind": row["room_kind"],
                "structural_score": round_or_none(row["structural_score"]),
                "cfd_score": round_or_none(row["cfd_score"]),
            }
            for row in tagged_rows
        ]

    case_index: dict[str, list[dict[str, Any]]] = {}
    for tag in TAG_ORDER:
        tagged_cases = [card for card in case_cards_payload if tag in set(card.get("dominant_tags", []))]
        case_index[tag] = [
            {
                "case": card["case"],
                "category": card["category"],
                "room_kind": card["room_kind"],
                "mean_structural_score": card["mean_structural_score"],
                "mean_cfd_score": card["mean_cfd_score"],
            }
            for card in tagged_cases
        ]

    return {
        "task_counts": {tag: len(task_index[tag]) for tag in TAG_ORDER},
        "task_examples": {tag: task_index[tag][:10] for tag in TAG_ORDER},
        "case_examples": {tag: case_index[tag][:10] for tag in TAG_ORDER},
        "descriptions": TAG_DESCRIPTIONS,
    }


def build_payload(rows: list[dict[str, Any]], evaluation_root: Path) -> dict[str, Any]:
    for row in rows:
        row["derived_tags"] = derive_task_tags(row)

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
            "derived_tags": row["derived_tags"],
        }
        for row in rows
    ]

    payload = {
        "ok": True,
        "project_root": str(PROJECT_ROOT),
        "evaluation_root": str(evaluation_root),
        "setting": rows[0].get("setting") if rows else None,
        "overall": summarize_rows(rows),
        "by_view": group_summary(rows, "view", VIEW_ORDER),
        "by_category": group_summary(rows, "category", CATEGORY_ORDER),
        "by_room_kind": group_summary(rows, "room_kind", ROOM_KIND_ORDER),
        "case_rankings": {
            "lowest_mean_cfd": cards_low[:10],
            "highest_mean_cfd": cards_high[:10],
        },
        "derived_tags": build_tag_summary(rows, cards),
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
            "laminar_fallback_tasks": [
                {
                    "case": row["case"],
                    "view": row["view"],
                    "mesh_size": row["mesh_size"],
                    "structural_score": round_or_none(row["structural_score"]),
                    "cfd_score": round_or_none(row["cfd_score"]),
                }
                for row in rows
                if row["preset"] == "laminar_fallback"
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
            "nonblocking_repair_sidecar_warning_tasks": [
                {
                    "case": row["case"],
                    "view": row["view"],
                    "preset": row["preset"],
                    "structural_score": round_or_none(row["structural_score"]),
                    "cfd_score": round_or_none(row["cfd_score"]),
                }
                for row in rows
                if row["nonblocking_repair_sidecar_warning"]
            ],
            "highest_cfd_tasks": sorted(task_cards, key=lambda card: card["cfd_score"] if card["cfd_score"] is not None else -1.0, reverse=True)[:10],
            "lowest_cfd_tasks": sorted(task_cards, key=lambda card: card["cfd_score"] if card["cfd_score"] is not None else 999.0)[:10],
        },
        "generated_from": {
            "task_count": len(rows),
            "evaluation_summary_pattern": f"{evaluation_root}/*/*/evaluation_summary.json",
        },
    }
    return payload


def write_markdown(payload: dict[str, Any], path: Path) -> None:
    overall = payload["overall"]
    by_view = payload["by_view"]
    by_category = payload["by_category"]
    by_room_kind = payload["by_room_kind"]
    tag_summary = payload["derived_tags"]

    setting = payload.get("setting") or "unknown_setting"
    title = f"# Benchmark Evaluation Aggregate Results ({setting})"

    lines: list[str] = []
    lines.append(title)
    lines.append("")
    lines.append("> Date: 2026-03-19")
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append(f"- Evaluation root: `{payload['evaluation_root']}`")
    lines.append(f"- Task count in this aggregate: **{overall['n']}**")
    lines.append(f"- Aggregate mean structural score: **{overall['mean_structural_score']}**")
    lines.append(f"- Aggregate mean CFD score: **{overall['mean_cfd_score']}**")
    lines.append(f"- Mean room bbox relative error: `Lx={overall['mean_room_bbox_rel_error_Lx']}`, `Ly={overall['mean_room_bbox_rel_error_Ly']}`, `Lz={overall['mean_room_bbox_rel_error_Lz']}`")
    lines.append(f"- Mean room volume relative error: **{overall['mean_room_volume_relative_error']}**")
    lines.append(f"- Room-kind match rate: **{overall['room_kind_match_rate']}**")
    lines.append(f"- Opening-wall match rate: **{overall['opening_wall_match_rate']}**")
    lines.append(f"- Repair was needed in **{overall['used_repair_count']}** tasks; mesh-size `0.25` fallback appeared in **{overall['mesh025_count']}** tasks.")
    lines.append(f"- Laminar fallback count: **{overall['laminar_fallback_count']}**")
    lines.append("")
    lines.append("## View-level results")
    lines.append("")
    for view in VIEW_ORDER:
        stats = by_view.get(view)
        if not stats:
            continue
        lines.append(
            f"- **{view}**: structural `{stats['mean_structural_score']}`, CFD `{stats['mean_cfd_score']}`, bbox-rel `(Lx={stats['mean_room_bbox_rel_error_Lx']}, Ly={stats['mean_room_bbox_rel_error_Ly']}, Lz={stats['mean_room_bbox_rel_error_Lz']})`, opening-wall match `{stats['opening_wall_match_rate']}`"
        )
    lines.append("")
    lines.append("## Category-level results")
    lines.append("")
    for category in CATEGORY_ORDER:
        stats = by_category.get(category)
        if not stats:
            continue
        lines.append(
            f"- **{category}**: structural `{stats['mean_structural_score']}`, CFD `{stats['mean_cfd_score']}`, bbox-rel `(Lx={stats['mean_room_bbox_rel_error_Lx']}, Ly={stats['mean_room_bbox_rel_error_Ly']}, Lz={stats['mean_room_bbox_rel_error_Lz']})`, opening-wall match `{stats['opening_wall_match_rate']}`"
        )
    lines.append("")
    lines.append("## Room-kind split")
    lines.append("")
    for room_kind in ROOM_KIND_ORDER:
        stats = by_room_kind.get(room_kind)
        if not stats:
            continue
        lines.append(
            f"- **{room_kind}**: structural `{stats['mean_structural_score']}`, CFD `{stats['mean_cfd_score']}`, bbox-rel `(Lx={stats['mean_room_bbox_rel_error_Lx']}, Ly={stats['mean_room_bbox_rel_error_Ly']}, Lz={stats['mean_room_bbox_rel_error_Lz']})`, opening-wall match `{stats['opening_wall_match_rate']}`"
        )
    lines.append("")
    lines.append("## Derived interpretation tags")
    lines.append("")
    for tag in TAG_ORDER:
        count = tag_summary["task_counts"].get(tag, 0)
        if count <= 0:
            continue
        lines.append(f"- **{tag}**: {count} tasks")
        lines.append(f"  - {tag_summary['descriptions'][tag]}")
    lines.append("")
    lines.append("## Lowest- and highest-performing cases by mean CFD")
    lines.append("")
    lines.append("### Lowest mean CFD")
    lines.append("")
    for row in payload["case_rankings"]["lowest_mean_cfd"][:5]:
        lines.append(
            f"- `{row['case']}` ({row['category']}, {row['room_kind']}): structural `{row['mean_structural_score']}`, CFD `{row['mean_cfd_score']}`, bbox-rel `(Lx={row['mean_room_bbox_rel_error_Lx']}, Ly={row['mean_room_bbox_rel_error_Ly']}, Lz={row['mean_room_bbox_rel_error_Lz']})`"
        )
    lines.append("")
    lines.append("### Highest mean CFD")
    lines.append("")
    for row in payload["case_rankings"]["highest_mean_cfd"][:5]:
        lines.append(
            f"- `{row['case']}` ({row['category']}, {row['room_kind']}): structural `{row['mean_structural_score']}`, CFD `{row['mean_cfd_score']}`, bbox-rel `(Lx={row['mean_room_bbox_rel_error_Lx']}, Ly={row['mean_room_bbox_rel_error_Ly']}, Lz={row['mean_room_bbox_rel_error_Lz']})`"
        )
    lines.append("")
    lines.append("## Artifact")
    lines.append("")
    lines.append(f"- JSON summary: `{path.with_suffix('.json').name}` (or CLI-provided target)")
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

    evaluation_root = args.evaluation_root.expanduser().resolve()
    rows = collect_rows(evaluation_root)
    payload = build_payload(rows, evaluation_root)

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
                "derived_tag_counts": payload["derived_tags"]["task_counts"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
