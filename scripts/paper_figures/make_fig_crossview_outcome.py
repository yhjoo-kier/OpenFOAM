#!/usr/bin/env python3
"""Fig 7: Representative cross-view outcome — 2 cases x 3 views."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
SCENES_DIR = PROJECT_ROOT / "benchmark" / "scenes"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_result_crossview_outcome.pdf"
PNG_OUT = OUT_DIR / "fig_result_crossview_outcome.png"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

# 2 cases x 3 views = 6 panels
CASES = [
    {"case": "bench_a2_03", "label": "A2-03 (rect. dense)", "scene_id": "a2_03"},
    {"case": "bench_a3_03", "label": "A3-03 (comp. simple)", "scene_id": "a3_03"},
]
VIEWS = ["floorplan", "perspective", "section"]
VIEW_LABELS = ["Floor plan", "Perspective", "Section"]

COLORS = {
    "room_fill": "#EEF3F8",
    "room_edge": "#334155",
    "ref_obstacle_fill": "#C8CED8",
    "ref_obstacle_edge": "#5B6573",
    "matched_fill": "#F5E6CC",
    "matched_edge": "#8B6A2F",
    "extra_fill": "#E8863C",
    "extra_edge": "#8B3A12",
    "gt_outline": "#7C3AED",
    "inlet": "#2563EB",
    "outlet": "#DC2626",
}


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def load_json(p: Path) -> Any:
    return json.loads(p.read_text(encoding="utf-8"))


def room_blocks(scene: dict) -> list[dict]:
    room = scene["room"]
    if "blocks" in room:
        return [{"x": float(b["origin"]["x"]), "y": float(b["origin"]["y"]),
                 "dx": float(b["size"]["dx"]), "dy": float(b["size"]["dy"])} for b in room["blocks"]]
    s = room["size"]
    return [{"x": 0.0, "y": 0.0, "dx": float(s["Lx"]), "dy": float(s["Ly"])}]


def obstacles(scene: dict) -> list[dict]:
    result = []
    for o in scene.get("obstacles", []):
        orig = o.get("origin") or o.get("min") or {"x": 0, "y": 0}
        size = o.get("size", {})
        result.append({"x": float(orig["x"]), "y": float(orig["y"]),
                       "dx": float(size.get("dx", 0)), "dy": float(size.get("dy", 0))})
    return result


def draw_geometry(ax, ref_scene, pred_scene, title, scores):
    """Draw predicted geometry with reference outline overlay."""
    ref_blocks = room_blocks(ref_scene)
    pred_blocks = room_blocks(pred_scene)
    ref_obs = obstacles(ref_scene)
    pred_obs = obstacles(pred_scene)

    # Draw predicted room blocks
    for b in pred_blocks:
        ax.add_patch(Rectangle((b["x"], b["y"]), b["dx"], b["dy"],
                                facecolor=COLORS["room_fill"], edgecolor=COLORS["room_edge"],
                                linewidth=1.2, zorder=1))

    # Draw reference room outline (dashed)
    for b in ref_blocks:
        ax.add_patch(Rectangle((b["x"], b["y"]), b["dx"], b["dy"],
                                facecolor="none", edgecolor=COLORS["gt_outline"],
                                linewidth=1.0, linestyle="--", zorder=4))

    # Draw predicted obstacles
    for o in pred_obs:
        ax.add_patch(Rectangle((o["x"], o["y"]), o["dx"], o["dy"],
                                facecolor=COLORS["matched_fill"], edgecolor=COLORS["matched_edge"],
                                linewidth=0.8, zorder=2))

    # Draw reference obstacles (dashed outline)
    for o in ref_obs:
        ax.add_patch(Rectangle((o["x"], o["y"]), o["dx"], o["dy"],
                                facecolor="none", edgecolor=COLORS["gt_outline"],
                                linewidth=0.8, linestyle="--", zorder=4))

    # Auto-scale
    all_blocks = pred_blocks + ref_blocks
    x_max = max(b["x"] + b["dx"] for b in all_blocks) if all_blocks else 10
    y_max = max(b["y"] + b["dy"] for b in all_blocks) if all_blocks else 10
    margin = max(x_max, y_max) * 0.08
    ax.set_xlim(-margin, x_max + margin)
    ax.set_ylim(-margin, y_max + margin)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=8.5, pad=4)

    # Score badge
    ax.text(0.97, 0.03, f"S={scores['structural']:.2f}\nC={scores['cfd']:.2f}",
            transform=ax.transAxes, fontsize=7, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85, edgecolor="#CBD5E1"))

    ax.set_xticks([])
    ax.set_yticks([])


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0))

    for row, case_cfg in enumerate(CASES):
        ref_scene = load_json(SCENES_DIR / f"{case_cfg['scene_id']}.json")

        for col, view in enumerate(VIEWS):
            ax = axes[row, col]
            eval_dir = EVAL_ROOT / case_cfg["case"] / view

            # Load predicted scene
            pred_scene_path = eval_dir / "scaled_scene.json"
            if not pred_scene_path.exists():
                pred_scene_path = eval_dir / "predicted_scene.json"

            # Load scores
            summary_path = eval_dir / "evaluation_summary.json"
            scores = {"structural": 0.0, "cfd": 0.0}
            if summary_path.exists():
                summary = load_json(summary_path)
                pred_sum = summary.get("prediction_summary") or {}
                cfd_agg = (summary.get("cfd_summary") or {}).get("aggregate_score") or {}
                scores["structural"] = pred_sum.get("structural_score", 0.0) or 0.0
                scores["cfd"] = cfd_agg.get("cfd_agreement_score", cfd_agg.get("cfd_score", 0.0)) or 0.0

            if pred_scene_path.exists():
                pred_scene = load_json(pred_scene_path)
                title = VIEW_LABELS[col] if row == 0 else ""
                draw_geometry(ax, ref_scene, pred_scene, title, scores)
            else:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                ax.axis("off")

        # Row label
        axes[row, 0].set_ylabel(case_cfg["label"], fontsize=9, fontweight="bold", labelpad=10)

    # Panel labels
    labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    for idx, ax in enumerate(axes.flat):
        ax.text(-0.08, 1.08, labels[idx], transform=ax.transAxes,
                fontsize=10, fontweight="bold", color="#34495E")

    # Legend
    legend_elements = [
        Patch(facecolor=COLORS["room_fill"], edgecolor=COLORS["room_edge"], label="Predicted room"),
        Patch(facecolor=COLORS["matched_fill"], edgecolor=COLORS["matched_edge"], label="Predicted obstacle"),
        Patch(facecolor="none", edgecolor=COLORS["gt_outline"], linestyle="--", label="Reference outline"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.01))

    fig.tight_layout(rect=[0.0, 0.04, 1.0, 1.0], h_pad=1.5)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
