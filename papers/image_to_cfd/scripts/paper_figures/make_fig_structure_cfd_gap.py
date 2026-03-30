#!/usr/bin/env python3
"""Fig 10: Structure-vs-CFD gap — two A4 cases with high structural but low CFD scores.

Simple 2x2 layout: reference geometry vs predicted geometry for each case.
Score badges show the structural-CFD gap. No CFD rendering needed.
"""
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
PDF_OUT = OUT_DIR / "fig_discuss_structure_cfd_gap.pdf"
PNG_OUT = OUT_DIR / "fig_discuss_structure_cfd_gap.png"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

CASES = [
    {"case": "bench_a4_02", "view": "birdseye", "scene_id": "a4_02",
     "label": "A4-02 / bird's eye"},
    {"case": "bench_a4_04", "view": "birdseye", "scene_id": "a4_04",
     "label": "A4-04 / bird's eye"},
]

COLORS = {
    "ref_room_fill": "#E8F0FE",
    "ref_room_edge": "#1565C0",
    "ref_obs_fill": "#BBDEFB",
    "ref_obs_edge": "#1565C0",
    "pred_room_fill": "#FFF3E0",
    "pred_room_edge": "#E65100",
    "pred_obs_fill": "#FFE0B2",
    "pred_obs_edge": "#E65100",
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


def draw_scene(ax, scene, room_fill, room_edge, obs_fill, obs_edge):
    blocks = room_blocks(scene)
    obs = obstacles(scene)
    for b in blocks:
        ax.add_patch(Rectangle((b["x"], b["y"]), b["dx"], b["dy"],
                                facecolor=room_fill, edgecolor=room_edge, linewidth=1.5, zorder=1))
    for o in obs:
        ax.add_patch(Rectangle((o["x"], o["y"]), o["dx"], o["dy"],
                                facecolor=obs_fill, edgecolor=obs_edge, linewidth=1.0, zorder=2))
    return blocks


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))

    panel_labels = ["(a)", "(b)", "(c)", "(d)"]
    col_titles = ["Reference geometry", "Predicted geometry"]

    for row, cfg in enumerate(CASES):
        ref_scene = load_json(SCENES_DIR / f"{cfg['scene_id']}.json")
        eval_dir = EVAL_ROOT / cfg["case"] / cfg["view"]
        pred_scene_path = eval_dir / "scaled_scene.json"
        if not pred_scene_path.exists():
            pred_scene_path = eval_dir / "predicted_scene.json"
        pred_scene = load_json(pred_scene_path)

        # Load scores
        summary = load_json(eval_dir / "evaluation_summary.json")
        pred_sum = summary.get("prediction_summary") or {}
        cfd_agg = (summary.get("cfd_summary") or {}).get("aggregate_score") or {}
        s_score = pred_sum.get("structural_score", 0.0) or 0.0
        c_score = cfd_agg.get("cfd_agreement_score", cfd_agg.get("cfd_score", 0.0)) or 0.0

        # Reference
        ax_ref = axes[row, 0]
        ref_blocks = draw_scene(ax_ref, ref_scene,
                                COLORS["ref_room_fill"], COLORS["ref_room_edge"],
                                COLORS["ref_obs_fill"], COLORS["ref_obs_edge"])

        # Predicted
        ax_pred = axes[row, 1]
        pred_blocks = draw_scene(ax_pred, pred_scene,
                                 COLORS["pred_room_fill"], COLORS["pred_room_edge"],
                                 COLORS["pred_obs_fill"], COLORS["pred_obs_edge"])

        # Shared extent
        all_blocks = ref_blocks + room_blocks(pred_scene)
        x_max = max(b["x"] + b["dx"] for b in all_blocks)
        y_max = max(b["y"] + b["dy"] for b in all_blocks)
        margin = max(x_max, y_max) * 0.08
        for ax in (ax_ref, ax_pred):
            ax.set_xlim(-margin, x_max + margin)
            ax.set_ylim(-margin, y_max + margin)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])

        # Column titles (first row only)
        if row == 0:
            ax_ref.set_title(col_titles[0], fontsize=10, fontweight="bold", pad=8)
            ax_pred.set_title(col_titles[1], fontsize=10, fontweight="bold", pad=8)

        # Row label
        ax_ref.set_ylabel(cfg["label"], fontsize=9.5, fontweight="bold", labelpad=10)

        # Score badge on predicted panel
        badge_color = "#C62828" if c_score < 0.4 else "#E65100"
        ax_pred.text(0.97, 0.03,
                     f"Structural: {s_score:.3f}\nCFD agr.: {c_score:.3f}\n"
                     f"Gap: {s_score - c_score:.3f}",
                     transform=ax_pred.transAxes, fontsize=8, ha="right", va="bottom",
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                               alpha=0.9, edgecolor=badge_color, linewidth=1.2),
                     color="#1a1a1a")

    # Panel labels
    for idx, ax in enumerate(axes.flat):
        ax.text(-0.05, 1.08, panel_labels[idx], transform=ax.transAxes,
                fontsize=11, fontweight="bold", color="#34495E")

    # Legend
    legend_elements = [
        Patch(facecolor=COLORS["ref_room_fill"], edgecolor=COLORS["ref_room_edge"], label="Reference room"),
        Patch(facecolor=COLORS["ref_obs_fill"], edgecolor=COLORS["ref_obs_edge"], label="Reference obstacle"),
        Patch(facecolor=COLORS["pred_room_fill"], edgecolor=COLORS["pred_room_edge"], label="Predicted room"),
        Patch(facecolor=COLORS["pred_obs_fill"], edgecolor=COLORS["pred_obs_edge"], label="Predicted obstacle"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=4,
               frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))

    fig.tight_layout(rect=[0.0, 0.04, 1.0, 1.0], h_pad=1.5)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
