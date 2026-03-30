#!/usr/bin/env python3
"""Fig 4: Evaluation pathway diagram — predicted vs reference paths."""
from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_method_eval_pathway.pdf"
PNG_OUT = OUT_DIR / "fig_method_eval_pathway.png"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def draw_box(ax, x, y, w, h, text, color="#E8F0FE", edge="#4285F4", fontsize=8.5, bold=False):
    box = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                          boxstyle="round,pad=0.02", facecolor=color,
                          edgecolor=edge, linewidth=1.2, zorder=3)
    ax.add_patch(box)
    fw = "bold" if bold else "normal"
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            fontweight=fw, color="#1a1a1a", zorder=4)


def draw_arrow(ax, x1, y1, x2, y2, color="#666666"):
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                             arrowstyle="-|>", mutation_scale=12,
                             color=color, linewidth=1.2, zorder=2)
    ax.add_patch(arrow)


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(7.2, 3.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.5)
    ax.axis("off")

    bw, bh = 1.4, 0.55  # box width, height

    # === Predicted path (top row, y=3.5) ===
    y_pred = 3.4
    pred_color = "#E8F4FD"
    pred_edge = "#2196F3"

    ax.text(0.3, y_pred + 0.5, "Predicted path", fontsize=10, fontweight="bold",
            color="#1565C0", va="center")

    boxes_pred = [
        (1.5, "2D Image"),
        (3.3, "VLM\n(Gemini)"),
        (5.1, "Scale\ncalibration"),
        (6.9, "Mesh +\nOpenFOAM"),
        (8.5, "Predicted\nCFD"),
    ]
    for x, label in boxes_pred:
        draw_box(ax, x, y_pred, bw, bh, label, color=pred_color, edge=pred_edge)
    for i in range(len(boxes_pred) - 1):
        draw_arrow(ax, boxes_pred[i][0] + bw / 2 + 0.02, y_pred,
                   boxes_pred[i + 1][0] - bw / 2 - 0.02, y_pred, color=pred_edge)

    # === Reference path (bottom row, y=1.5) ===
    y_ref = 1.5
    ref_color = "#FFF3E0"
    ref_edge = "#E65100"

    ax.text(0.3, y_ref + 0.5, "Reference path", fontsize=10, fontweight="bold",
            color="#BF360C", va="center")

    boxes_ref = [
        (1.5, "Rule-based\ngenerator"),
        (3.3, "Ground-truth\nscene JSON"),
        (5.1, "2D rendering\n(for VLM input)"),
        (6.9, "Mesh +\nOpenFOAM"),
        (8.5, "Reference\nCFD"),
    ]
    for x, label in boxes_ref:
        draw_box(ax, x, y_ref, bw, bh, label, color=ref_color, edge=ref_edge)
    for i in range(len(boxes_ref) - 1):
        draw_arrow(ax, boxes_ref[i][0] + bw / 2 + 0.02, y_ref,
                   boxes_ref[i + 1][0] - bw / 2 - 0.02, y_ref, color=ref_edge)

    # Cross-path arrow: 2D rendering → 2D Image (provides input)
    draw_arrow(ax, 5.1, y_ref + bh / 2 + 0.02, 1.5, y_pred - bh / 2 - 0.02, color="#9E9E9E")
    ax.text(2.8, 2.45, "input to VLM", fontsize=7.5, color="#757575",
            ha="center", fontstyle="italic", rotation=30)

    # Cross-path arrow: scale calibration uses one ref dimension
    draw_arrow(ax, 3.3, y_ref + bh / 2 + 0.02, 5.1, y_pred - bh / 2 - 0.02, color="#9E9E9E")
    ax.text(4.5, 2.45, "1 reference\ndimension", fontsize=7.5, color="#757575",
            ha="center", fontstyle="italic")

    # === Metric comparison (center-right) ===
    y_metric = 2.45
    draw_box(ax, 9.3, y_metric, 1.3, 0.8, "Metric\ncomparison",
             color="#E8F5E9", edge="#2E7D32", fontsize=9, bold=True)

    # Arrows from both CFD to metric comparison
    draw_arrow(ax, 8.5, y_pred - bh / 2 - 0.02, 9.3, y_metric + 0.4 + 0.02, color="#2196F3")
    draw_arrow(ax, 8.5, y_ref + bh / 2 + 0.02, 9.3, y_metric - 0.4 - 0.02, color="#E65100")

    # "No VLM" annotation on reference path
    ax.text(4.2, y_ref - 0.55, "No VLM involvement — no data leakage",
            fontsize=8, color="#BF360C", fontstyle="italic", ha="center")

    fig.tight_layout(pad=0.3)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
