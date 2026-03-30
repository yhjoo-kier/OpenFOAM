#!/usr/bin/env python3
"""Fig 1: Overall framework diagram — image-to-CFD pipeline with scale calibration."""
from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_method_framework.pdf"
PNG_OUT = OUT_DIR / "fig_method_framework.png"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def draw_box(ax, x, y, w, h, text, color, edge, fontsize=8.5, bold=False, text_color="#1a1a1a"):
    box = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                          boxstyle="round,pad=0.03", facecolor=color,
                          edgecolor=edge, linewidth=1.3, zorder=3)
    ax.add_patch(box)
    fw = "bold" if bold else "normal"
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            fontweight=fw, color=text_color, zorder=4)


def draw_arrow(ax, x1, y1, x2, y2, color="#555555", lw=1.3):
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                             arrowstyle="-|>", mutation_scale=13,
                             color=color, linewidth=lw, zorder=2)
    ax.add_patch(arrow)


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5.5)
    ax.axis("off")

    bw, bh = 1.6, 0.65
    arrow_c = "#555555"

    # Stage numbers + colors
    stages = [
        # (x, y, w, h, text, bg, edge, stage_num)
        (1.2, 4.2, 1.5, 0.7, "Indoor\nimage", "#FFF9C4", "#F9A825", "Input"),
        (3.5, 4.2, 1.6, 0.7, "VLM\n(Gemini 3.1 Pro)", "#E3F2FD", "#1976D2", "Stage 1"),
        (6.0, 4.2, 1.6, 0.7, "3D scene\nJSON", "#E8F5E9", "#388E3C", ""),
        (8.5, 4.2, 1.8, 0.7, "Post-hoc scale\ncalibration", "#FCE4EC", "#C62828", "Stage 2"),
        (11.0, 4.2, 1.4, 0.7, "Calibrated\nscene", "#E8F5E9", "#388E3C", ""),
    ]

    for x, y, w, h, text, bg, edge, stage in stages:
        draw_box(ax, x, y, w, h, text, bg, edge, fontsize=8)
        if stage:
            ax.text(x, y + h / 2 + 0.18, stage, ha="center", va="bottom",
                    fontsize=7, color="#666666", fontstyle="italic")

    # Arrows for top row
    pairs_top = [(1.2, 3.5), (3.5, 6.0), (6.0, 8.5), (8.5, 11.0)]
    for x1, x2 in pairs_top:
        draw_arrow(ax, x1 + 0.75 + 0.05, 4.2, x2 - 0.8 - 0.05, 4.2)

    # "1 reference dimension" annotation on scale calibration
    ax.text(8.5, 3.35, "requires 1 reference\ndimension (longest span)",
            ha="center", fontsize=7, color="#C62828", fontstyle="italic")

    # Second row: meshing + CFD
    row2 = [
        (2.5, 2.2, 1.8, 0.7, "Gmsh\ntetrahedral mesh", "#F3E5F5", "#7B1FA2", "Stage 3"),
        (5.5, 2.2, 2.0, 0.7, "OpenFOAM\nsimpleFoam (RANS)", "#E3F2FD", "#1565C0", "Stage 4"),
        (8.8, 2.2, 1.8, 0.7, "Robustness\nladder", "#FFF3E0", "#E65100", ""),
        (11.0, 2.2, 1.4, 0.7, "Converged\nCFD solution", "#E8F5E9", "#2E7D32", "Output"),
    ]
    for x, y, w, h, text, bg, edge, stage in row2:
        draw_box(ax, x, y, w, h, text, bg, edge, fontsize=8)
        if stage:
            ax.text(x, y + h / 2 + 0.18, stage, ha="center", va="bottom",
                    fontsize=7, color="#666666", fontstyle="italic")

    # Arrow: calibrated scene → mesh
    draw_arrow(ax, 11.0, 4.2 - 0.35 - 0.05, 2.5, 2.2 + 0.35 + 0.05, color="#388E3C")

    # Arrows for bottom row
    pairs_bot = [(2.5, 5.5), (5.5, 8.8), (8.8, 11.0)]
    for x1, x2 in pairs_bot:
        w1 = next(w for xx, _, w, *_ in row2 if xx == x1)
        w2 = next(w for xx, _, w, *_ in row2 if xx == x2)
        draw_arrow(ax, x1 + w1 / 2 + 0.05, 2.2, x2 - w2 / 2 - 0.05, 2.2)

    # Robustness annotation
    ax.text(8.8, 1.35, "mesh ladder: 0.18 → 0.25 → 0.35 m\npreset escalation: robust → laminar",
            ha="center", fontsize=6.5, color="#E65100", fontstyle="italic")

    # Bottom row: evaluation (optional)
    draw_box(ax, 6.0, 0.5, 2.2, 0.6, "Evaluation: structural +\nCFD agreement scores",
             "#E8F5E9", "#2E7D32", fontsize=7.5, bold=True)
    draw_arrow(ax, 11.0, 2.2 - 0.35 - 0.05, 6.0 + 1.1 + 0.05, 0.5, color="#2E7D32")

    fig.tight_layout(pad=0.3)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
