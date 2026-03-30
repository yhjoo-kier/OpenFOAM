#!/usr/bin/env python3
"""Figure 1: Overall framework with mandatory scale calibration.

End-to-end pipeline from 2D image to CFD results, with evaluation pathway.
Publication-quality, double-column width, clean academic style.

Layout: 2-row design
  Top row:    2D Image -> VLM -> Scale Cal. -> 3D Scene -> OF Setup -> CFD -> Results
  Bottom row: Reference -> Comparison -> Scores
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_method_framework.pdf"
PNG_OUT = OUT_DIR / "fig_method_framework.png"

# ---------------------------------------------------------------------------
# Font selection
# ---------------------------------------------------------------------------
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


FONT = pick_font()

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
C_INPUT = "#DBEAFE"        # light blue
C_VLM = "#EDE9FE"          # light purple
C_SCALE = "#FEF3C7"        # light amber
C_SCALE_BORDER = "#D97706" # amber-600 for emphasis
C_PROCESS = "#D1FAE5"      # light green
C_OUTPUT = "#F1F5F9"       # light gray
C_EVAL = "#FEE2E2"         # light red for evaluation
C_ARROW = "#64748B"        # slate-500
C_TEXT = "#1E293B"          # slate-800
C_SUBTEXT = "#64748B"      # slate-500
C_EVAL_ARROW = "#94A3B8"   # slate-400

# ---------------------------------------------------------------------------
# Pipeline stage definitions
# ---------------------------------------------------------------------------
STAGES = [
    ("2D Input\nImage",       "photograph /\narch. drawing",     C_INPUT,   "#93C5FD"),
    ("VLM\nAbstraction",      "Gemini 2.0\nFlash",              C_VLM,     "#C4B5FD"),
    ("Scale\nCalibration",    "post-hoc\nlongest-span",          C_SCALE,   C_SCALE_BORDER),
    ("3D Scene\nJSON",        "room, obstacles,\nopenings",      C_PROCESS, "#6EE7B7"),
    ("OpenFOAM\nSetup",       "blockMesh +\nboundary cond.",     C_PROCESS, "#6EE7B7"),
    ("CFD\nSimulation",       "steady RANS\n(simpleFoam)",       C_PROCESS, "#6EE7B7"),
    ("Flow\nResults",         "velocity /\npressure fields",     C_OUTPUT,  "#CBD5E1"),
]

EVAL_STAGES = [
    ("Reference\nResults",    "rule-based\nground truth",        C_OUTPUT,  "#CBD5E1"),
    ("Metric\nComparison",    "structural +\nCFD scores",        C_EVAL,    "#FCA5A5"),
    ("Evaluation\nScores",    "per-view /\nper-category",        C_EVAL,    "#FCA5A5"),
]


# ---------------------------------------------------------------------------
# Drawing helpers -- all in axes-fraction coordinates [0,1]
# ---------------------------------------------------------------------------
def draw_box(ax, cx, cy, w, h, label, facecolor, edgecolor,
             linewidth=1.4, fontsize=11, bold=False, extra_lw=0.0,
             transform=None):
    """Draw a rounded box with centered label."""
    t = transform or ax.transAxes
    box = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.012",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth + extra_lw,
        transform=t,
        zorder=2,
    )
    ax.add_patch(box)
    weight = "bold" if bold else "medium"
    ax.text(cx, cy, label,
            ha="center", va="center",
            fontsize=fontsize, fontweight=weight,
            color=C_TEXT, fontfamily=FONT,
            linespacing=1.20, zorder=3,
            transform=t)


def draw_detail(ax, cx, cy, text, fontsize=7.5, transform=None):
    """Draw italic detail text (below a box)."""
    t = transform or ax.transAxes
    ax.text(cx, cy, text,
            ha="center", va="top",
            fontsize=fontsize, fontweight="normal",
            color=C_SUBTEXT, fontfamily=FONT,
            linespacing=1.10, zorder=3,
            style="italic", transform=t)


def h_arrow(ax, x0, y0, x1, y1, color=C_ARROW, lw=1.4, ms=12, transform=None):
    t = transform or ax.transAxes
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="-|>", color=color,
        linewidth=lw, mutation_scale=ms, zorder=1,
        transform=t))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    t = ax.transAxes

    # -------------------------------------------------------------------
    # Layout constants (axes fraction coordinates)
    # -------------------------------------------------------------------
    # Top row
    N_TOP = len(STAGES)
    TOP_Y = 0.72                # vertical center of top row boxes
    BOX_W = 0.088               # box width
    BOX_H = 0.17                # box height
    SCALE_W = 0.100             # wider for Scale Calibration
    TOP_LEFT = 0.065            # left edge of first box center
    TOP_RIGHT = 0.935           # right edge of last box center
    top_gap = (TOP_RIGHT - TOP_LEFT) / (N_TOP - 1)

    top_x = [TOP_LEFT + i * top_gap for i in range(N_TOP)]

    # Bottom row
    N_EVAL = len(EVAL_STAGES)
    EVAL_Y = 0.22               # vertical center of eval row
    EVAL_W = 0.105
    # Eval boxes centered under right portion of top row
    eval_left = top_x[3] - 0.02
    eval_right = top_x[6] + 0.02
    eval_gap = (eval_right - eval_left) / (N_EVAL - 1)
    eval_x = [eval_left + i * eval_gap for i in range(N_EVAL)]

    # -------------------------------------------------------------------
    # Draw top row
    # -------------------------------------------------------------------
    for i, (label, detail, fc, ec) in enumerate(STAGES):
        cx = top_x[i]
        is_scale = (i == 2)
        w = SCALE_W if is_scale else BOX_W
        draw_box(ax, cx, TOP_Y, w, BOX_H,
                 label, fc, ec,
                 linewidth=2.4 if is_scale else 1.4,
                 fontsize=9.5,
                 bold=is_scale,
                 extra_lw=1.0 if is_scale else 0.0,
                 transform=t)
        # Detail below box
        draw_detail(ax, cx, TOP_Y - BOX_H / 2 - 0.015, detail, fontsize=7, transform=t)
        # "Key contribution" above Scale Calibration
        if is_scale:
            ax.text(cx, TOP_Y + BOX_H / 2 + 0.025, "Key contribution",
                    ha="center", va="bottom",
                    fontsize=8, fontweight="bold",
                    color=C_SCALE_BORDER, fontfamily=FONT, zorder=3,
                    transform=t)

    # Arrows between top-row stages
    for i in range(N_TOP - 1):
        w_cur = SCALE_W if i == 2 else BOX_W
        w_nxt = SCALE_W if (i + 1) == 2 else BOX_W
        x0 = top_x[i] + w_cur / 2 + 0.008
        x1 = top_x[i + 1] - w_nxt / 2 - 0.008
        h_arrow(ax, x0, TOP_Y, x1, TOP_Y, transform=t)

    # -------------------------------------------------------------------
    # Draw bottom row (evaluation)
    # -------------------------------------------------------------------
    for i, (label, detail, fc, ec) in enumerate(EVAL_STAGES):
        cx = eval_x[i]
        draw_box(ax, cx, EVAL_Y, EVAL_W, BOX_H,
                 label, fc, ec,
                 fontsize=10, transform=t)
        draw_detail(ax, cx, EVAL_Y - BOX_H / 2 - 0.015, detail, fontsize=7, transform=t)

    # Arrows between eval stages
    for i in range(N_EVAL - 1):
        x0 = eval_x[i] + EVAL_W / 2 + 0.008
        x1 = eval_x[i + 1] - EVAL_W / 2 - 0.008
        h_arrow(ax, x0, EVAL_Y, x1, EVAL_Y, color=C_EVAL_ARROW, transform=t)

    # --- Connector: Flow Results -> Metric Comparison (L-shaped dashed) ---
    flow_cx = top_x[6]
    comp_cx = eval_x[1]
    bend_y = EVAL_Y + BOX_H / 2 + 0.06

    # Vertical dashed from Flow Results down
    ax.plot([flow_cx, flow_cx],
            [TOP_Y - BOX_H / 2 - 0.01, bend_y],
            color=C_EVAL_ARROW, linewidth=1.0,
            linestyle=(0, (4, 3)), zorder=0,
            transform=t, clip_on=False)
    # Horizontal dashed leftward to above Comparison
    ax.plot([flow_cx, comp_cx + EVAL_W / 2 + 0.02],
            [bend_y, bend_y],
            color=C_EVAL_ARROW, linewidth=1.0,
            linestyle=(0, (4, 3)), zorder=0,
            transform=t, clip_on=False)
    # Arrow into Comparison from above-right
    h_arrow(ax,
            comp_cx + EVAL_W / 2 + 0.02, bend_y,
            comp_cx + EVAL_W / 2 + 0.005, EVAL_Y + BOX_H / 2 + 0.008,
            color=C_EVAL_ARROW, lw=1.0, ms=10, transform=t)

    # "predicted" label on vertical dashed
    mid_y = (TOP_Y - BOX_H / 2 + bend_y) / 2
    ax.text(flow_cx - 0.018, mid_y,
            "predicted",
            ha="right", va="center", rotation=90,
            fontsize=7.5, color=C_SUBTEXT,
            fontfamily=FONT, style="italic", zorder=3,
            transform=t)

    # "Evaluation pathway" label to the left of bottom row
    ax.text(eval_x[0] - EVAL_W / 2 - 0.02, EVAL_Y,
            "Evaluation\npathway",
            ha="right", va="center",
            fontsize=10, fontweight="medium",
            color=C_SUBTEXT, fontfamily=FONT,
            linespacing=1.25, zorder=3,
            transform=t)

    # -------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------
    fig.savefig(str(PDF_OUT), format="pdf", bbox_inches="tight", dpi=600)
    fig.savefig(str(PNG_OUT), format="png", bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"Saved: {PDF_OUT}")
    print(f"Saved: {PNG_OUT}")


if __name__ == "__main__":
    main()
