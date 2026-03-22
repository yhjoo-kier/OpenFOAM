#!/usr/bin/env python3
"""Figure 4 – Reference benchmark and evaluation pathway.

Flow diagram showing how the evaluation pipeline works:
  2D image → VLM → 3D scene JSON → Scale Calibration → CFD pipeline → predicted results
  Reference case → CFD pipeline → reference results
  Both converge at Comparison → Structural Score + CFD Score.

Publication-quality, double-column width, clean academic style.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_method_eval_pathway.pdf"
PNG_OUT = OUT_DIR / "fig_method_eval_pathway.png"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

# Colour palette -----------------------------------------------------------
C_INPUT = "#D6EAF8"
C_INPUT_EDGE = "#5B9BD5"
C_PROCESS = "#D5F5E3"
C_PROCESS_EDGE = "#58B87A"
C_OUTPUT = "#FDEBD0"
C_OUTPUT_EDGE = "#E8A838"
C_COMPARE = "#E8DAEF"
C_COMPARE_EDGE = "#9B72B0"
C_SCALE = "#FFF3CD"
C_SCALE_EDGE = "#D4A843"
C_ARROW = "#4A4A4A"
C_TEXT = "#2C3E50"
C_LABEL = "#7F8C8D"


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


def _box(ax, cx, cy, w, h, text, fc, ec, fontsize=10.0, bold=False):
    """Draw a rounded box centred at (cx, cy) with ONLY its main label."""
    rect = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.05",
        facecolor=fc, edgecolor=ec, linewidth=1.3,
        transform=ax.transData, zorder=2,
    )
    ax.add_patch(rect)
    ax.text(cx, cy, text, ha="center", va="center",
            fontsize=fontsize, color=C_TEXT,
            fontweight="bold" if bold else "normal", zorder=3,
            linespacing=1.15)
    return rect


def _subtitle(ax, cx, cy, text, fontsize=7.5):
    """Small italic annotation BELOW a box."""
    ax.text(cx, cy, text, ha="center", va="top",
            fontsize=fontsize, color=C_LABEL, style="italic", zorder=3,
            linespacing=1.1)


def _arrow(ax, x0, y0, x1, y1, style="-|>", lw=1.3, color=None):
    """Draw an arrow with visible arrowhead."""
    arrow = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle=style, mutation_scale=14,
        lw=lw, color=color or C_ARROW, zorder=1,
        shrinkA=0, shrinkB=0,
    )
    ax.add_patch(arrow)
    return arrow


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    selected_font = pick_font()
    plt.rcParams.update({
        "font.family": selected_font,
        "font.size": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.set_xlim(0, 10.5)
    ax.set_ylim(0, 5)
    ax.set_aspect("auto")
    ax.axis("off")

    # ── Dimensions ──
    bw = 1.05    # standard box width
    bh = 0.72    # standard box height
    gap = 0.10   # arrow gap from box edge

    # ── Rows ──
    y_top = 3.70
    y_bot = 1.30
    y_mid = 2.50

    # ── Columns for predicted path (6 boxes) ──
    # Spread across x=0.6 to x=7.2
    x1 = 0.70    # 2D Input Image
    x2 = 1.95    # VLM
    x3 = 3.20    # 3D Scene JSON
    x_sc = 4.45  # Scale Calibration
    x4 = 5.85    # CFD Pipeline (merged OF setup + CFD)
    x5 = 7.10    # Predicted Results

    # ── Right-side elements ──
    x_comp = 8.45   # Comparison box
    x_score = 9.65  # Score boxes

    sub_dy = 0.10

    # ══════════════════════════════════════════════════════════════════════
    #  PREDICTED PATHWAY  (top row)
    # ══════════════════════════════════════════════════════════════════════
    ax.text(x1 - bw / 2, y_top + bh / 2 + 0.20, "Predicted path",
            fontsize=10, fontweight="bold", color="#2F5C85",
            ha="left", va="bottom")

    _box(ax, x1, y_top, bw, bh, "2D Input\nImage",
         C_INPUT, C_INPUT_EDGE, bold=True, fontsize=9.5)

    _box(ax, x2, y_top, bw, bh, "VLM",
         C_PROCESS, C_PROCESS_EDGE, bold=True, fontsize=11)
    _subtitle(ax, x2, y_top - bh / 2 - sub_dy, "Gemini 2.0 Flash", fontsize=7.5)

    _box(ax, x3, y_top, bw, bh, "3D Scene\nJSON",
         C_OUTPUT, C_OUTPUT_EDGE, bold=True, fontsize=9.5)

    _box(ax, x_sc, y_top, bw, bh, "Scale\nCalibration",
         C_SCALE, C_SCALE_EDGE, fontsize=9.5)
    _subtitle(ax, x_sc, y_top - bh / 2 - sub_dy, "post-hoc scaling", fontsize=7.5)

    # Merged box for OF setup + CFD
    cfd_w = bw + 0.10
    _box(ax, x4, y_top, cfd_w, bh, "CFD\nPipeline",
         C_PROCESS, C_PROCESS_EDGE, fontsize=9.5)
    _subtitle(ax, x4, y_top - bh / 2 - sub_dy,
              "mesh + simulation", fontsize=7.5)

    _box(ax, x5, y_top, bw, bh, "Predicted\nResults",
         C_OUTPUT, C_OUTPUT_EDGE, bold=True, fontsize=9.5)

    # Horizontal arrows – predicted path
    pred_xs = [x1, x2, x3, x_sc, x4, x5]
    pred_ws = [bw, bw, bw, bw, cfd_w, bw]
    for i in range(len(pred_xs) - 1):
        _arrow(ax, pred_xs[i] + pred_ws[i] / 2 + gap, y_top,
               pred_xs[i + 1] - pred_ws[i + 1] / 2 - gap, y_top)

    # ══════════════════════════════════════════════════════════════════════
    #  REFERENCE PATHWAY  (bottom row)
    # ══════════════════════════════════════════════════════════════════════
    ax.text(x1 - bw / 2, y_bot + bh / 2 + 0.20, "Reference path",
            fontsize=10, fontweight="bold", color="#7D3C98",
            ha="left", va="bottom")

    ref_bw = bw + 0.15
    _box(ax, x1, y_bot, ref_bw, bh,
         "Reference\nBenchmark", C_INPUT, C_INPUT_EDGE, bold=True, fontsize=9.5)
    _subtitle(ax, x1, y_bot - bh / 2 - sub_dy, "rule-based", fontsize=7.5)

    _box(ax, x4, y_bot, cfd_w, bh, "CFD\nPipeline",
         C_PROCESS, C_PROCESS_EDGE, fontsize=9.5)
    _subtitle(ax, x4, y_bot - bh / 2 - sub_dy,
              "mesh + simulation", fontsize=7.5)

    _box(ax, x5, y_bot, bw, bh, "Reference\nResults",
         C_OUTPUT, C_OUTPUT_EDGE, bold=True, fontsize=9.5)

    # Long arrow: Reference Benchmark → CFD Pipeline
    _arrow(ax, x1 + ref_bw / 2 + gap, y_bot,
           x4 - cfd_w / 2 - gap, y_bot)
    # Short arrow: CFD Pipeline → Reference Results
    _arrow(ax, x4 + cfd_w / 2 + gap, y_bot,
           x5 - bw / 2 - gap, y_bot)

    # ══════════════════════════════════════════════════════════════════════
    #  COMPARISON  (between results and scores)
    # ══════════════════════════════════════════════════════════════════════
    comp_w = bw + 0.10
    comp_h = 1.20
    _box(ax, x_comp, y_mid, comp_w, comp_h,
         "Comparison", C_COMPARE, C_COMPARE_EDGE, bold=True, fontsize=10.5)

    # Diagonal arrows: predicted results → comparison
    _arrow(ax, x5 + bw / 2 + gap, y_top - bh * 0.15,
           x_comp - comp_w / 2 - gap, y_mid + comp_h * 0.22)
    # Diagonal arrows: reference results → comparison
    _arrow(ax, x5 + bw / 2 + gap, y_bot + bh * 0.15,
           x_comp - comp_w / 2 - gap, y_mid - comp_h * 0.22)

    # ══════════════════════════════════════════════════════════════════════
    #  SCORE OUTPUTS  (stacked right of comparison)
    # ══════════════════════════════════════════════════════════════════════
    score_bw = bw
    score_bh = bh - 0.05
    score_y_top = y_mid + 0.55
    score_y_bot = y_mid - 0.55

    _box(ax, x_score, score_y_top, score_bw, score_bh,
         "Structural\nScore", C_COMPARE, C_COMPARE_EDGE, fontsize=8.5, bold=True)

    _box(ax, x_score, score_y_bot, score_bw, score_bh,
         "CFD\nScore", C_COMPARE, C_COMPARE_EDGE, fontsize=8.5, bold=True)

    # Arrows: comparison → score boxes
    _arrow(ax, x_comp + comp_w / 2 + gap, y_mid + 0.15,
           x_score - score_bw / 2 - gap, score_y_top)
    _arrow(ax, x_comp + comp_w / 2 + gap, y_mid - 0.15,
           x_score - score_bw / 2 - gap, score_y_bot)

    # ── Subtle horizontal lane separator ──
    ax.plot([0.1, x5 + bw / 2 + 0.15], [y_mid, y_mid],
            ls=":", color="#D5D8DC", lw=0.7, zorder=0)

    # ══════════════════════════════════════════════════════════════════════
    #  Save
    # ══════════════════════════════════════════════════════════════════════
    fig.tight_layout(pad=0.3)
    fig.savefig(str(PDF_OUT), bbox_inches="tight", pad_inches=0.10)
    fig.savefig(str(PNG_OUT), dpi=600, bbox_inches="tight", pad_inches=0.10)
    plt.close(fig)
    print(f"Saved: {PDF_OUT}")
    print(f"Saved: {PNG_OUT}")


if __name__ == "__main__":
    main()
