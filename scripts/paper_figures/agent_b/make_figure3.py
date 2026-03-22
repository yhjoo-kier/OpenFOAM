#!/usr/bin/env python3
"""Figure 3 -- Multi-view rendering protocol.

Shows a single benchmark case (bench_a4_03) rendered into 5 different
2D view types in a 2-row layout (3 top, 2 bottom centred).
Publication-quality, double-column width, 600 dpi.

Note: Switched from bench_a2_03 to bench_a4_03 because the former had
a nearly-blank section view that failed external QC.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.image as mpimg
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CASE = "bench_a4_03"
EVAL_DIR = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span" / CASE
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_bench_multiview.pdf"
PNG_OUT = OUT_DIR / "fig_bench_multiview.png"

# Top row: perspective, birdseye, floorplan
# Bottom row: wireframe, section (centred)
TOP_VIEWS = ["perspective", "birdseye", "floorplan"]
BOT_VIEWS = ["wireframe", "section"]

PANEL_LABELS = {
    "perspective": "(a) Perspective",
    "birdseye": "(b) Bird\u2019s-eye",
    "floorplan": "(c) Floor plan",
    "wireframe": "(d) Wireframe",
    "section": "(e) Section",
}

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]
PANEL_COLOR = "#34495E"
BORDER_COLOR = "#9E9E9E"


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load images
    images: dict[str, np.ndarray] = {}
    for view in TOP_VIEWS + BOT_VIEWS:
        img_path = EVAL_DIR / view / "input.png"
        if not img_path.exists():
            raise FileNotFoundError(f"Missing input image: {img_path}")
        images[view] = mpimg.imread(str(img_path))

    selected_font = pick_font()
    plt.rcParams.update({
        "font.family": selected_font,
        "font.size": 10.0,
        "axes.titlesize": 12.0,
        "axes.labelsize": 10.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    # Layout: 2 rows x 6 columns grid
    # Top row spans columns 0-1, 2-3, 4-5 (3 panels)
    # Bottom row spans columns 1-2, 3-4 (2 panels, centred)
    fig_width = 7.2
    fig_height = 5.2  # tall enough for 2 rows of images + labels
    fig = plt.figure(figsize=(fig_width, fig_height))

    gs = GridSpec(
        2, 6,
        figure=fig,
        wspace=0.08,
        hspace=0.22,
        left=0.01,
        right=0.99,
        bottom=0.06,
        top=0.95,
    )

    # Top row: 3 panels each spanning 2 columns
    top_axes = [
        fig.add_subplot(gs[0, 0:2]),
        fig.add_subplot(gs[0, 2:4]),
        fig.add_subplot(gs[0, 4:6]),
    ]
    # Bottom row: 2 panels centred (span cols 1-2 and 3-4)
    bot_axes = [
        fig.add_subplot(gs[1, 1:3]),
        fig.add_subplot(gs[1, 3:5]),
    ]

    all_pairs = list(zip(top_axes, TOP_VIEWS)) + list(zip(bot_axes, BOT_VIEWS))

    for ax, view in all_pairs:
        img = images[view]
        ax.imshow(img, aspect="equal", interpolation="lanczos")
        ax.set_axis_off()

        # Panel label above
        ax.set_title(
            PANEL_LABELS[view],
            fontsize=12,
            fontweight="bold",
            color=PANEL_COLOR,
            pad=6,
        )

        # Simple thin rectangular border
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.6)
            spine.set_edgecolor(BORDER_COLOR)

    # Inlet / Outlet legend at bottom
    inlet_handle = mlines.Line2D(
        [], [], color="#0000FF", linewidth=3, label="Inlet",
    )
    outlet_handle = mlines.Line2D(
        [], [], color="#FF0000", linewidth=3, label="Outlet",
    )
    fig.legend(
        handles=[inlet_handle, outlet_handle],
        loc="lower center",
        ncol=2,
        fontsize=11,
        frameon=False,
        handlelength=2.0,
        columnspacing=3.0,
        bbox_to_anchor=(0.5, 0.0),
    )

    fig.savefig(str(PDF_OUT), dpi=600, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(str(PNG_OUT), dpi=600, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    print(f"Saved PDF: {PDF_OUT}")
    print(f"Saved PNG: {PNG_OUT}")
    print(f"PDF size : {PDF_OUT.stat().st_size / 1024:.1f} KB")
    print(f"PNG size : {PNG_OUT.stat().st_size / 1024:.1f} KB")


if __name__ == "__main__":
    main()
