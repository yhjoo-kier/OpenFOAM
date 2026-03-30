#!/usr/bin/env python3
"""Figure 2: Benchmark dataset design and difficulty matrix.

2x2 grid of representative floorplan renderings (one per architecture level)
with a summary statistics table below. Publication-quality, double-column width.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.lines as mlines
import numpy as np
from PIL import Image

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
AGG_PATH = (
    PROJECT_ROOT
    / "benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json"
)
EVAL_ROOT = (
    PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
)
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_bench_design.pdf"
PNG_OUT = OUT_DIR / "fig_bench_design.png"

# ---------------------------------------------------------------------------
# Representative cases (one per architecture level, floorplan view)
# ---------------------------------------------------------------------------
PANELS = [
    {
        "case": "bench_a1_01",
        "view": "floorplan",
        "arch": "A1",
        "label": "(a)",
        "desc": "Simple rectangular",
    },
    {
        "case": "bench_a2_03",
        "view": "floorplan",
        "arch": "A2",
        "label": "(b)",
        "desc": "Rectangular + obstacles",
    },
    {
        "case": "bench_a3_03",
        "view": "floorplan",
        "arch": "A3",
        "label": "(c)",
        "desc": "Composite (L-shaped)",
    },
    {
        "case": "bench_a4_03",
        "view": "floorplan",
        "arch": "A4",
        "label": "(d)",
        "desc": "Dense composite",
    },
]

# ---------------------------------------------------------------------------
# Style constants (matching other paper figures)
# ---------------------------------------------------------------------------
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

COLORS = {
    "panel": "#34495E",
    "grid": "#D9E2EC",
    "header_bg": "#2F5C85",
    "header_fg": "#FFFFFF",
    "row_even": "#F7F9FB",
    "row_odd": "#FFFFFF",
    "text": "#1E293B",
    "border": "#94A3B8",
    "accent": "#D97A2B",
}


def autocrop(img_path: str, border: int = 10) -> np.ndarray:
    """Load an image and crop whitespace borders, keeping a small margin."""
    pil_img = Image.open(img_path).convert("RGB")
    arr = np.array(pil_img)
    # Detect non-white rows/cols (threshold at 250 to handle near-white)
    gray = arr.mean(axis=2)
    rows = np.where(gray.min(axis=1) < 250)[0]
    cols = np.where(gray.min(axis=0) < 250)[0]
    if len(rows) == 0 or len(cols) == 0:
        return arr
    r0 = max(rows[0] - border, 0)
    r1 = min(rows[-1] + border + 1, arr.shape[0])
    c0 = max(cols[0] - border, 0)
    c1 = min(cols[-1] + border + 1, arr.shape[1])
    return arr[r0:r1, c0:c1]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


def main() -> None:
    # Load aggregate data
    payload = json.loads(AGG_PATH.read_text(encoding="utf-8"))
    by_cat = payload["by_category"]

    selected_font = pick_font()
    plt.rcParams.update(
        {
            "font.family": selected_font,
            "font.size": 9.0,
            "axes.titlesize": 10.0,
            "axes.labelsize": 9.0,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    # Figure: double-column width (~7.2 in), height to fit 2x2 images + table
    fig = plt.figure(figsize=(7.2, 7.0), dpi=150)

    # Top region: 2x2 image grid; Middle: legend; Bottom: summary table
    gs_main = gridspec.GridSpec(
        3, 1, figure=fig, height_ratios=[4.2, 0.25, 1.6],
        hspace=0.18, left=0.04, right=0.96, top=0.97, bottom=0.03,
    )

    gs_images = gridspec.GridSpecFromSubplotSpec(
        2, 2, subplot_spec=gs_main[0], wspace=0.12, hspace=0.20,
    )

    # ------------------------------------------------------------------
    # Draw the 2x2 image panels
    # ------------------------------------------------------------------
    for idx, panel in enumerate(PANELS):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs_images[row, col])

        img_path = EVAL_ROOT / panel["case"] / panel["view"] / "input.png"
        if img_path.exists():
            img = autocrop(str(img_path), border=8)
            ax.imshow(img, aspect="equal")
        else:
            ax.text(
                0.5, 0.5, "Image not found",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=9, color="#999",
            )

        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(COLORS["border"])
            spine.set_linewidth(0.8)

        # Panel label and description below image
        ax.set_title(
            f"{panel['label']}  {panel['arch']} \u2014 {panel['desc']}",
            fontsize=9.5,
            fontweight="bold",
            color=COLORS["panel"],
            pad=6,
        )

    # ------------------------------------------------------------------
    # Inlet / Outlet legend (between panels and table)
    # ------------------------------------------------------------------
    ax_legend = fig.add_subplot(gs_main[1])
    ax_legend.axis("off")
    inlet_handle = mlines.Line2D(
        [], [], color="#0000FF", linewidth=3, label="Inlet",
    )
    outlet_handle = mlines.Line2D(
        [], [], color="#FF0000", linewidth=3, label="Outlet",
    )
    ax_legend.legend(
        handles=[inlet_handle, outlet_handle],
        loc="center",
        ncol=2,
        fontsize=11,
        frameon=False,
        handlelength=2.0,
        columnspacing=3.0,
    )

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    ax_table = fig.add_subplot(gs_main[2])
    ax_table.axis("off")

    arch_levels = ["A1", "A2", "A3", "A4"]
    arch_descs = [
        "Simple rectangular",
        "Rectangular + obstacles",
        "Composite (L-shaped)",
        "Dense composite",
    ]

    # Build table data
    col_labels = [
        "Level",
        "Description",
        "Cases",
        "Views",
        "Total\nevals",
        "Struct.\nscore",
        "CFD\nscore",
        "Room\nmatch",
        "Opening\nmatch",
    ]

    cell_text = []
    for i, arch in enumerate(arch_levels):
        d = by_cat[arch]
        n_cases = d["n"] // 5  # 5 views per case
        cell_text.append([
            arch,
            arch_descs[i],
            str(n_cases),
            "5",
            str(d["n"]),
            f"{d['mean_structural_score']:.3f}",
            f"{d['mean_cfd_score']:.3f}",
            f"{d['room_kind_match_rate']:.0%}",
            f"{d['opening_wall_match_rate']:.0%}",
        ])

    # Totals row
    overall = payload["overall"]
    cell_text.append([
        "All",
        "Overall",
        "20",
        "5",
        str(overall["n"]),
        f"{overall['mean_structural_score']:.3f}",
        f"{overall['mean_cfd_score']:.3f}",
        f"{overall['room_kind_match_rate']:.0%}",
        f"{overall['opening_wall_match_rate']:.0%}",
    ])

    table = ax_table.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )

    table.auto_set_font_size(False)
    table.set_fontsize(8.0)
    table.scale(1.0, 1.45)

    # Style the table
    for (row_idx, col_idx), cell in table.get_celld().items():
        cell.set_edgecolor(COLORS["border"])
        cell.set_linewidth(0.5)

        if row_idx == 0:
            # Header row
            cell.set_facecolor(COLORS["header_bg"])
            cell.set_text_props(
                color=COLORS["header_fg"], fontweight="bold", fontsize=7.8,
            )
        elif row_idx == len(cell_text):
            # Totals row (last data row)
            cell.set_facecolor("#E8EDF2")
            cell.set_text_props(
                color=COLORS["text"], fontweight="bold", fontsize=7.8,
            )
        else:
            # Alternating row colors
            bg = COLORS["row_even"] if row_idx % 2 == 0 else COLORS["row_odd"]
            cell.set_facecolor(bg)
            cell.set_text_props(color=COLORS["text"], fontsize=7.8)

        # Left-align description column
        if col_idx == 1:
            cell.set_text_props(ha="left")
            cell.PAD = 0.05

    # Adjust column widths
    col_widths = [0.06, 0.22, 0.07, 0.06, 0.07, 0.1, 0.1, 0.1, 0.1]
    for col_idx, w in enumerate(col_widths):
        for row_idx in range(len(cell_text) + 1):
            table[(row_idx, col_idx)].set_width(w)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        str(PDF_OUT), format="pdf", bbox_inches="tight", pad_inches=0.05,
    )
    fig.savefig(
        str(PNG_OUT), format="png", dpi=600, bbox_inches="tight", pad_inches=0.05,
    )
    plt.close(fig)

    print(f"PDF -> {PDF_OUT}")
    print(f"PNG -> {PNG_OUT}")


if __name__ == "__main__":
    main()
