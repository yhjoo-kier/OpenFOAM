#!/usr/bin/env python3
"""Fig 2: Benchmark dataset design — 2x2 category matrix with representative renders."""
from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_bench_design.pdf"
PNG_OUT = OUT_DIR / "fig_bench_design.png"

# Representative cases: one per category
PANELS = [
    {"case": "bench_a1_01", "view": "floorplan", "label": "(a)", "title": "A1: Rectangular, simple",
     "desc": "Single rect. room\n0 obstacles"},
    {"case": "bench_a2_03", "view": "floorplan", "label": "(b)", "title": "A2: Rectangular, dense",
     "desc": "Single rect. room\n3–5 obstacles"},
    {"case": "bench_a3_03", "view": "floorplan", "label": "(c)", "title": "A3: Composite, simple",
     "desc": "L-shaped (2 blocks)\n0–1 obstacle"},
    {"case": "bench_a4_03", "view": "floorplan", "label": "(d)", "title": "A4: Composite, dense",
     "desc": "L-shaped (2 blocks)\n3–5 obstacles"},
]

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0))

    for idx, panel in enumerate(PANELS):
        row, col = divmod(idx, 2)
        ax = axes[row, col]

        img_path = EVAL_ROOT / panel["case"] / panel["view"] / "input.png"
        if img_path.exists():
            img = mpimg.imread(str(img_path))
            ax.imshow(img, aspect="equal")
        else:
            ax.text(0.5, 0.5, "Image not found", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="red")

        ax.set_title(panel["title"], fontsize=10.5, fontweight="bold", pad=8)
        ax.axis("off")

        # Panel label
        ax.text(-0.05, 1.05, panel["label"], transform=ax.transAxes,
                fontsize=11, fontweight="bold", color="#34495E", va="bottom")

        # Description annotation (bottom-right)
        ax.text(0.97, 0.03, panel["desc"], transform=ax.transAxes,
                fontsize=8.5, ha="right", va="bottom", color="#334155",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#CBD5E1"))

    # Add row/column headers as text
    fig.text(0.28, 0.98, "Simple (0–1 obstacles)", ha="center", fontsize=10, color="#555555", fontstyle="italic")
    fig.text(0.74, 0.98, "Dense (3–5 obstacles)", ha="center", fontsize=10, color="#555555", fontstyle="italic")
    fig.text(0.02, 0.72, "Rectangular", ha="center", fontsize=10, color="#555555", fontstyle="italic", rotation=90)
    fig.text(0.02, 0.30, "Composite", ha="center", fontsize=10, color="#555555", fontstyle="italic", rotation=90)

    fig.tight_layout(rect=[0.04, 0.0, 1.0, 0.96], h_pad=1.5, w_pad=1.0)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
