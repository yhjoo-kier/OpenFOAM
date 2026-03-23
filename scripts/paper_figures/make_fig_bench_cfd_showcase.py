#!/usr/bin/env python3
"""New Figure: Benchmark CFD showcase — representative cases per category with 3D geometry + CFD."""
from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_bench_cfd_showcase.pdf"
PNG_OUT = OUT_DIR / "fig_bench_cfd_showcase.png"

# One reference case per category — show what the benchmark looks like
CASES = [
    {"dir": "phase2_ref_a1_01", "label": "A1: Rectangular, simple"},
    {"dir": "phase2_ref_a2_03", "label": "A2: Rectangular, dense"},
    {"dir": "phase2_ref_a3_03", "label": "A3: Composite, simple"},
    {"dir": "phase2_ref_a4_03", "label": "A4: Composite, dense"},
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
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for idx, case in enumerate(CASES):
        row, col = divmod(idx, 2)
        ax = axes[row, col]

        # Use individual panel PNGs side by side (no baked-in annotations)
        geo_path = RESULTS / case["dir"] / "panel_geometry_3d.png"
        flow_path = RESULTS / case["dir"] / "panel_flow_3d.png"
        if flow_path.exists() and geo_path.exists():
            import numpy as np
            from PIL import Image as PILImage
            geo_pil = PILImage.open(str(geo_path))
            flow_pil = PILImage.open(str(flow_path))
            # Resize to match heights
            target_h = min(geo_pil.height, flow_pil.height)
            geo_pil = geo_pil.resize((int(geo_pil.width * target_h / geo_pil.height), target_h))
            flow_pil = flow_pil.resize((int(flow_pil.width * target_h / flow_pil.height), target_h))
            combined = np.concatenate([np.array(geo_pil), np.array(flow_pil)], axis=1)
            ax.imshow(combined, interpolation="lanczos")
        else:
            ax.text(0.5, 0.5, "Render not found", transform=ax.transAxes,
                    ha="center", fontsize=10, color="red")

        ax.set_title(case["label"], fontsize=10, fontweight="bold", pad=8)
        ax.axis("off")

        ax.text(-0.02, 1.08, panel_labels[idx], transform=ax.transAxes,
                fontsize=11, fontweight="bold", color="#34495E")

    fig.tight_layout(h_pad=2.0, w_pad=0.8)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
