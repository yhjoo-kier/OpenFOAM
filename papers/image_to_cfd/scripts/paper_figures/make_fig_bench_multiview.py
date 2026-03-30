#!/usr/bin/env python3
"""Fig 3: Multi-view rendering protocol — 5 views of one representative case."""
from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RENDER_ROOT = PROJECT_ROOT / "benchmark" / "renderings"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_bench_multiview.pdf"
PNG_OUT = OUT_DIR / "fig_bench_multiview.png"

CASE = "bench_a4_03"
VIEWS = [
    ("perspective", "(a) Perspective"),
    ("birdseye", "(b) Bird's eye"),
    ("floorplan", "(c) Floor plan"),
    ("wireframe", "(d) Wireframe"),
    ("section", "(e) Section"),
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

    fig, axes = plt.subplots(1, 5, figsize=(7.2, 2.2))

    for idx, (view, title) in enumerate(VIEWS):
        ax = axes[idx]
        img_path = RENDER_ROOT / CASE / view / f"{CASE}_{view}.png"
        if img_path.exists():
            img = mpimg.imread(str(img_path))
            ax.imshow(img, aspect="equal")
        ax.set_title(title, fontsize=8.5, pad=4)
        ax.axis("off")

    fig.suptitle(f"Case {CASE.replace('bench_', '').upper()} (composite dense)",
                 fontsize=10, fontweight="bold", y=0.04)

    fig.tight_layout(rect=[0.0, 0.10, 1.0, 1.0], w_pad=0.5)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight", pad_inches=0.05)
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
