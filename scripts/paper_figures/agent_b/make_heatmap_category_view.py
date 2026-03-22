#!/usr/bin/env python3
"""Category x View interaction heatmap (structural & CFD scores)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
STATS_PATH = PROJECT_ROOT / "benchmark/manifests/evaluation_statistics.json"
OUT_DIR = PROJECT_ROOT / "results/paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_result_heatmap_category_view.pdf"
PNG_OUT = OUT_DIR / "fig_result_heatmap_category_view.png"

CATEGORIES = ["A1", "A2", "A3", "A4"]
VIEW_ORDER = ["perspective", "birdseye", "floorplan", "wireframe", "section"]
VIEW_LABELS = ["Persp.", "Bird's-eye", "Floor plan", "Wireframe", "Section"]

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]
PANEL_COLOR = "#34495E"


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


def load_cross_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (structural_mean, cfd_mean, counts) each shape (4, 5)."""
    with open(STATS_PATH) as f:
        data = json.load(f)

    records = data["records"]

    cat_idx = {c: i for i, c in enumerate(CATEGORIES)}
    view_idx = {v: i for i, v in enumerate(VIEW_ORDER)}

    structural_acc = np.zeros((4, 5))
    cfd_acc = np.zeros((4, 5))
    counts = np.zeros((4, 5), dtype=int)

    for rec in records:
        ci = cat_idx.get(rec["category"])
        vi = view_idx.get(rec["view"])
        if ci is None or vi is None:
            continue
        structural_acc[ci, vi] += rec["structural_score"]
        cfd_acc[ci, vi] += rec["cfd_score"]
        counts[ci, vi] += 1

    with np.errstate(divide="ignore", invalid="ignore"):
        structural_mean = np.where(counts > 0, structural_acc / counts, np.nan)
        cfd_mean = np.where(counts > 0, cfd_acc / counts, np.nan)

    return structural_mean, cfd_mean, counts


def make_figure() -> None:
    font_family = pick_font()
    plt.rcParams.update({
        "font.family": font_family,
        "font.size": 11,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
    })

    structural, cfd, counts = load_cross_table()

    fig, (ax_s, ax_c) = plt.subplots(1, 2, figsize=(7.2, 3.5))

    cmap = plt.cm.RdYlGn
    vmin, vmax = 0.3, 1.0

    for ax, matrix, title, label in [
        (ax_s, structural, "Structural score", "(a)"),
        (ax_c, cfd, "CFD score", "(b)"),
    ]:
        im = ax.imshow(
            matrix,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
            interpolation="nearest",
        )

        # Annotate cells
        for i in range(4):
            for j in range(5):
                val = matrix[i, j]
                n = counts[i, j]
                if np.isnan(val):
                    txt = "n/a"
                    color = "gray"
                else:
                    txt = f"{val:.2f}"
                    # dark text on light cells, white on dark
                    color = "white" if val < 0.50 else "black"
                suffix = f"\n(n={n})" if n < 5 else ""
                ax.text(
                    j, i, txt + suffix,
                    ha="center", va="center",
                    fontsize=10, fontweight="medium",
                    color=color,
                )

        ax.set_xticks(range(5))
        ax.set_xticklabels(VIEW_LABELS, rotation=30, ha="right")
        ax.set_yticks(range(4))
        ax.set_yticklabels(CATEGORIES)
        ax.set_title(title, pad=8)

        # Panel label
        ax.text(
            -0.22, 1.08, label,
            transform=ax.transAxes,
            ha="left", va="bottom",
            fontsize=12, fontweight="bold",
            color=PANEL_COLOR,
        )

    # Shared colorbar
    fig.subplots_adjust(right=0.88, wspace=0.35, bottom=0.22, top=0.88)
    cbar_ax = fig.add_axes([0.90, 0.22, 0.02, 0.66])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Score", fontsize=11)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, dpi=600, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {PDF_OUT}")
    print(f"Saved: {PNG_OUT}")

    # Print cross-table for verification
    print("\n--- Structural score (category x view) ---")
    print(f"{'':>6s}", *[f"{v:>12s}" for v in VIEW_LABELS])
    for i, cat in enumerate(CATEGORIES):
        vals = "".join(f"{structural[i,j]:12.3f}" for j in range(5))
        print(f"{cat:>6s}{vals}")

    print("\n--- CFD score (category x view) ---")
    print(f"{'':>6s}", *[f"{v:>12s}" for v in VIEW_LABELS])
    for i, cat in enumerate(CATEGORIES):
        vals = "".join(f"{cfd[i,j]:12.3f}" for j in range(5))
        print(f"{cat:>6s}{vals}")

    print(f"\n--- Sample counts ---")
    print(f"{'':>6s}", *[f"{v:>12s}" for v in VIEW_LABELS])
    for i, cat in enumerate(CATEGORIES):
        vals = "".join(f"{counts[i,j]:12d}" for j in range(5))
        print(f"{cat:>6s}{vals}")


if __name__ == "__main__":
    make_figure()
