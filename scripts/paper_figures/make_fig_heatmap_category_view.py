#!/usr/bin/env python3
"""Fig 13: Category x View interaction heatmaps (Phase 2)."""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_result_heatmap_category_view.pdf"
PNG_OUT = OUT_DIR / "fig_result_heatmap_category_view.png"

CATEGORIES = ["A1", "A2", "A3", "A4"]
VIEWS = ["perspective", "birdseye", "floorplan", "wireframe", "section"]
VIEW_LABELS = ["Persp.", "Bird's eye", "Floor plan", "Wireframe", "Section"]
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    # Collect per-cell data
    cells: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for sp in sorted(EVAL_ROOT.glob("*/*/evaluation_summary.json")):
        d = json.loads(sp.read_text())
        task = d.get("task", {})
        pred = d.get("prediction_summary") or {}
        cfd_agg = (d.get("cfd_summary") or {}).get("aggregate_score") or {}
        case = task.get("case_name", "")
        cat = case.split("_")[1].upper() if "_" in case else "?"
        view = task.get("view", "")
        struct = pred.get("structural_score")
        cfd = cfd_agg.get("cfd_agreement_score")
        if cat in CATEGORIES and view in VIEWS:
            cells[(cat, view)].append({"structural": struct, "cfd": cfd})

    # Build matrices
    struct_mat = np.full((len(CATEGORIES), len(VIEWS)), np.nan)
    cfd_mat = np.full((len(CATEGORIES), len(VIEWS)), np.nan)
    count_mat = np.zeros((len(CATEGORIES), len(VIEWS)), dtype=int)

    for i, cat in enumerate(CATEGORIES):
        for j, view in enumerate(VIEWS):
            vals = cells[(cat, view)]
            count_mat[i, j] = len(vals)
            s_vals = [v["structural"] for v in vals if v["structural"] is not None]
            c_vals = [v["cfd"] for v in vals if v["cfd"] is not None]
            if s_vals:
                struct_mat[i, j] = np.mean(s_vals)
            if c_vals:
                cfd_mat[i, j] = np.mean(c_vals)

    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2))

    for ax, mat, title, cmap, vrange in [
        (ax1, struct_mat, "Structural score", "Blues", (0.4, 1.0)),
        (ax2, cfd_mat, "CFD agreement score", "Oranges", (0.2, 0.7)),
    ]:
        im = ax.imshow(mat, cmap=cmap, vmin=vrange[0], vmax=vrange[1], aspect="auto")
        ax.set_xticks(range(len(VIEWS)))
        ax.set_xticklabels(VIEW_LABELS, rotation=35, ha="right", fontsize=8.5)
        ax.set_yticks(range(len(CATEGORIES)))
        ax.set_yticklabels(CATEGORIES, fontsize=9.5)
        ax.set_title(title, fontsize=10.5, pad=8)

        # Annotate cells
        for i in range(len(CATEGORIES)):
            for j in range(len(VIEWS)):
                val = mat[i, j]
                if not np.isnan(val):
                    color = "white" if val > (vrange[0] + vrange[1]) / 2 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=8.5, color=color, fontweight="bold")

        cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03)
        cbar.ax.tick_params(labelsize=8)

    ax1.text(-0.22, 1.05, "(a)", transform=ax1.transAxes, fontsize=11, fontweight="bold", color="#34495E")
    ax2.text(-0.15, 1.05, "(b)", transform=ax2.transAxes, fontsize=11, fontweight="bold", color="#34495E")

    fig.tight_layout(pad=1.0)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
