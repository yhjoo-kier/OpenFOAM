#!/usr/bin/env python3
"""Fig 12: Scatter plot — structural score vs CFD agreement score (Phase 2)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_discuss_scatter_struct_cfd.pdf"
PNG_OUT = OUT_DIR / "fig_discuss_scatter_struct_cfd.png"

CATEGORY_COLORS = {
    "A1": "#1f77b4",
    "A2": "#ff7f0e",
    "A3": "#2ca02c",
    "A4": "#d62728",
}
VIEW_MARKERS = {
    "perspective": "o",
    "birdseye": "s",
    "floorplan": "D",
    "wireframe": "^",
    "section": "v",
}
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    rows = []
    for summary_path in sorted(EVAL_ROOT.glob("*/*/evaluation_summary.json")):
        d = json.loads(summary_path.read_text())
        task = d.get("task", {})
        pred = d.get("prediction_summary") or {}
        cfd_agg = (d.get("cfd_summary") or {}).get("aggregate_score") or {}
        case_name = task.get("case_name", "")
        category = case_name.split("_")[1].upper() if "_" in case_name else "?"
        struct = pred.get("structural_score")
        cfd = cfd_agg.get("cfd_agreement_score")
        if struct is not None and cfd is not None:
            rows.append({
                "category": category,
                "view": task.get("view", ""),
                "structural": struct,
                "cfd": cfd,
            })

    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 10,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(4.0, 3.8))

    # Add small jitter to reduce overplotting
    rng = np.random.default_rng(42)
    jitter_x = rng.uniform(-0.012, 0.012, len(rows))
    jitter_y = rng.uniform(-0.008, 0.008, len(rows))

    # Plot by category and view
    legend_handles_cat = {}
    legend_handles_view = {}
    for idx, row in enumerate(rows):
        cat = row["category"]
        view = row["view"]
        color = CATEGORY_COLORS.get(cat, "#888888")
        marker = VIEW_MARKERS.get(view, "o")
        ax.scatter(
            row["structural"] + jitter_x[idx],
            row["cfd"] + jitter_y[idx],
            c=color, marker=marker, s=44, alpha=0.70,
            edgecolors="#333333", linewidths=0.4, zorder=3,
        )
        if cat not in legend_handles_cat:
            legend_handles_cat[cat] = ax.scatter([], [], c=color, marker="o", s=48, label=cat, edgecolors="#333333", linewidths=0.4)
        if view not in legend_handles_view:
            legend_handles_view[view] = ax.scatter([], [], c="#888888", marker=marker, s=48, label=view.capitalize(), edgecolors="#333333", linewidths=0.4)

    # Correlation annotation
    s_arr = np.array([r["structural"] for r in rows])
    c_arr = np.array([r["cfd"] for r in rows])
    r_val = float(np.corrcoef(s_arr, c_arr)[0, 1])
    ax.text(0.97, 0.03, f"r = {r_val:.2f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9, color="#555555", fontstyle="italic")

    ax.set_xlabel("Structural score")
    ax.set_ylabel("CFD agreement score")
    ax.set_xlim(0.35, 1.06)
    ax.set_ylim(0.0, 0.85)
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Two legends: category (color) and view (marker)
    leg1 = ax.legend(
        handles=list(legend_handles_cat.values()),
        loc="upper left", fontsize=8.5, framealpha=0.8,
        title="Category", title_fontsize=8.5,
        handletextpad=0.3, borderpad=0.4,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=list(legend_handles_view.values()),
        loc="lower right", fontsize=7.5, framealpha=0.85,
        title="View", title_fontsize=7.5,
        handletextpad=0.3, borderpad=0.3, labelspacing=0.3,
    )

    fig.tight_layout(pad=0.5)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)
    print(f"Font: {font}")
    print(f"Points: {len(rows)}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
