#!/usr/bin/env python3
"""Fig NEW: CFD agreement score component breakdown by view type (Phase 2)."""
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
PDF_OUT = OUT_DIR / "fig_result_cfd_component_breakdown.pdf"
PNG_OUT = OUT_DIR / "fig_result_cfd_component_breakdown.png"

VIEWS = ["perspective", "birdseye", "floorplan", "wireframe", "section"]
VIEW_LABELS = ["Perspective", "Bird's eye", "Floor plan", "Wireframe", "Section"]
COMPONENTS = ["overlap_ratio_vs_union", "velocity_magnitude_similarity",
              "velocity_direction_similarity", "pressure_similarity"]
COMP_LABELS = ["Overlap", "Vel. magnitude", "Vel. direction", "Pressure"]
COMP_COLORS = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759"]

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    # Collect per-view component means
    by_view: dict[str, list[dict]] = defaultdict(list)
    for sp in sorted(EVAL_ROOT.glob("*/*/evaluation_summary.json")):
        d = json.loads(sp.read_text())
        task = d.get("task", {})
        cfd_agg = (d.get("cfd_summary") or {}).get("aggregate_score") or {}
        comps = cfd_agg.get("components") or {}
        view = task.get("view", "")
        if view in VIEWS and comps:
            by_view[view].append(comps)

    # Build component matrices: views x components
    means = np.zeros((len(VIEWS), len(COMPONENTS)))
    sds = np.zeros((len(VIEWS), len(COMPONENTS)))
    for i, view in enumerate(VIEWS):
        for j, comp in enumerate(COMPONENTS):
            vals = [c.get(comp) for c in by_view[view] if c.get(comp) is not None]
            if vals:
                means[i, j] = np.mean(vals)
                sds[i, j] = np.std(vals, ddof=1) if len(vals) > 1 else 0.0

    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 10,
        "axes.labelsize": 10.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(7.0, 3.8))

    x = np.arange(len(VIEWS))
    n_comp = len(COMPONENTS)
    bar_w = 0.18
    offsets = np.linspace(-(n_comp - 1) * bar_w / 2, (n_comp - 1) * bar_w / 2, n_comp)

    for j, (comp, label, color) in enumerate(zip(COMPONENTS, COMP_LABELS, COMP_COLORS)):
        ax.bar(x + offsets[j], means[:, j], width=bar_w, color=color,
               label=label, zorder=3,
               yerr=sds[:, j], error_kw=dict(ecolor="#4B5563", elinewidth=0.7, capsize=2, capthick=0.7))

    ax.set_xticks(x)
    ax.set_xticklabels(VIEW_LABELS, fontsize=9.5)
    ax.set_ylabel("Component score")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.3, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=4,
              frameon=False, fontsize=9, handlelength=1.2, columnspacing=1.0)

    # Add overall mean line
    overall_mean = float(np.mean(means))
    ax.axhline(overall_mean, color="#666666", linestyle="--", linewidth=1.0, zorder=2, alpha=0.7)
    ax.text(0.97, 0.97, f"Overall mean: {overall_mean:.3f}",
            transform=ax.transAxes, fontsize=8.5, color="#444444",
            ha="right", va="top", fontstyle="italic")

    fig.tight_layout(pad=1.0)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
