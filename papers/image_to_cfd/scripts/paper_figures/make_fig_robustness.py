#!/usr/bin/env python3
"""Fig 11: Robustness and convergence summary (Phase 2)."""
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
PDF_OUT = OUT_DIR / "fig_result_robustness.pdf"
PNG_OUT = OUT_DIR / "fig_result_robustness.png"

VIEWS = ["perspective", "birdseye", "floorplan", "wireframe", "section"]
VIEW_LABELS = ["Perspective", "Bird's eye", "Floor plan", "Wireframe", "Section"]

# Solver divergence cases (excluded from converged stats)
SOLVER_FAILED = {
    ("bench_a3_04", "perspective"),
    ("bench_a4_02", "perspective"),
    ("bench_a4_02", "wireframe"),
}

PRESET_COLORS = {
    "nominal": "#4e79a7",
    "robust": "#59a14f",
    "conservative": "#f28e2b",
    "ultra_robust": "#e15759",
    "laminar_fallback": "#b07aa1",
    "failed": "#333333",
}
PRESET_ORDER = ["nominal", "robust", "conservative", "ultra_robust", "laminar_fallback", "failed"]
PRESET_LABELS = {
    "nominal": "Nominal",
    "robust": "Robust",
    "conservative": "Conservative",
    "ultra_robust": "Ultra-robust",
    "laminar_fallback": "Laminar fallback",
    "failed": "Solver divergence",
}

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    # Collect per-task robustness data
    by_view: dict[str, dict[str, int]] = {v: defaultdict(int) for v in VIEWS}
    total_by_preset: dict[str, int] = defaultdict(int)

    for sp in sorted(EVAL_ROOT.glob("*/*/evaluation_summary.json")):
        d = json.loads(sp.read_text())
        task = d.get("task", {})
        pipeline = d.get("pipeline_summary") or {}
        case_name = task.get("case_name", "")
        view = task.get("view", "")

        if (case_name, view) in SOLVER_FAILED:
            preset = "failed"
        else:
            preset = pipeline.get("successful_preset", "nominal") or "nominal"
            # Normalize preset names
            if preset not in PRESET_ORDER:
                preset = "nominal"

        if view in VIEWS:
            by_view[view][preset] += 1
            total_by_preset[preset] += 1

    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 10,
        "axes.labelsize": 10.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.5),
                                     gridspec_kw={"width_ratios": [1.6, 1.0]})

    # Panel (a): Stacked bar by view
    x = np.arange(len(VIEWS))
    bottoms = np.zeros(len(VIEWS))
    for preset in PRESET_ORDER:
        vals = np.array([by_view[v][preset] for v in VIEWS], dtype=float)
        if vals.sum() == 0:
            continue
        ax1.bar(x, vals, bottom=bottoms, width=0.6,
                color=PRESET_COLORS[preset], label=PRESET_LABELS[preset], zorder=3)
        # Annotate non-zero segments
        for i, val in enumerate(vals):
            if val > 0:
                ax1.text(x[i], bottoms[i] + val / 2, str(int(val)),
                        ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        bottoms += vals

    ax1.set_xticks(x)
    ax1.set_xticklabels(VIEW_LABELS, fontsize=8.5, rotation=30, ha="right")
    ax1.set_ylabel("Number of cases")
    ax1.set_ylim(0, 22)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(axis="y", alpha=0.3, linewidth=0.6, zorder=0)
    ax1.set_axisbelow(True)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=3,
               frameon=False, fontsize=8, handlelength=1.0, columnspacing=0.8)

    # Panel (b): Overall horizontal bar summary (cleaner than donut for small segments)
    presets_used = [(p, total_by_preset[p]) for p in PRESET_ORDER if total_by_preset[p] > 0]
    bar_labels = [PRESET_LABELS[p] for p, _ in presets_used]
    bar_vals = [c for _, c in presets_used]
    bar_colors = [PRESET_COLORS[p] for p, _ in presets_used]

    y2 = np.arange(len(presets_used))
    ax2.barh(y2, bar_vals, height=0.55, color=bar_colors, zorder=3, edgecolor="white", linewidth=0.8)
    for i, val in enumerate(bar_vals):
        ax2.text(val + 0.8, y2[i], str(val), ha="left", va="center", fontsize=9, fontweight="bold")
    ax2.set_yticks(y2)
    ax2.set_yticklabels(bar_labels, fontsize=8.5)
    ax2.invert_yaxis()
    ax2.set_xlabel("Cases", fontsize=9.5)
    ax2.set_xlim(0, 95)
    ax2.set_title("Overall (n = 100)", fontsize=10, pad=8)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="x", alpha=0.3, linewidth=0.6, zorder=0)
    ax2.set_axisbelow(True)

    # Panel labels
    ax1.text(-0.12, 1.08, "(a)", transform=ax1.transAxes, fontsize=11,
             fontweight="bold", color="#34495E")
    ax2.text(-0.15, 1.05, "(b)", transform=ax2.transAxes, fontsize=11,
             fontweight="bold", color="#34495E")

    fig.tight_layout(pad=1.0)
    fig.subplots_adjust(bottom=0.20)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)
    print(f"Font: {font}")
    print(f"Presets: {dict(total_by_preset)}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
