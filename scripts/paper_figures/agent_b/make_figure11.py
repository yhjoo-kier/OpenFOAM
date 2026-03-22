#!/usr/bin/env python3
"""Figure 11 – Robustness and convergence summary.

Three-panel layout:
  (a) Match rates (room kind, opening wall) by view type — grouped horizontal bars
  (b) Solver fallback counts by view type — stacked horizontal bars
  (c) Overall summary annotation box (table layout)

Double-column width, publication-quality output.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# ── paths ────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[3]
AGG_PATH = (
    PROJECT_ROOT
    / "benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json"
)
OUT_DIR = PROJECT_ROOT / "results/paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_result_robustness.pdf"
PNG_OUT = OUT_DIR / "fig_result_robustness.png"

# ── constants ────────────────────────────────────────────────────────────
VIEW_ORDER = ["perspective", "birdseye", "floorplan", "wireframe", "section"]
VIEW_LABELS = {
    "perspective": "Perspective",
    "birdseye": "Bird's-eye",
    "floorplan": "Floor plan",
    "wireframe": "Wireframe",
    "section": "Section",
}

COLORS = {
    "room_kind": "#2F5C85",
    "opening_wall": "#D97A2B",
    "ultra_robust": "#1E3A5F",   # dark blue — clearly distinct
    "laminar": "#E67E22",        # orange
    "conservative": "#2CA58D",   # teal
    "mesh025": "#E74C3C",        # coral/red
    "grid": "#D9E2EC",
    "panel": "#34495E",
    "annotation_bg": "#FFFFFF",          # clean white background
    "annotation_border": "#CED4DA",
}

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

FALLBACK_KEYS = [
    ("ultra_robust_count", "Ultra-robust"),
    ("laminar_fallback_count", "Laminar"),
    ("conservative_count", "Conservative"),
    ("mesh025_count", "Mesh 0.25"),
]
FALLBACK_COLORS = [
    COLORS["ultra_robust"],
    COLORS["laminar"],
    COLORS["conservative"],
    COLORS["mesh025"],
]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


def main() -> None:
    payload = json.loads(AGG_PATH.read_text(encoding="utf-8"))
    overall = payload["overall"]
    by_view = payload["by_view"]

    # ── extract arrays ───────────────────────────────────────────────────
    labels = [VIEW_LABELS[v] for v in VIEW_ORDER]
    n_views = len(VIEW_ORDER)
    y = np.arange(n_views)

    room_kind = np.array(
        [by_view[v]["room_kind_match_rate"] for v in VIEW_ORDER], dtype=float
    )
    opening_wall = np.array(
        [by_view[v]["opening_wall_match_rate"] for v in VIEW_ORDER], dtype=float
    )

    fallback_data = {}
    for key, _ in FALLBACK_KEYS:
        fallback_data[key] = np.array(
            [by_view[v][key] for v in VIEW_ORDER], dtype=float
        )

    # ── font / rcParams ──────────────────────────────────────────────────
    selected_font = pick_font()
    plt.rcParams.update(
        {
            "font.family": selected_font,
            "font.size": 10.0,
            "axes.titlesize": 11.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 10.0,
            "ytick.labelsize": 10.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    # ── layout ───────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(8.5, 4.2), constrained_layout=False)
    gs = fig.add_gridspec(
        1, 3, width_ratios=[1.0, 1.0, 1.0], wspace=0.55
    )
    ax_match = fig.add_subplot(gs[0, 0])
    ax_fallback = fig.add_subplot(gs[0, 1])
    ax_summary = fig.add_subplot(gs[0, 2])

    bar_h = 0.32

    # ── panel (a): match rates ───────────────────────────────────────────
    ax_match.barh(
        y - bar_h / 2, room_kind, height=bar_h,
        color=COLORS["room_kind"], zorder=3, label="Room kind",
    )
    ax_match.barh(
        y + bar_h / 2, opening_wall, height=bar_h,
        color=COLORS["opening_wall"], zorder=3, label="Opening wall",
    )

    for i in range(n_views):
        rk = room_kind[i]
        ax_match.text(
            rk + 0.03, y[i] - bar_h / 2, f"{rk:.0%}",
            va="center", ha="left", fontsize=10.5,
            color=COLORS["room_kind"], fontweight="bold",
        )
        ow = opening_wall[i]
        ax_match.text(
            ow + 0.03, y[i] + bar_h / 2, f"{ow:.0%}",
            va="center", ha="left", fontsize=10.5,
            color=COLORS["opening_wall"], fontweight="bold",
        )

    ax_match.set_yticks(y, labels)
    ax_match.invert_yaxis()
    ax_match.set_xlim(0, 1.35)
    ax_match.set_xlabel("Match rate")
    ax_match.set_title("Agreement rates", pad=8, fontsize=11)
    ax_match.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.14),
        ncol=2, frameon=False, fontsize=10.0,
        handlelength=1.5, handleheight=1.4, columnspacing=1.5,
    )

    # ── panel (b): fallback counts (stacked) ─────────────────────────────
    left = np.zeros(n_views)
    for idx, (key, lbl) in enumerate(FALLBACK_KEYS):
        vals = fallback_data[key]
        ax_fallback.barh(
            y, vals, left=left, height=0.52,
            color=FALLBACK_COLORS[idx], zorder=3, label=lbl,
        )
        left += vals

    totals = left
    for i in range(n_views):
        if totals[i] > 0:
            ax_fallback.text(
                totals[i] + 0.2, y[i], f"{int(totals[i])}/20",
                va="center", ha="left", fontsize=9.0, color="#555555",
            )

    ax_fallback.set_yticks(y, labels)
    ax_fallback.invert_yaxis()
    x_max_fb = max(totals.max() + 2, 7)
    ax_fallback.set_xlim(0, x_max_fb)
    ax_fallback.set_xlabel("Count (out of 20)")
    # Title centered, with enough pad so (b) label doesn't collide
    ax_fallback.set_title("Solver fallback triggers", pad=8, fontsize=11)
    ax_fallback.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.14),
        ncol=2, frameon=False, fontsize=10.0,
        handlelength=1.5, handleheight=1.4, columnspacing=1.5,
    )

    # ── panel (c): overall summary — two-column table ────────────────────
    ax_summary.set_xlim(0, 1)
    ax_summary.set_ylim(0, 1)
    ax_summary.axis("off")

    total_fallbacks = (
        overall["ultra_robust_count"]
        + overall["laminar_fallback_count"]
        + overall["conservative_count"]
        + overall["mesh025_count"]
    )
    nominal_count = overall["n"] - total_fallbacks

    # Short labels so they never overlap right-aligned values
    summary_lines = [
        ("Evaluations", f"{overall['n']}", False),
        ("Room kind", f"{overall['room_kind_match_rate']:.0%}", False),
        ("Opening wall", f"{overall['opening_wall_match_rate']:.0%}", False),
        (None, None, None),  # separator
        ("Nominal", f"{nominal_count}/{overall['n']}", False),
        ("Ultra-robust", f"{overall['ultra_robust_count']}", False),
        ("Laminar f/b", f"{overall['laminar_fallback_count']}", False),
        ("Conservative", f"{overall['conservative_count']}", False),
        ("Mesh 0.25", f"{overall['mesh025_count']}", False),
        (None, None, None),  # separator
        ("Geom. repair", f"{overall['used_repair_count']}", False),
        ("Sidecar warn.", f"{overall['nonblocking_repair_sidecar_warning_count']}", False),
    ]

    # Clean white background with subtle border
    rect = mpatches.FancyBboxPatch(
        (0.02, 0.02), 0.96, 0.96,
        boxstyle="round,pad=0.03",
        facecolor="#FFFFFF",
        edgecolor="#D0D0D0",
        linewidth=0.8,
        transform=ax_summary.transAxes, zorder=0,
    )
    ax_summary.add_patch(rect)

    ax_summary.text(
        0.50, 0.96, "Overall summary",
        transform=ax_summary.transAxes,
        ha="center", va="top",
        fontsize=11.5, fontweight="bold", color=COLORS["panel"],
    )

    def draw_hsep(yy: float) -> None:
        ax_summary.plot(
            [0.08, 0.92], [yy, yy],
            transform=ax_summary.transAxes,
            color=COLORS["annotation_border"], linewidth=0.5, zorder=1,
        )

    x_label = 0.08
    x_value = 0.92
    y_start = 0.86
    dy = 0.065
    row = 0
    for label, value, _ in summary_lines:
        if label is None:
            yy = y_start - (row - 0.2) * dy
            draw_hsep(yy)
            row += 0.5
            continue
        yy = y_start - row * dy
        ax_summary.text(
            x_label, yy, label,
            transform=ax_summary.transAxes,
            ha="left", va="center",
            fontsize=10.5, color="#444444",
        )
        ax_summary.text(
            x_value, yy, value,
            transform=ax_summary.transAxes,
            ha="right", va="center",
            fontsize=10.5, fontweight="bold", color=COLORS["panel"],
        )
        row += 1

    # ── panel labels — placed far left, above title area ─────────────────
    for ax, lbl in [(ax_match, "(a)"), (ax_fallback, "(b)"), (ax_summary, "(c)")]:
        x_pos = -0.06 if ax is ax_summary else -0.22
        ax.text(
            x_pos, 1.06, lbl,
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=12.0, fontweight="bold", color=COLORS["panel"],
        )

    # ── shared styling ───────────────────────────────────────────────────
    for ax in (ax_match, ax_fallback):
        ax.margins(y=0.18)
        ax.grid(axis="x", color=COLORS["grid"], linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#6B7280")
        ax.spines["bottom"].set_color("#6B7280")

    # ── save ─────────────────────────────────────────────────────────────
    fig.subplots_adjust(left=0.12, right=0.98, top=0.90, bottom=0.22, wspace=0.52)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)
    plt.close(fig)

    print(f"Font: {selected_font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
