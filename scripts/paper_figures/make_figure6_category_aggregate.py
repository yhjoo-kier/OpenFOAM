#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

PROJECT_ROOT = Path(__file__).resolve().parents[2]
AGG_PATH = PROJECT_ROOT / "benchmark/manifests/evaluation_aggregate_summary.json"
OUT_DIR = PROJECT_ROOT / "results/paper_figures"
PDF_OUT = OUT_DIR / "figure6_category_aggregate_performance.pdf"
PNG_OUT = OUT_DIR / "figure6_category_aggregate_performance.png"

CATEGORY_ORDER = ["A1", "A2", "A3", "A4"]
CATEGORY_LABELS = {
    "A1": "A1  Rect. simple",
    "A2": "A2  Rect. dense",
    "A3": "A3  Comp. simple",
    "A4": "A4  Comp. dense",
}
COLORS = {
    "structural": "#2F5C85",
    "cfd": "#E08D2D",
    "room": "#3F7F6C",
    "opening": "#8D5A97",
    "grid": "#D9E2EC",
}


def main() -> None:
    payload = json.loads(AGG_PATH.read_text(encoding="utf-8"))
    by_category = payload["by_category"]

    structural = np.array([by_category[c]["mean_structural_score"] for c in CATEGORY_ORDER], dtype=float)
    cfd = np.array([by_category[c]["mean_cfd_score"] for c in CATEGORY_ORDER], dtype=float)
    room = np.array([by_category[c]["room_kind_match_rate"] for c in CATEGORY_ORDER], dtype=float)
    opening = np.array([by_category[c]["opening_wall_match_rate"] for c in CATEGORY_ORDER], dtype=float)
    labels = [CATEGORY_LABELS[c] for c in CATEGORY_ORDER]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.2,
            "axes.titlesize": 11.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 9.4,
            "ytick.labelsize": 9.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(6.95, 3.70), constrained_layout=False)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.16)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)

    y = np.arange(len(CATEGORY_ORDER))
    bar_h = 0.32

    ax0.barh(y - bar_h / 2, structural, height=bar_h, color=COLORS["structural"], zorder=3)
    ax0.barh(y + bar_h / 2, cfd, height=bar_h, color=COLORS["cfd"], zorder=3)

    ax1.barh(y - bar_h / 2, room, height=bar_h, color=COLORS["room"], zorder=3)
    ax1.barh(y + bar_h / 2, opening, height=bar_h, color=COLORS["opening"], zorder=3)

    for ax in (ax0, ax1):
        ax.invert_yaxis()
        ax.grid(axis="x", color=COLORS["grid"], linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.8)
        ax.spines["bottom"].set_linewidth(0.8)
        ax.spines["left"].set_color("#6B7280")
        ax.spines["bottom"].set_color("#6B7280")

    ax0.set_yticks(y, labels)
    ax1.tick_params(axis="y", which="both", left=False, labelleft=False)
    for yi in y:
        ax0.axhline(yi, color="#F3F6F9", linewidth=0.6, zorder=1)
        ax1.axhline(yi, color="#F3F6F9", linewidth=0.6, zorder=1)

    ax0.set_xlim(0, 0.94)
    ax1.set_xlim(0, 1.22)
    ax0.set_xlabel("Mean score")
    ax1.set_xlabel("Match rate (%)")
    ax1.xaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))
    ax0.set_title("Structural / CFD", loc="left", pad=6, fontweight="bold", fontsize=10.2)
    ax1.set_title("Room-kind / opening-wall", loc="left", pad=6, fontweight="bold", fontsize=10.2)
    ax0.text(-0.13, 1.03, "(a)", transform=ax0.transAxes, fontsize=10.0, fontweight="bold")
    ax1.text(-0.13, 1.03, "(b)", transform=ax1.transAxes, fontsize=10.0, fontweight="bold")

    for i, (s, c) in enumerate(zip(structural, cfd)):
        ax0.text(s + 0.014, i - bar_h / 2, f"{s:.2f}", va="center", ha="left", color="#334155", fontsize=8.0, fontweight="medium")
        ax0.text(c + 0.014, i + bar_h / 2, f"{c:.2f}", va="center", ha="left", color="#334155", fontsize=8.0, fontweight="medium")

    for i, (r, o) in enumerate(zip(room, opening)):
        ax1.text(r + 0.018, i - bar_h / 2, f"{r * 100:.0f}%", va="center", ha="left", color="#334155", fontsize=8.0, fontweight="medium")
        ax1.text(o + 0.018, i + bar_h / 2, f"{o * 100:.0f}%", va="center", ha="left", color="#334155", fontsize=8.0, fontweight="medium")

    fig.legend(
        [
            plt.Line2D([0], [0], color=COLORS["structural"], lw=7),
            plt.Line2D([0], [0], color=COLORS["cfd"], lw=7),
            plt.Line2D([0], [0], color=COLORS["room"], lw=7),
            plt.Line2D([0], [0], color=COLORS["opening"], lw=7),
        ],
        ["Structural", "CFD", "Room-kind", "Opening-wall"],
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.02),
        loc="lower center",
        columnspacing=1.3,
        handlelength=1.6,
        fontsize=8.5,
    )

    fig.subplots_adjust(left=0.24, right=0.985, top=0.80, bottom=0.26, wspace=0.15)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
