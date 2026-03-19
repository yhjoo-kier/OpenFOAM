#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
AGG_PATH = PROJECT_ROOT / "benchmark/manifests/evaluation_aggregate_summary.json"
OUT_DIR = PROJECT_ROOT / "results/paper_figures"
PDF_OUT = OUT_DIR / "figure5_view_aggregate_performance.pdf"
PNG_OUT = OUT_DIR / "figure5_view_aggregate_performance.png"

VIEW_ORDER = ["perspective", "birdseye", "floorplan", "wireframe", "section"]
VIEW_LABELS = {
    "perspective": "Perspective",
    "birdseye": "Bird's eye",
    "floorplan": "Floor plan",
    "wireframe": "Wireframe",
    "section": "Section",
}

COLORS = {
    "structural": "#2F5C85",
    "cfd": "#D97A2B",
    "room": "#3B7A57",
    "opening": "#A23E48",
    "grid": "#D9E2EC",
}


def main() -> None:
    payload = json.loads(AGG_PATH.read_text(encoding="utf-8"))
    by_view = payload["by_view"]

    structural = np.array([by_view[v]["mean_structural_score"] for v in VIEW_ORDER], dtype=float)
    cfd = np.array([by_view[v]["mean_cfd_score"] for v in VIEW_ORDER], dtype=float)
    room = np.array([by_view[v]["room_kind_match_rate"] for v in VIEW_ORDER], dtype=float)
    opening = np.array([by_view[v]["opening_wall_match_rate"] for v in VIEW_ORDER], dtype=float)
    labels = [VIEW_LABELS[v] for v in VIEW_ORDER]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.8,
            "axes.titlesize": 12.2,
            "axes.labelsize": 10.8,
            "xtick.labelsize": 10.2,
            "ytick.labelsize": 10.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(7.65, 4.70), constrained_layout=False)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1.0], wspace=0.34)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])

    y = np.arange(len(VIEW_ORDER))
    bar_h = 0.34

    ax0.barh(y - bar_h / 2, structural, height=bar_h, color=COLORS["structural"], zorder=3)
    ax0.barh(y + bar_h / 2, cfd, height=bar_h, color=COLORS["cfd"], zorder=3)

    ax1.barh(y - bar_h / 2, room, height=bar_h, color=COLORS["room"], zorder=3)
    ax1.barh(y + bar_h / 2, opening, height=bar_h, color=COLORS["opening"], zorder=3)

    for ax in (ax0, ax1):
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.grid(axis="x", color=COLORS["grid"], linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax0.set_xlim(0, 0.90)
    ax1.set_xlim(0, 1.12)
    ax0.set_xlabel("Mean score")
    ax1.set_xlabel("Match rate")
    ax0.set_title("Performance fidelity\n(blue: structural, orange: CFD)", pad=10)
    ax1.set_title("Topology / opening agreement\n(green: room-kind, red: opening-wall)", pad=10)

    for i, (s, c) in enumerate(zip(structural, cfd)):
        ax0.text(s + 0.010, i - bar_h / 2, f"{s:.3f}", va="center", ha="left", color=COLORS["structural"], fontsize=9.2, fontweight="medium")
        ax0.text(c + 0.010, i + bar_h / 2, f"{c:.3f}", va="center", ha="left", color=COLORS["cfd"], fontsize=9.2, fontweight="medium")

    for i, (r, o) in enumerate(zip(room, opening)):
        ax1.text(r + 0.015, i - bar_h / 2, f"{r * 100:.0f}%", va="center", ha="left", color=COLORS["room"], fontsize=9.2, fontweight="medium")
        ax1.text(max(o - 0.020, 0.06), i + bar_h / 2, f"{o * 100:.0f}%", va="center", ha="right", color="white", fontsize=9.2, fontweight="bold")

    fig.suptitle(
        "Figure 5. Aggregate performance across input views (frozen-20 benchmark, 100 tasks)",
        x=0.5,
        y=0.955,
        fontsize=11.6,
        fontweight="bold",
    )

    fig.subplots_adjust(left=0.21, right=0.985, top=0.78, bottom=0.12, wspace=0.34)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=300, bbox_inches="tight")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
