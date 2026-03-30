#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
AGG_PATH = PROJECT_ROOT / "benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json"
OUT_DIR = PROJECT_ROOT / "results/paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_result_view_aggregate.pdf"
PNG_OUT = OUT_DIR / "fig_result_view_aggregate.png"
META_OUT = OUT_DIR / "fig_result_view_aggregate_meta.json"

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
    "panel": "#34495E",
}

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.18,
        1.05,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.0,
        fontweight="bold",
        color=COLORS["panel"],
    )


def main() -> None:
    payload = json.loads(AGG_PATH.read_text(encoding="utf-8"))
    by_view = payload["by_view"]

    structural = np.array([by_view[v]["mean_structural_score"] for v in VIEW_ORDER], dtype=float)
    cfd = np.array([by_view[v]["mean_cfd_score"] for v in VIEW_ORDER], dtype=float)
    room = np.array([by_view[v]["room_kind_match_rate"] for v in VIEW_ORDER], dtype=float)
    opening = np.array([by_view[v]["opening_wall_match_rate"] for v in VIEW_ORDER], dtype=float)
    labels = [VIEW_LABELS[v] for v in VIEW_ORDER]

    selected_font = pick_font()
    plt.rcParams.update(
        {
            "font.family": selected_font,
            "font.size": 10.2,
            "axes.titlesize": 10.8,
            "axes.labelsize": 10.2,
            "xtick.labelsize": 9.6,
            "ytick.labelsize": 10.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(7.20, 3.95), constrained_layout=False)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.34)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])

    y = np.arange(len(VIEW_ORDER))
    bar_h = 0.32

    structural_bars = ax0.barh(y - bar_h / 2, structural, height=bar_h, color=COLORS["structural"], zorder=3, label="Structural")
    cfd_bars = ax0.barh(y + bar_h / 2, cfd, height=bar_h, color=COLORS["cfd"], zorder=3, label="CFD")

    room_bars = ax1.barh(y - bar_h / 2, room, height=bar_h, color=COLORS["room"], zorder=3, label="Room kind")
    opening_bars = ax1.barh(y + bar_h / 2, opening, height=bar_h, color=COLORS["opening"], zorder=3, label="Opening wall")

    for ax in (ax0, ax1):
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.margins(y=0.20)
        ax.grid(axis="x", color=COLORS["grid"], linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#6B7280")
        ax.spines["bottom"].set_color("#6B7280")

    ax0.set_xlim(0, 0.94)
    ax1.set_xlim(0, 1.18)
    ax0.set_xlabel("Mean score")
    ax1.set_xlabel("Match rate")
    ax0.set_title("Score metrics", pad=18)
    ax1.set_title("Agreement rates", pad=18)
    ax0.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        frameon=False,
        fontsize=9.4,
        handlelength=1.5,
        columnspacing=1.0,
        borderaxespad=0.0,
    )
    ax1.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        frameon=False,
        fontsize=9.4,
        handlelength=1.5,
        columnspacing=1.0,
        borderaxespad=0.0,
    )
    add_panel_label(ax0, "(a)")
    add_panel_label(ax1, "(b)")

    fig.subplots_adjust(left=0.15, right=0.985, top=0.82, bottom=0.16, wspace=0.38)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)

    meta = {
        "source_artifact": str(AGG_PATH),
        "setting": payload.get("setting"),
        "font_family_requested": FONT_CANDIDATES,
        "font_family_selected": selected_font,
        "intended_width": "double-column",
        "panel_layout": "1x2",
        "subfigure_labels": ["(a)", "(b)"],
        "png_dpi": 600,
        "pdf_vector": True,
        "internal_caption_text": False,
    }
    META_OUT.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"Selected font: {selected_font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {META_OUT}")


if __name__ == "__main__":
    main()
