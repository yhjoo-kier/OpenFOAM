#!/usr/bin/env python3
"""Plot summary heatmaps for paint-booth panel-frame sweep."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def grid(rows, key):
    supplies = sorted({float(r["supply_velocity_mps"]) for r in rows})
    fs = sorted({float(r["filter_forchheimer"]) for r in rows})
    z = np.full((len(fs), len(supplies)), np.nan)
    for r in rows:
        i = fs.index(float(r["filter_forchheimer"]))
        j = supplies.index(float(r["supply_velocity_mps"]))
        z[i, j] = float(r[key])
    return supplies, fs, z


def annotate(ax, data, fmt=".2f"):
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if np.isfinite(data[i, j]):
                ax.text(j, i, format(data[i, j], fmt), ha="center", va="center", color="white", fontsize=9, fontweight="bold")


def main():
    root = Path("/home/yhjoo/projects/OpenFOAM")
    rows = json.loads((root / "cases/paint_booth_panel_frame_sweep/sweep_summary.json").read_text())
    figdir = root / "docs/figures"
    figdir.mkdir(parents=True, exist_ok=True)

    panels = [
        ("filter_dp_pa_rho1p2", "Filter Δp proxy [Pa], ρ=1.2", ".0f", "magma"),
        ("filter_below_Uz_mean", "Filter-below |mean Uz| [m/s]", ".2f", "viridis"),
        ("work_zone_Uz_mean", "Work-zone |mean Uz| [m/s]", ".2f", "viridis"),
        ("near_car_Uz_mean", "Near-car |mean Uz| [m/s]", ".2f", "viridis"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for ax, (key, title, fmt, cmap) in zip(axes.ravel(), panels):
        supplies, fs, data = grid(rows, key)
        if "Uz" in key:
            data = np.abs(data)
        im = ax.imshow(data, origin="lower", aspect="auto", cmap=cmap)
        annotate(ax, data, fmt)
        ax.set_title(title)
        ax.set_xticks(range(len(supplies)), [f"{s:.2f}" for s in supplies])
        ax.set_yticks(range(len(fs)), [f"{f:.0f}" for f in fs])
        ax.set_xlabel("Supply velocity [m/s]")
        ax.set_ylabel("Filter Forchheimer f")
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle("Paint-booth panel-frame sweep summary", fontsize=16)
    out = figdir / "26-05-02_paint_booth_panel_frame_sweep_heatmaps.png"
    fig.savefig(out, dpi=180)
    print(out)


if __name__ == "__main__":
    main()
