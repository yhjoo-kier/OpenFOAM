#!/usr/bin/env python3
"""Fig: Nielsen validation — CFD vs Experimental scatter plot."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CASE_DIR = PROJECT_ROOT / "cases" / "validation_nielsen_1990"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_method_validation_scatter.pdf"
PNG_OUT = OUT_DIR / "fig_method_validation_scatter.png"

L, H, D, U0 = 9.0, 3.0, 0.1, 0.455
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

PROFILE_COLORS = {"x_over_L_0.33": "#2563EB", "x_over_L_0.50": "#16A34A", "x_over_L_0.67": "#DC2626"}
PROFILE_LABELS = {"x_over_L_0.33": "$x/L$ = 1/3", "x_over_L_0.50": "$x/L$ = 1/2", "x_over_L_0.67": "$x/L$ = 2/3"}


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font, "font.size": 10,
        "axes.labelsize": 10.5, "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    exp = json.load(open(CASE_DIR / "experimental_data.json"))

    # Load VTK
    vtk_path = None
    for d in sorted((CASE_DIR / "VTK").iterdir()):
        internal = d / "internal.vtu"
        if internal.exists():
            vtk_path = internal
    mesh = pv.read(vtk_path)
    if "U" in mesh.cell_data and "U" not in mesh.point_data:
        mesh = mesh.cell_data_to_point_data()

    fig, ax = plt.subplots(figsize=(4.0, 4.0))

    all_exp, all_cfd = [], []
    for pname, data in exp.items():
        y_H = np.array(data["y_over_H"])
        u_exp = np.array(data["u_over_U0"])
        points = np.column_stack([np.full_like(y_H, float(pname.split("_")[-1]) * L),
                                   np.full_like(y_H, D / 2), y_H * H])
        sampled = pv.PolyData(points).sample(mesh)
        u_cfd = np.array(sampled.point_data["U"])[:, 0] / U0

        ax.scatter(u_exp, u_cfd, c=PROFILE_COLORS[pname], s=40, alpha=0.8,
                   edgecolors="#333333", linewidths=0.4, label=PROFILE_LABELS[pname], zorder=3)
        all_exp.extend(u_exp.tolist())
        all_cfd.extend(u_cfd.tolist())

    # 1:1 line
    lim = [-0.25, 0.45]
    ax.plot(lim, lim, "k--", linewidth=1.0, alpha=0.5, label="1:1 line", zorder=1)

    # ±25% band
    x_band = np.linspace(lim[0], lim[1], 100)
    ax.fill_between(x_band, x_band * 0.75, x_band * 1.25, alpha=0.08, color="#888888", zorder=0)

    # Correlation
    r = float(np.corrcoef(all_exp, all_cfd)[0, 1])
    ax.text(0.03, 0.97, f"$R$ = {r:.2f}\n$n$ = {len(all_exp)}",
            transform=ax.transAxes, fontsize=9.5, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#CBD5E1"))

    ax.set_xlabel("Experimental $u/U_0$")
    ax.set_ylabel("CFD $u/U_0$")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25, linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", fontsize=8.5, framealpha=0.85)

    fig.tight_layout(pad=0.5)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}, R={r:.3f}, n={len(all_exp)}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
