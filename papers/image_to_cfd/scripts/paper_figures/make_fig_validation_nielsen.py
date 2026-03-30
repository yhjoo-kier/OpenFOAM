#!/usr/bin/env python3
"""Fig: Nielsen (1990) solver validation — CFD vs experimental velocity profiles."""
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
PDF_OUT = OUT_DIR / "fig_method_validation_nielsen.pdf"
PNG_OUT = OUT_DIR / "fig_method_validation_nielsen.png"

L, H, D, U0 = 9.0, 3.0, 0.1, 0.455
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

PROFILES = [
    {"name": "x_over_L_0.33", "x_frac": 0.33, "label": "$x/L$ = 1/3"},
    {"name": "x_over_L_0.50", "x_frac": 0.50, "label": "$x/L$ = 1/2"},
    {"name": "x_over_L_0.67", "x_frac": 0.67, "label": "$x/L$ = 2/3"},
]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 10,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    # Load experimental data
    exp = json.load(open(CASE_DIR / "experimental_data.json"))

    # Load CFD (refined mesh)
    vtk_dirs = sorted((CASE_DIR / "VTK").iterdir())
    vtk_path = None
    for d in vtk_dirs:
        internal = d / "internal.vtu"
        if internal.exists():
            vtk_path = internal
    if not vtk_path:
        raise FileNotFoundError("No VTK internal.vtu found")

    mesh = pv.read(vtk_path)
    if "U" in mesh.cell_data and "U" not in mesh.point_data:
        mesh = mesh.cell_data_to_point_data()

    # Load validation results for annotations
    val_results = json.load(open(CASE_DIR / "validation_results_refined.json"))

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.2), sharey=True)
    panel_labels = ["(a)", "(b)", "(c)"]

    for idx, prof in enumerate(PROFILES):
        ax = axes[idx]
        data = exp[prof["name"]]
        y_H = np.array(data["y_over_H"])
        u_exp = np.array(data["u_over_U0"])

        # Sample CFD
        points = np.column_stack([
            np.full_like(y_H, prof["x_frac"] * L),
            np.full_like(y_H, D / 2),
            y_H * H,
        ])
        sampled = pv.PolyData(points).sample(mesh)
        u_cfd = np.array(sampled.point_data["U"])[:, 0] / U0

        # Also sample a dense CFD profile for smooth line
        z_dense = np.linspace(0.01 * H, 0.99 * H, 100)
        pts_dense = np.column_stack([
            np.full_like(z_dense, prof["x_frac"] * L),
            np.full_like(z_dense, D / 2),
            z_dense,
        ])
        sampled_dense = pv.PolyData(pts_dense).sample(mesh)
        u_cfd_dense = np.array(sampled_dense.point_data["U"])[:, 0] / U0
        z_H_dense = z_dense / H

        # Plot
        ax.plot(u_cfd_dense, z_H_dense, color="#2563EB", linewidth=1.5, label="CFD (k-ω SST)", zorder=3)
        ax.scatter(u_exp, y_H, color="#DC2626", s=28, marker="o", zorder=4,
                   edgecolors="white", linewidths=0.4, label="Exp. (Nielsen 1990)")

        ax.set_xlabel("$u / U_0$")
        if idx == 0:
            ax.set_ylabel("$z / H$")
        ax.set_title(prof["label"], fontsize=10, pad=6)
        ax.set_xlim(-0.25, 0.45)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.25, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.axvline(0, color="#888888", linewidth=0.5, linestyle="--", alpha=0.5)

        # Annotation: R and q
        prof_result = val_results["profiles"][prof["name"]]
        r_val = prof_result["R"]
        q_val = prof_result["q"]
        ax.text(0.97, 0.03, f"$R$ = {r_val:.2f}\n$q$ = {q_val:.2f}",
                transform=ax.transAxes, fontsize=8.5, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="#CBD5E1"))

        # Panel label
        ax.text(-0.15, 1.08, panel_labels[idx], transform=ax.transAxes,
                fontsize=11, fontweight="bold", color="#34495E")

    # Shared legend
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
                   frameon=False, fontsize=9, handletextpad=0.4)

    fig.tight_layout(rect=[0.0, 0.08, 1.0, 1.0], w_pad=1.0)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
