#!/usr/bin/env python3
"""Fig: Nielsen validation — 2D velocity contour with measurement locations."""
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
PDF_OUT = OUT_DIR / "fig_method_validation_contour.pdf"
PNG_OUT = OUT_DIR / "fig_method_validation_contour.png"

L, H, D, U0 = 9.0, 3.0, 0.1, 0.455
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


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

    # Load VTK and slice at y=D/2
    vtk_path = None
    for d in sorted((CASE_DIR / "VTK").iterdir()):
        internal = d / "internal.vtu"
        if internal.exists():
            vtk_path = internal
    mesh = pv.read(vtk_path)
    if "U" in mesh.cell_data and "U" not in mesh.point_data:
        mesh = mesh.cell_data_to_point_data()

    # Compute velocity magnitude
    U = np.array(mesh.point_data["U"])
    mesh.point_data["Umag"] = np.linalg.norm(U, axis=1)

    # Sample on a 2D grid (x, z) at y=D/2
    nx, nz = 360, 120
    x_vals = np.linspace(0.01, L - 0.01, nx)
    z_vals = np.linspace(0.01, H - 0.01, nz)
    xx, zz = np.meshgrid(x_vals, z_vals)
    yy = np.full_like(xx, D / 2)
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    sampled = pv.PolyData(points).sample(mesh)
    umag = np.array(sampled.point_data["Umag"]).reshape(nz, nx)
    ux = np.array(sampled.point_data["U"])[:, 0].reshape(nz, nx)
    uz = np.array(sampled.point_data["U"])[:, 2].reshape(nz, nx)

    # Load experimental data for marker positions
    exp = json.load(open(CASE_DIR / "experimental_data.json"))

    fig, ax = plt.subplots(figsize=(7.2, 2.8))

    # Velocity magnitude contour
    levels = np.linspace(0, 0.5, 25)
    cf = ax.contourf(xx, zz, umag, levels=levels, cmap="viridis", extend="max")
    cbar = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.02)
    cbar.set_label("$|U|$ [m/s]", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Streamlines (subsample for clarity)
    skip = 6
    ax.streamplot(x_vals, z_vals, ux, uz, color="white", linewidth=0.5,
                  density=1.5, arrowsize=0.8, arrowstyle="->")

    # Measurement locations
    for pname, data in exp.items():
        x_frac = float(pname.split("_")[-1])
        y_H = np.array(data["y_over_H"])
        ax.scatter(np.full_like(y_H, x_frac * L), y_H * H,
                   c="red", s=12, marker="|", linewidths=1.0, zorder=5)

    # Profile location labels
    for x_frac, label in [(0.33, "A"), (0.50, "B"), (0.67, "C")]:
        ax.text(x_frac * L, H + 0.12, label, ha="center", fontsize=9,
                fontweight="bold", color="#DC2626")

    # Inlet/outlet annotations
    ax.annotate("Inlet", xy=(0, H - 0.084), fontsize=8, color="white",
                ha="left", va="center", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="white", lw=1.5),
                xytext=(0.5, H - 0.3))
    ax.annotate("Outlet", xy=(L, 0.084), fontsize=8, color="white",
                ha="right", va="center", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="white", lw=1.5),
                xytext=(L - 0.5, 0.4))

    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$z$ [m]")
    ax.set_xlim(0, L)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")

    fig.tight_layout(pad=0.5)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
