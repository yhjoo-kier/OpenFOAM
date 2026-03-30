#!/usr/bin/env python3
"""Render BFS flow field visualizations for CFD Visual QA benchmark."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CASE_BASE = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark" / "cases"
IMAGE_DIR = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark" / "images"
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

# BFS geometry
H = 1.0
L_UP = 5.0
L_DOWN = 30.0

# Visualization configs
VIZ_TYPES = {
    "V1": {"field": "Umag", "label": "$|U|$ [m/s]", "cmap": "viridis", "type": "contour"},
    "V3": {"field": "p", "label": "$p$ [m²/s²]", "cmap": "RdBu_r", "type": "contour"},
    "V4": {"field": "U", "label": "$U$ [m/s]", "cmap": "viridis", "type": "vector"},
    "V5": {"field": "U", "label": "", "cmap": "viridis", "type": "streamline"},
    "V6": {"field": "Umag", "label": "$|U|$ [m/s]", "cmap": "viridis", "type": "composite"},
}


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def load_case(case_dir: Path) -> pv.DataSet | None:
    vtk_dir = case_dir / "VTK"
    if not vtk_dir.exists():
        return None
    vtk_path = None
    for d in sorted(vtk_dir.iterdir()):
        internal = d / "internal.vtu"
        if internal.exists():
            vtk_path = internal
    if not vtk_path:
        return None
    mesh = pv.read(vtk_path)
    if "U" in mesh.cell_data and "U" not in mesh.point_data:
        mesh = mesh.cell_data_to_point_data()
    U = np.array(mesh.point_data["U"])
    mesh.point_data["Umag"] = np.linalg.norm(U, axis=1)
    return mesh


def sample_2d(mesh: pv.DataSet, nx: int = 400, ny: int = 160) -> dict:
    """Sample on a 2D grid at z=W/2."""
    x_vals = np.linspace(-L_UP + 0.01, L_DOWN - 0.01, nx)
    y_vals = np.linspace(0.01, 2 * H - 0.01, ny)
    xx, yy = np.meshgrid(x_vals, y_vals)
    zz = np.full_like(xx, 0.05)  # mid-plane
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    sampled = pv.PolyData(points).sample(mesh)

    umag = np.array(sampled.point_data["Umag"]).reshape(ny, nx)
    U = np.array(sampled.point_data["U"])
    ux = U[:, 0].reshape(ny, nx)
    uy = U[:, 1].reshape(ny, nx)
    p = np.array(sampled.point_data["p"]).reshape(ny, nx)

    # Mask upstream lower region (inside step)
    for j, y in enumerate(y_vals):
        for i, x in enumerate(x_vals):
            if x < 0 and y < H:
                umag[j, i] = np.nan
                ux[j, i] = np.nan
                uy[j, i] = np.nan
                p[j, i] = np.nan

    return {
        "xx": xx, "yy": yy, "x_vals": x_vals, "y_vals": y_vals,
        "umag": umag, "ux": ux, "uy": uy, "p": p,
    }


def render_contour(ax, data, field_key, cmap, label, fig):
    field = data[field_key]
    vmin = np.nanmin(field)
    vmax = np.nanmax(field)
    levels = np.linspace(vmin, vmax, 30)
    cf = ax.contourf(data["xx"], data["yy"], field, levels=levels, cmap=cmap, extend="both")
    cbar = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.02, aspect=25)
    cbar.set_label(label, fontsize=9)
    cbar.ax.tick_params(labelsize=7)


def render_vectors(ax, data, cmap, fig):
    skip_x, skip_y = 12, 4
    ux_s = data["ux"][::skip_y, ::skip_x]
    uy_s = data["uy"][::skip_y, ::skip_x]
    xx_s = data["xx"][::skip_y, ::skip_x]
    yy_s = data["yy"][::skip_y, ::skip_x]
    mag_s = np.sqrt(ux_s**2 + uy_s**2)
    ax.quiver(xx_s, yy_s, ux_s, uy_s, mag_s, cmap=cmap, scale=25, width=0.002, alpha=0.8)


def render_streamlines(ax, data):
    x = data["x_vals"]
    y = data["y_vals"]
    ux = np.nan_to_num(data["ux"], nan=0.0)
    uy = np.nan_to_num(data["uy"], nan=0.0)
    speed = np.sqrt(ux**2 + uy**2)
    ax.streamplot(x, y, ux, uy, color=speed, cmap="viridis",
                  density=2.0, linewidth=0.7, arrowsize=0.8)


def draw_step(ax):
    """Draw step geometry outline."""
    step_x = [-L_UP, -L_UP, 0, 0, L_DOWN, L_DOWN, -L_UP]
    step_y = [2*H, H, H, 0, 0, 2*H, 2*H]
    ax.plot(step_x, step_y, "k-", linewidth=1.2)


def render_single(case_name: str, viz_id: str, data: dict, meta: dict, font: str) -> Path:
    viz = VIZ_TYPES[viz_id]
    plt.rcParams.update({
        "font.family": font, "font.size": 10,
        "axes.labelsize": 10, "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(10, 3.5))

    if viz["type"] == "contour":
        field_key = "umag" if viz["field"] == "Umag" else viz["field"]
        render_contour(ax, data, field_key, viz["cmap"], viz["label"], fig)
    elif viz["type"] == "vector":
        render_contour(ax, data, "umag", viz["cmap"], viz["label"], fig)
        render_vectors(ax, data, "gray_r", fig)
    elif viz["type"] == "streamline":
        render_streamlines(ax, data)
    elif viz["type"] == "composite":
        render_contour(ax, data, "umag", viz["cmap"], viz["label"], fig)
        x = data["x_vals"]
        y = data["y_vals"]
        ux = np.nan_to_num(data["ux"], nan=0.0)
        uy = np.nan_to_num(data["uy"], nan=0.0)
        ax.streamplot(x, y, ux, uy, color="white", density=1.5,
                      linewidth=0.5, arrowsize=0.7, arrowstyle="->")

    draw_step(ax)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_xlim(-L_UP, L_DOWN)
    ax.set_ylim(0, 2 * H)
    ax.set_aspect("equal")

    fig.tight_layout(pad=0.5)

    out_dir = IMAGE_DIR / case_name
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{case_name}_{viz_id}.png"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return png_path


def render_case(case_name: str) -> list[Path]:
    case_dir = CASE_BASE / case_name
    meta_path = case_dir / "case_meta.json"
    if not meta_path.exists():
        print(f"  Skip {case_name}: no metadata")
        return []

    meta = json.load(open(meta_path))
    mesh = load_case(case_dir)
    if mesh is None:
        print(f"  Skip {case_name}: no VTK data")
        return []

    print(f"  Sampling 2D grid for {case_name}...")
    data = sample_2d(mesh)

    paths = []
    for viz_id in VIZ_TYPES:
        print(f"    Rendering {viz_id}...")
        p = render_single(case_name, viz_id, data, meta, pick_font())
        paths.append(p)
        print(f"    Wrote {p}")

    return paths


def main():
    print("=" * 60)
    print("CFD Visual QA — Rendering BFS Visualizations")
    print("=" * 60)

    font = pick_font()
    print(f"Font: {font}")

    cases = sorted(CASE_BASE.iterdir())
    all_paths = []
    for case_dir in cases:
        if not case_dir.is_dir():
            continue
        name = case_dir.name
        if not name.startswith("S4_"):
            continue
        print(f"\nRendering: {name}")
        paths = render_case(name)
        all_paths.extend(paths)

    print(f"\nTotal images rendered: {len(all_paths)}")


if __name__ == "__main__":
    main()
