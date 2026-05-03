#!/usr/bin/env python3
"""Diagnostic visualizations for paint-booth downdraft uniformity.

This renders fields that reveal small deviations hidden by the uniform mean downdraft:
  1. horizontal relative Uz deviation stack below the filter;
  2. reverse-flow mask on y=0 mid-plane;
  3. work-zone clipped Uz on y=0 mid-plane;
  4. vorticity magnitude on y=0 mid-plane.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True


def latest_vtk_dir(case_dir: Path) -> Path:
    dirs = [p for p in (case_dir / "VTK").iterdir() if p.is_dir()]
    if not dirs:
        raise FileNotFoundError(f"No VTK time directories under {case_dir / 'VTK'}")
    def key(p: Path):
        try:
            return (1, float(p.name))
        except ValueError:
            return (0, p.name)
    return sorted(dirs, key=key)[-1]


def load_internal(case_dir: Path) -> pv.DataSet:
    vtk_dir = latest_vtk_dir(case_dir)
    files = list(vtk_dir.glob("internal.vtu")) + list(vtk_dir.glob("*internal*.vtu"))
    if not files:
        raise FileNotFoundError(f"No internal.vtu found in {vtk_dir}")
    mesh = pv.read(files[0])
    # Work in point data for clean sampling/slicing.
    if "U" in mesh.cell_data and "U" not in mesh.point_data:
        mesh = mesh.cell_data_to_point_data(pass_cell_data=True)
    U = np.asarray(mesh.point_data["U"])
    mesh.point_data["U_mag"] = np.linalg.norm(U, axis=1)
    mesh.point_data["Uz"] = U[:, 2]
    return mesh


def regular_xy_sample(mesh: pv.DataSet, z: float, nx: int = 260, ny: int = 160):
    x = np.linspace(-1.5, 6.5, nx)
    y = np.linspace(-2.0, 2.0, ny)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    zz = np.full_like(xx, z)
    grid = pv.StructuredGrid(xx, yy, zz)
    sampled = grid.sample(mesh)
    valid = np.asarray(sampled.point_data.get("vtkValidPointMask", np.ones(xx.size, dtype=bool))).reshape(yy.shape).astype(bool)
    U = np.asarray(sampled.point_data["U"]).reshape((*yy.shape, 3))
    Uz = U[:, :, 2]
    Uz = np.where(valid, Uz, np.nan)
    return x, y, Uz


def regular_xz_sample(mesh: pv.DataSet, y0: float = 0.0, nx: int = 280, nz: int = 180):
    x = np.linspace(-1.5, 6.5, nx)
    z = np.linspace(0.0, 3.8, nz)
    xx, zz = np.meshgrid(x, z, indexing="xy")
    yy = np.full_like(xx, y0)
    grid = pv.StructuredGrid(xx, yy, zz)
    sampled = grid.sample(mesh)
    valid = np.asarray(sampled.point_data.get("vtkValidPointMask", np.ones(xx.size, dtype=bool))).reshape(zz.shape).astype(bool)
    U = np.asarray(sampled.point_data["U"]).reshape((*zz.shape, 3))
    Uz = np.where(valid, U[:, :, 2], np.nan)
    Umag = np.where(valid, np.linalg.norm(U, axis=2), np.nan)
    return x, z, Uz, Umag


def draw_car_side(ax):
    # Approximate side-view envelope of the procedural smooth car shell.
    xs = [0.0, 0.25, 1.0, 2.0, 3.3, 4.2, 4.5, 4.5, 0.0, 0.0]
    zs = [0.35, 0.55, 0.95, 1.55, 1.35, 0.75, 0.45, 0.18, 0.18, 0.35]
    ax.fill(xs, zs, color="black", alpha=0.72, zorder=5)


def draw_car_top(ax):
    # Approximate top-view footprint.
    t = np.linspace(0, 2*np.pi, 200)
    x = 2.25 + 2.25*np.cos(t)
    y = 0.0 + 0.90*np.sin(t)
    ax.plot(x, y, color="black", lw=1.0, alpha=0.75)
    ax.fill(x, y, color="black", alpha=0.06)


def render_horizontal_deviation_stack(mesh: pv.DataSet, out: Path):
    zs = [2.70, 2.00, 1.20, 0.70]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.0), constrained_layout=True)
    for ax, z in zip(axes.ravel(), zs):
        x, y, Uz = regular_xy_sample(mesh, z)
        mean = np.nanmean(Uz)
        rel = (Uz - mean) / max(abs(mean), 1e-9) * 100.0
        rel = np.clip(rel, -60, 60)
        im = ax.imshow(rel, origin="lower", extent=[x.min(), x.max(), y.min(), y.max()], cmap="coolwarm", vmin=-60, vmax=60, aspect="equal")
        ax.contour(x, y, np.where(np.isfinite(Uz), Uz, np.nan), levels=[0.0], colors="red", linewidths=0.8)
        draw_car_top(ax)
        ax.set_title(f"z = {z:.2f} m, mean Uz = {mean:.3f} m/s")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_xlim(-1.5, 6.5)
        ax.set_ylim(-2.0, 2.0)
        ax.grid(color="k", alpha=0.12, lw=0.4)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.92, pad=0.015)
    cbar.set_label("Relative Uz deviation from plane mean [%]\nblue: stronger downdraft, red: weaker/upward")
    fig.suptitle("L2 horizontal slice stack: relative vertical-velocity deviation below filter", fontsize=15)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def render_reverse_mask(mesh: pv.DataSet, out: Path):
    x, z, Uz, _ = regular_xz_sample(mesh)
    fig, ax = plt.subplots(figsize=(12.5, 6.3), constrained_layout=True)
    base = np.ma.masked_invalid(np.clip(Uz, -0.45, 0.12))
    im = ax.imshow(base, origin="lower", extent=[x.min(), x.max(), z.min(), z.max()], cmap="coolwarm", vmin=-0.45, vmax=0.12, aspect="auto")
    reverse = np.ma.masked_where(~(Uz > 0.0), Uz)
    ax.imshow(reverse, origin="lower", extent=[x.min(), x.max(), z.min(), z.max()], cmap="Reds", vmin=0.0, vmax=0.18, alpha=0.78, aspect="auto")
    ax.contour(x, z, Uz, levels=[0.0], colors="black", linewidths=1.0)
    draw_car_side(ax)
    ax.axhspan(2.95, 3.05, color="gold", alpha=0.22, label="filter layer")
    ax.set_title("L2 y=0 mid-plane: reverse-flow mask (Uz > 0) over clipped Uz")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_xlim(-1.5, 6.5)
    ax.set_ylim(0, 3.8)
    ax.grid(color="k", alpha=0.12, lw=0.4)
    cbar = fig.colorbar(im, ax=ax, pad=0.015)
    cbar.set_label("Uz [m/s], clipped for work-zone contrast")
    ax.legend(loc="upper right")
    fig.savefig(out, dpi=180)
    plt.close(fig)


def render_clipped_uz(mesh: pv.DataSet, out: Path):
    x, z, Uz, _ = regular_xz_sample(mesh)
    fig, ax = plt.subplots(figsize=(12.5, 6.3), constrained_layout=True)
    im = ax.imshow(np.clip(Uz, -0.35, 0.10), origin="lower", extent=[x.min(), x.max(), z.min(), z.max()], cmap="coolwarm", vmin=-0.35, vmax=0.10, aspect="auto")
    ax.contour(x, z, Uz, levels=[-0.25, -0.15, -0.05, 0.0], colors=["#08306b", "#2171b5", "#fdae6b", "black"], linewidths=[0.7, 0.7, 0.7, 1.0])
    draw_car_side(ax)
    ax.axhspan(2.95, 3.05, color="gold", alpha=0.20)
    ax.set_title("L2 y=0 mid-plane: work-zone clipped vertical velocity Uz")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_xlim(-1.5, 6.5)
    ax.set_ylim(0, 3.8)
    ax.grid(color="k", alpha=0.12, lw=0.4)
    cbar = fig.colorbar(im, ax=ax, pad=0.015)
    cbar.set_label("Uz [m/s] clipped to [-0.35, 0.10]")
    fig.savefig(out, dpi=180)
    plt.close(fig)


def render_vorticity(mesh: pv.DataSet, out: Path):
    # Compute vorticity on point data. Polyhedral cells can warn but still produce a usable diagnostic field.
    derived = mesh.compute_derivative(scalars="U", vorticity=True, preference="point")
    if "vorticity" not in derived.point_data:
        # Some VTK versions name it differently.
        keys = list(derived.point_data.keys())
        raise KeyError(f"No vorticity field after derivative; point_data keys={keys}")
    vort = np.asarray(derived.point_data["vorticity"])
    derived.point_data["vort_mag"] = np.linalg.norm(vort, axis=1)
    # Sample the derived field on x-z plane.
    x = np.linspace(-1.5, 6.5, 280)
    z = np.linspace(0.0, 3.8, 180)
    xx, zz = np.meshgrid(x, z, indexing="xy")
    yy = np.zeros_like(xx)
    grid = pv.StructuredGrid(xx, yy, zz)
    sampled = grid.sample(derived)
    valid = np.asarray(sampled.point_data.get("vtkValidPointMask", np.ones(xx.size, dtype=bool))).reshape(zz.shape).astype(bool)
    w = np.asarray(sampled.point_data["vort_mag"]).reshape(zz.shape)
    w = np.where(valid, w, np.nan)
    vmax = np.nanpercentile(w, 98)
    fig, ax = plt.subplots(figsize=(12.5, 6.3), constrained_layout=True)
    im = ax.imshow(np.clip(w, 0, vmax), origin="lower", extent=[x.min(), x.max(), z.min(), z.max()], cmap="magma", vmin=0, vmax=vmax, aspect="auto")
    draw_car_side(ax)
    ax.axhspan(2.95, 3.05, color="cyan", alpha=0.12)
    ax.set_title("L2 y=0 mid-plane: vorticity magnitude |curl(U)|")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_xlim(-1.5, 6.5)
    ax.set_ylim(0, 3.8)
    ax.grid(color="w", alpha=0.16, lw=0.4)
    cbar = fig.colorbar(im, ax=ax, pad=0.015)
    cbar.set_label("|curl(U)| [1/s], clipped at p98")
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", type=Path, default=Path("cases/paint_booth_panel_frame_komega_l2_fine_yplus"))
    ap.add_argument("--out-dir", type=Path, default=Path("docs/figures"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    mesh = load_internal(args.case)
    prefix = "26-05-02_paint_booth_l2_diagnostic"
    outputs = {
        "uz_deviation_stack": args.out_dir / f"{prefix}_uz_deviation_horizontal_stack.png",
        "reverse_mask": args.out_dir / f"{prefix}_reverse_flow_mask_midplane.png",
        "clipped_uz": args.out_dir / f"{prefix}_clipped_uz_midplane.png",
        "vorticity": args.out_dir / f"{prefix}_vorticity_midplane.png",
    }
    render_horizontal_deviation_stack(mesh, outputs["uz_deviation_stack"])
    render_reverse_mask(mesh, outputs["reverse_mask"])
    render_clipped_uz(mesh, outputs["clipped_uz"])
    render_vorticity(mesh, outputs["vorticity"])
    for name, path in outputs.items():
        print(f"{name}: {path}")
    os._exit(0)


if __name__ == "__main__":
    main()
