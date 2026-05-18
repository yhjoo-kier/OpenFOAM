#!/usr/bin/env python3
"""QC-safe diagnostic visualizations using cell-center slab binning.

The first diagnostic renderer used point sampling on an unstructured/polyhedral VTK
mesh and produced strong row-wise artifacts. This version avoids that pathway:
  - uses original cell-centered U values;
  - selects finite-thickness slabs around each plotting plane;
  - bins cells into coarse 2-D bins by cell center;
  - masks low-count bins and solid car footprint explicitly;
  - uses robust/explicit scales and significant reverse threshold.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.path import Path as MplPath
import numpy as np
import pyvista as pv


def latest_vtk_dir(case_dir: Path) -> Path:
    dirs = [p for p in (case_dir / "VTK").iterdir() if p.is_dir()]
    if not dirs:
        raise FileNotFoundError(f"No VTK time directories under {case_dir/'VTK'}")
    def key(p: Path):
        try:
            return (1, float(p.name))
        except ValueError:
            return (0, p.name)
    return sorted(dirs, key=key)[-1]


def load_cell_data(case_dir: Path):
    vtk_dir = latest_vtk_dir(case_dir)
    files = list(vtk_dir.glob("internal.vtu")) + list(vtk_dir.glob("*internal*.vtu"))
    mesh = pv.read(files[0])
    centers = np.asarray(mesh.cell_centers().points)
    if "U" in mesh.cell_data:
        U = np.asarray(mesh.cell_data["U"])
    else:
        # Fallback, but keep values at cells after interpolation.
        cell_mesh = mesh.point_data_to_cell_data(pass_point_data=False)
        U = np.asarray(cell_mesh.cell_data["U"])
    return centers, U


def weighted_bin2d(x, y, values, x_edges, y_edges, min_count=1):
    # np.histogram2d returns shape (len(x_edges)-1, len(y_edges)-1) when called x,y.
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(values)
    x = x[mask]; y = y[mask]; values = values[mask]
    sums, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=values)
    counts, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sums / counts
    mean[counts < min_count] = np.nan
    return mean.T, counts.T  # rows=y, cols=x for pcolormesh


def fill_nan_neighbor(a, iterations=2):
    out = a.copy()
    for _ in range(iterations):
        nan = ~np.isfinite(out)
        if not nan.any():
            break
        padded = np.pad(out, 1, mode="constant", constant_values=np.nan)
        vals = []
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                vals.append(padded[1+di:1+di+out.shape[0], 1+dj:1+dj+out.shape[1]])
        neigh = np.nanmean(np.stack(vals), axis=0)
        out[nan & np.isfinite(neigh)] = neigh[nan & np.isfinite(neigh)]
    return out


def car_side_polygon():
    return np.array([
        [0.00, 0.18], [0.00, 0.35], [0.25, 0.55], [1.00, 0.95],
        [2.00, 1.55], [3.30, 1.35], [4.20, 0.75], [4.50, 0.45],
        [4.50, 0.18], [0.00, 0.18]
    ])


def draw_car_side(ax, face="0.15", edge="black", alpha=0.95):
    poly = car_side_polygon()
    ax.fill(poly[:,0], poly[:,1], color=face, alpha=alpha, zorder=10)
    ax.plot(poly[:,0], poly[:,1], color=edge, lw=1.2, zorder=11)


def draw_car_top(ax):
    # Better than previous ellipse: rectangular-ish footprint matching smooth shell bounds.
    x0, x1 = 0.0, 4.5
    y0, y1 = -0.9, 0.9
    ax.add_patch(plt.Rectangle((x0, y0), x1-x0, y1-y0, fill=False, edgecolor="black", lw=1.2, zorder=5))
    ax.text(0.08, 0.76, "car footprint", fontsize=7, color="black", zorder=6)


def mask_car_side_array(xc, zc, arr):
    X, Z = np.meshgrid(xc, zc)
    pts = np.c_[X.ravel(), Z.ravel()]
    inside = MplPath(car_side_polygon()).contains_points(pts).reshape(arr.shape)
    out = arr.copy()
    out[inside] = np.nan
    return out, inside


def render_horizontal_deviation(centers, U, out: Path):
    Uz = U[:,2]
    z_levels = [2.70, 2.00, 1.20, 0.70]
    x_edges = np.linspace(-1.2, 6.2, 75)  # coarser, within filter/work-zone footprint neighborhood
    y_edges = np.linspace(-1.8, 1.8, 49)
    x_cent = 0.5*(x_edges[:-1]+x_edges[1:])
    y_cent = 0.5*(y_edges[:-1]+y_edges[1:])
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.2), constrained_layout=True)
    im = None
    for ax, z0 in zip(axes.ravel(), z_levels):
        slab_half_thickness = 0.08
        slab = np.abs(centers[:,2] - z0) <= slab_half_thickness
        # Limit stats to central work-zone footprint, not side-wall bins.
        stat_mask = slab & (centers[:,0] >= -1.0) & (centers[:,0] <= 6.0) & (centers[:,1] >= -1.6) & (centers[:,1] <= 1.6)
        if not np.any(stat_mask):
            slab_half_thickness = 0.12
            slab = np.abs(centers[:,2] - z0) <= slab_half_thickness
            stat_mask = slab & (centers[:,0] >= -1.0) & (centers[:,0] <= 6.0) & (centers[:,1] >= -1.6) & (centers[:,1] <= 1.6)
        mean = float(np.nanmean(Uz[stat_mask]))
        binned, counts = weighted_bin2d(centers[slab,0], centers[slab,1], Uz[slab], x_edges, y_edges, min_count=1)
        # Neighbour fill only for visualization continuity, not for stats.
        binned_filled = fill_nan_neighbor(binned, iterations=1)
        rel = 100.0 * (binned_filled - mean) / max(abs(mean), 1e-9)
        # sign convention: negative rel = stronger downward when mean is negative.
        rel = np.clip(rel, -40, 40)
        im = ax.pcolormesh(x_edges, y_edges, rel, shading="auto", cmap="coolwarm", norm=TwoSlopeNorm(vmin=-40, vcenter=0, vmax=40))
        ax.contour(x_cent, y_cent, rel, levels=[-20,-10,0,10,20], colors="k", linewidths=[0.35,0.35,0.7,0.35,0.35], alpha=0.55)
        # Show central filter footprint and car footprint explicitly.
        ax.add_patch(plt.Rectangle((-1.0, -1.6), 7.0, 3.2, fill=False, edgecolor="goldenrod", lw=1.4, ls="--", zorder=6))
        draw_car_top(ax)
        ax.set_title(f"z={z0:.2f} m; slab ±{slab_half_thickness:.2f} m; <Uz>={mean:.3f} m/s")
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        ax.set_xlim(-1.2, 6.2); ax.set_ylim(-1.8, 1.8); ax.set_aspect("equal")
        ax.grid(color="k", alpha=0.08, lw=0.3)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.92, pad=0.015)
    cbar.set_label("100×(Uz-<Uz>)/|<Uz>| [%]\nblue=stronger downward, red=weaker/upward; clipped ±40%")
    fig.suptitle("QC-safe cell-center slab binning: horizontal Uz deviation below filter", fontsize=15)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def render_midplane(centers, U, out_reverse: Path, out_clipped: Path, out_shear: Path):
    Uz = U[:,2]
    Umag = np.linalg.norm(U, axis=1)
    slab = np.abs(centers[:,1]) <= 0.04
    x_edges = np.linspace(-1.5, 6.5, 130)
    z_edges = np.linspace(0.0, 3.8, 96)
    x_cent = 0.5*(x_edges[:-1]+x_edges[1:])
    z_cent = 0.5*(z_edges[:-1]+z_edges[1:])
    uz_bin, counts = weighted_bin2d(centers[slab,0], centers[slab,2], Uz[slab], x_edges, z_edges, min_count=1)
    uz_fill = fill_nan_neighbor(uz_bin, iterations=1)
    uz_masked, car_inside = mask_car_side_array(x_cent, z_cent, uz_fill)

    # Reverse plot: significant positive Uz only, not threshold-noise at 0.
    fig, ax = plt.subplots(figsize=(12.5, 6.3), constrained_layout=True)
    im = ax.pcolormesh(x_edges, z_edges, np.clip(uz_masked, -0.32, 0.12), shading="auto", cmap="coolwarm", norm=TwoSlopeNorm(vmin=-0.32, vcenter=0, vmax=0.12))
    reverse = np.ma.masked_where(~(uz_masked > 0.02), uz_masked)
    ax.pcolormesh(x_edges, z_edges, reverse, shading="auto", cmap="spring", vmin=0.02, vmax=max(0.12, np.nanpercentile(uz_masked,99)), alpha=0.62)
    ax.contour(x_cent, z_cent, uz_masked, levels=[0.0], colors="black", linewidths=1.3)
    draw_car_side(ax)
    ax.axhspan(2.95, 3.05, facecolor="gold", alpha=0.20, edgecolor="goldenrod", lw=1.0, label="filter layer")
    rev_frac = np.nanmean(uz_masked > 0.02) * 100.0
    ax.set_title(f"QC-safe reverse-flow map: Uz > +0.02 m/s, y-slab ±0.04 m; area bins ≈ {rev_frac:.1f}%")
    ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.set_xlim(-1.5,6.5); ax.set_ylim(0,3.8)
    ax.grid(color="k", alpha=0.08, lw=0.3); ax.legend(loc="upper right")
    cb = fig.colorbar(im, ax=ax, pad=0.015); cb.set_label("Uz [m/s], clipped; magenta overlay = significant upward/reverse")
    fig.savefig(out_reverse, dpi=180); plt.close(fig)

    # Clipped Uz plot with solid mask.
    fig, ax = plt.subplots(figsize=(12.5, 6.3), constrained_layout=True)
    im = ax.pcolormesh(x_edges, z_edges, np.clip(uz_masked, -0.35, 0.10), shading="auto", cmap="coolwarm", norm=TwoSlopeNorm(vmin=-0.35, vcenter=0, vmax=0.10))
    ax.contour(x_cent, z_cent, uz_masked, levels=[-0.25,-0.15,-0.05,0.0], colors=["#08306b","#2171b5","#fdae6b","black"], linewidths=[0.7,0.7,0.8,1.2])
    draw_car_side(ax)
    ax.axhspan(2.95, 3.05, color="gold", alpha=0.18)
    ax.set_title("QC-safe work-zone Uz: cell-center y-slab binning, solid car masked")
    ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.set_xlim(-1.5,6.5); ax.set_ylim(0,3.8)
    ax.grid(color="k", alpha=0.08, lw=0.3)
    cb = fig.colorbar(im, ax=ax, pad=0.015); cb.set_label("Uz [m/s], clipped to [-0.35, 0.10]")
    fig.savefig(out_clipped, dpi=180); plt.close(fig)

    # QC-safe shear proxy: gradient magnitude of binned Uz, not raw VTK vorticity.
    uz_for_grad = fill_nan_neighbor(uz_masked, iterations=4)
    dx = np.mean(np.diff(x_cent)); dz = np.mean(np.diff(z_cent))
    dudz, dudx = np.gradient(uz_for_grad, dz, dx)
    shear = np.sqrt(dudx**2 + dudz**2)
    shear[car_inside] = np.nan
    vmax = np.nanpercentile(shear, 97)
    fig, ax = plt.subplots(figsize=(12.5, 6.3), constrained_layout=True)
    cmap = plt.get_cmap("magma").copy(); cmap.set_bad("#e0e0e0")
    im = ax.pcolormesh(x_edges, z_edges, np.clip(shear, 0, vmax), shading="auto", cmap=cmap, vmin=0, vmax=vmax)
    draw_car_side(ax, face="0.08", edge="white", alpha=0.95)
    ax.axhspan(2.95, 3.05, color="cyan", alpha=0.12)
    ax.set_title("QC-safe shear proxy |∇Uz| on y=0 slab (raw vorticity plot failed QC)")
    ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.set_xlim(-1.5,6.5); ax.set_ylim(0,3.8)
    ax.grid(color="w", alpha=0.10, lw=0.3)
    cb = fig.colorbar(im, ax=ax, pad=0.015); cb.set_label("|∇Uz| [1/s], clipped at p97")
    fig.savefig(out_shear, dpi=180); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", type=Path, default=Path("cases/paint_booth_panel_frame_komega_l2_fine_yplus"))
    ap.add_argument("--out-dir", type=Path, default=Path("docs/figures"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    centers, U = load_cell_data(args.case)
    prefix = "26-05-02_paint_booth_l2_qc_safe"
    render_horizontal_deviation(centers, U, args.out_dir / f"{prefix}_uz_deviation_horizontal_stack.png")
    render_midplane(
        centers, U,
        args.out_dir / f"{prefix}_reverse_flow_mask_midplane.png",
        args.out_dir / f"{prefix}_clipped_uz_midplane.png",
        args.out_dir / f"{prefix}_shear_proxy_midplane.png",
    )
    for p in sorted(args.out_dir.glob(f"{prefix}_*.png")):
        print(p)
    os._exit(0)


if __name__ == "__main__":
    main()
