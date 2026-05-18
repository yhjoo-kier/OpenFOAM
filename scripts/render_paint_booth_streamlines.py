#!/usr/bin/env python3
"""Render QC-conscious streamlines for the L2 paint-booth CFD case.

The intent is qualitative flow-path visualization, not streamline-based
quantification. Seeds are restricted to physically meaningful regions:
  - just below the filter panel;
  - around the car roof/wake region.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

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


def load_mesh_and_car(case_dir: Path):
    vtk_dir = latest_vtk_dir(case_dir)
    internal = list(vtk_dir.glob("internal.vtu")) + list(vtk_dir.glob("*internal*.vtu"))
    if not internal:
        raise FileNotFoundError(f"No internal .vtu found under {vtk_dir}")
    mesh = pv.read(internal[0])
    if "U" in mesh.cell_data and "U" not in mesh.point_data:
        mesh = mesh.cell_data_to_point_data(pass_cell_data=True)
    U = np.asarray(mesh.point_data["U"])
    mesh.point_data["U_mag"] = np.linalg.norm(U, axis=1)
    mesh.point_data["Uz"] = U[:, 2]

    boundary = vtk_dir / "boundary"
    car_files = list(boundary.glob("carBody.vtp")) + list(boundary.glob("carBody*.vtp"))
    car = pv.read(car_files[0]) if car_files else None
    return mesh, car


def add_context(pl: pv.Plotter, car: pv.PolyData | None, *, car_opacity=0.88):
    domain = pv.Box(bounds=(-1.5, 6.5, -2.0, 2.0, 0.0, 3.8))
    pl.add_mesh(domain, style="wireframe", color="#333333", line_width=1.0, opacity=0.35)
    filter_panel = pv.Plane(center=(2.5, 0.0, 3.0), direction=(0, 0, 1), i_size=7.0, j_size=3.2)
    pl.add_mesh(filter_panel, color="#e6b800", opacity=0.20)
    supply = pv.Plane(center=(2.5, 0.0, 3.805), direction=(0, 0, 1), i_size=1.5, j_size=1.2)
    pl.add_mesh(supply, color="#1f78b4", opacity=0.55)
    if car is not None:
        pl.add_mesh(car, color="#3b3b3b", opacity=car_opacity, smooth_shading=True)
        pl.add_mesh(car, style="wireframe", color="black", opacity=0.10, line_width=0.2)


def seed_grid(x_values, y_values, z: float) -> pv.PolyData:
    pts = np.array([[x, y, z] for x in x_values for y in y_values], dtype=float)
    return pv.PolyData(pts)


def make_streamlines(mesh: pv.DataSet, seed: pv.PolyData, *, max_length: float, terminal_speed: float = 1e-4):
    # streamlines_from_source expects point vector field.
    sl = mesh.streamlines_from_source(
        seed,
        vectors="U",
        integrator_type=45,
        integration_direction="both",
        max_time=max_length,
        max_steps=4000,
        initial_step_length=0.015,
        min_step_length=0.003,
        max_step_length=0.06,
        terminal_speed=terminal_speed,
        surface_streamlines=False,
    )
    if sl.n_points and "U" in sl.point_data:
        U = np.asarray(sl.point_data["U"])
        sl.point_data["U_mag"] = np.linalg.norm(U, axis=1)
        sl.point_data["Uz"] = U[:, 2]
    return sl


def add_streamlines(pl: pv.Plotter, sl: pv.PolyData, *, scalar="U_mag", clim=(0, 1.2), radius=0.008):
    if sl.n_points == 0:
        return
    tube = sl.tube(radius=radius, n_sides=8)
    pl.add_mesh(tube, scalars=scalar, cmap="turbo", clim=clim, opacity=0.92,
                scalar_bar_args={"title": "|U| [m/s]", "vertical": True})


def render_filter_overview(mesh, car, out: Path):
    # Seeds just below filter, inside central panel but coarser to avoid clutter.
    seed = seed_grid(np.linspace(-0.6, 5.6, 9), np.linspace(-1.35, 1.35, 7), 2.72)
    sl = make_streamlines(mesh, seed, max_length=4.8)
    pl = pv.Plotter(off_screen=True, window_size=(1900, 1300))
    pl.set_background("white")
    add_context(pl, car)
    add_streamlines(pl, sl, clim=(0, 1.25), radius=0.009)
    pl.add_mesh(seed, color="magenta", point_size=12, render_points_as_spheres=True)
    pl.add_text("L2 streamlines seeded below filter panel", font_size=18, color="black")
    pl.add_text("magenta dots=seed points at z=2.72 m; yellow=filter, blue=supply", position="lower_left", font_size=11, color="black")
    pl.camera_position = [(8.2, -7.4, 4.6), (2.35, 0.0, 1.85), (0, 0, 1)]
    pl.camera.zoom(1.04)
    pl.show(screenshot=str(out), auto_close=True)


def render_car_closeup(mesh, car, out: Path):
    # Seeds around roof and wake, not inside car. Use multiple y levels for 3D paths.
    pts = []
    for x in np.linspace(0.35, 4.35, 9):
        for y in [-0.95, -0.55, 0.0, 0.55, 0.95]:
            pts.append([x, y, 1.82])
    for x in np.linspace(0.2, 4.5, 8):
        for y in [-1.15, 1.15]:
            pts.append([x, y, 1.05])
    seed = pv.PolyData(np.array(pts, dtype=float))
    sl = make_streamlines(mesh, seed, max_length=3.2)
    pl = pv.Plotter(off_screen=True, window_size=(1900, 1300))
    pl.set_background("white")
    # Add cropped near-car box instead of full domain for readability.
    box = pv.Box(bounds=(-0.4, 5.0, -1.55, 1.55, 0.0, 2.4))
    pl.add_mesh(box, style="wireframe", color="#555555", opacity=0.35, line_width=1.1)
    # Add filter-height guide line/plane for orientation.
    filter_guide = pv.Plane(center=(2.5, 0.0, 2.95), direction=(0, 0, 1), i_size=5.4, j_size=3.1)
    pl.add_mesh(filter_guide, color="#e6b800", opacity=0.12)
    if car is not None:
        pl.add_mesh(car, color="#3c3c3c", opacity=0.93, smooth_shading=True)
    add_streamlines(pl, sl, clim=(0, 0.85), radius=0.007)
    pl.add_mesh(seed, color="magenta", point_size=10, render_points_as_spheres=True)
    pl.add_text("L2 near-car local streamlines", font_size=18, color="black")
    pl.add_text("short seeded paths around roof/side/wake; qualitative only", position="lower_left", font_size=11, color="black")
    pl.camera_position = [(3.2, -6.2, 2.85), (2.2, 0.0, 1.05), (0, 0, 1)]
    pl.camera.zoom(1.05)
    pl.show(screenshot=str(out), auto_close=True)


def render_side_projection(mesh, car, out: Path):
    # y=0-ish seeded line to make path curvature readable in side view.
    seed = seed_grid(np.linspace(-0.5, 5.5, 13), [0.0], 2.72)
    seed2 = seed_grid(np.linspace(0.2, 4.4, 10), [0.0], 1.75)
    seed_all = pv.PolyData(np.vstack([np.asarray(seed.points), np.asarray(seed2.points)]))
    sl = make_streamlines(mesh, seed_all, max_length=4.2)
    pl = pv.Plotter(off_screen=True, window_size=(1800, 1150))
    pl.set_background("white")
    # Add center-plane scalar field faintly for context.
    slc = mesh.slice(normal="y", origin=(2.5, 0.0, 1.9))
    pl.add_mesh(slc, scalars="Uz", cmap="coolwarm", clim=(-0.35, 0.10), opacity=0.38,
                scalar_bar_args={"title": "Uz [m/s]", "vertical": True})
    if car is not None:
        pl.add_mesh(car, color="#404040", opacity=0.84, smooth_shading=True)
    add_streamlines(pl, sl, clim=(0, 1.0), radius=0.007)
    pl.add_mesh(seed_all, color="magenta", point_size=8, render_points_as_spheres=True)
    # filter band guide
    panel = pv.Plane(center=(2.5, 0.0, 3.0), direction=(0, 0, 1), i_size=7.0, j_size=0.04)
    pl.add_mesh(panel, color="#e6b800", opacity=0.40)
    pl.add_text("L2 side-view streamlines over clipped Uz context", font_size=18, color="black")
    pl.add_text("center-line seeds show downdraft deflection around car", position="lower_left", font_size=11, color="black")
    pl.camera_position = [(2.5, -10.5, 1.9), (2.5, 0.0, 1.9), (0, 0, 1)]
    pl.camera.zoom(1.18)
    pl.show(screenshot=str(out), auto_close=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", type=Path, default=Path("cases/paint_booth_panel_frame_komega_l2_fine_yplus"))
    ap.add_argument("--out-dir", type=Path, default=Path("docs/figures"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    mesh, car = load_mesh_and_car(args.case)
    prefix = "26-05-02_paint_booth_l2_streamlines"
    outputs = {
        "filter_overview": args.out_dir / f"{prefix}_filter_seed_overview.png",
        "near_car": args.out_dir / f"{prefix}_near_car_closeup.png",
        "side_projection": args.out_dir / f"{prefix}_side_projection.png",
    }
    render_filter_overview(mesh, car, outputs["filter_overview"])
    render_car_closeup(mesh, car, outputs["near_car"])
    render_side_projection(mesh, car, outputs["side_projection"])
    for k, p in outputs.items():
        print(f"{k}: {p}")
    os._exit(0)


if __name__ == "__main__":
    main()
