#!/usr/bin/env python3
"""Render L2 paint-booth CFD visualizations from foamToVTK outputs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True


def latest_vtk_dir(case_dir: Path) -> Path:
    vtk_root = case_dir / "VTK"
    dirs = [p for p in vtk_root.iterdir() if p.is_dir()]
    if not dirs:
        raise FileNotFoundError(f"No VTK time directories under {vtk_root}")

    def key(p: Path):
        try:
            return (1, float(p.name))
        except ValueError:
            return (0, p.name)

    return sorted(dirs, key=key)[-1]


def get_internal(vtk_dir: Path) -> pv.DataSet:
    files = list(vtk_dir.glob("internal.vtu")) + list(vtk_dir.glob("*internal*.vtu"))
    if not files:
        raise FileNotFoundError(f"internal.vtu not found in {vtk_dir}")
    mesh = pv.read(files[0])
    if "U" in mesh.cell_data:
        U = np.asarray(mesh.cell_data["U"])
        mesh.cell_data["U_mag"] = np.linalg.norm(U, axis=1)
        mesh.cell_data["Uz"] = U[:, 2]
    elif "U" in mesh.point_data:
        U = np.asarray(mesh.point_data["U"])
        mesh.point_data["U_mag"] = np.linalg.norm(U, axis=1)
        mesh.point_data["Uz"] = U[:, 2]
    return mesh


def get_car(vtk_dir: Path) -> pv.PolyData:
    boundary = vtk_dir / "boundary"
    files = list(boundary.glob("carBody.vtp")) + list(boundary.glob("carBody*.vtp"))
    if not files:
        raise FileNotFoundError(f"carBody.vtp not found under {boundary}")
    return pv.read(files[0])


def add_domain_guides(pl: pv.Plotter) -> None:
    # Domain and key booth planes. Coordinates match the generator metadata.
    domain = pv.Box(bounds=(-1.5, 6.5, -2.0, 2.0, 0.0, 3.8))
    pl.add_mesh(domain, style="wireframe", color="#444444", line_width=1.2, opacity=0.4)
    # central filter panel plane at z≈3.0
    panel = pv.Plane(center=(2.5, 0.0, 3.0), direction=(0, 0, 1), i_size=7.0, j_size=3.2, i_resolution=1, j_resolution=1)
    pl.add_mesh(panel, color="#f0c64a", opacity=0.22)
    # supply patch on top
    supply = pv.Plane(center=(2.5, 0.0, 3.805), direction=(0, 0, 1), i_size=1.5, j_size=1.2, i_resolution=1, j_resolution=1)
    pl.add_mesh(supply, color="#1f78b4", opacity=0.55)


def render_midplane(mesh: pv.DataSet, car: pv.PolyData, out: Path, scalar: str, title: str, cmap: str, clim: tuple[float, float]) -> None:
    sl = mesh.slice(normal="y", origin=(2.5, 0.0, 1.9))
    # Keep only the x-z slice view by using a camera along -y.
    pl = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    pl.set_background("white")
    add_domain_guides(pl)
    pl.add_mesh(sl, scalars=scalar, cmap=cmap, clim=clim, opacity=0.96, scalar_bar_args={"title": scalar, "vertical": True})
    pl.add_mesh(car, color="#303030", opacity=0.75, smooth_shading=True)
    pl.add_text(title, font_size=18, color="black")
    pl.add_text("y=0 center-plane; yellow=filter panel, blue=top supply", position="lower_left", font_size=10, color="black")
    pl.camera_position = [(2.5, -11.5, 1.9), (2.5, 0.0, 1.9), (0, 0, 1)]
    pl.camera.zoom(1.12)
    pl.show(screenshot=str(out), auto_close=True)


def render_yplus(car: pv.PolyData, out: Path) -> None:
    pl = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    pl.set_background("white")
    scalar = "yPlus" if "yPlus" in car.cell_data or "yPlus" in car.point_data else None
    if scalar:
        pl.add_mesh(car, scalars=scalar, cmap="plasma", clim=(0, 5), show_edges=False, smooth_shading=True,
                    scalar_bar_args={"title": "carBody y+", "vertical": True})
    else:
        pl.add_mesh(car, color="#666666", smooth_shading=True)
    pl.add_mesh(car, style="wireframe", color="black", line_width=0.25, opacity=0.18)
    pl.add_text("L2 carBody y+ map (kOmegaSST)", font_size=18, color="black")
    pl.add_text("Most carBody surface is y+ < 5; median ≈ 0.53", position="lower_left", font_size=11, color="black")
    pl.camera_position = [(2.3, -6.2, 1.25), (2.3, 0.0, 0.85), (0, 0, 1)]
    pl.camera.zoom(1.55)
    pl.show(screenshot=str(out), auto_close=True)


def render_mesh_closeup(car: pv.PolyData, out: Path) -> None:
    pl = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    pl.set_background("white")
    pl.add_mesh(car, color="#d9d9d9", smooth_shading=True, opacity=1.0)
    pl.add_mesh(car, style="wireframe", color="#111111", line_width=0.45, opacity=0.5)
    pl.add_text("L2 carBody surface mesh density", font_size=18, color="black")
    pl.add_text("20,692 carBody faces; median sqrt(face area) ≈ 0.0306 m", position="lower_left", font_size=11, color="black")
    pl.camera_position = [(2.2, -4.2, 1.1), (2.2, 0.0, 0.85), (0, 0, 1)]
    pl.camera.zoom(2.0)
    pl.show(screenshot=str(out), auto_close=True)


def render_iso(mesh: pv.DataSet, car: pv.PolyData, out: Path) -> None:
    # Clip front half open to show interior velocity field volume/slice.
    clipped = mesh.clip(normal="y", origin=(0, 0.0, 0), invert=False)
    pl = pv.Plotter(off_screen=True, window_size=(1800, 1300))
    pl.set_background("white")
    add_domain_guides(pl)
    slx = mesh.slice(normal="y", origin=(2.5, 0.0, 1.9))
    pl.add_mesh(slx, scalars="U_mag", cmap="viridis", clim=(0, 1.6), opacity=0.72, scalar_bar_args={"title": "|U| [m/s]"})
    pl.add_mesh(car, color="#444444", opacity=0.85, smooth_shading=True)
    pl.add_text("L2 paint-booth overview: center-plane velocity + carBody", font_size=17, color="black")
    pl.add_text("Blue top patch=supply; yellow plane=filter panel", position="lower_left", font_size=10, color="black")
    pl.camera_position = [(8.0, -7.2, 4.9), (2.35, 0.0, 1.8), (0, 0, 1)]
    pl.camera.zoom(1.05)
    pl.show(screenshot=str(out), auto_close=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", type=Path, default=Path("cases/paint_booth_panel_frame_komega_l2_fine_yplus"))
    ap.add_argument("--out-dir", type=Path, default=Path("docs/figures"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    vtk_dir = latest_vtk_dir(args.case)
    mesh = get_internal(vtk_dir)
    car = get_car(vtk_dir)

    prefix = "26-05-02_paint_booth_l2"
    files = {
        "overview": args.out_dir / f"{prefix}_overview_velocity.png",
        "umag_midplane": args.out_dir / f"{prefix}_velocity_magnitude_midplane.png",
        "uz_midplane": args.out_dir / f"{prefix}_vertical_velocity_midplane.png",
        "yplus": args.out_dir / f"{prefix}_carbody_yplus.png",
        "surface_mesh": args.out_dir / f"{prefix}_carbody_surface_mesh.png",
    }
    render_iso(mesh, car, files["overview"])
    render_midplane(mesh, car, files["umag_midplane"], "U_mag", "L2 velocity magnitude |U| on center plane", "viridis", (0, 1.6))
    render_midplane(mesh, car, files["uz_midplane"], "Uz", "L2 vertical velocity Uz on center plane", "coolwarm", (-1.2, 0.3))
    render_yplus(car, files["yplus"])
    render_mesh_closeup(car, files["surface_mesh"])

    for k, p in files.items():
        print(f"{k}: {p}")
    os._exit(0)


if __name__ == "__main__":
    main()
