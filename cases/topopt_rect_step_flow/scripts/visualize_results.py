#!/usr/bin/env python3
"""
Visualize TopOpt channel-flow results (velocity/pressure) using PyVista.

Steps:
  1) foamToVTK (latestTime recommended)
  2) Read VTK/internal.vtu
  3) Slice at mid-plane z = (zmin+zmax)/2
  4) Save images to images/
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

try:
    import pyvista as pv
except ImportError as e:
    print(f"Required package not found: {e}")
    print("Install with: pip install pyvista numpy")
    sys.exit(1)

pv.OFF_SCREEN = True


CASE_DIR = Path(__file__).resolve().parent.parent
IMAGES_DIR = CASE_DIR / "images"
IMAGES_DIR.mkdir(exist_ok=True)


def _run_foam_to_vtk() -> None:
    cmd = f"cd {CASE_DIR} && foamToVTK -latestTime > /dev/null 2>&1"
    os.system(cmd)


def _find_latest_internal_vtu() -> Path | None:
    vtk_dir = CASE_DIR / "VTK"
    if not vtk_dir.exists():
        return None

    best = None
    best_time = None

    for item in vtk_dir.iterdir():
        if not item.is_dir():
            continue
        # Directory name typically: <caseName>_<time>
        try:
            t = float(item.name.split("_")[-1])
        except ValueError:
            continue
        candidate = item / "internal.vtu"
        if candidate.exists() and (best_time is None or t > best_time):
            best = candidate
            best_time = t

    return best


def _slice_midplane(mesh: pv.DataSet) -> pv.DataSet:
    bounds = mesh.bounds
    z_center = 0.5 * (bounds[4] + bounds[5])
    return mesh.slice(normal="z", origin=(0.0, 0.0, z_center))


def main() -> int:
    print("Converting to VTK (latestTime)...")
    _run_foam_to_vtk()

    vtu = _find_latest_internal_vtu()
    if vtu is None:
        print("No VTK internal.vtu found. Did foamToVTK succeed?")
        return 2

    print(f"Reading: {vtu}")
    mesh = pv.read(str(vtu))

    # Prefer point data for plotting
    if mesh.n_arrays > 0:
        mesh = mesh.cell_data_to_point_data()

    mesh_slice = _slice_midplane(mesh)

    # Geometry image (edges)
    print("Saving geometry slice...")
    plotter = pv.Plotter(off_screen=True, window_size=(1600, 500))
    plotter.add_mesh(mesh_slice, color="white", show_edges=True, opacity=1.0)
    plotter.view_xy()
    plotter.camera.zoom(1.2)
    plotter.set_background("black")
    geom_path = IMAGES_DIR / "geometry.png"
    plotter.screenshot(str(geom_path))
    plotter.close()
    print(f"Saved: {geom_path}")

    # Velocity magnitude
    if "U" in mesh_slice.array_names:
        print("Saving velocity magnitude slice...")
        U = mesh_slice["U"]
        mesh_slice["U_mag"] = np.linalg.norm(U, axis=1)

        plotter = pv.Plotter(off_screen=True, window_size=(1600, 500))
        plotter.add_mesh(
            mesh_slice,
            scalars="U_mag",
            cmap="plasma",
            scalar_bar_args={"title": "Velocity [m/s]"},
        )
        plotter.view_xy()
        plotter.camera.zoom(1.2)
        vel_path = IMAGES_DIR / "velocity_field.png"
        plotter.screenshot(str(vel_path))
        plotter.close()
        print(f"Saved: {vel_path}")
    else:
        print("Field 'U' not found in VTK output.")

    # Pressure
    if "p" in mesh_slice.array_names:
        print("Saving pressure slice...")
        plotter = pv.Plotter(off_screen=True, window_size=(1600, 500))
        plotter.add_mesh(
            mesh_slice,
            scalars="p",
            cmap="coolwarm",
            scalar_bar_args={"title": "Pressure [Pa]"},
        )
        plotter.view_xy()
        plotter.camera.zoom(1.2)
        p_path = IMAGES_DIR / "pressure_field.png"
        plotter.screenshot(str(p_path))
        plotter.close()
        print(f"Saved: {p_path}")
    else:
        print("Field 'p' not found in VTK output.")

    print("Visualization complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


