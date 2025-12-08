#!/usr/bin/env python3
"""
Visualize channelFlow simulation results.
Generates velocity and pressure field plots.
"""

import os
import sys

try:
    import pyvista as pv
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
except ImportError as e:
    print(f"Required package not found: {e}")
    print("Install with: pip install pyvista numpy matplotlib")
    sys.exit(1)

# Use non-interactive backend
import matplotlib
matplotlib.use('Agg')
pv.OFF_SCREEN = True


def get_case_dir():
    """Get the case directory (parent of scripts/)."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(script_dir)


def main():
    case_dir = get_case_dir()
    images_dir = os.path.join(case_dir, "images")
    os.makedirs(images_dir, exist_ok=True)

    # Convert to VTK first
    print("Converting to VTK format...")
    os.system(f"cd {case_dir} && foamToVTK > /dev/null 2>&1")

    vtk_dir = os.path.join(case_dir, "VTK")
    if not os.path.exists(vtk_dir):
        print("VTK conversion failed")
        return

    # Find the latest internal.vtu file
    vtk_path = None
    max_time = -1

    for item in os.listdir(vtk_dir):
        item_path = os.path.join(vtk_dir, item)
        if os.path.isdir(item_path):
            try:
                time_str = item.split('_')[-1]
                t = float(time_str)
                internal_vtu = os.path.join(item_path, "internal.vtu")
                if os.path.exists(internal_vtu) and t > max_time:
                    max_time = t
                    vtk_path = internal_vtu
            except (ValueError, IndexError):
                continue

    if vtk_path is None:
        print("No VTK internal.vtu file found")
        return

    print(f"Reading VTK file: {vtk_path}")

    try:
        mesh = pv.read(vtk_path)
    except Exception as e:
        print(f"Error reading VTK file: {e}")
        return

    # Convert cell data to point data if needed
    if mesh.n_arrays > 0:
        mesh = mesh.cell_data_to_point_data()

    # Create velocity magnitude plot
    print("Generating velocity field plot...")

    if 'U' in mesh.array_names:
        U = mesh['U']
        U_mag = np.linalg.norm(U, axis=1)
        mesh['U_magnitude'] = U_mag

        # Use pyvista for plotting
        plotter = pv.Plotter(off_screen=True)
        plotter.add_mesh(mesh, scalars='U_magnitude', cmap='jet',
                        scalar_bar_args={'title': 'Velocity [m/s]'})
        plotter.view_xy()
        plotter.camera.zoom(1.2)

        velocity_path = os.path.join(images_dir, 'velocity_field.png')
        plotter.screenshot(velocity_path)
        plotter.close()
        print(f"Saved: {velocity_path}")

    # Create pressure plot
    print("Generating pressure field plot...")

    if 'p' in mesh.array_names:
        plotter = pv.Plotter(off_screen=True)
        plotter.add_mesh(mesh, scalars='p', cmap='coolwarm',
                        scalar_bar_args={'title': 'Pressure [Pa]'})
        plotter.view_xy()
        plotter.camera.zoom(1.2)

        pressure_path = os.path.join(images_dir, 'pressure_field.png')
        plotter.screenshot(pressure_path)
        plotter.close()
        print(f"Saved: {pressure_path}")

    print("Visualization complete!")


if __name__ == "__main__":
    main()
