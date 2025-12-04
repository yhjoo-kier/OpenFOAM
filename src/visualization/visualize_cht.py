#!/usr/bin/env python3
"""
CHT Results Visualization using PyVista
Visualizes temperature and velocity fields from OpenFOAM CHT simulation
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import os

pv.OFF_SCREEN = True


def load_openfoam_case(case_dir: str, time_dir: str = None, region: str = None):
    """
    Load OpenFOAM results as PyVista mesh.

    Args:
        case_dir: Path to OpenFOAM case
        time_dir: Time directory (e.g., '100'). If None, uses latest.
        region: Region name for multi-region cases
    """
    # Find VTK directory
    vtk_base = os.path.join(case_dir, "VTK")

    if region:
        vtk_base = os.path.join(vtk_base, region)

    if not os.path.exists(vtk_base):
        raise FileNotFoundError(f"VTK directory not found: {vtk_base}")

    # Find time directories
    time_dirs = []
    for d in os.listdir(vtk_base):
        full_path = os.path.join(vtk_base, d)
        if os.path.isdir(full_path) and d.startswith(region if region else ""):
            time_dirs.append(d)

    if not time_dirs:
        # Look for .vtu files directly
        vtu_files = [f for f in os.listdir(vtk_base) if f.endswith('.vtu')]
        if vtu_files:
            # Use the latest one
            vtu_files.sort()
            return pv.read(os.path.join(vtk_base, vtu_files[-1]))

    time_dirs.sort(key=lambda x: float(x.split('_')[-1]) if '_' in x else 0)

    if time_dir:
        target = f"{region}_{time_dir}" if region else time_dir
        if target in time_dirs:
            vtk_dir = os.path.join(vtk_base, target)
        else:
            vtk_dir = os.path.join(vtk_base, time_dirs[-1])
    else:
        vtk_dir = os.path.join(vtk_base, time_dirs[-1])

    # Find internal.vtu
    internal_vtu = os.path.join(vtk_dir, "internal.vtu")
    if os.path.exists(internal_vtu):
        return pv.read(internal_vtu)

    # Try to find any vtu file
    for f in os.listdir(vtk_dir):
        if f.endswith('.vtu'):
            return pv.read(os.path.join(vtk_dir, f))

    raise FileNotFoundError(f"No VTU files found in {vtk_dir}")


def visualize_temperature_field(mesh, output_file: str, title: str = "Temperature Field"):
    """Create temperature contour plot."""
    print(f"Creating temperature field visualization...")

    if 'T' not in mesh.array_names:
        print("Warning: Temperature field 'T' not found")
        return

    # Convert cell data to point data if needed
    if mesh.n_arrays > 0 and 'T' in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()

    # Create slice through center
    center = mesh.center
    slice_mesh = mesh.slice(normal='y', origin=center)

    fig, ax = plt.subplots(figsize=(14, 8))

    points = slice_mesh.points
    x = points[:, 0] * 1000  # Convert to mm
    z = points[:, 2] * 1000
    T = slice_mesh['T']

    scatter = ax.scatter(x, z, c=T, cmap='hot', s=3, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax, label='Temperature [K]')

    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Z [mm]')
    ax.set_title(title)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def visualize_velocity_field(mesh, output_file: str, title: str = "Velocity Field"):
    """Create velocity magnitude contour plot."""
    print(f"Creating velocity field visualization...")

    if 'U' not in mesh.array_names:
        print("Warning: Velocity field 'U' not found")
        return

    # Convert cell data to point data if needed
    if 'U' in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()

    # Calculate velocity magnitude
    U = mesh['U']
    U_mag = np.linalg.norm(U, axis=1)
    mesh['U_magnitude'] = U_mag

    # Create slice through center
    center = mesh.center
    slice_mesh = mesh.slice(normal='y', origin=center)

    fig, ax = plt.subplots(figsize=(14, 8))

    points = slice_mesh.points
    x = points[:, 0] * 1000
    z = points[:, 2] * 1000
    U_mag_slice = slice_mesh['U_magnitude']

    scatter = ax.scatter(x, z, c=U_mag_slice, cmap='jet', s=3, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax, label='Velocity Magnitude [m/s]')

    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Z [mm]')
    ax.set_title(title)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def visualize_combined_cht(fluid_mesh, solid_mesh, output_dir: str):
    """Create combined visualization for CHT results."""

    os.makedirs(output_dir, exist_ok=True)

    # Temperature distribution - combined view
    print("Creating combined temperature visualization...")

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Fluid temperature
    if 'T' in fluid_mesh.cell_data:
        fluid_mesh = fluid_mesh.cell_data_to_point_data()

    center_y = fluid_mesh.center[1]
    fluid_slice = fluid_mesh.slice(normal='y', origin=[0, center_y, 0])

    points = fluid_slice.points
    x = points[:, 0] * 1000
    z = points[:, 2] * 1000
    T = fluid_slice['T']

    scatter1 = axes[0].scatter(x, z, c=T, cmap='coolwarm', s=3, edgecolors='none')
    plt.colorbar(scatter1, ax=axes[0], label='Temperature [K]')
    axes[0].set_xlabel('X [mm]')
    axes[0].set_ylabel('Z [mm]')
    axes[0].set_title('Fluid Temperature')
    axes[0].set_aspect('equal')

    # Solid temperature
    if 'T' in solid_mesh.cell_data:
        solid_mesh = solid_mesh.cell_data_to_point_data()

    solid_slice = solid_mesh.slice(normal='y', origin=[0, center_y, 0])

    points = solid_slice.points
    x = points[:, 0] * 1000
    z = points[:, 2] * 1000
    T = solid_slice['T']

    scatter2 = axes[1].scatter(x, z, c=T, cmap='hot', s=5, edgecolors='none')
    plt.colorbar(scatter2, ax=axes[1], label='Temperature [K]')
    axes[1].set_xlabel('X [mm]')
    axes[1].set_ylabel('Z [mm]')
    axes[1].set_title('Solid (Heatsink) Temperature')
    axes[1].set_aspect('equal')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'temperature_combined.png'), dpi=150, bbox_inches='tight')
    plt.close()

    # Velocity field
    visualize_velocity_field(
        fluid_mesh,
        os.path.join(output_dir, 'velocity_field.png'),
        'Air Velocity Around Heatsink'
    )

    # Individual temperature plots
    visualize_temperature_field(
        fluid_mesh,
        os.path.join(output_dir, 'temperature_fluid.png'),
        'Air Temperature Distribution'
    )

    visualize_temperature_field(
        solid_mesh,
        os.path.join(output_dir, 'temperature_solid.png'),
        'Heatsink Temperature Distribution'
    )

    print(f"All visualizations saved to: {output_dir}")


def create_heatsink_3d_view(solid_mesh, output_file: str):
    """Create 3D view of heatsink temperature distribution."""
    print("Creating 3D heatsink visualization...")

    if 'T' in solid_mesh.cell_data:
        solid_mesh = solid_mesh.cell_data_to_point_data()

    # Create plotter
    plotter = pv.Plotter(off_screen=True, window_size=[1200, 900])

    # Add mesh with temperature coloring
    plotter.add_mesh(
        solid_mesh,
        scalars='T',
        cmap='hot',
        scalar_bar_args={'title': 'Temperature [K]'}
    )

    # Set camera position
    plotter.camera_position = 'iso'
    plotter.add_axes()

    # Save screenshot
    plotter.screenshot(output_file)
    plotter.close()
    print(f"Saved: {output_file}")


def visualize_cht_results(case_dir: str, output_dir: str, time_dir: str = None):
    """
    Main function to visualize CHT results.

    Args:
        case_dir: OpenFOAM case directory
        output_dir: Output directory for images
        time_dir: Specific time to visualize (default: latest)
    """
    print("=" * 60)
    print("CHT Results Visualization")
    print("=" * 60)

    os.makedirs(output_dir, exist_ok=True)

    try:
        # Load fluid region
        print("Loading fluid region...")
        fluid_mesh = load_openfoam_case(case_dir, time_dir, region="fluid")
        print(f"  Fluid mesh: {fluid_mesh.n_cells} cells, {fluid_mesh.n_points} points")

        # Load solid region
        print("Loading solid region...")
        solid_mesh = load_openfoam_case(case_dir, time_dir, region="solid")
        print(f"  Solid mesh: {solid_mesh.n_cells} cells, {solid_mesh.n_points} points")

        # Create visualizations
        visualize_combined_cht(fluid_mesh, solid_mesh, output_dir)

        # 3D view
        create_heatsink_3d_view(solid_mesh, os.path.join(output_dir, 'heatsink_3d.png'))

    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Make sure foamToVTK has been run on the case")
        return False

    print("=" * 60)
    print("Visualization complete!")
    print("=" * 60)

    return True


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Visualize CHT results')
    parser.add_argument('case_dir', help='OpenFOAM case directory')
    parser.add_argument('--output', '-o', default='./results', help='Output directory')
    parser.add_argument('--time', '-t', default=None, help='Time directory to visualize')

    args = parser.parse_args()

    visualize_cht_results(args.case_dir, args.output, args.time)


if __name__ == "__main__":
    main()
