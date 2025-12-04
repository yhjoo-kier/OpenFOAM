#!/usr/bin/env python3
"""
Heatsink Temperature Visualization (matplotlib-based)
Includes domain geometry, cross-section contours, and temperature profiles.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.tri import Triangulation
import numpy as np
import pyvista as pv
import os

pv.OFF_SCREEN = True


def find_latest_vtk(case_dir: str) -> str:
    """Find the latest VTK timestep directory."""
    vtk_base = os.path.join(case_dir, "VTK")
    vtk_dirs = []

    for d in os.listdir(vtk_base):
        full_path = os.path.join(vtk_base, d)
        if os.path.isdir(full_path):
            # Extract timestep number from directory name
            parts = d.split('_')
            try:
                timestep = int(parts[-1])
                vtk_dirs.append((timestep, d))
            except ValueError:
                continue

    if not vtk_dirs:
        raise FileNotFoundError("No VTK directories found!")

    # Sort by timestep and return latest
    vtk_dirs.sort(key=lambda x: x[0])
    return os.path.join(vtk_base, vtk_dirs[-1][1], "internal.vtu")


def create_triangulation(points_2d, scalar_data):
    """Create triangulation for contour plotting."""
    from scipy.spatial import Delaunay

    # Remove duplicate points
    unique_points, indices = np.unique(points_2d, axis=0, return_inverse=True)
    unique_scalars = np.zeros(len(unique_points))
    np.add.at(unique_scalars, indices, scalar_data)
    counts = np.bincount(indices).astype(float)
    unique_scalars /= counts

    if len(unique_points) < 4:
        return None, None, None

    try:
        tri = Delaunay(unique_points)
        return unique_points[:, 0], unique_points[:, 1], Triangulation(
            unique_points[:, 0], unique_points[:, 1], tri.simplices
        ), unique_scalars
    except Exception:
        return None, None, None, None


def visualize_heatsink_results(case_dir: str, output_dir: str):
    """Visualize heatsink thermal analysis results."""

    os.makedirs(output_dir, exist_ok=True)

    # Find and load VTK data
    latest_vtk = find_latest_vtk(case_dir)
    print(f"Loading: {latest_vtk}")

    mesh = pv.read(latest_vtk)
    print(f"Mesh: {mesh.n_cells} cells, {mesh.n_points} points")
    print(f"Arrays: {mesh.array_names}")
    print(f"Bounds: {mesh.bounds}")

    # Convert cell data to point data
    if 'T' in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()

    T = mesh['T']
    T_min, T_max = T.min(), T.max()
    print(f"Temperature range: {T_min:.2f} K to {T_max:.2f} K")
    print(f"Temperature rise: {T_max - T_min:.2f} K")

    # Get mesh bounds in mm
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    xmin_mm, xmax_mm = xmin * 1000, xmax * 1000
    ymin_mm, ymax_mm = ymin * 1000, ymax * 1000
    zmin_mm, zmax_mm = zmin * 1000, zmax * 1000

    # ==================================================
    # 1. Domain Geometry Visualization (3D wireframe)
    # ==================================================
    print("\nCreating domain geometry visualization...")

    fig = plt.figure(figsize=(14, 5))

    # 3D isometric view of domain
    ax1 = fig.add_subplot(131, projection='3d')

    # Draw box edges
    vertices = np.array([
        [xmin_mm, ymin_mm, zmin_mm],
        [xmax_mm, ymin_mm, zmin_mm],
        [xmax_mm, ymax_mm, zmin_mm],
        [xmin_mm, ymax_mm, zmin_mm],
        [xmin_mm, ymin_mm, zmax_mm],
        [xmax_mm, ymin_mm, zmax_mm],
        [xmax_mm, ymax_mm, zmax_mm],
        [xmin_mm, ymax_mm, zmax_mm],
    ])

    # Define edges
    edges = [
        [0, 1], [1, 2], [2, 3], [3, 0],  # bottom
        [4, 5], [5, 6], [6, 7], [7, 4],  # top
        [0, 4], [1, 5], [2, 6], [3, 7],  # vertical
    ]

    for edge in edges:
        pts = vertices[edge]
        ax1.plot3D(pts[:, 0], pts[:, 1], pts[:, 2], 'b-', linewidth=2)

    # Add face labels
    ax1.text(25, 25, zmin_mm - 2, 'Bottom\n(Heat Source)', ha='center', fontsize=8, color='red')
    ax1.text(25, 25, zmax_mm + 2, 'Top\n(Fixed 300K)', ha='center', fontsize=8, color='blue')

    ax1.set_xlabel('X [mm]')
    ax1.set_ylabel('Y [mm]')
    ax1.set_zlabel('Z [mm]')
    ax1.set_title('Domain Geometry\n(50 x 50 x 25 mm)')

    # 2D projections
    ax2 = fig.add_subplot(132)
    rect = plt.Rectangle((xmin_mm, zmin_mm), xmax_mm - xmin_mm, zmax_mm - zmin_mm,
                          fill=False, edgecolor='blue', linewidth=2)
    ax2.add_patch(rect)
    ax2.set_xlim(xmin_mm - 5, xmax_mm + 5)
    ax2.set_ylim(zmin_mm - 5, zmax_mm + 5)
    ax2.set_xlabel('X [mm]')
    ax2.set_ylabel('Z [mm]')
    ax2.set_title('Side View (XZ)')
    ax2.set_aspect('equal')
    ax2.axhline(y=zmin_mm, color='red', linestyle='--', label='Heat source')
    ax2.axhline(y=zmax_mm, color='blue', linestyle='--', label='Fixed T=300K')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    ax3 = fig.add_subplot(133)
    rect = plt.Rectangle((xmin_mm, ymin_mm), xmax_mm - xmin_mm, ymax_mm - ymin_mm,
                          fill=False, edgecolor='blue', linewidth=2)
    ax3.add_patch(rect)
    ax3.set_xlim(xmin_mm - 5, xmax_mm + 5)
    ax3.set_ylim(ymin_mm - 5, ymax_mm + 5)
    ax3.set_xlabel('X [mm]')
    ax3.set_ylabel('Y [mm]')
    ax3.set_title('Top View (XY)')
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)

    plt.suptitle('Heat Sink Domain Geometry', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatsink_geometry.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.join(output_dir, 'heatsink_geometry.png')}")

    # ==================================================
    # 2. Cross-section Temperature Contours (improved)
    # ==================================================
    print("\nCreating cross-section contour views...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    norm = Normalize(vmin=T_min, vmax=T_max)

    # XZ slice (y = center) - Side view
    center_y = (ymin + ymax) / 2
    slice_xz = mesh.slice(normal='y', origin=[0, center_y, 0])
    if slice_xz.n_points > 0:
        points = slice_xz.points
        x = points[:, 0] * 1000
        z = points[:, 2] * 1000
        T_slice = slice_xz['T']

        # Create grid for contour
        xi = np.linspace(x.min(), x.max(), 100)
        zi = np.linspace(z.min(), z.max(), 60)
        Xi, Zi = np.meshgrid(xi, zi)

        from scipy.interpolate import griddata
        Ti = griddata((x, z), T_slice, (Xi, Zi), method='linear')

        contour = axes[0].contourf(Xi, Zi, Ti, levels=30, cmap='hot', norm=norm)
        axes[0].contour(Xi, Zi, Ti, levels=10, colors='k', linewidths=0.3, alpha=0.5)
        plt.colorbar(contour, ax=axes[0], label='Temperature [K]')
        axes[0].set_xlabel('X [mm]')
        axes[0].set_ylabel('Z [mm]')
        axes[0].set_title(f'XZ Cross-section (Y={center_y*1000:.1f}mm)\nSide View')
    axes[0].set_aspect('equal')

    # XY slice (z = near bottom) - Base temperature
    z_slice = zmin + (zmax - zmin) * 0.1
    slice_xy = mesh.slice(normal='z', origin=[0, 0, z_slice])
    if slice_xy.n_points > 0:
        points = slice_xy.points
        x = points[:, 0] * 1000
        y = points[:, 1] * 1000
        T_slice = slice_xy['T']

        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        Xi, Yi = np.meshgrid(xi, yi)

        from scipy.interpolate import griddata
        Ti = griddata((x, y), T_slice, (Xi, Yi), method='linear')

        contour = axes[1].contourf(Xi, Yi, Ti, levels=30, cmap='hot', norm=norm)
        axes[1].contour(Xi, Yi, Ti, levels=10, colors='k', linewidths=0.3, alpha=0.5)
        plt.colorbar(contour, ax=axes[1], label='Temperature [K]')
        axes[1].set_xlabel('X [mm]')
        axes[1].set_ylabel('Y [mm]')
        axes[1].set_title(f'XY Cross-section (Z={z_slice*1000:.1f}mm)\nNear Bottom')
    axes[1].set_aspect('equal')

    # YZ slice (x = center) - Front view
    center_x = (xmin + xmax) / 2
    slice_yz = mesh.slice(normal='x', origin=[center_x, 0, 0])
    if slice_yz.n_points > 0:
        points = slice_yz.points
        y = points[:, 1] * 1000
        z = points[:, 2] * 1000
        T_slice = slice_yz['T']

        yi = np.linspace(y.min(), y.max(), 100)
        zi = np.linspace(z.min(), z.max(), 60)
        Yi, Zi = np.meshgrid(yi, zi)

        from scipy.interpolate import griddata
        Ti = griddata((y, z), T_slice, (Yi, Zi), method='linear')

        contour = axes[2].contourf(Yi, Zi, Ti, levels=30, cmap='hot', norm=norm)
        axes[2].contour(Yi, Zi, Ti, levels=10, colors='k', linewidths=0.3, alpha=0.5)
        plt.colorbar(contour, ax=axes[2], label='Temperature [K]')
        axes[2].set_xlabel('Y [mm]')
        axes[2].set_ylabel('Z [mm]')
        axes[2].set_title(f'YZ Cross-section (X={center_x*1000:.1f}mm)\nFront View')
    axes[2].set_aspect('equal')

    plt.suptitle(f'Heat Sink Temperature Distribution\n'
                 f'Heat Flux: 500,000 W/m² (bottom), Fixed T: 300 K (top)\n'
                 f'Temperature Range: {T_min:.1f} K - {T_max:.1f} K',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatsink_cross_sections.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.join(output_dir, 'heatsink_cross_sections.png')}")

    # ==================================================
    # 3. Temperature Profile Along Height
    # ==================================================
    print("Creating temperature profile...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Profile along z at center
    n_sample = 100
    line = pv.Line([center_x, center_y, zmin], [center_x, center_y, zmax], resolution=n_sample-1)
    sampled = line.sample(mesh)

    z_coords = sampled.points[:, 2] * 1000
    T_line = sampled['T']

    sort_idx = np.argsort(z_coords)
    z_sorted = z_coords[sort_idx]
    T_sorted = T_line[sort_idx]

    axes[0].plot(z_sorted, T_sorted, 'r-', linewidth=2, marker='o', markersize=3, label='Temperature')
    axes[0].axhline(y=300, color='b', linestyle='--', alpha=0.7, label='Top BC (300 K)')
    axes[0].fill_between([zmin_mm, zmax_mm], 300, T_max, alpha=0.1, color='red')

    axes[0].set_xlabel('Height Z [mm]')
    axes[0].set_ylabel('Temperature [K]')
    axes[0].set_title('Temperature Profile Along Height (Center)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xlim([zmin_mm, zmax_mm])

    # Analytical comparison for 1D heat conduction
    # T(z) = T_top + (q/k) * (L - z) where L is height, q is heat flux, k is conductivity
    L = zmax - zmin  # height in meters
    q = 500000  # W/m² (high heat flux for visible gradient)
    k = 205   # W/(m·K) for aluminum

    z_analytical = np.linspace(zmin, zmax, 100) * 1000  # mm
    T_analytical = 300 + (q / k) * (L - (z_analytical / 1000 - zmin))

    axes[0].plot(z_analytical, T_analytical, 'g--', linewidth=2, alpha=0.7,
                 label=f'Analytical (1D): ΔT = q·L/k = {q*L/k:.1f} K')
    axes[0].legend()

    # Temperature distribution histogram
    axes[1].hist(T, bins=50, color='orange', edgecolor='black', alpha=0.7)
    axes[1].axvline(x=T_min, color='blue', linestyle='--', linewidth=2, label=f'Min: {T_min:.1f} K')
    axes[1].axvline(x=T_max, color='red', linestyle='--', linewidth=2, label=f'Max: {T_max:.1f} K')
    axes[1].axvline(x=T.mean(), color='green', linestyle='-', linewidth=2, label=f'Mean: {T.mean():.1f} K')

    axes[1].set_xlabel('Temperature [K]')
    axes[1].set_ylabel('Frequency (# of points)')
    axes[1].set_title('Temperature Distribution Histogram')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatsink_temperature_profile.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.join(output_dir, 'heatsink_temperature_profile.png')}")

    # ==================================================
    # 4. Bottom Surface Temperature (Heat Source)
    # ==================================================
    print("Creating bottom surface temperature map...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Bottom slice
    slice_bottom = mesh.slice(normal='z', origin=[0, 0, zmin + 0.0005])
    if slice_bottom.n_points > 0:
        points = slice_bottom.points
        x = points[:, 0] * 1000
        y = points[:, 1] * 1000
        T_bottom = slice_bottom['T']

        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        Xi, Yi = np.meshgrid(xi, yi)

        from scipy.interpolate import griddata
        Ti = griddata((x, y), T_bottom, (Xi, Yi), method='linear')

        contour = axes[0].contourf(Xi, Yi, Ti, levels=30, cmap='hot')
        axes[0].contour(Xi, Yi, Ti, levels=10, colors='k', linewidths=0.3, alpha=0.5)
        cbar = plt.colorbar(contour, ax=axes[0], label='Temperature [K]')

        axes[0].set_xlabel('X [mm]')
        axes[0].set_ylabel('Y [mm]')
        axes[0].set_title(f'Bottom Surface Temperature (Heat Source)\nMax: {T_bottom.max():.1f} K')
        axes[0].set_aspect('equal')

    # Top slice
    slice_top = mesh.slice(normal='z', origin=[0, 0, zmax - 0.0005])
    if slice_top.n_points > 0:
        points = slice_top.points
        x = points[:, 0] * 1000
        y = points[:, 1] * 1000
        T_top = slice_top['T']

        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        Xi, Yi = np.meshgrid(xi, yi)

        from scipy.interpolate import griddata
        Ti = griddata((x, y), T_top, (Xi, Yi), method='linear')

        contour = axes[1].contourf(Xi, Yi, Ti, levels=30, cmap='hot')
        axes[1].contour(Xi, Yi, Ti, levels=10, colors='k', linewidths=0.3, alpha=0.5)
        cbar = plt.colorbar(contour, ax=axes[1], label='Temperature [K]')

        axes[1].set_xlabel('X [mm]')
        axes[1].set_ylabel('Y [mm]')
        axes[1].set_title(f'Top Surface Temperature (Fixed BC)\nValue: {T_top.mean():.1f} K')
        axes[1].set_aspect('equal')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatsink_surface_temperature.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.join(output_dir, 'heatsink_surface_temperature.png')}")

    # ==================================================
    # 5. Summary Report
    # ==================================================
    # Calculate analytical solution for comparison
    L = zmax - zmin
    q = 500000  # W/m² (high heat flux)
    k = 205
    T_max_analytical = 300 + q * L / k

    print("\n" + "=" * 60)
    print("Heatsink Thermal Analysis Summary")
    print("=" * 60)
    print(f"Geometry:")
    print(f"  Domain: {xmax_mm:.0f} x {ymax_mm:.0f} x {zmax_mm:.0f} mm")
    print(f"  Mesh: {mesh.n_cells} cells, {mesh.n_points} points")
    print(f"\nBoundary Conditions:")
    print(f"  Bottom: Heat flux q = {q:,} W/m²")
    print(f"  Top: Fixed temperature T = 300 K")
    print(f"  Sides: Adiabatic (zero gradient)")
    print(f"\nMaterial: Aluminum")
    print(f"  Thermal conductivity k = 205 W/(m·K)")
    print(f"\nResults:")
    print(f"  Minimum temperature: {T_min:.2f} K ({T_min-273.15:.2f} °C)")
    print(f"  Maximum temperature: {T_max:.2f} K ({T_max-273.15:.2f} °C)")
    print(f"  Mean temperature: {T.mean():.2f} K ({T.mean()-273.15:.2f} °C)")
    print(f"  Temperature rise: {T_max - T_min:.2f} K")
    print(f"\nAnalytical 1D Solution:")
    print(f"  T_max = T_top + q·L/k = 300 + {q}×{L:.4f}/205 = {T_max_analytical:.2f} K")
    print(f"  Simulation T_max: {T_max:.2f} K")
    print(f"  Error: {abs(T_max - T_max_analytical):.2f} K ({abs(T_max - T_max_analytical)/T_max_analytical*100:.2f}%)")
    print(f"\nOutput directory: {output_dir}")
    print("=" * 60)

    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Visualize heatsink thermal analysis results')
    parser.add_argument('case_dir', help='OpenFOAM case directory')
    parser.add_argument('--output', '-o', default='./results', help='Output directory for images')

    args = parser.parse_args()

    visualize_heatsink_results(args.case_dir, args.output)
