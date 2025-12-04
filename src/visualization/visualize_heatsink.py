#!/usr/bin/env python3
"""
Heatsink Temperature Visualization (matplotlib only - no 3D rendering)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import os

pv.OFF_SCREEN = True


def visualize_heatsink_results(case_dir: str, output_dir: str):
    """Visualize heatsink thermal analysis results."""

    os.makedirs(output_dir, exist_ok=True)

    # Find latest VTK directory
    vtk_base = os.path.join(case_dir, "VTK")
    vtk_dirs = sorted([d for d in os.listdir(vtk_base) if os.path.isdir(os.path.join(vtk_base, d))],
                      key=lambda x: int(x.split('_')[-1]) if '_' in x else 0)

    if not vtk_dirs:
        print("No VTK directories found!")
        return

    latest_vtk = os.path.join(vtk_base, vtk_dirs[-1], "internal.vtu")
    print(f"Loading: {latest_vtk}")

    mesh = pv.read(latest_vtk)
    print(f"Mesh: {mesh.n_cells} cells, {mesh.n_points} points")
    print(f"Arrays: {mesh.array_names}")

    # Convert cell data to point data
    if 'T' in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()

    T = mesh['T']
    T_min, T_max = T.min(), T.max()
    print(f"Temperature range: {T_min:.2f} K to {T_max:.2f} K")
    print(f"Temperature rise: {T_max - 300:.2f} K")

    # 1. Cross-section views (main visualization)
    print("\nCreating cross-section views...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # XZ slice (y = center)
    center_y = mesh.center[1]
    slice_xz = mesh.slice(normal='y', origin=[0, center_y, 0])
    points = slice_xz.points
    x, z = points[:, 0] * 1000, points[:, 2] * 1000
    T_slice = slice_xz['T']

    scatter = axes[0].scatter(x, z, c=T_slice, cmap='hot', s=15, edgecolors='none',
                               vmin=T_min, vmax=T_max)
    plt.colorbar(scatter, ax=axes[0], label='Temperature [K]')
    axes[0].set_xlabel('X [mm]')
    axes[0].set_ylabel('Z [mm]')
    axes[0].set_title('XZ Cross-section (Side View)')
    axes[0].set_aspect('equal')

    # XY slice (z = base center)
    slice_xy = mesh.slice(normal='z', origin=[0, 0, 0.0025])
    points = slice_xy.points
    x, y = points[:, 0] * 1000, points[:, 1] * 1000
    T_slice = slice_xy['T']

    scatter = axes[1].scatter(x, y, c=T_slice, cmap='hot', s=15, edgecolors='none',
                               vmin=T_min, vmax=T_max)
    plt.colorbar(scatter, ax=axes[1], label='Temperature [K]')
    axes[1].set_xlabel('X [mm]')
    axes[1].set_ylabel('Y [mm]')
    axes[1].set_title('XY Cross-section (Base, Z=2.5mm)')
    axes[1].set_aspect('equal')

    # YZ slice (x = center)
    center_x = mesh.center[0]
    slice_yz = mesh.slice(normal='x', origin=[center_x, 0, 0])
    points = slice_yz.points
    y, z = points[:, 1] * 1000, points[:, 2] * 1000
    T_slice = slice_yz['T']

    scatter = axes[2].scatter(y, z, c=T_slice, cmap='hot', s=15, edgecolors='none',
                               vmin=T_min, vmax=T_max)
    plt.colorbar(scatter, ax=axes[2], label='Temperature [K]')
    axes[2].set_xlabel('Y [mm]')
    axes[2].set_ylabel('Z [mm]')
    axes[2].set_title('YZ Cross-section (Front View)')
    axes[2].set_aspect('equal')

    plt.suptitle(f'Heatsink Temperature Distribution\n(Heat Flux: 5000 W/m², Ambient: 300 K)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatsink_cross_sections.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.join(output_dir, 'heatsink_cross_sections.png')}")

    # 2. Temperature profile along height
    print("Creating temperature profile...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Profile along z at center
    n_points = 50
    line = pv.Line([0.025, 0.025, 0], [0.025, 0.025, 0.025], resolution=n_points-1)
    sampled = line.sample(mesh)

    z_coords = sampled.points[:, 2] * 1000
    T_line = sampled['T']

    sort_idx = np.argsort(z_coords)
    axes[0].plot(z_coords[sort_idx], T_line[sort_idx], 'r-', linewidth=2, marker='o', markersize=4)

    axes[0].axhline(y=300, color='b', linestyle='--', alpha=0.5, label='Ambient (300 K)')
    axes[0].axvline(x=5, color='gray', linestyle=':', alpha=0.5, label='Base/Fin interface')

    axes[0].set_xlabel('Height Z [mm]')
    axes[0].set_ylabel('Temperature [K]')
    axes[0].set_title('Temperature Profile Along Height (Center)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xlim([0, 25])

    # Temperature distribution histogram
    axes[1].hist(T, bins=50, color='orange', edgecolor='black', alpha=0.7)
    axes[1].axvline(x=T_min, color='blue', linestyle='--', label=f'Min: {T_min:.1f} K')
    axes[1].axvline(x=T_max, color='red', linestyle='--', label=f'Max: {T_max:.1f} K')
    axes[1].axvline(x=T.mean(), color='green', linestyle='-', label=f'Mean: {T.mean():.1f} K')

    axes[1].set_xlabel('Temperature [K]')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Temperature Distribution')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatsink_temperature_profile.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.join(output_dir, 'heatsink_temperature_profile.png')}")

    # 3. Bottom surface temperature (heat source)
    print("Creating bottom surface temperature map...")

    fig, ax = plt.subplots(figsize=(10, 8))

    # Get bottom surface temperature
    slice_bottom = mesh.slice(normal='z', origin=[0, 0, 0.001])
    points = slice_bottom.points
    x, y = points[:, 0] * 1000, points[:, 1] * 1000
    T_bottom = slice_bottom['T']

    scatter = ax.scatter(x, y, c=T_bottom, cmap='hot', s=30, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax, label='Temperature [K]')

    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    ax.set_title(f'Bottom Surface Temperature (Heat Source)\nMax: {T_bottom.max():.1f} K')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatsink_bottom_temperature.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.join(output_dir, 'heatsink_bottom_temperature.png')}")

    # 4. Summary report
    print("\n" + "="*60)
    print("Heatsink Thermal Analysis Summary")
    print("="*60)
    print(f"Geometry:")
    print(f"  Base: 50 x 50 x 5 mm")
    print(f"  Fin:  50 x 10 x 20 mm (centered)")
    print(f"Boundary Conditions:")
    print(f"  Bottom: Heat flux = 5000 W/m²")
    print(f"  Top: Fixed temperature = 300 K")
    print(f"  Sides: Adiabatic (zero gradient)")
    print(f"Material: Aluminum (k = 205 W/m·K)")
    print(f"\nResults:")
    print(f"  Minimum temperature: {T_min:.2f} K ({T_min-273.15:.2f} °C)")
    print(f"  Maximum temperature: {T_max:.2f} K ({T_max-273.15:.2f} °C)")
    print(f"  Mean temperature: {T.mean():.2f} K ({T.mean()-273.15:.2f} °C)")
    print(f"  Temperature rise (max): {T_max - 300:.2f} K")
    print(f"\nOutput directory: {output_dir}")
    print("="*60)

    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('case_dir', help='Case directory')
    parser.add_argument('--output', '-o', default='./results', help='Output directory')

    args = parser.parse_args()

    visualize_heatsink_results(args.case_dir, args.output)
