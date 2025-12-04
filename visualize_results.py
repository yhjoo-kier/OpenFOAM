#!/usr/bin/env python3
"""
OpenFOAM 2D Channel Flow Results Visualization
Generates velocity field, pressure field, and centerline velocity plots
"""

import matplotlib
matplotlib.use('Agg')  # Headless backend
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# Disable interactive plotting for headless environment
pv.OFF_SCREEN = True

# Path to VTK results (final time step)
VTK_PATH = '/home/user/OpenFOAM/channelFlow/VTK/channelFlow_1000/internal.vtu'
OUTPUT_DIR = '/mnt/user-data/outputs'

def load_mesh():
    """Load the VTK mesh and convert cell data to point data"""
    mesh = pv.read(VTK_PATH)
    print(f"Mesh loaded: {mesh.n_cells} cells, {mesh.n_points} points")
    print(f"Available arrays: {mesh.array_names}")

    # Convert cell data to point data for smooth visualization
    mesh = mesh.cell_data_to_point_data()
    return mesh

def plot_velocity_field(mesh):
    """Plot velocity magnitude contour"""
    print("Creating velocity field plot...")

    # Calculate velocity magnitude
    U = mesh['U']
    U_mag = np.linalg.norm(U, axis=1)
    mesh['U_magnitude'] = U_mag

    # Create a slice at z=0.005 (middle of the 2D domain)
    slice_mesh = mesh.slice(normal='z', origin=[0.5, 0.05, 0.005])

    # Create figure with matplotlib
    fig, ax = plt.subplots(figsize=(14, 4))

    # Extract coordinates and values
    points = slice_mesh.points
    x = points[:, 0]
    y = points[:, 1]
    U_mag_slice = slice_mesh['U_magnitude']

    # Create scatter plot (works better with unstructured data)
    scatter = ax.scatter(x, y, c=U_mag_slice, cmap='jet', s=5, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax, label='Velocity Magnitude [m/s]')

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title('2D Channel Flow - Velocity Magnitude (t = 1s)')
    ax.set_aspect('equal')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 0.1])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/velocity_field.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/velocity_field.png")

def plot_pressure_field(mesh):
    """Plot pressure distribution contour"""
    print("Creating pressure field plot...")

    # Create a slice at z=0.005 (middle of the 2D domain)
    slice_mesh = mesh.slice(normal='z', origin=[0.5, 0.05, 0.005])

    # Create figure with matplotlib
    fig, ax = plt.subplots(figsize=(14, 4))

    # Extract coordinates and values
    points = slice_mesh.points
    x = points[:, 0]
    y = points[:, 1]
    p = slice_mesh['p']

    # Create scatter plot
    scatter = ax.scatter(x, y, c=p, cmap='coolwarm', s=5, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax, label='Pressure [m²/s²]')

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title('2D Channel Flow - Pressure Distribution (t = 1s)')
    ax.set_aspect('equal')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 0.1])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/pressure_field.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/pressure_field.png")

def plot_centerline_velocity(mesh):
    """Plot velocity profile along channel centerline"""
    print("Creating centerline velocity plot...")

    # Create a line along the channel centerline (y = 0.05, z = 0.005)
    n_points = 200
    start = [0.0, 0.05, 0.005]
    end = [1.0, 0.05, 0.005]

    line = pv.Line(start, end, resolution=n_points-1)

    # Sample the mesh along the line
    sampled = line.sample(mesh)

    # Extract data
    x_coords = sampled.points[:, 0]
    U = sampled['U']
    Ux = U[:, 0]  # x-component of velocity

    # Sort by x coordinate
    sort_idx = np.argsort(x_coords)
    x_sorted = x_coords[sort_idx]
    Ux_sorted = Ux[sort_idx]

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(x_sorted, Ux_sorted, 'b-', linewidth=2, label='Ux (x-component)')
    ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Inlet velocity (1 m/s)')

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Velocity Ux [m/s]')
    ax.set_title('Channel Centerline Velocity Profile (y = 0.05m, t = 1s)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/centerline_velocity.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/centerline_velocity.png")

def plot_velocity_profile_at_outlet(mesh):
    """Additional: Plot velocity profile at outlet (x = 0.9)"""
    print("Creating outlet velocity profile plot...")

    # Create a line at x = 0.9 from y = 0 to y = 0.1
    n_points = 100
    start = [0.9, 0.0, 0.005]
    end = [0.9, 0.1, 0.005]

    line = pv.Line(start, end, resolution=n_points-1)

    # Sample the mesh along the line
    sampled = line.sample(mesh)

    # Extract data
    y_coords = sampled.points[:, 1]
    U = sampled['U']
    Ux = U[:, 0]  # x-component of velocity

    # Sort by y coordinate
    sort_idx = np.argsort(y_coords)
    y_sorted = y_coords[sort_idx]
    Ux_sorted = Ux[sort_idx]

    # Theoretical parabolic profile for developed flow
    h = 0.1  # channel height
    y_theory = np.linspace(0, h, 100)
    # Parabolic velocity profile: U(y) = Umax * (1 - (2y/h - 1)^2)
    # For Umean = 1 m/s, Umax = 1.5 m/s
    Umax = 1.5
    Ux_theory = Umax * (1 - (2*y_theory/h - 1)**2)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(Ux_sorted, y_sorted*1000, 'b-', linewidth=2, label='CFD Result')
    ax.plot(Ux_theory, y_theory*1000, 'r--', linewidth=2, label='Theoretical Parabolic')

    ax.set_xlabel('Velocity Ux [m/s]')
    ax.set_ylabel('Y [mm]')
    ax.set_title('Velocity Profile at x = 0.9m (Near Outlet)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 100])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/outlet_velocity_profile.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/outlet_velocity_profile.png")

def main():
    print("=" * 60)
    print("OpenFOAM Channel Flow Results Visualization")
    print("=" * 60)

    # Load mesh
    mesh = load_mesh()

    # Generate all plots
    plot_velocity_field(mesh)
    plot_pressure_field(mesh)
    plot_centerline_velocity(mesh)
    plot_velocity_profile_at_outlet(mesh)

    print("=" * 60)
    print("All visualizations completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
