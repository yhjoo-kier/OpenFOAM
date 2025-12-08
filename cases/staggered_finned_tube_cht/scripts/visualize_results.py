#!/usr/bin/env python3
"""Visualization script for staggered finned tube CHT results using matplotlib."""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Circle
import pyvista as pv

# Get case directory (parent of scripts/)
script_dir = os.path.dirname(os.path.abspath(__file__))
case_dir = os.path.dirname(script_dir)
os.chdir(case_dir)

# Create images directory
images_dir = os.path.join(case_dir, "images")
os.makedirs(images_dir, exist_ok=True)

vtk_dir = os.path.join(case_dir, "VTK")

# Load fluid mesh
fluid_mesh = pv.read(os.path.join(vtk_dir, "fluid/staggered_finned_tube_cht_1000/internal.vtu"))

# Load solid meshes
solid_mesh = pv.read(os.path.join(vtk_dir, "solid/staggered_finned_tube_cht_1000/internal.vtu"))
domain0_mesh = pv.read(os.path.join(vtk_dir, "domain0/staggered_finned_tube_cht_1000/internal.vtu"))
domain1_mesh = pv.read(os.path.join(vtk_dir, "domain1/staggered_finned_tube_cht_1000/internal.vtu"))

print(f"Fluid mesh: {fluid_mesh.n_cells} cells")
print(f"Solid mesh: {solid_mesh.n_cells} cells")
print(f"Domain0 mesh: {domain0_mesh.n_cells} cells")
print(f"Domain1 mesh: {domain1_mesh.n_cells} cells")

# Get mesh bounds
bounds = fluid_mesh.bounds
print(f"Domain bounds: x=[{bounds[0]:.4f}, {bounds[1]:.4f}], y=[{bounds[2]:.4f}, {bounds[3]:.4f}], z=[{bounds[4]:.4f}, {bounds[5]:.4f}]")

# Create z=0 slice for 2D view
z_center = (bounds[4] + bounds[5]) / 2

# Extract cell centers and data for fluid
fluid_centers = fluid_mesh.cell_centers()
fluid_x = fluid_centers.points[:, 0]
fluid_y = fluid_centers.points[:, 1]
fluid_z = fluid_centers.points[:, 2]

# Filter cells near z=0 slice (within tolerance)
z_tol = 0.001
fluid_mask = np.abs(fluid_z - z_center) < z_tol

# Get field data
U = fluid_mesh['U']
U_mag = np.linalg.norm(U, axis=1)
T_fluid = fluid_mesh['T']

# For visualization, use all cells projected to xy plane (thin slice)
# Extract a proper slice
slice_normal = [0, 0, 1]
slice_origin = [0, 0, z_center]

# Use cell centers for scatter plot approach
# Get all cells, sort by z to find middle layer
z_unique = np.unique(np.round(fluid_z, 5))
z_mid = z_unique[len(z_unique)//2]
fluid_mask = np.abs(fluid_z - z_mid) < z_tol

fluid_x_slice = fluid_x[fluid_mask]
fluid_y_slice = fluid_y[fluid_mask]
U_mag_slice = U_mag[fluid_mask]
T_fluid_slice = T_fluid[fluid_mask]

print(f"Slice has {np.sum(fluid_mask)} cells")

# Also get solid centers
solid_centers = solid_mesh.cell_centers()
domain0_centers = domain0_mesh.cell_centers()
domain1_centers = domain1_mesh.cell_centers()

# Combine solid regions
solid_all_x = np.concatenate([solid_centers.points[:, 0], domain0_centers.points[:, 0], domain1_centers.points[:, 0]])
solid_all_y = np.concatenate([solid_centers.points[:, 1], domain0_centers.points[:, 1], domain1_centers.points[:, 1]])
solid_all_z = np.concatenate([solid_centers.points[:, 2], domain0_centers.points[:, 2], domain1_centers.points[:, 2]])
solid_all_T = np.concatenate([solid_mesh['T'], domain0_mesh['T'], domain1_mesh['T']])

# Filter solid cells near slice
solid_mask = np.abs(solid_all_z - z_mid) < z_tol
solid_x_slice = solid_all_x[solid_mask]
solid_y_slice = solid_all_y[solid_mask]
solid_T_slice = solid_all_T[solid_mask]

print(f"Solid slice has {np.sum(solid_mask)} cells")

# ============== Figure 1: Velocity Field ==============
fig, ax = plt.subplots(figsize=(14, 6), dpi=150)

# Plot velocity magnitude as scatter
scatter = ax.scatter(fluid_x_slice, fluid_y_slice, c=U_mag_slice,
                     cmap='jet', s=1, marker='s')
cbar = plt.colorbar(scatter, ax=ax, label='Velocity [m/s]')

# Plot solid region in gray
if len(solid_x_slice) > 0:
    ax.scatter(solid_x_slice, solid_y_slice, c='gray', s=1, marker='s')

ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Velocity Field - Staggered Finned Tube Heat Exchanger (z=0 slice)')
ax.set_aspect('equal')
ax.set_xlim([bounds[0], bounds[1]])
ax.set_ylim([bounds[2], bounds[3]])

vel_path = os.path.join(images_dir, "velocity_field.png")
plt.savefig(vel_path, dpi=150, bbox_inches='tight')
print(f"Saved: {vel_path}")
plt.close()

# ============== Figure 2: Temperature Field ==============
fig, ax = plt.subplots(figsize=(14, 6), dpi=150)

# Get temperature range
T_min = min(T_fluid_slice.min(), solid_T_slice.min()) if len(solid_T_slice) > 0 else T_fluid_slice.min()
T_max = max(T_fluid_slice.max(), solid_T_slice.max()) if len(solid_T_slice) > 0 else T_fluid_slice.max()

# Plot fluid temperature
scatter = ax.scatter(fluid_x_slice, fluid_y_slice, c=T_fluid_slice,
                     cmap='coolwarm', s=1, marker='s', vmin=T_min, vmax=T_max)
cbar = plt.colorbar(scatter, ax=ax, label='Temperature [K]')

# Plot solid temperature
if len(solid_x_slice) > 0:
    ax.scatter(solid_x_slice, solid_y_slice, c=solid_T_slice,
               cmap='coolwarm', s=1, marker='s', vmin=T_min, vmax=T_max)

ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Temperature Field - Staggered Finned Tube Heat Exchanger (z=0 slice)')
ax.set_aspect('equal')
ax.set_xlim([bounds[0], bounds[1]])
ax.set_ylim([bounds[2], bounds[3]])

temp_path = os.path.join(images_dir, "temperature_field.png")
plt.savefig(temp_path, dpi=150, bbox_inches='tight')
print(f"Saved: {temp_path}")
plt.close()

# ============== Figure 3: Overview with Both Fields ==============
fig, axes = plt.subplots(2, 1, figsize=(14, 10), dpi=150)

# Velocity
ax = axes[0]
scatter = ax.scatter(fluid_x_slice, fluid_y_slice, c=U_mag_slice,
                     cmap='jet', s=1, marker='s')
plt.colorbar(scatter, ax=ax, label='Velocity [m/s]')
if len(solid_x_slice) > 0:
    ax.scatter(solid_x_slice, solid_y_slice, c='gray', s=1, marker='s')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Velocity Magnitude')
ax.set_aspect('equal')
ax.set_xlim([bounds[0], bounds[1]])
ax.set_ylim([bounds[2], bounds[3]])

# Temperature
ax = axes[1]
scatter = ax.scatter(fluid_x_slice, fluid_y_slice, c=T_fluid_slice,
                     cmap='coolwarm', s=1, marker='s', vmin=T_min, vmax=T_max)
plt.colorbar(scatter, ax=ax, label='Temperature [K]')
if len(solid_x_slice) > 0:
    ax.scatter(solid_x_slice, solid_y_slice, c=solid_T_slice,
               cmap='coolwarm', s=1, marker='s', vmin=T_min, vmax=T_max)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Temperature')
ax.set_aspect('equal')
ax.set_xlim([bounds[0], bounds[1]])
ax.set_ylim([bounds[2], bounds[3]])

plt.tight_layout()
overview_path = os.path.join(images_dir, "overview.png")
plt.savefig(overview_path, dpi=150, bbox_inches='tight')
print(f"Saved: {overview_path}")
plt.close()

# ============== Figure 4: Geometry check (all cells xy projection) ==============
fig, ax = plt.subplots(figsize=(14, 6), dpi=150)

# Plot all fluid cells
ax.scatter(fluid_x, fluid_y, c='lightblue', s=0.1, alpha=0.5, label='Fluid')

# Plot all solid cells
ax.scatter(solid_all_x, solid_all_y, c='silver', s=0.1, alpha=0.5, label='Solid')

# Mark expected tube positions with circles
SL = 0.06  # Longitudinal pitch
ST = 0.06  # Transverse pitch
D = 0.025  # Tube diameter

tube_centers = [
    (0.0, 0.0),       # Row 1: at y=0 (bottom boundary)
    (0.0, ST),        # Row 1: at y=S_T (top boundary)
    (SL, ST / 2.0),   # Row 2: offset by half pitch
]

for cx, cy in tube_centers:
    circle = Circle((cx, cy), D/2, fill=False, edgecolor='red', linewidth=2)
    ax.add_patch(circle)
    ax.plot(cx, cy, 'r+', markersize=10)

ax.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='y boundaries')
ax.axhline(y=ST, color='k', linestyle='--', alpha=0.5)
ax.axvline(x=-SL/2, color='g', linestyle='--', alpha=0.5, label='x boundaries (periodic)')
ax.axvline(x=3*SL/2, color='g', linestyle='--', alpha=0.5)

ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Geometry Overview - Staggered Finned Tube REV\n(Red circles: expected tube positions)')
ax.set_aspect('equal')
ax.legend(loc='upper right')

geom_path = os.path.join(images_dir, "geometry_check.png")
plt.savefig(geom_path, dpi=150, bbox_inches='tight')
print(f"Saved: {geom_path}")
plt.close()

print("\nVisualization complete!")
print(f"Generated files in {images_dir}:")
print("  - velocity_field.png")
print("  - temperature_field.png")
print("  - overview.png")
print("  - geometry_check.png")
