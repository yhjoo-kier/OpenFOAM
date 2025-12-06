#!/usr/bin/env python3
"""
Visualize y-z cross-section temperature distribution with contour plots
"""

import pyvista as pv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Load the VTK file
vtm_file = "VTK/heatsink_flow_3000.vtm"
print(f"Loading {vtm_file}...")

mesh = pv.read(vtm_file)

# Get internal mesh
if isinstance(mesh, pv.MultiBlock):
    internal = None
    for block in mesh:
        if block is not None and hasattr(block, 'n_cells') and block.n_cells > 0:
            internal = block
            break
else:
    internal = mesh

print(f"Mesh loaded: {internal.n_cells} cells")

# Channel dimensions
channel_length = 0.12
channel_height = 0.04
channel_width = 0.06

# Create y-z slices at different x positions
x_positions = [0.03, 0.06, 0.09]  # 30mm, 60mm, 90mm (flow direction)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('Y-Z Cross-Section Temperature Distribution (Flow direction: x)\n' +
             f'Inlet: 25°C, Wall: 80°C, Velocity: 0.1 m/s', fontsize=12)

for idx, x_pos in enumerate(x_positions):
    # Create a plane at x = x_pos
    plane_origin = [x_pos, channel_height/2, channel_width/2]
    plane_normal = [1, 0, 0]  # Normal pointing in x direction

    slice_mesh = internal.slice(normal=plane_normal, origin=plane_origin)

    if slice_mesh.n_points > 0:
        # Convert cell data to point data for proper visualization
        slice_mesh = slice_mesh.cell_data_to_point_data()

        # Extract coordinates and temperature
        points = slice_mesh.points
        T = slice_mesh['T'] - 273.15  # Convert to Celsius

        # Get y and z coordinates (in mm)
        y = points[:, 1] * 1000
        z = points[:, 2] * 1000

        ax = axes[idx]

        # Create regular grid for contour plot
        yi = np.linspace(y.min(), y.max(), 100)
        zi = np.linspace(z.min(), z.max(), 100)
        Yi, Zi = np.meshgrid(yi, zi)

        # Interpolate temperature onto regular grid
        Ti = griddata((y, z), T, (Yi, Zi), method='cubic')

        # Create filled contour plot
        levels = np.linspace(25, 80, 23)  # Temperature levels
        contour = ax.contourf(Zi, Yi, Ti, levels=levels, cmap='hot', extend='both')

        # Add contour lines
        ax.contour(Zi, Yi, Ti, levels=levels[::2], colors='black', linewidths=0.3, alpha=0.5)

        ax.set_xlabel('Z (mm)')
        ax.set_ylabel('Y (mm)')
        ax.set_title(f'x = {x_pos*1000:.0f} mm')
        ax.set_aspect('equal')

        # Add colorbar
        cbar = plt.colorbar(contour, ax=ax)
        cbar.set_label('Temperature (°C)')

plt.tight_layout()
plt.savefig('yz_temperature_contour.png', dpi=150, bbox_inches='tight')
print("Saved: yz_temperature_contour.png")

# Also create centerline temperature plot
print("\nCreating centerline temperature plot...")
fig2, ax2 = plt.subplots(figsize=(10, 6))

# Sample along centerline (y=channel_height/2, z=channel_width/2)
x_samples = np.linspace(0.001, 0.119, 50)
centerline_T = []

for x in x_samples:
    # Find nearest cell to centerline
    sample_point = [x, channel_height/2, channel_width/2]
    closest_idx = internal.find_closest_cell(sample_point)
    if closest_idx >= 0:
        centerline_T.append(internal['T'][closest_idx] - 273.15)
    else:
        centerline_T.append(np.nan)

ax2.plot(x_samples*1000, centerline_T, 'b-', linewidth=2)
ax2.axhline(y=25, color='g', linestyle='--', label='Inlet temp (25°C)')
ax2.axhline(y=80, color='r', linestyle='--', label='Wall temp (80°C)')
ax2.axhspan(30, 70, alpha=0.2, color='green', label='Target outlet range')
ax2.set_xlabel('X position (mm)')
ax2.set_ylabel('Temperature (°C)')
ax2.set_title('Centerline Temperature Profile\n(y = 20mm, z = 30mm)')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 120)
ax2.set_ylim(20, 85)

plt.savefig('centerline_temperature.png', dpi=150, bbox_inches='tight')
print("Saved: centerline_temperature.png")

# Wall temperature profile at bottom
print("\nCreating bottom wall temperature plot...")
fig3, ax3 = plt.subplots(figsize=(10, 6))

# Sample along bottom wall (y near 0)
y_bottom = 0.001  # 1mm from bottom
x_samples_wall = np.linspace(0.001, 0.119, 100)
wall_T = []

for x in x_samples_wall:
    sample_point = [x, y_bottom, channel_width/2]
    closest_idx = internal.find_closest_cell(sample_point)
    if closest_idx >= 0:
        wall_T.append(internal['T'][closest_idx] - 273.15)
    else:
        wall_T.append(np.nan)

ax3.plot(x_samples_wall*1000, wall_T, 'r-', linewidth=2, label='Near-wall temperature')
ax3.axhline(y=80, color='k', linestyle='--', label='Wall BC (80°C)')
ax3.axhline(y=25, color='b', linestyle='--', label='Inlet (25°C)')
ax3.set_xlabel('X position (mm)')
ax3.set_ylabel('Temperature (°C)')
ax3.set_title('Near-Wall Temperature Profile (y = 1mm from heated wall)')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 120)
ax3.set_ylim(20, 85)

plt.savefig('wall_temperature_profile.png', dpi=150, bbox_inches='tight')
print("Saved: wall_temperature_profile.png")

print("\nVisualization complete!")
