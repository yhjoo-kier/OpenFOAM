#!/usr/bin/env python3
"""
Analyze heatsink cooling simulation results
"""

import pyvista as pv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load the VTK file
vtm_file = "VTK/heatsink_flow_5000.vtm"
print(f"Loading {vtm_file}...")

try:
    mesh = pv.read(vtm_file)
    print(f"Mesh type: {type(mesh)}")

    # For MultiBlock dataset, iterate through blocks
    if isinstance(mesh, pv.MultiBlock):
        print(f"Number of blocks: {len(mesh)}")
        for i, block in enumerate(mesh):
            if block is not None:
                print(f"  Block {i}: {block}")
                if hasattr(block, 'array_names'):
                    print(f"    Arrays: {block.array_names}")

        # Get the internal mesh (usually the first non-None block)
        internal = None
        for block in mesh:
            if block is not None and hasattr(block, 'n_cells') and block.n_cells > 0:
                internal = block
                break

        if internal is None:
            # Try to get composite as single mesh
            internal = mesh.combine()
    else:
        internal = mesh

    print(f"\nInternal mesh: {internal}")
    print(f"  Number of cells: {internal.n_cells}")
    print(f"  Number of points: {internal.n_points}")
    print(f"  Arrays: {internal.array_names}")

    # Check if T field exists
    if 'T' in internal.array_names:
        T = internal['T']
        print(f"\nTemperature field statistics:")
        print(f"  Min: {T.min():.2f} K ({T.min() - 273.15:.2f} °C)")
        print(f"  Max: {T.max():.2f} K ({T.max() - 273.15:.2f} °C)")
        print(f"  Mean: {T.mean():.2f} K ({T.mean() - 273.15:.2f} °C)")
    else:
        print("T field not found in mesh")
        print(f"Available arrays: {internal.array_names}")

    # Calculate outlet temperature (x = channel_length = 0.12)
    channel_length = 0.12
    outlet_tolerance = 0.005  # 5mm tolerance

    # Get cell centers
    centers = internal.cell_centers()
    x_coords = centers.points[:, 0]

    # Find cells near outlet
    outlet_mask = x_coords > (channel_length - outlet_tolerance)

    if 'T' in internal.array_names:
        T_outlet = internal['T'][outlet_mask]
        if len(T_outlet) > 0:
            avg_outlet_T = T_outlet.mean()
            print(f"\nOutlet temperature (x > {channel_length - outlet_tolerance:.3f} m):")
            print(f"  Number of cells: {len(T_outlet)}")
            print(f"  Average: {avg_outlet_T:.2f} K ({avg_outlet_T - 273.15:.2f} °C)")
            print(f"  Min: {T_outlet.min():.2f} K ({T_outlet.min() - 273.15:.2f} °C)")
            print(f"  Max: {T_outlet.max():.2f} K ({T_outlet.max() - 273.15:.2f} °C)")

            # Check if in desired range (30-70°C)
            outlet_temp_C = avg_outlet_T - 273.15
            if 30 <= outlet_temp_C <= 70:
                print(f"\n✓ Outlet temperature ({outlet_temp_C:.1f}°C) is within target range (30-70°C)")
            else:
                print(f"\n✗ Outlet temperature ({outlet_temp_C:.1f}°C) is OUTSIDE target range (30-70°C)")
                if outlet_temp_C < 30:
                    print("  Suggestion: Decrease inlet velocity to increase heat transfer")
                else:
                    print("  Suggestion: Increase inlet velocity to reduce outlet temperature")

except Exception as e:
    print(f"Error loading VTK file: {e}")
    import traceback
    traceback.print_exc()
