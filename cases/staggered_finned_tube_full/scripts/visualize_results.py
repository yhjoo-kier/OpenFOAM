"""
Visualize results for full staggered finned-tube heat exchanger.

Generates images showing temperature and velocity fields through the
N-REV domain, demonstrating flow development and heat transfer.
"""
import os
from pathlib import Path

import numpy as np

try:
    import pyvista as pv
    pv.OFF_SCREEN = True
except ImportError:
    print("PyVista not available. Install with: pip install pyvista")
    exit(1)

CASE_DIR = Path(__file__).resolve().parent.parent
IMAGES_DIR = CASE_DIR / "images"
IMAGES_DIR.mkdir(exist_ok=True)


def find_latest_vtk():
    """Find the latest VTK result directory."""
    vtk_dir = CASE_DIR / "VTK" / "fluid"
    if not vtk_dir.exists():
        return None
    
    dirs = [d for d in vtk_dir.iterdir() if d.is_dir()]
    if not dirs:
        return None
    
    # Extract time from directory names
    times = []
    for d in dirs:
        try:
            time = int(d.name.split('_')[-1])
            times.append((time, d))
        except ValueError:
            continue
    
    if not times:
        return None
    
    times.sort(key=lambda x: x[0], reverse=True)
    return times[0][1]


def load_fluid_mesh(vtk_dir):
    """Load fluid VTK mesh."""
    internal_file = vtk_dir / "internal.vtu"
    if internal_file.exists():
        return pv.read(str(internal_file))
    return None


def load_solid_mesh(vtk_dir):
    """Load solid VTK mesh."""
    solid_dir = vtk_dir.parent.parent / "solid"
    if not solid_dir.exists():
        return None
    
    # Find matching time directory
    time_suffix = vtk_dir.name.split('_')[-1]
    solid_time_dir = None
    for d in solid_dir.iterdir():
        if d.name.endswith(time_suffix):
            solid_time_dir = d
            break
    
    if solid_time_dir is None:
        return None
    
    internal_file = solid_time_dir / "internal.vtu"
    if internal_file.exists():
        return pv.read(str(internal_file))
    return None


def create_temperature_plot(fluid_mesh, solid_mesh, output_path):
    """Create temperature field visualization."""
    plotter = pv.Plotter(off_screen=True, window_size=(1920, 800))
    
    # Slice at z=0 (center plane)
    if fluid_mesh is not None and 'T' in fluid_mesh.array_names:
        # Create slice
        bounds = fluid_mesh.bounds
        z_center = (bounds[4] + bounds[5]) / 2
        slice_fluid = fluid_mesh.slice(normal='z', origin=(0, 0, z_center))
        
        plotter.add_mesh(
            slice_fluid,
            scalars='T',
            cmap='coolwarm',
            show_scalar_bar=True,
            scalar_bar_args={
                'title': 'Temperature [K]',
                'vertical': True,
                'position_x': 0.85,
            }
        )
    
    if solid_mesh is not None and 'T' in solid_mesh.array_names:
        bounds = solid_mesh.bounds
        z_center = (bounds[4] + bounds[5]) / 2
        slice_solid = solid_mesh.slice(normal='z', origin=(0, 0, z_center))
        
        plotter.add_mesh(
            slice_solid,
            scalars='T',
            cmap='coolwarm',
            show_scalar_bar=False,
        )
    
    plotter.view_xy()
    plotter.camera.zoom(1.2)
    plotter.add_text(
        "Temperature Field - Full Domain (MPI Parallel)",
        position='upper_edge',
        font_size=14,
        color='black'
    )
    
    plotter.screenshot(str(output_path))
    plotter.close()
    print(f"Saved: {output_path}")


def create_velocity_plot(fluid_mesh, output_path):
    """Create velocity field visualization."""
    if fluid_mesh is None or 'U' not in fluid_mesh.array_names:
        print("No velocity data available")
        return
    
    plotter = pv.Plotter(off_screen=True, window_size=(1920, 800))
    
    # Compute velocity magnitude
    U = fluid_mesh['U']
    U_mag = np.linalg.norm(U, axis=1)
    fluid_mesh['U_mag'] = U_mag
    
    # Slice at z=0
    bounds = fluid_mesh.bounds
    z_center = (bounds[4] + bounds[5]) / 2
    slice_mesh = fluid_mesh.slice(normal='z', origin=(0, 0, z_center))
    
    plotter.add_mesh(
        slice_mesh,
        scalars='U_mag',
        cmap='viridis',
        show_scalar_bar=True,
        scalar_bar_args={
            'title': 'Velocity Magnitude [m/s]',
            'vertical': True,
            'position_x': 0.85,
        }
    )
    
    plotter.view_xy()
    plotter.camera.zoom(1.2)
    plotter.add_text(
        "Velocity Field - Full Domain (MPI Parallel)",
        position='upper_edge',
        font_size=14,
        color='black'
    )
    
    plotter.screenshot(str(output_path))
    plotter.close()
    print(f"Saved: {output_path}")


def create_streamwise_temperature_profile(fluid_mesh, output_path):
    """Create streamwise temperature profile plot using matplotlib."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib not available for profile plot")
        return
    
    if fluid_mesh is None or 'T' not in fluid_mesh.array_names:
        return
    
    # Sample along centerline
    bounds = fluid_mesh.bounds
    x_min, x_max = bounds[0], bounds[1]
    y_center = (bounds[2] + bounds[3]) / 2
    z_center = (bounds[4] + bounds[5]) / 2
    
    # Create sample points
    n_points = 100
    x_vals = np.linspace(x_min, x_max, n_points)
    points = np.column_stack([
        x_vals,
        np.full(n_points, y_center),
        np.full(n_points, z_center)
    ])
    
    # Sample temperature
    sample = fluid_mesh.sample_over_line(
        (x_min, y_center, z_center),
        (x_max, y_center, z_center),
        resolution=n_points
    )
    
    if 'T' not in sample.array_names:
        return
    
    T = sample['T']
    x = sample.points[:, 0]
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(x * 1000, T, 'b-', linewidth=2)
    ax.set_xlabel('Streamwise Position [mm]', fontsize=12)
    ax.set_ylabel('Temperature [K]', fontsize=12)
    ax.set_title('Centerline Temperature Development', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([x.min() * 1000, x.max() * 1000])
    
    plt.tight_layout()
    plt.savefig(str(output_path), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def create_geometry_view(fluid_mesh, solid_mesh, output_path):
    """Create 3D geometry view."""
    plotter = pv.Plotter(off_screen=True, window_size=(1600, 800))
    
    if fluid_mesh is not None:
        # Show fluid domain outline
        outline = fluid_mesh.outline()
        plotter.add_mesh(outline, color='blue', line_width=2)
        
        # Show fluid surface (semi-transparent)
        if 'T' in fluid_mesh.array_names:
            plotter.add_mesh(
                fluid_mesh,
                scalars='T',
                cmap='coolwarm',
                opacity=0.3,
                show_scalar_bar=False,
            )
    
    if solid_mesh is not None:
        plotter.add_mesh(
            solid_mesh,
            color='gray',
            opacity=0.8,
            show_edges=False,
        )
    
    plotter.add_text(
        "Full Domain Geometry - 10 REV Units",
        position='upper_edge',
        font_size=14,
        color='black'
    )
    
    plotter.view_isometric()
    plotter.camera.zoom(1.0)
    plotter.screenshot(str(output_path))
    plotter.close()
    print(f"Saved: {output_path}")


def main():
    print("Finding VTK results...")
    vtk_dir = find_latest_vtk()
    
    if vtk_dir is None:
        print("No VTK results found. Run foamToVTK first.")
        return
    
    print(f"Using VTK directory: {vtk_dir}")
    
    print("Loading meshes...")
    fluid_mesh = load_fluid_mesh(vtk_dir)
    solid_mesh = load_solid_mesh(vtk_dir)
    
    if fluid_mesh is None:
        print("Could not load fluid mesh")
        return
    
    print(f"Fluid mesh: {fluid_mesh.n_cells} cells")
    if solid_mesh is not None:
        print(f"Solid mesh: {solid_mesh.n_cells} cells")
    
    print("\nGenerating visualizations...")
    
    # Temperature field
    create_temperature_plot(
        fluid_mesh, solid_mesh,
        IMAGES_DIR / "temperature_field.png"
    )
    
    # Velocity field
    create_velocity_plot(
        fluid_mesh,
        IMAGES_DIR / "velocity_field.png"
    )
    
    # Streamwise temperature profile
    create_streamwise_temperature_profile(
        fluid_mesh,
        IMAGES_DIR / "temperature_profile.png"
    )
    
    # 3D geometry view
    create_geometry_view(
        fluid_mesh, solid_mesh,
        IMAGES_DIR / "geometry_3d.png"
    )
    
    print(f"\nAll images saved to: {IMAGES_DIR}")


if __name__ == "__main__":
    main()

