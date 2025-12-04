#!/usr/bin/env python3
"""
Heatsink Geometry Generator using Gmsh
Creates a heatsink with fins and surrounding air domain for CHT analysis
"""

import gmsh
import sys
import os


def create_heatsink_geometry(
    # Base dimensions (mm -> converted to m)
    base_length: float = 50e-3,      # x direction
    base_width: float = 50e-3,       # y direction
    base_height: float = 5e-3,       # z direction

    # Fin dimensions
    fin_height: float = 25e-3,       # z direction (above base)
    fin_thickness: float = 2e-3,     # y direction
    num_fins: int = 5,

    # Air domain
    air_height: float = 40e-3,       # z direction (above fins)
    air_inlet_length: float = 20e-3, # x direction (before heatsink)
    air_outlet_length: float = 30e-3,# x direction (after heatsink)

    # Mesh parameters
    mesh_size_solid: float = 1.5e-3,
    mesh_size_fluid: float = 2.0e-3,
    mesh_size_interface: float = 1.0e-3,

    output_file: str = "heatsink.msh"
):
    """
    Create heatsink geometry with air domain for conjugate heat transfer.

    The geometry consists of:
    - Solid region: Heatsink (base + fins)
    - Fluid region: Air channel around and above the heatsink
    """

    gmsh.initialize()
    gmsh.model.add("heatsink_cht")

    # Calculate derived dimensions
    total_air_length = air_inlet_length + base_length + air_outlet_length
    total_height = base_height + fin_height + air_height
    fin_spacing = (base_width - num_fins * fin_thickness) / (num_fins + 1)

    print(f"Creating heatsink geometry:")
    print(f"  Base: {base_length*1000:.1f} x {base_width*1000:.1f} x {base_height*1000:.1f} mm")
    print(f"  Fins: {num_fins} fins, {fin_thickness*1000:.1f} mm thick, {fin_height*1000:.1f} mm tall")
    print(f"  Air domain: {total_air_length*1000:.1f} x {base_width*1000:.1f} x {total_height*1000:.1f} mm")

    # =========================================================================
    # Create Solid Region (Heatsink)
    # =========================================================================

    # Base
    base = gmsh.model.occ.addBox(
        air_inlet_length, 0, 0,
        base_length, base_width, base_height
    )

    # Fins
    fins = []
    for i in range(num_fins):
        y_pos = fin_spacing + i * (fin_thickness + fin_spacing)
        fin = gmsh.model.occ.addBox(
            air_inlet_length, y_pos, base_height,
            base_length, fin_thickness, fin_height
        )
        fins.append(fin)

    # Fuse base and fins to create solid heatsink
    all_solid_parts = [(3, base)] + [(3, f) for f in fins]
    solid_fused, _ = gmsh.model.occ.fuse([all_solid_parts[0]], all_solid_parts[1:])
    solid_volume = solid_fused[0][1]

    # =========================================================================
    # Create Fluid Region (Air)
    # =========================================================================

    # Full air box (will subtract heatsink)
    air_box = gmsh.model.occ.addBox(
        0, 0, 0,
        total_air_length, base_width, total_height
    )

    # Cut heatsink from air to create fluid region
    fluid_result, _ = gmsh.model.occ.cut([(3, air_box)], [(3, solid_volume)], removeObject=True, removeTool=False)
    fluid_volume = fluid_result[0][1]

    gmsh.model.occ.synchronize()

    # =========================================================================
    # Define Physical Groups (for boundary conditions)
    # =========================================================================

    # Get all surfaces
    solid_surfaces = gmsh.model.getBoundary([(3, solid_volume)], oriented=False)
    fluid_surfaces = gmsh.model.getBoundary([(3, fluid_volume)], oriented=False)

    # Identify surfaces by position
    eps = 1e-6

    # Helper function to get surface center
    def get_surface_center(surf_tag):
        bbox = gmsh.model.occ.getBoundingBox(2, surf_tag)
        return [(bbox[0] + bbox[3])/2, (bbox[1] + bbox[4])/2, (bbox[2] + bbox[5])/2]

    def get_surface_bounds(surf_tag):
        return gmsh.model.occ.getBoundingBox(2, surf_tag)

    # Classify fluid surfaces
    inlet_surfs = []
    outlet_surfs = []
    top_surfs = []
    front_back_surfs = []
    bottom_fluid_surfs = []
    interface_surfs = []

    for surf in fluid_surfaces:
        tag = abs(surf[1])
        bbox = get_surface_bounds(tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox

        # Inlet (x = 0)
        if abs(xmin) < eps and abs(xmax) < eps:
            inlet_surfs.append(tag)
        # Outlet (x = total_air_length)
        elif abs(xmin - total_air_length) < eps and abs(xmax - total_air_length) < eps:
            outlet_surfs.append(tag)
        # Top (z = total_height)
        elif abs(zmin - total_height) < eps and abs(zmax - total_height) < eps:
            top_surfs.append(tag)
        # Front (y = 0)
        elif abs(ymin) < eps and abs(ymax) < eps:
            front_back_surfs.append(tag)
        # Back (y = base_width)
        elif abs(ymin - base_width) < eps and abs(ymax - base_width) < eps:
            front_back_surfs.append(tag)
        # Bottom fluid (z = 0, outside heatsink region)
        elif abs(zmin) < eps and abs(zmax) < eps:
            if xmax < air_inlet_length + eps or xmin > air_inlet_length + base_length - eps:
                bottom_fluid_surfs.append(tag)
        # Interface surfaces (shared with solid)
        else:
            interface_surfs.append(tag)

    # Classify solid surfaces
    heatsink_bottom_surfs = []
    solid_interface_surfs = []

    for surf in solid_surfaces:
        tag = abs(surf[1])
        bbox = get_surface_bounds(tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox

        # Heatsink bottom (heat source)
        if abs(zmin) < eps and abs(zmax) < eps:
            heatsink_bottom_surfs.append(tag)
        else:
            solid_interface_surfs.append(tag)

    # Create physical groups
    # Fluid boundaries
    if inlet_surfs:
        gmsh.model.addPhysicalGroup(2, inlet_surfs, name="inlet")
    if outlet_surfs:
        gmsh.model.addPhysicalGroup(2, outlet_surfs, name="outlet")
    if top_surfs:
        gmsh.model.addPhysicalGroup(2, top_surfs, name="top")
    if front_back_surfs:
        gmsh.model.addPhysicalGroup(2, front_back_surfs, name="frontAndBack")
    if bottom_fluid_surfs:
        gmsh.model.addPhysicalGroup(2, bottom_fluid_surfs, name="bottomFluid")

    # Solid boundaries
    if heatsink_bottom_surfs:
        gmsh.model.addPhysicalGroup(2, heatsink_bottom_surfs, name="heatsinkBottom")

    # Interface (fluid-solid coupling)
    if interface_surfs:
        gmsh.model.addPhysicalGroup(2, interface_surfs, name="fluid_to_solid")
    if solid_interface_surfs:
        gmsh.model.addPhysicalGroup(2, solid_interface_surfs, name="solid_to_fluid")

    # Volume physical groups
    gmsh.model.addPhysicalGroup(3, [fluid_volume], name="fluid")
    gmsh.model.addPhysicalGroup(3, [solid_volume], name="solid")

    # =========================================================================
    # Mesh Generation
    # =========================================================================

    # Set mesh sizes
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size_fluid)

    # Finer mesh on interface
    interface_points = []
    for surf in interface_surfs + solid_interface_surfs:
        pts = gmsh.model.getBoundary([(2, surf)], combined=False, oriented=False)
        for pt in pts:
            interface_points.append(pt)

    if interface_points:
        # Get unique points
        unique_pts = list(set([p[1] for p in interface_points]))
        for pt in unique_pts:
            gmsh.model.mesh.setSize([(0, pt)], mesh_size_interface)

    # Set mesh algorithm
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT

    # Generate mesh
    print("Generating 3D mesh...")
    gmsh.model.mesh.generate(3)

    # Optimize mesh
    print("Optimizing mesh...")
    gmsh.model.mesh.optimize("Relocate3D")

    # Save mesh
    gmsh.write(output_file)
    print(f"Mesh saved to: {output_file}")

    # Print mesh statistics
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=3)
    total_3d_elements = sum(len(tags) for tags in elem_tags)
    print(f"Mesh statistics: {len(node_tags)} nodes, {total_3d_elements} 3D elements")

    gmsh.finalize()

    return output_file


def main():
    """Main function to create heatsink geometry."""
    import argparse

    parser = argparse.ArgumentParser(description='Create heatsink geometry for CHT analysis')
    parser.add_argument('--output', '-o', default='heatsink.msh', help='Output mesh file')
    parser.add_argument('--num-fins', type=int, default=5, help='Number of fins')
    parser.add_argument('--fin-height', type=float, default=25, help='Fin height in mm')
    parser.add_argument('--mesh-size', type=float, default=1.5, help='Base mesh size in mm')

    args = parser.parse_args()

    create_heatsink_geometry(
        num_fins=args.num_fins,
        fin_height=args.fin_height * 1e-3,
        mesh_size_solid=args.mesh_size * 1e-3,
        mesh_size_fluid=args.mesh_size * 1.5e-3,
        mesh_size_interface=args.mesh_size * 0.7e-3,
        output_file=args.output
    )


if __name__ == "__main__":
    main()
