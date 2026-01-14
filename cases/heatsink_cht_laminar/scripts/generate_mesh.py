#!/usr/bin/env python3
"""
Generate mesh for heatsink CHT laminar case using Gmsh.

Geometry:
- Main heatsink domain: 80mm x 83mm x 17mm (including 2mm base plate)
- Cylindrical pin fins: 6mm diameter, 15mm height, 8x8 array, 10mm spacing
- Base plate: 2mm thickness with heat flux BC on bottom
- Inlet/Outlet tubes: 5mm diameter, 25mm length

Coordinate system:
- x: 0 to 80mm (main domain), -25mm to 0 (inlet tube), 80mm to 105mm (outlet tube)
- y: 0 to 83mm
- z: 0 to 2mm (base plate), 2mm to 17mm (fluid + fins)
"""

import argparse
import os

import gmsh


DEFAULT_PARAMS = {
    # Domain dimensions (m)
    "L_dom": 0.080,      # x-direction length
    "W_dom": 0.083,      # y-direction length
    "H_fluid": 0.015,    # fluid/fin height
    "H_base": 0.002,     # base plate thickness

    # Pin fin parameters
    "D_fin": 0.006,      # fin diameter
    "fin_pitch_x": 0.010,  # x-direction pitch
    "fin_pitch_y": 0.010,  # y-direction pitch
    "fin_origin_x": 0.005,  # first fin center x
    "fin_origin_y": 0.0065, # first fin center y
    "n_fin_x": 8,        # number of fins in x
    "n_fin_y": 8,        # number of fins in y

    # Inlet/Outlet tube parameters
    "D_tube": 0.005,     # tube diameter
    "L_tube": 0.025,     # tube length
    "inlet_y": 0.005,    # inlet tube center y
    "outlet_y": 0.078,   # outlet tube center y
    "tube_z": 0.0075,    # tube center z (relative to fluid bottom, i.e., z=H_base)

    # Mesh parameters
    "mesh_size": 0.0015,
    "interface_size": 0.0008,
    "tube_size": 0.0006,
}


def _surface_sets_for_volume(volume_tag):
    """Get set of surface tags for a volume."""
    boundary = gmsh.model.getBoundary([(3, volume_tag)], oriented=False, recursive=False)
    return {tag for dim, tag in boundary if dim == 2}


def _select_surfaces_by_coord(surfaces, axis, target, tol=1e-5):
    """Select surfaces that lie on a plane at given coordinate."""
    selected = []
    for tag in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        coord_min, coord_max = {
            "x": (xmin, xmax),
            "y": (ymin, ymax),
            "z": (zmin, zmax),
        }[axis]
        if abs(coord_min - target) < tol and abs(coord_max - target) < tol:
            selected.append(tag)
    return selected


def _select_cylindrical_surfaces(surfaces, axis, radius, tol=1e-5):
    """Select cylindrical surfaces with given radius oriented along axis."""
    selected = []
    for tag in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)

        if axis == "z":
            dx = xmax - xmin
            dy = ymax - ymin
            # Cylindrical surface: bounding box width ~ 2*radius in both x and y
            if abs(dx - 2*radius) < tol and abs(dy - 2*radius) < tol:
                selected.append(tag)
        elif axis == "x":
            dy = ymax - ymin
            dz = zmax - zmin
            if abs(dy - 2*radius) < tol and abs(dz - 2*radius) < tol:
                selected.append(tag)
    return selected


def build_geometry(params):
    """Build the heatsink geometry with separate fluid and solid regions."""

    L_dom = params["L_dom"]
    W_dom = params["W_dom"]
    H_fluid = params["H_fluid"]
    H_base = params["H_base"]
    D_fin = params["D_fin"]
    R_fin = D_fin / 2.0
    D_tube = params["D_tube"]
    R_tube = D_tube / 2.0
    L_tube = params["L_tube"]

    occ = gmsh.model.occ

    # ========== Create base plate (solid) ==========
    base_plate = occ.addBox(0, 0, 0, L_dom, W_dom, H_base)

    # ========== Create pin fins (solid) ==========
    fin_volumes = []
    fin_centers = []

    for i in range(params["n_fin_x"]):
        for j in range(params["n_fin_y"]):
            cx = params["fin_origin_x"] + i * params["fin_pitch_x"]
            cy = params["fin_origin_y"] + j * params["fin_pitch_y"]
            fin_centers.append((cx, cy))

            # Create cylinder from base plate top to fluid top
            cyl = occ.addCylinder(cx, cy, H_base, 0, 0, H_fluid, R_fin)
            fin_volumes.append(cyl)

    # Fuse all solid parts (base plate + fins)
    all_solid_parts = [(3, base_plate)] + [(3, v) for v in fin_volumes]
    if len(all_solid_parts) > 1:
        fused_solid, _ = occ.fuse([all_solid_parts[0]], all_solid_parts[1:])
        solid_tags = [tag for dim, tag in fused_solid if dim == 3]
    else:
        solid_tags = [base_plate]

    occ.synchronize()

    # ========== Create fluid domain ==========
    # Main fluid box (above base plate)
    fluid_main = occ.addBox(0, 0, H_base, L_dom, W_dom, H_fluid)

    # Inlet tube (extending in -x direction)
    inlet_center_z = H_base + params["tube_z"]
    inlet_tube = occ.addCylinder(-L_tube, params["inlet_y"], inlet_center_z,
                                  L_tube, 0, 0, R_tube)

    # Outlet tube (extending in +x direction)
    outlet_tube = occ.addCylinder(L_dom, params["outlet_y"], inlet_center_z,
                                   L_tube, 0, 0, R_tube)

    # Fuse fluid parts
    fluid_fused, _ = occ.fuse([(3, fluid_main)], [(3, inlet_tube), (3, outlet_tube)])
    fluid_tags = [tag for dim, tag in fluid_fused if dim == 3]

    occ.synchronize()

    # ========== Cut fins from fluid domain ==========
    # Re-create fin cylinders for cutting (since they were fused into solid)
    fin_cut_volumes = []
    for cx, cy in fin_centers:
        cyl = occ.addCylinder(cx, cy, H_base, 0, 0, H_fluid, R_fin)
        fin_cut_volumes.append((3, cyl))

    occ.synchronize()

    # Cut fins from fluid
    if fin_cut_volumes:
        fluid_cut, _ = occ.cut([(3, tag) for tag in fluid_tags], fin_cut_volumes)
        fluid_tags = [tag for dim, tag in fluid_cut if dim == 3]

    occ.synchronize()

    # ========== Fragment for conformal mesh ==========
    # Get current solid and fluid volumes
    all_volumes = gmsh.model.getEntities(dim=3)

    # Identify solid and fluid by checking if they overlap with base plate region
    solid_volumes = []
    fluid_volumes = []

    for dim, tag in all_volumes:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        # Solid is in z: 0 to H_base+H_fluid, but includes base plate
        # Fluid is in z: H_base to H_base+H_fluid and tubes
        if zmin < H_base / 2:  # Contains base plate region
            solid_volumes.append(tag)
        else:
            fluid_volumes.append(tag)

    # Fragment all volumes for conformal interface
    all_vol_tuples = [(3, t) for t in solid_volumes + fluid_volumes]
    if len(all_vol_tuples) > 1:
        frag_result, _ = occ.fragment(all_vol_tuples[:1], all_vol_tuples[1:])

    occ.synchronize()

    return solid_volumes, fluid_volumes, fin_centers, params


def assign_physical_groups(params):
    """Assign physical groups to volumes and surfaces."""

    L_dom = params["L_dom"]
    W_dom = params["W_dom"]
    H_fluid = params["H_fluid"]
    H_base = params["H_base"]
    D_tube = params["D_tube"]
    R_tube = D_tube / 2.0
    L_tube = params["L_tube"]

    # Get all volumes
    all_volumes = gmsh.model.getEntities(dim=3)

    # Classify volumes by their z-extent
    solid_volumes = []
    fluid_volumes = []

    for dim, tag in all_volumes:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        if zmin < H_base * 0.5:  # Contains base plate
            solid_volumes.append(tag)
        else:
            fluid_volumes.append(tag)

    if not solid_volumes or not fluid_volumes:
        raise RuntimeError(f"Failed to identify regions. Solid: {solid_volumes}, Fluid: {fluid_volumes}")

    # Create physical groups for volumes
    gmsh.model.addPhysicalGroup(3, solid_volumes, name="solid")
    gmsh.model.addPhysicalGroup(3, fluid_volumes, name="fluid")

    # Get all surfaces
    fluid_surfaces = set()
    for fv in fluid_volumes:
        fluid_surfaces.update(_surface_sets_for_volume(fv))

    solid_surfaces = set()
    for sv in solid_volumes:
        solid_surfaces.update(_surface_sets_for_volume(sv))

    # Interface surfaces (shared between fluid and solid)
    interface_surfaces = fluid_surfaces.intersection(solid_surfaces)

    # Identify boundary surfaces
    inlet_center_z = H_base + params["tube_z"]

    # Inlet: circular surface at x = -L_tube
    inlet = _select_surfaces_by_coord(fluid_surfaces, "x", -L_tube)

    # Outlet: circular surface at x = L_dom + L_tube
    outlet = _select_surfaces_by_coord(fluid_surfaces, "x", L_dom + L_tube)

    # Base bottom: z = 0
    base_bottom = _select_surfaces_by_coord(solid_surfaces, "z", 0.0)

    # Fluid top: z = H_base + H_fluid
    fluid_top = _select_surfaces_by_coord(fluid_surfaces, "z", H_base + H_fluid)

    # Side walls (y = 0 and y = W_dom)
    fluid_y0 = _select_surfaces_by_coord(fluid_surfaces - interface_surfaces, "y", 0.0)
    fluid_ymax = _select_surfaces_by_coord(fluid_surfaces - interface_surfaces, "y", W_dom)
    solid_y0 = _select_surfaces_by_coord(solid_surfaces - interface_surfaces, "y", 0.0)
    solid_ymax = _select_surfaces_by_coord(solid_surfaces - interface_surfaces, "y", W_dom)

    # Front/back walls (x = 0 and x = L_dom) for main domain only
    fluid_x0 = _select_surfaces_by_coord(fluid_surfaces - interface_surfaces, "x", 0.0)
    fluid_xmax = _select_surfaces_by_coord(fluid_surfaces - interface_surfaces, "x", L_dom)
    solid_x0 = _select_surfaces_by_coord(solid_surfaces - interface_surfaces, "x", 0.0)
    solid_xmax = _select_surfaces_by_coord(solid_surfaces - interface_surfaces, "x", L_dom)

    # Tube outer walls (cylindrical surfaces on tubes)
    # These are fluid surfaces that are not on planes and not interface
    tube_walls = []
    for surf in fluid_surfaces - interface_surfaces:
        if surf in inlet or surf in outlet:
            continue
        if surf in fluid_top or surf in fluid_y0 or surf in fluid_ymax:
            continue
        if surf in fluid_x0 or surf in fluid_xmax:
            continue
        # Check if it's a cylindrical surface (tube wall)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, surf)
        # Tube walls: either in inlet region (x < 0) or outlet region (x > L_dom)
        if xmax < 0.001 or xmin > L_dom - 0.001:
            tube_walls.append(surf)

    # Create physical groups for surfaces
    if inlet:
        gmsh.model.addPhysicalGroup(2, inlet, name="inlet")
    if outlet:
        gmsh.model.addPhysicalGroup(2, outlet, name="outlet")
    if base_bottom:
        gmsh.model.addPhysicalGroup(2, base_bottom, name="base_bottom")
    if fluid_top:
        gmsh.model.addPhysicalGroup(2, fluid_top, name="fluid_top")
    if interface_surfaces:
        gmsh.model.addPhysicalGroup(2, list(interface_surfaces), name="interface_fluid_solid")

    # Combine side walls
    side_walls_fluid = list(set(fluid_y0 + fluid_ymax + fluid_x0 + fluid_xmax) - set(tube_walls))
    side_walls_solid = list(set(solid_y0 + solid_ymax + solid_x0 + solid_xmax))

    if side_walls_fluid:
        gmsh.model.addPhysicalGroup(2, side_walls_fluid, name="walls_fluid")
    if side_walls_solid:
        gmsh.model.addPhysicalGroup(2, side_walls_solid, name="walls_solid")
    if tube_walls:
        gmsh.model.addPhysicalGroup(2, tube_walls, name="tube_walls")

    return {
        "solid_volumes": solid_volumes,
        "fluid_volumes": fluid_volumes,
        "interface": list(interface_surfaces),
    }


def set_mesh_sizes(params, region_info):
    """Configure mesh size fields."""

    mesh_size = params["mesh_size"]
    interface_size = params["interface_size"]
    tube_size = params["tube_size"]

    # Global mesh size
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", tube_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 24)

    # Use Delaunay algorithm
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    # Refinement near interface
    if region_info["interface"]:
        field_id = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_id, "SurfacesList", region_info["interface"])

        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", field_id)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", interface_size)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", mesh_size)
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0005)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", 0.003)
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold)


def build_model(params):
    """Build complete mesh model."""

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("heatsink_cht_laminar")

    # Build geometry
    solid_volumes, fluid_volumes, fin_centers, params = build_geometry(params)

    # Assign physical groups
    region_info = assign_physical_groups(params)

    # Set mesh sizes
    set_mesh_sizes(params, region_info)

    # Generate mesh
    gmsh.model.mesh.generate(3)

    # Optimize mesh
    try:
        gmsh.model.mesh.optimize("Netgen")
    except Exception as e:
        print(f"Mesh optimization warning: {e}")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate heatsink CHT laminar mesh with Gmsh"
    )
    parser.add_argument(
        "--output",
        default="mesh.msh",
        help="Output mesh filename (default: mesh.msh)"
    )
    parser.add_argument(
        "--mesh-size",
        type=float,
        default=DEFAULT_PARAMS["mesh_size"],
        help=f"Global mesh size in meters (default: {DEFAULT_PARAMS['mesh_size']})"
    )
    parser.add_argument(
        "--interface-size",
        type=float,
        default=DEFAULT_PARAMS["interface_size"],
        help=f"Interface refinement size (default: {DEFAULT_PARAMS['interface_size']})"
    )
    parser.add_argument(
        "--coarse",
        action="store_true",
        help="Use coarse mesh for quick testing"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    params = dict(DEFAULT_PARAMS)

    if args.coarse:
        params["mesh_size"] = 0.003
        params["interface_size"] = 0.002
        params["tube_size"] = 0.0015
    else:
        params["mesh_size"] = args.mesh_size
        params["interface_size"] = args.interface_size

    try:
        build_model(params)
        gmsh.write(args.output)
        print(f"Mesh written to {os.path.abspath(args.output)}")

        # Print mesh statistics
        nodes = gmsh.model.mesh.getNodes()
        elements = gmsh.model.mesh.getElements()
        print(f"Nodes: {len(nodes[0])}")
        print(f"Element types: {len(elements[0])}")

    except Exception as exc:
        gmsh.finalize()
        raise SystemExit(f"Mesh generation failed: {exc}")

    gmsh.finalize()


if __name__ == "__main__":
    main()
