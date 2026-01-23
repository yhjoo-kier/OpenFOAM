"""Mesh generator with parametric refinement for grid study."""

import subprocess
import sys
from pathlib import Path
from typing import Optional

from .config import MeshLevel


def generate_mesh(
    output_file: Path,
    mesh_level: MeshLevel,
    base_lc_fluid: float = 0.0015,
    base_lc_solid: float = 0.002,
) -> dict:
    """
    Generate mesh with specified refinement level.

    Args:
        output_file: Path to output .msh file
        mesh_level: Mesh refinement configuration
        base_lc_fluid: Base characteristic length for fluid [m]
        base_lc_solid: Base characteristic length for solid [m]

    Returns:
        Dictionary with mesh statistics
    """
    try:
        import gmsh
    except ImportError:
        raise ImportError("gmsh Python API is required. Install with: pip install gmsh")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # Suppress terminal output
    gmsh.model.add("heatsink_water")

    # Apply mesh factor
    lc_fluid = base_lc_fluid * mesh_level.mesh_factor
    lc_solid = base_lc_solid * mesh_level.mesh_factor

    # Geometry dimensions (same as original)
    L_dom = 0.080  # 80mm
    W_dom = 0.083  # 83mm
    H_base = 0.002  # 2mm
    H_fin = 0.015  # 15mm

    D_fin = 0.006  # 6mm
    R_fin = D_fin / 2
    N_rows = 8
    N_cols = 8
    pitch = 0.010  # 10mm

    x0_fin = 0.005
    y0_fin = 0.0065

    D_tube = 0.005  # 5mm
    R_tube = D_tube / 2
    L_tube = 0.025  # 25mm
    z_tube = 0.0075  # 7.5mm (from bottom of base)

    y_inlet_center = 0.005
    y_outlet_center = 0.078

    factory = gmsh.model.occ

    # 1. Create Base Plate (Solid)
    base_tag = factory.addBox(0, 0, 0, L_dom, W_dom, H_base)

    # 2. Create Fins (Solid)
    fin_tags = []
    for i in range(N_rows):
        for j in range(N_cols):
            xc = x0_fin + i * pitch
            yc = y0_fin + j * pitch
            ft = factory.addCylinder(xc, yc, H_base, 0, 0, H_fin, R_fin)
            fin_tags.append(ft)

    # Fuse base and fins
    solid_tag, _ = factory.fuse([(3, base_tag)], [(3, t) for t in fin_tags])
    solid_vol_tag = solid_tag[0][1]

    # 3. Create Fluid Domain
    fluid_block = factory.addBox(0, 0, H_base, L_dom, W_dom, H_fin)
    inlet_tube = factory.addCylinder(0, y_inlet_center, z_tube, -L_tube, 0, 0, R_tube)
    outlet_tube = factory.addCylinder(L_dom, y_outlet_center, z_tube, L_tube, 0, 0, R_tube)

    fluid_parts, _ = factory.fuse([(3, fluid_block)], [(3, inlet_tube), (3, outlet_tube)])
    fluid_vol_pre = fluid_parts[0][1]

    # 4. Interface Handling (Boolean Fragments)
    ov, ov_map = factory.fragment([(3, fluid_vol_pre)], [(3, solid_vol_tag)])
    factory.synchronize()

    # Identify volumes
    solid_final_tag = ov_map[1][0][1]
    fluid_final_tag = -1
    for t in ov_map[0]:
        if t[1] != solid_final_tag:
            fluid_final_tag = t[1]
            break

    if fluid_final_tag == -1:
        gmsh.finalize()
        raise RuntimeError("Could not identify fluid volume")

    # 5. Physical Groups
    gmsh.model.addPhysicalGroup(3, [solid_final_tag], name="solid")
    gmsh.model.addPhysicalGroup(3, [fluid_final_tag], name="fluid")

    # 6. Identify and tag surfaces
    surfaces = gmsh.model.getEntities(dim=2)
    inlet_surfs = []
    outlet_surfs = []
    heat_source_surfs = []
    interface_surfs = []  # fluid-solid interface for boundary layer

    eps = 1e-4

    for surf in surfaces:
        tag = surf[1]
        bb = gmsh.model.getBoundingBox(2, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bb
        com = gmsh.model.occ.getCenterOfMass(2, tag)
        cx, cy, cz = com

        # Inlet (x approx -L_tube)
        if abs(cx - (-L_tube)) < eps:
            inlet_surfs.append(tag)
            continue

        # Outlet (x approx L_dom + L_tube)
        if abs(cx - (L_dom + L_tube)) < eps:
            outlet_surfs.append(tag)
            continue

        # Heat Source (z approx 0)
        if abs(zmin - 0) < eps and abs(zmax - 0) < eps:
            heat_source_surfs.append(tag)
            continue

        # Interface surfaces (between z=H_base and z=H_base+H_fin, touching both volumes)
        # These are fin surfaces and base top surface
        if zmin >= H_base - eps and zmax <= H_base + H_fin + eps:
            # Check if surface is adjacent to fluid volume
            # Heuristic: surfaces at z=H_base (base top) or cylindrical (fins)
            if abs(zmin - H_base) < eps and abs(zmax - H_base) < eps:
                # Base top surface (excluding fin holes)
                interface_surfs.append(tag)
            elif zmin >= H_base - eps:
                # Could be fin surface
                interface_surfs.append(tag)

    gmsh.model.addPhysicalGroup(2, inlet_surfs, name="inlet")
    gmsh.model.addPhysicalGroup(2, outlet_surfs, name="outlet")
    gmsh.model.addPhysicalGroup(2, heat_source_surfs, name="heat_source")

    # 7. Mesh Size Fields
    # Global mesh size
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_fluid)

    # Near-wall refinement using Distance + Threshold fields
    # This creates a refined region near the fluid-solid interface
    if mesh_level.bl_num_layers > 0 and interface_surfs:
        # Distance field from interface surfaces
        dist_field = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist_field, "SurfacesList", interface_surfs)

        # Threshold field for gradual size transition
        # Size = bl_first_height near surface, transitions to lc_fluid at distance
        thresh_field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
        gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", mesh_level.bl_first_height)
        gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", lc_fluid)
        gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", 0)
        # Transition distance based on BL parameters
        bl_total_thickness = mesh_level.bl_first_height * (
            (mesh_level.bl_growth_ratio ** mesh_level.bl_num_layers - 1)
            / (mesh_level.bl_growth_ratio - 1)
        )
        gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", bl_total_thickness * 2)

        # Set as background mesh field
        gmsh.model.mesh.field.setAsBackgroundMesh(thresh_field)

    # 8. Generate mesh
    gmsh.model.mesh.generate(3)

    # Get mesh statistics
    nodes = gmsh.model.mesh.getNodes()
    num_nodes = len(nodes[0])

    # Count 3D elements
    elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(dim=3)
    num_cells = sum(len(tags) for tags in elem_tags)

    # Write mesh in Gmsh format 2.2 for OpenFOAM compatibility
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(str(output_file))
    gmsh.finalize()

    return {
        "num_nodes": num_nodes,
        "num_cells": num_cells,
        "mesh_file": str(output_file),
        "mesh_level": mesh_level.name,
        "lc_fluid": lc_fluid,
        "lc_solid": lc_solid,
    }


def convert_to_openfoam(msh_file: Path, case_dir: Path) -> bool:
    """
    Convert Gmsh mesh to OpenFOAM format.

    Args:
        msh_file: Path to .msh file
        case_dir: OpenFOAM case directory

    Returns:
        True if successful
    """
    cmd = ["gmshToFoam", str(msh_file), "-case", str(case_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=case_dir)

    if result.returncode != 0:
        print(f"gmshToFoam failed: {result.stderr}")
        return False

    return True


def split_mesh_regions(case_dir: Path) -> bool:
    """
    Split mesh into fluid and solid regions.

    Args:
        case_dir: OpenFOAM case directory

    Returns:
        True if successful
    """
    cmd = ["splitMeshRegions", "-cellZones", "-overwrite", "-case", str(case_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=case_dir)

    if result.returncode != 0:
        print(f"splitMeshRegions failed: {result.stderr}")
        return False

    return True
