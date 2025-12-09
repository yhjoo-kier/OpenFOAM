"""
Generate mesh for full staggered finned-tube heat exchanger with N REV repeats.

This script creates a mesh with N repetitions of the basic REV unit in the
streamwise (x) direction. Unlike the cyclic case, this enables MPI parallelization
and captures actual flow development.

Includes inlet and outlet extension regions for:
- Inlet: uniform flow development before first tube row
- Outlet: wake recovery after last tube row (prevents backflow issues)
"""
import argparse
import os

import gmsh


DEFAULT_PARAMS = {
    "SL": 0.060,  # longitudinal pitch (m)
    "ST": 0.060,  # transverse pitch (m)
    "R_tube_out": 0.012,
    "R_tube_in": 0.010,
    "R_fin": 0.022,
    "t_fin": 0.001,
    "P_fin": 0.050,
    "mesh_size": 0.004,  # Slightly coarser for larger domain
    "interface_size": 0.0025,
    "N_repeat": 10,  # Number of REV repetitions in x direction
    "inlet_extension": 1.0,   # Inlet extension in units of S_L
    "outlet_extension": 5.0,  # Outlet extension in units of S_L (for wake recovery)
}


def _add_geometry(params):
    """
    Create staggered finned-tube geometry with N REV repetitions and extension regions.
    
    Domain layout:
    
    |<- inlet ext ->|<------ N REV heat exchanger core ------>|<-- outlet ext -->|
    |    (empty)    |  tubes with fins in staggered pattern   |     (empty)      |
    
    Each REV unit spans 2*S_L in x direction.
    Total domain: x in [0, total_x], y in [0, S_T], z in [-P_fin/2, P_fin/2]
    
    Tube positions in each REV (offset by inlet_ext + i*2*S_L):
      - Row 1: tubes at y = 0, S_T (half tubes at boundaries)
      - Row 2: tube at y = S_T/2 (full tube, offset)
    """
    SL = params["SL"]
    ST = params["ST"]
    P_fin = params["P_fin"]
    R_tube_out = params["R_tube_out"]
    R_tube_in = params["R_tube_in"]
    R_fin = params["R_fin"]
    t_fin = params["t_fin"]
    N = params["N_repeat"]
    
    # Extension regions
    inlet_ext = params["inlet_extension"] * SL   # e.g., 1*SL = 0.06m
    outlet_ext = params["outlet_extension"] * SL  # e.g., 5*SL = 0.30m
    
    occ = gmsh.model.occ
    
    # Domain sizes
    REV_length = 2.0 * SL
    core_length = N * REV_length  # Heat exchanger core
    domain_x = inlet_ext + core_length + outlet_ext  # Total domain length
    domain_y = ST
    domain_z = P_fin
    
    print(f"Domain layout:")
    print(f"  Inlet extension:  {inlet_ext*1000:.1f} mm")
    print(f"  HX core (N={N}):  {core_length*1000:.1f} mm")
    print(f"  Outlet extension: {outlet_ext*1000:.1f} mm")
    print(f"  Total length:     {domain_x*1000:.1f} mm")
    
    # Create main fluid box (entire domain including extensions)
    fluid_box = occ.addBox(0.0, 0.0, -P_fin / 2.0, domain_x, domain_y, domain_z)
    
    # Create tube units for each REV repetition
    # Tubes start after inlet extension
    trimmed_solids = []
    
    for i in range(N):
        # Offset tubes by inlet extension
        x_offset = inlet_ext + i * REV_length
        
        # Tube positions within each REV (shifted by S_L/2 from REV start)
        # Row 1: x = x_offset + S_L/2
        # Row 2: x = x_offset + 3*S_L/2
        tube_centers = [
            (x_offset + SL / 2.0, 0.0),       # Row 1: half tube at y=0
            (x_offset + SL / 2.0, ST),        # Row 1: half tube at y=S_T
            (x_offset + 3.0 * SL / 2.0, ST / 2.0),  # Row 2: full tube at center
        ]
        
        for cx, cy in tube_centers:
            # Create tube outer cylinder
            tube_out = occ.addCylinder(cx, cy, -P_fin / 2.0, 0.0, 0.0, P_fin, R_tube_out)
            # Create fin disk
            fin = occ.addCylinder(cx, cy, -t_fin / 2.0, 0.0, 0.0, t_fin, R_fin)
            # Create inner bore
            tube_in = occ.addCylinder(cx, cy, -P_fin / 2.0, 0.0, 0.0, P_fin, R_tube_in)
            
            # Fuse tube outer with fin
            solid_unit, _ = occ.fuse([(3, tube_out)], [(3, fin)])
            solid_tag = solid_unit[0][1]
            
            # Subtract inner bore
            solid_cut, _ = occ.cut([(3, solid_tag)], [(3, tube_in)])
            solid_tag = solid_cut[0][1]
            
            # Create domain box copy for intersection (to trim tubes at boundaries)
            domain_copy = occ.addBox(0.0, 0.0, -P_fin / 2.0, domain_x, domain_y, domain_z)
            
            # Intersect tube with domain to trim parts outside
            trimmed, _ = occ.intersect([(3, solid_tag)], [(3, domain_copy)], 
                                        removeObject=True, removeTool=True)
            if trimmed:
                trimmed_solids.extend(trimmed)
    
    if not trimmed_solids:
        raise RuntimeError("No solid volumes created after trimming tubes to domain.")
    
    # Fragment fluid box with all trimmed solids for conformal mesh
    frag_entities, _ = occ.fragment([(3, fluid_box)], trimmed_solids)
    occ.synchronize()
    
    # Classify volumes: largest is fluid
    volumes = gmsh.model.getEntities(dim=3)
    if len(volumes) < 2:
        raise RuntimeError("Expected at least two volumes (fluid and solid) after fragment.")
    
    volumes_with_mass = []
    for dim, tag in volumes:
        mass = gmsh.model.occ.getMass(dim, tag)
        volumes_with_mass.append((tag, mass))
    volumes_with_mass.sort(key=lambda item: item[1], reverse=True)
    
    fluid_volume = volumes_with_mass[0][0]
    fluid_surfaces = _surface_sets_for_volume(fluid_volume)
    
    # Find all solid volumes sharing interface with fluid
    solid_volumes = []
    for dim, tag in volumes:
        if tag == fluid_volume:
            continue
        surfaces = _surface_sets_for_volume(tag)
        if fluid_surfaces.intersection(surfaces):
            solid_volumes.append(tag)
    
    if not solid_volumes:
        raise RuntimeError("Failed to identify solid volumes sharing interface with fluid.")
    
    # Remove stray fragments
    for dim, tag in volumes:
        if tag != fluid_volume and tag not in solid_volumes:
            gmsh.model.occ.remove([(dim, tag)], recursive=True)
    gmsh.model.occ.synchronize()
    
    return fluid_volume, solid_volumes, domain_x, domain_y, domain_z


def _surface_sets_for_volume(volume_tag):
    boundary = gmsh.model.getBoundary([(3, volume_tag)], oriented=False, recursive=False)
    return {tag for dim, tag in boundary if dim == 2}


def _select_surfaces_by_coord(surfaces, axis, target, tol=1e-5):
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


def build_model(params):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("staggered_finned_tube_full")
    
    fluid_volume, solid_volumes, domain_x, domain_y, domain_z = _add_geometry(params)
    
    fluid_surfaces = _surface_sets_for_volume(fluid_volume)
    
    # Collect all solid surfaces
    all_solid_surfaces = set()
    for sv in solid_volumes:
        all_solid_surfaces.update(_surface_sets_for_volume(sv))
    interface_surfaces = fluid_surfaces.intersection(all_solid_surfaces)
    
    P_fin = params["P_fin"]
    
    # Identify boundary surfaces (NO periodicity - actual inlet/outlet)
    inlet = _select_surfaces_by_coord(fluid_surfaces, "x", 0.0)
    outlet = _select_surfaces_by_coord(fluid_surfaces, "x", domain_x)
    bottom = _select_surfaces_by_coord(fluid_surfaces, "y", 0.0)
    top = _select_surfaces_by_coord(fluid_surfaces, "y", domain_y)
    back = _select_surfaces_by_coord(fluid_surfaces, "z", -P_fin / 2.0)
    front = _select_surfaces_by_coord(fluid_surfaces, "z", P_fin / 2.0)
    
    if not all([inlet, outlet, bottom, top, back, front]):
        raise RuntimeError("Failed to locate all boundary surfaces.")
    
    # Physical groups (NO periodicity constraints)
    gmsh.model.addPhysicalGroup(3, [fluid_volume], name="fluid")
    gmsh.model.addPhysicalGroup(3, solid_volumes, name="solid")
    
    gmsh.model.addPhysicalGroup(2, inlet, name="inlet")
    gmsh.model.addPhysicalGroup(2, outlet, name="outlet")
    gmsh.model.addPhysicalGroup(2, top, name="top")
    gmsh.model.addPhysicalGroup(2, bottom, name="bottom")
    gmsh.model.addPhysicalGroup(2, front, name="front")
    gmsh.model.addPhysicalGroup(2, back, name="back")
    gmsh.model.addPhysicalGroup(2, list(interface_surfaces), name="interface_fluid_solid")
    
    gmsh.model.occ.synchronize()
    
    # Mesh settings
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", params["mesh_size"] * 0.5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", params["mesh_size"])
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 16)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT algorithm
    
    # Refine near fluid-solid interface
    if interface_surfaces:
        field_id = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_id, "FacesList", list(interface_surfaces))
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", field_id)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", params["interface_size"])
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", params["mesh_size"])
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.002)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", 0.01)
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
    
    print(f"Generating mesh for {params['N_repeat']} REV units...")
    print(f"Domain size: {domain_x:.3f} x {domain_y:.3f} x {domain_z:.3f} m")
    
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    
    # Report mesh statistics
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    print(f"Total nodes: {len(node_tags)}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate full staggered finned-tube mesh with N REV repetitions."
    )
    parser.add_argument("--output", default="mesh.msh", help="Output mesh filename")
    parser.add_argument("--mesh-size", type=float, default=DEFAULT_PARAMS["mesh_size"],
                        help="Target mesh size (m)")
    parser.add_argument("--interface-size", type=float, default=DEFAULT_PARAMS["interface_size"],
                        help="Interface refinement size (m)")
    parser.add_argument("-N", "--n-repeat", type=int, default=DEFAULT_PARAMS["N_repeat"],
                        help="Number of REV repetitions in x direction")
    parser.add_argument("--inlet-ext", type=float, default=DEFAULT_PARAMS["inlet_extension"],
                        help="Inlet extension length in units of S_L (default: 1.0)")
    parser.add_argument("--outlet-ext", type=float, default=DEFAULT_PARAMS["outlet_extension"],
                        help="Outlet extension length in units of S_L (default: 5.0)")
    return parser.parse_args()


def main():
    args = parse_args()
    params = dict(DEFAULT_PARAMS)
    params["mesh_size"] = args.mesh_size
    params["interface_size"] = args.interface_size
    params["N_repeat"] = args.n_repeat
    params["inlet_extension"] = args.inlet_ext
    params["outlet_extension"] = args.outlet_ext
    
    try:
        build_model(params)
        gmsh.write(args.output)
    except Exception as exc:
        gmsh.finalize()
        raise SystemExit(f"Mesh generation failed: {exc}")
    gmsh.finalize()
    print(f"Mesh written to {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()

