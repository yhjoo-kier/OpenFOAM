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
    "mesh_size": 0.0035,
    "interface_size": 0.002,
}


def _add_geometry(params):
    """
    Create staggered finned-tube REV geometry with clean inlet/outlet boundaries.

    Domain: x in [-S_L/2, 3*S_L/2], y in [0, S_T], z in [-P_fin/2, P_fin/2]

    To ensure inlet and outlet have identical flat surfaces (no tube intersection),
    we shift the domain so tubes are fully inside. The periodic length is still 2*S_L.

    Tube positions (viewed from z-axis):
      - (0, 0): full tube at left (inside domain)
      - (0, S_T): full tube at left-top (half inside)
      - (S_L, S_T/2): full tube at center
      - (S_L, 0): half tube at center-bottom (new, for symmetry)
      - (S_L, S_T): half tube at center-top (new, for symmetry)

    This creates a shifted REV with clean inlet/outlet at x=-S_L/2 and x=3*S_L/2.
    """
    SL = params["SL"]
    ST = params["ST"]
    P_fin = params["P_fin"]
    R_tube_out = params["R_tube_out"]
    R_tube_in = params["R_tube_in"]
    R_fin = params["R_fin"]
    t_fin = params["t_fin"]

    occ = gmsh.model.occ

    # Shift domain by S_L/2 so inlet/outlet don't intersect with tubes
    # Domain: x in [-S_L/2, 3*S_L/2], y in [0, S_T], z in [-P_fin/2, P_fin/2]
    x_start = -SL / 2.0
    domain_x = 2.0 * SL  # periodic length stays the same
    domain_y = ST
    domain_z = P_fin
    fluid_box = occ.addBox(x_start, 0.0, -P_fin / 2.0, domain_x, domain_y, domain_z)

    # Tube positions for proper staggered arrangement
    # Row 1 (x=0): tubes at y = n×S_T (at y-boundaries)
    # Row 2 (x=S_L): tubes at y = (n+0.5)×S_T (offset by half pitch)
    tube_centers = [
        (0.0, 0.0),       # Row 1: half tube at y=0 boundary
        (0.0, ST),        # Row 1: half tube at y=S_T boundary
        (SL, ST / 2.0),   # Row 2: full tube at center (offset)
    ]

    # Create each tube unit separately and intersect with domain box individually
    # This ensures all tube parts inside the domain are captured correctly
    trimmed_solids = []

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

        # Create a fresh copy of domain box for intersection
        domain_copy = occ.addBox(x_start, 0.0, -P_fin / 2.0, domain_x, domain_y, domain_z)

        # Intersect this tube unit with domain box to trim parts outside REV
        trimmed, _ = occ.intersect([(3, solid_tag)], [(3, domain_copy)], removeObject=True, removeTool=True)
        if trimmed:
            trimmed_solids.extend(trimmed)

    if not trimmed_solids:
        raise RuntimeError("No solid volumes created after trimming tubes to domain.")

    # Fragment fluid box with all trimmed solids for conformal mesh
    frag_entities, _ = occ.fragment([(3, fluid_box)], trimmed_solids)
    occ.synchronize()

    # Classify volumes: largest is fluid, solid volumes touch the fluid
    volumes = gmsh.model.getEntities(dim=3)
    if len(volumes) < 2:
        raise RuntimeError("Expected at least two volumes (fluid and solid) after fragment operation.")

    volumes_with_mass = []
    for dim, tag in volumes:
        mass = gmsh.model.occ.getMass(dim, tag)
        volumes_with_mass.append((tag, mass))
    volumes_with_mass.sort(key=lambda item: item[1], reverse=True)

    fluid_volume = volumes_with_mass[0][0]
    fluid_surfaces = _surface_sets_for_volume(fluid_volume)

    # Find all solid volumes that share interface with fluid
    solid_volumes = []
    for dim, tag in volumes:
        if tag == fluid_volume:
            continue
        surfaces = _surface_sets_for_volume(tag)
        if fluid_surfaces.intersection(surfaces):
            solid_volumes.append(tag)

    if not solid_volumes:
        raise RuntimeError("Failed to identify solid volumes sharing an interface with the fluid region.")

    # Remove stray fragments (e.g., inner bores outside domain)
    for dim, tag in volumes:
        if tag != fluid_volume and tag not in solid_volumes:
            gmsh.model.occ.remove([(dim, tag)], recursive=True)
    gmsh.model.occ.synchronize()

    return fluid_volume, solid_volumes


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


def _pair_by_centroid(surfaces, axis):
    def centroid(tag):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        cx, cy, cz = (xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0
        if axis == "x":
            return (cy, cz)
        if axis == "y":
            return (cx, cz)
        return (cx, cy)

    return sorted(surfaces, key=centroid)


def _set_periodicity(min_surfaces, max_surfaces, translation, axis):
    """Set periodicity between matching surfaces. Skip if counts don't match."""
    if len(min_surfaces) != len(max_surfaces):
        print(f"Warning: Periodic surface counts do not match ({len(min_surfaces)} vs {len(max_surfaces)}) for axis {axis}. Skipping Gmsh periodicity.")
        return False
    min_sorted = _pair_by_centroid(min_surfaces, axis=axis)
    max_sorted = _pair_by_centroid(max_surfaces, axis=axis)
    affine = [
        1.0, 0.0, 0.0, translation[0],
        0.0, 1.0, 0.0, translation[1],
        0.0, 0.0, 1.0, translation[2],
        0.0, 0.0, 0.0, 1.0,
    ]
    try:
        gmsh.model.mesh.setPeriodic(2, max_sorted, min_sorted, affine)
        return True
    except Exception as e:
        print(f"Warning: Could not set periodicity for axis {axis}: {e}")
        return False


def build_model(params):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("staggered_finned_tube")

    fluid_volume, solid_volumes = _add_geometry(params)

    fluid_surfaces = _surface_sets_for_volume(fluid_volume)

    # Collect all solid surfaces from multiple solid volumes
    all_solid_surfaces = set()
    for sv in solid_volumes:
        all_solid_surfaces.update(_surface_sets_for_volume(sv))
    interface_surfaces = fluid_surfaces.intersection(all_solid_surfaces)

    SL = params["SL"]
    ST = params["ST"]
    P_fin = params["P_fin"]

    # Domain shifted: x in [-S_L/2, 3*S_L/2], y in [0, S_T], z in [-P_fin/2, P_fin/2]
    x_start = -SL / 2.0
    x_end = x_start + 2.0 * SL  # = 3*S_L/2
    inlet = _select_surfaces_by_coord(fluid_surfaces, "x", x_start)
    outlet = _select_surfaces_by_coord(fluid_surfaces, "x", x_end)
    bottom = _select_surfaces_by_coord(fluid_surfaces, "y", 0.0)
    top = _select_surfaces_by_coord(fluid_surfaces, "y", ST)
    back = _select_surfaces_by_coord(fluid_surfaces, "z", -P_fin / 2.0)
    front = _select_surfaces_by_coord(fluid_surfaces, "z", P_fin / 2.0)

    if not all([inlet, outlet, bottom, top, back, front]):
        raise RuntimeError("Failed to locate all periodic boundary surfaces.")

    # Periodicity: translation vectors match domain dimensions
    _set_periodicity(inlet, outlet, (2.0 * SL, 0.0, 0.0), axis="x")
    _set_periodicity(bottom, top, (0.0, ST, 0.0), axis="y")  # Full S_T
    _set_periodicity(back, front, (0.0, 0.0, P_fin), axis="z")

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

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", params["mesh_size"] * 0.5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", params["mesh_size"])
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    # Refine near the fluid-solid interface.
    if interface_surfaces:
        field_id = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_id, "FacesList", list(interface_surfaces))
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", field_id)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", params["interface_size"])
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", params["mesh_size"])
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.002)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", 0.008)
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate staggered finned-tube REV mesh with periodicity.")
    parser.add_argument("--output", default="mesh.msh", help="Output mesh filename (msh)")
    parser.add_argument("--mesh-size", type=float, default=DEFAULT_PARAMS["mesh_size"], help="Target mesh size (m)")
    parser.add_argument("--interface-size", type=float, default=DEFAULT_PARAMS["interface_size"], help="Interface refinement size (m)")
    return parser.parse_args()


def main():
    args = parse_args()
    params = dict(DEFAULT_PARAMS)
    params["mesh_size"] = args.mesh_size
    params["interface_size"] = args.interface_size

    try:
        build_model(params)
        gmsh.write(args.output)
    except Exception as exc:  # noqa: BLE001
        gmsh.finalize()
        raise SystemExit(f"Mesh generation failed: {exc}")
    gmsh.finalize()
    print(f"Mesh written to {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
