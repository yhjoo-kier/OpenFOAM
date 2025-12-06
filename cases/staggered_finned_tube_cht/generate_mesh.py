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
    "P_fin": 0.005,
    "mesh_size": 0.0035,
    "interface_size": 0.002,
}


def _add_geometry(params):
    SL = params["SL"]
    ST = params["ST"]
    P_fin = params["P_fin"]
    R_tube_out = params["R_tube_out"]
    R_tube_in = params["R_tube_in"]
    R_fin = params["R_fin"]
    t_fin = params["t_fin"]

    occ = gmsh.model.occ
    fluid_box = occ.addBox(-SL / 2.0, -ST / 2.0, -P_fin / 2.0, SL, ST, P_fin)
    tube_outer = occ.addCylinder(0.0, 0.0, -P_fin / 2.0, 0.0, 0.0, P_fin, R_tube_out)
    tube_inner = occ.addCylinder(0.0, 0.0, -P_fin / 2.0, 0.0, 0.0, P_fin, R_tube_in)
    fin = occ.addCylinder(0.0, 0.0, -t_fin / 2.0, 0.0, 0.0, t_fin, R_fin)

    # Merge the tube wall and fin, then subtract the inner bore.
    union_entities, _ = occ.fuse([(3, tube_outer)], [(3, fin)])
    solid_tag = union_entities[0][1]
    cut_entities, _ = occ.cut([(3, solid_tag)], [(3, tube_inner)])
    solid_tag = cut_entities[0][1]

    # Fragment fluid and solid to enforce a conformal interface.
    frag_entities, _ = occ.fragment([(3, fluid_box)], [(3, solid_tag)])
    occ.synchronize()

    # Identify the resulting fluid and solid volumes by bounding box extents.
    volumes = gmsh.model.getEntities(dim=3)
    fluid_volume = None
    solid_volume = None
    for dim, tag in volumes:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
        if abs(dx - SL) < 1e-4 and abs(dy - ST) < 1e-4:
            fluid_volume = tag
        else:
            solid_volume = tag
    if fluid_volume is None or solid_volume is None:
        raise RuntimeError("Failed to classify fluid and solid volumes after fragment operation.")

    return fluid_volume, solid_volume


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
    if len(min_surfaces) != len(max_surfaces):
        raise RuntimeError("Periodic surface counts do not match.")
    min_sorted = _pair_by_centroid(min_surfaces, axis=axis)
    max_sorted = _pair_by_centroid(max_surfaces, axis=axis)
    affine = [
        1.0, 0.0, 0.0, translation[0],
        0.0, 1.0, 0.0, translation[1],
        0.0, 0.0, 1.0, translation[2],
        0.0, 0.0, 0.0, 1.0,
    ]
    gmsh.model.mesh.setPeriodic(2, max_sorted, min_sorted, affine)


def build_model(params):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("staggered_finned_tube")

    fluid_volume, solid_volume = _add_geometry(params)

    fluid_surfaces = _surface_sets_for_volume(fluid_volume)
    solid_surfaces = _surface_sets_for_volume(solid_volume)
    interface_surfaces = fluid_surfaces.intersection(solid_surfaces)

    SL = params["SL"]
    ST = params["ST"]
    P_fin = params["P_fin"]

    inlet = _select_surfaces_by_coord(fluid_surfaces, "x", -SL / 2.0)
    outlet = _select_surfaces_by_coord(fluid_surfaces, "x", SL / 2.0)
    bottom = _select_surfaces_by_coord(fluid_surfaces, "y", -ST / 2.0)
    top = _select_surfaces_by_coord(fluid_surfaces, "y", ST / 2.0)
    back = _select_surfaces_by_coord(fluid_surfaces, "z", -P_fin / 2.0)
    front = _select_surfaces_by_coord(fluid_surfaces, "z", P_fin / 2.0)

    if not all([inlet, outlet, bottom, top, back, front]):
        raise RuntimeError("Failed to locate all periodic boundary surfaces.")

    _set_periodicity(inlet, outlet, (SL, 0.0, 0.0), axis="x")
    _set_periodicity(bottom, top, (0.0, ST, 0.0), axis="y")
    _set_periodicity(back, front, (0.0, 0.0, P_fin), axis="z")

    gmsh.model.addPhysicalGroup(3, [fluid_volume], name="fluid")
    gmsh.model.addPhysicalGroup(3, [solid_volume], name="solid")

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
