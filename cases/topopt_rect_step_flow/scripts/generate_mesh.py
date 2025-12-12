#!/usr/bin/env python3
"""
Generate 3D mesh for topology-optimized obstacle in a rectangular channel.

Pipeline:
  1) Import STEP obstacle (units: mm in the original export)
  2) Scale geometry by 1e-3 (mm -> m)
  3) Build fluid domain as: fluid = Box(0..L, 0..H, 0..T) \\ obstacle
  4) Tag boundaries with Physical Groups for OpenFOAM:
       inlet, outlet, bottomWall, topWall, front, back, obstacle, fluid
  5) Write mesh.msh (Msh2.2 recommended for gmshToFoam)

Usage:
  python scripts/generate_mesh.py --output mesh.msh
"""

from __future__ import annotations

import argparse
from pathlib import Path

import gmsh


def _case_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def _add_physical_group(dim: int, entity_tags: list[int], name: str) -> int | None:
    if not entity_tags:
        return None
    phys = gmsh.model.addPhysicalGroup(dim, entity_tags)
    gmsh.model.setPhysicalName(dim, phys, name)
    return phys


def _select_surfaces_by_plane(
    surface_tags: list[int],
    *,
    axis: str,
    value: float,
    tol: float,
) -> list[int]:
    """Select planar surfaces by bounding-box test."""
    if axis not in {"x", "y", "z"}:
        raise ValueError(f"axis must be one of x/y/z, got: {axis}")

    selected: list[int] = []
    for s in surface_tags:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, s)
        if axis == "x":
            cmin, cmax = xmin, xmax
        elif axis == "y":
            cmin, cmax = ymin, ymax
        else:
            cmin, cmax = zmin, zmax

        if abs(cmin - value) < tol and abs(cmax - value) < tol:
            selected.append(s)
    return selected


def build_model(
    *,
    step_path: Path,
    scale: float,
    L: float,
    H: float,
    T: float,
    mesh_size: float,
    interface_size: float | None,
    msh_version: float,
) -> None:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", msh_version)
    gmsh.model.add("topopt_rect_step_flow")

    occ = gmsh.model.occ

    # --- Import obstacle STEP (in mm) and scale to meters ---
    if not step_path.exists():
        raise FileNotFoundError(f"STEP not found: {step_path}")

    obstacle_dimtags = occ.importShapes(str(step_path))
    occ.synchronize()

    # Unit conversion: mm -> m
    occ.dilate(obstacle_dimtags, 0.0, 0.0, 0.0, scale, scale, scale)
    occ.synchronize()

    # Keep only volume entities as cut tools
    obstacle_vols = [(dim, tag) for (dim, tag) in obstacle_dimtags if dim == 3]
    if not obstacle_vols:
        # Some STEP importers may return only surfaces; try to collect volumes anyway.
        obstacle_vols = gmsh.model.getEntities(dim=3)
    if not obstacle_vols:
        raise RuntimeError("No volume entities found in imported STEP.")

    # --- Build fluid box (dimensions provided in original STEP units) ---
    Lm, Hm, Tm = L * scale, H * scale, T * scale
    box_tag = occ.addBox(0.0, 0.0, 0.0, Lm, Hm, Tm)
    occ.synchronize()

    # --- Boolean cut: fluid = box \\ obstacle ---
    cut_out, _ = occ.cut([(3, box_tag)], obstacle_vols, removeObject=True, removeTool=True)
    occ.synchronize()

    fluid_vols = [(dim, tag) for (dim, tag) in cut_out if dim == 3]
    if not fluid_vols:
        # Fallback: pick remaining volumes
        fluid_vols = gmsh.model.getEntities(dim=3)
    if not fluid_vols:
        raise RuntimeError("No fluid volume produced after boolean cut.")

    # Get all fluid volume tags
    fluid_tags = [tag for (_, tag) in fluid_vols]

    # --- Boundary classification (by domain planes) ---
    # Collect surfaces from ALL fluid volumes (handle fragmented domains)
    surfaces = set()
    for ftag in fluid_tags:
        for (dim, tag) in gmsh.model.getBoundary([(3, ftag)], oriented=False):
            if dim == 2:
                surfaces.add(tag)
    surfaces = list(surfaces)
    if not surfaces:
        raise RuntimeError("Failed to extract boundary surfaces from fluid volume.")

    # OCC geometry has ~1e-7 tolerance, use larger tolerance for plane detection
    tol = max(1e-6, 1e-4 * min(Lm, Hm, Tm))

    inlet = _select_surfaces_by_plane(surfaces, axis="x", value=0.0, tol=tol)
    outlet = _select_surfaces_by_plane(surfaces, axis="x", value=Lm, tol=tol)
    bottom = _select_surfaces_by_plane(surfaces, axis="y", value=0.0, tol=tol)
    top = _select_surfaces_by_plane(surfaces, axis="y", value=Hm, tol=tol)
    front = _select_surfaces_by_plane(surfaces, axis="z", value=0.0, tol=tol)
    back = _select_surfaces_by_plane(surfaces, axis="z", value=Tm, tol=tol)

    known = set(inlet + outlet + bottom + top + front + back)
    obstacle = [s for s in surfaces if s not in known]

    # --- Physical groups ---
    _add_physical_group(3, fluid_tags, "fluid")
    _add_physical_group(2, inlet, "inlet")
    _add_physical_group(2, outlet, "outlet")
    _add_physical_group(2, bottom, "bottomWall")
    _add_physical_group(2, top, "topWall")
    _add_physical_group(2, front, "front")
    _add_physical_group(2, back, "back")
    _add_physical_group(2, obstacle, "obstacle")

    # --- Mesh sizing ---
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size * 0.5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT (robust tetra meshing)

    # Optional refinement near obstacle
    if interface_size is not None and obstacle:
        fd = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fd, "FacesList", obstacle)
        ft = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(ft, "InField", fd)
        gmsh.model.mesh.field.setNumber(ft, "SizeMin", interface_size)
        gmsh.model.mesh.field.setNumber(ft, "SizeMax", mesh_size)
        gmsh.model.mesh.field.setNumber(ft, "DistMin", 2.0 * interface_size)
        gmsh.model.mesh.field.setNumber(ft, "DistMax", 10.0 * mesh_size)
        gmsh.model.mesh.field.setAsBackgroundMesh(ft)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate Gmsh mesh from TopOpt STEP obstacle.")
    p.add_argument("--step", type=str, default="", help="Input STEP path (default: case/geometry/solid.step)")
    p.add_argument("--output", type=str, default="mesh.msh", help="Output .msh filename (default: mesh.msh)")

    # Unit conversion (STEP uses mm)
    p.add_argument("--scale", type=float, default=1.0e-3, help="Uniform scale applied to STEP (default: 1e-3)")

    # Channel dimensions in original STEP units (mm)
    p.add_argument("--L", type=float, default=2.2, help="Channel length in STEP units (default: 2.2)")
    p.add_argument("--H", type=float, default=0.41, help="Channel height in STEP units (default: 0.41)")
    p.add_argument("--T", type=float, default=10.0, help="Channel thickness in STEP units (default: 10.0)")

    # Mesh sizes in meters (after scaling)
    p.add_argument("--mesh-size", type=float, default=5.0e-5, help="Global mesh size [m] (default: 5e-5)")
    p.add_argument(
        "--interface-size",
        type=float,
        default=2.0e-5,
        help="Refinement size near obstacle [m] (default: 2e-5, set 0 to disable)",
    )
    p.add_argument("--msh-version", type=float, default=2.2, help="Msh file version (default: 2.2)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    case_dir = _case_dir()

    step_path = Path(args.step) if args.step else (case_dir / "geometry" / "solid.step")
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = case_dir / output_path

    interface_size = None if float(args.interface_size) <= 0.0 else float(args.interface_size)

    try:
        build_model(
            step_path=step_path,
            scale=float(args.scale),
            L=float(args.L),
            H=float(args.H),
            T=float(args.T),
            mesh_size=float(args.mesh_size),
            interface_size=interface_size,
            msh_version=float(args.msh_version),
        )
        gmsh.write(str(output_path))
        print(f"[OK] Wrote mesh: {output_path}")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()


