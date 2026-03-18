#!/usr/bin/env python3
"""Generate a Gmsh model from indoor_cfd_scene_v1 JSON.

Supported room geometries:
- legacy rectangular room via room.size
- composite room via room.blocks (1-2 joined rectangular blocks)

Current scope:
- box obstacles
- wall-aligned rectangular inlet/outlet openings
- fluid volume = room union - openings - obstacles
- physical groups for fluid, inlet, outlet, roomWalls, obstacleWalls
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import gmsh

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from validate_indoor_scene import validate_scene  # noqa: E402

EPS = 1e-7
OPENING_DEPTH_FACTOR = 1e-3


class SceneToGmshError(RuntimeError):
    pass


def get_room_blocks(scene: dict) -> list[dict]:
    room = scene["room"]
    if "blocks" in room:
        return room["blocks"]
    size = room["size"]
    return [{
        "origin": {"x": 0.0, "y": 0.0, "z": 0.0},
        "size": {"dx": size["Lx"], "dy": size["Ly"], "dz": size["Lz"]},
    }]


def overall_room_bounds(room_blocks: list[dict]) -> dict[str, float]:
    xmax = max(block["origin"]["x"] + block["size"]["dx"] for block in room_blocks)
    ymax = max(block["origin"]["y"] + block["size"]["dy"] for block in room_blocks)
    zmax = max(block["origin"]["z"] + block["size"]["dz"] for block in room_blocks)
    return {"Lx": xmax, "Ly": ymax, "Lz": zmax}


def wall_supporting_blocks(room_blocks: list[dict], wall: str, opening: dict) -> list[dict]:
    u = opening["center"]["u"]
    v = opening["center"]["v"]
    du = opening["size"]["du"]
    dv = opening["size"]["dv"]
    supported: list[dict] = []
    for block in room_blocks:
        ox, oy, oz = block["origin"]["x"], block["origin"]["y"], block["origin"]["z"]
        dx, dy, dz = block["size"]["dx"], block["size"]["dy"], block["size"]["dz"]
        if wall in {"west", "east"}:
            ok = (
                u - du / 2 >= oy - EPS and u + du / 2 <= oy + dy + EPS
                and v - dv / 2 >= oz - EPS and v + dv / 2 <= oz + dz + EPS
            )
        elif wall in {"south", "north"}:
            ok = (
                u - du / 2 >= ox - EPS and u + du / 2 <= ox + dx + EPS
                and v - dv / 2 >= oz - EPS and v + dv / 2 <= oz + dz + EPS
            )
        else:
            ok = (
                u - du / 2 >= ox - EPS and u + du / 2 <= ox + dx + EPS
                and v - dv / 2 >= oy - EPS and v + dv / 2 <= oy + dy + EPS
            )
        if ok:
            supported.append(block)
    return supported


def opening_global_bounds(room_blocks: list[dict], opening: dict) -> dict[str, float]:
    """Convert wall-local opening definition to global bounding box.

    For composite rooms, an opening is attached to the exposed face of the block(s)
    whose local wall coordinates contain the opening rectangle.
    """
    room_bounds = overall_room_bounds(room_blocks)
    wall = opening["wall"]
    u = opening["center"]["u"]
    v = opening["center"]["v"]
    du = opening["size"]["du"]
    dv = opening["size"]["dv"]
    depth = max(room_bounds["Lx"], room_bounds["Ly"], room_bounds["Lz"]) * OPENING_DEPTH_FACTOR

    supporting = wall_supporting_blocks(room_blocks, wall, opening)
    if not supporting:
        raise SceneToGmshError(f"Opening {opening.get('id', '<unknown>')} is not supported by any room block on wall {wall}")

    if wall == "west":
        x_plane = min(block["origin"]["x"] for block in supporting)
        return {
            "xmin": x_plane,
            "xmax": x_plane + depth,
            "ymin": u - du / 2,
            "ymax": u + du / 2,
            "zmin": v - dv / 2,
            "zmax": v + dv / 2,
        }
    if wall == "east":
        x_plane = max(block["origin"]["x"] + block["size"]["dx"] for block in supporting)
        return {
            "xmin": x_plane - depth,
            "xmax": x_plane,
            "ymin": u - du / 2,
            "ymax": u + du / 2,
            "zmin": v - dv / 2,
            "zmax": v + dv / 2,
        }
    if wall == "south":
        y_plane = min(block["origin"]["y"] for block in supporting)
        return {
            "xmin": u - du / 2,
            "xmax": u + du / 2,
            "ymin": y_plane,
            "ymax": y_plane + depth,
            "zmin": v - dv / 2,
            "zmax": v + dv / 2,
        }
    if wall == "north":
        y_plane = max(block["origin"]["y"] + block["size"]["dy"] for block in supporting)
        return {
            "xmin": u - du / 2,
            "xmax": u + du / 2,
            "ymin": y_plane - depth,
            "ymax": y_plane,
            "zmin": v - dv / 2,
            "zmax": v + dv / 2,
        }
    if wall == "floor":
        z_plane = min(block["origin"]["z"] for block in supporting)
        return {
            "xmin": u - du / 2,
            "xmax": u + du / 2,
            "ymin": v - dv / 2,
            "ymax": v + dv / 2,
            "zmin": z_plane,
            "zmax": z_plane + depth,
        }
    if wall == "ceiling":
        z_plane = max(block["origin"]["z"] + block["size"]["dz"] for block in supporting)
        return {
            "xmin": u - du / 2,
            "xmax": u + du / 2,
            "ymin": v - dv / 2,
            "ymax": v + dv / 2,
            "zmin": z_plane - depth,
            "zmax": z_plane,
        }
    raise SceneToGmshError(f"Unsupported wall: {wall}")


def add_box_from_bounds(factory, bounds: dict[str, float]) -> int:
    return factory.addBox(
        bounds["xmin"],
        bounds["ymin"],
        bounds["zmin"],
        bounds["xmax"] - bounds["xmin"],
        bounds["ymax"] - bounds["ymin"],
        bounds["zmax"] - bounds["zmin"],
    )


def load_scene(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        scene = json.load(f)
    report = validate_scene(scene)
    if not report.ok:
        raise SceneToGmshError(
            "Scene validation failed:\n- " + "\n- ".join(report.errors)
        )
    return scene


def classify_surface(room_bounds: dict[str, float], opening_boxes: list[dict], obstacle_boxes: list[dict], bbox: tuple[float, ...]) -> tuple[str | None, str | None]:
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    xmid = 0.5 * (xmin + xmax)
    ymid = 0.5 * (ymin + ymax)
    zmid = 0.5 * (zmin + zmax)
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    for op in opening_boxes:
        b = op["bounds"]
        wall = op["wall"]
        if wall in {"west", "east"}:
            matches = (
                abs(dy - (b["ymax"] - b["ymin"])) < 1e-4
                and abs(dz - (b["zmax"] - b["zmin"])) < 1e-4
                and ymid >= b["ymin"] - EPS and ymid <= b["ymax"] + EPS
                and zmid >= b["zmin"] - EPS and zmid <= b["zmax"] + EPS
                and (
                    abs(xmin - b["xmin"]) < 1e-4 or abs(xmax - b["xmax"]) < 1e-4
                    or abs(xmid - 0.5 * (b["xmin"] + b["xmax"])) < 1e-4
                )
            )
        elif wall in {"south", "north"}:
            matches = (
                abs(dx - (b["xmax"] - b["xmin"])) < 1e-4
                and abs(dz - (b["zmax"] - b["zmin"])) < 1e-4
                and xmid >= b["xmin"] - EPS and xmid <= b["xmax"] + EPS
                and zmid >= b["zmin"] - EPS and zmid <= b["zmax"] + EPS
                and (
                    abs(ymin - b["ymin"]) < 1e-4 or abs(ymax - b["ymax"]) < 1e-4
                    or abs(ymid - 0.5 * (b["ymin"] + b["ymax"])) < 1e-4
                )
            )
        else:
            matches = (
                abs(dx - (b["xmax"] - b["xmin"])) < 1e-4
                and abs(dy - (b["ymax"] - b["ymin"])) < 1e-4
                and xmid >= b["xmin"] - EPS and xmid <= b["xmax"] + EPS
                and ymid >= b["ymin"] - EPS and ymid <= b["ymax"] + EPS
                and (
                    abs(zmin - b["zmin"]) < 1e-4 or abs(zmax - b["zmax"]) < 1e-4
                    or abs(zmid - 0.5 * (b["zmin"] + b["zmax"])) < 1e-4
                )
            )
        if matches:
            return op["type"], op["id"]

    for obs in obstacle_boxes:
        b = obs["bounds"]
        on_x_face = (
            abs(xmin - b["xmin"]) < 1e-4 and abs(xmax - b["xmin"]) < 1e-4
        ) or (
            abs(xmin - b["xmax"]) < 1e-4 and abs(xmax - b["xmax"]) < 1e-4
        )
        on_y_face = (
            abs(ymin - b["ymin"]) < 1e-4 and abs(ymax - b["ymin"]) < 1e-4
        ) or (
            abs(ymin - b["ymax"]) < 1e-4 and abs(ymax - b["ymax"]) < 1e-4
        )
        on_z_face = (
            abs(zmin - b["zmin"]) < 1e-4 and abs(zmax - b["zmin"]) < 1e-4
        ) or (
            abs(zmin - b["zmax"]) < 1e-4 and abs(zmax - b["zmax"]) < 1e-4
        )
        within_x = xmin >= b["xmin"] - EPS and xmax <= b["xmax"] + EPS
        within_y = ymin >= b["ymin"] - EPS and ymax <= b["ymax"] + EPS
        within_z = zmin >= b["zmin"] - EPS and zmax <= b["zmax"] + EPS
        if (on_x_face and within_y and within_z) or (on_y_face and within_x and within_z) or (on_z_face and within_x and within_y):
            return "obstacleWalls", obs["id"]

    if abs(xmin) < EPS and abs(xmax) < EPS:
        return "roomWalls", "west"
    if abs(xmin - room_bounds["Lx"]) < EPS and abs(xmax - room_bounds["Lx"]) < EPS:
        return "roomWalls", "east"
    if abs(ymin) < EPS and abs(ymax) < EPS:
        return "roomWalls", "south"
    if abs(ymin - room_bounds["Ly"]) < EPS and abs(ymax - room_bounds["Ly"]) < EPS:
        return "roomWalls", "north"
    if abs(zmin) < EPS and abs(zmax) < EPS:
        return "roomWalls", "floor"
    if abs(zmin - room_bounds["Lz"]) < EPS and abs(zmax - room_bounds["Lz"]) < EPS:
        return "roomWalls", "ceiling"

    return "roomWalls", None


def build_model(scene: dict, mesh_size: float, model_name: str) -> dict[str, list[int] | int]:
    room_blocks = get_room_blocks(scene)
    room_bounds = overall_room_bounds(room_blocks)
    factory = gmsh.model.occ
    gmsh.model.add(model_name)

    block_entities = []
    for block in room_blocks:
        tag = factory.addBox(
            block["origin"]["x"],
            block["origin"]["y"],
            block["origin"]["z"],
            block["size"]["dx"],
            block["size"]["dy"],
            block["size"]["dz"],
        )
        block_entities.append((3, tag))

    room_entities = block_entities
    if len(block_entities) > 1:
        fused, _ = factory.fuse([block_entities[0]], block_entities[1:], removeObject=True, removeTool=True)
        room_entities = [(dim, tag) for dim, tag in fused if dim == 3]
    if not room_entities:
        raise SceneToGmshError("Failed to build room volume from room.blocks")

    opening_tools = []
    opening_boxes = []
    for opening in scene["openings"]:
        bounds = opening_global_bounds(room_blocks, opening)
        tag = add_box_from_bounds(factory, bounds)
        opening_tools.append((3, tag))
        opening_boxes.append({
            "id": opening["id"],
            "type": opening["type"],
            "wall": opening["wall"],
            "bounds": bounds,
        })

    if opening_tools:
        cut_room, _ = factory.cut(room_entities, opening_tools, removeObject=True, removeTool=False)
        room_entities = [(dim, tag) for dim, tag in cut_room if dim == 3]
        if not room_entities:
            raise SceneToGmshError("Failed to cut room with opening tools")

    obstacle_tools = []
    obstacle_boxes = []
    for obstacle in scene["obstacles"]:
        min_corner = obstacle["min"]
        size = obstacle["size"]
        tag = factory.addBox(
            min_corner["x"],
            min_corner["y"],
            min_corner["z"],
            size["dx"],
            size["dy"],
            size["dz"],
        )
        obstacle_tools.append((3, tag))
        obstacle_boxes.append({
            "id": obstacle["id"],
            "bounds": {
                "xmin": min_corner["x"],
                "xmax": min_corner["x"] + size["dx"],
                "ymin": min_corner["y"],
                "ymax": min_corner["y"] + size["dy"],
                "zmin": min_corner["z"],
                "zmax": min_corner["z"] + size["dz"],
            },
        })

    if obstacle_tools:
        cut_out, _ = factory.cut(room_entities, obstacle_tools, removeObject=True, removeTool=False)
        if not cut_out:
            raise SceneToGmshError("Boolean cut failed; no fluid volume produced")
        fluid_entities = [(dim, tag) for dim, tag in cut_out if dim == 3]
    else:
        fluid_entities = list(room_entities)
    if len(fluid_entities) > 1:
        fused, _ = factory.fuse([fluid_entities[0]], fluid_entities[1:], removeObject=True, removeTool=True)
        fluid_entities = [(dim, tag) for dim, tag in fused if dim == 3]

    factory.synchronize()

    fluid_volumes = [tag for dim, tag in fluid_entities]
    if not fluid_volumes:
        raise SceneToGmshError("No 3D fluid volume found after boolean cut")

    boundary = gmsh.model.getBoundary([(3, tag) for tag in fluid_volumes], oriented=False, recursive=False)
    surface_tags = sorted({tag for dim, tag in boundary if dim == 2})

    inlet_surfs: list[int] = []
    outlet_surfs: list[int] = []
    room_wall_surfs: list[int] = []
    obstacle_wall_surfs: list[int] = []

    for tag in surface_tags:
        bbox = gmsh.model.getBoundingBox(2, tag)
        group, _ = classify_surface(room_bounds, opening_boxes, obstacle_boxes, bbox)
        if group == "inlet":
            inlet_surfs.append(tag)
        elif group == "outlet":
            outlet_surfs.append(tag)
        elif group == "obstacleWalls":
            obstacle_wall_surfs.append(tag)
        else:
            room_wall_surfs.append(tag)

    gmsh.model.addPhysicalGroup(3, fluid_volumes, name="fluid")
    if inlet_surfs:
        gmsh.model.addPhysicalGroup(2, inlet_surfs, name="inlet")
    if outlet_surfs:
        gmsh.model.addPhysicalGroup(2, outlet_surfs, name="outlet")
    if room_wall_surfs:
        gmsh.model.addPhysicalGroup(2, room_wall_surfs, name="roomWalls")
    if obstacle_wall_surfs:
        gmsh.model.addPhysicalGroup(2, obstacle_wall_surfs, name="obstacleWalls")

    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    return {
        "fluid_volumes": fluid_volumes,
        "inlet_surfs": inlet_surfs,
        "outlet_surfs": outlet_surfs,
        "room_wall_surfs": room_wall_surfs,
        "obstacle_wall_surfs": obstacle_wall_surfs,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Gmsh geometry/mesh from indoor CFD scene JSON")
    parser.add_argument("scene_json", type=Path, help="Validated indoor scene JSON path")
    parser.add_argument("-o", "--output", type=Path, default=Path("generated/indoor_scene.msh"), help="Output .msh path")
    parser.add_argument("--geo", type=Path, default=None, help="Optional .geo_unrolled output path")
    parser.add_argument("--mesh-size", type=float, default=0.25, help="Global target mesh size in meters")
    parser.add_argument("--no-mesh", action="store_true", help="Only build CAD and write .geo_unrolled")
    parser.add_argument("--binary", action="store_true", help="Write binary .msh when meshing")
    args = parser.parse_args()

    scene = load_scene(args.scene_json)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.geo is not None:
        args.geo.parent.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    try:
        summary = build_model(scene, args.mesh_size, model_name=args.scene_json.stem)
        if args.geo is not None:
            gmsh.write(str(args.geo))
        if not args.no_mesh:
            gmsh.option.setNumber("Mesh.Binary", 1 if args.binary else 0)
            gmsh.model.mesh.generate(3)
            gmsh.write(str(args.output))
    finally:
        gmsh.finalize()

    print(json.dumps({
        "ok": True,
        "scene": str(args.scene_json),
        "output": None if args.no_mesh else str(args.output),
        "geo": str(args.geo) if args.geo is not None else None,
        "mesh_size": args.mesh_size,
        "num_fluid_volumes": len(summary["fluid_volumes"]),
        "num_inlet_surfaces": len(summary["inlet_surfs"]),
        "num_outlet_surfaces": len(summary["outlet_surfs"]),
        "num_room_wall_surfaces": len(summary["room_wall_surfs"]),
        "num_obstacle_wall_surfaces": len(summary["obstacle_wall_surfs"]),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
