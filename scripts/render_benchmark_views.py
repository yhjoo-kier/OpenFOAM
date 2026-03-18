#!/usr/bin/env python3
"""Render benchmark input views from scene JSON for the image-to-CFD paper.

Outputs canonical 2D image assets under:
  benchmark/renderings/<case_name>/

Initial supported views:
- perspective : shaded 3D geometry view
- wireframe   : structure-only 3D view
- birdseye    : elevated near-top 3D view
- floorplan   : 2D top plan view
- section     : canonical vertical section view aligned to the opening pair when possible

This renderer is geometry-only on purpose: these images are intended as benchmark
inputs for the VLM scene-generation stage, distinct from the reference CFD outputs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SCENES_DIR = PROJECT_ROOT / "benchmark" / "scenes"
DEFAULT_RENDER_ROOT = PROJECT_ROOT / "benchmark" / "renderings"


def load_scene(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def room_blocks(scene: dict) -> list[dict]:
    room = scene["room"]
    if "blocks" in room:
        return room["blocks"]
    size = room["size"]
    return [{
        "origin": {"x": 0.0, "y": 0.0, "z": 0.0},
        "size": {"dx": size["Lx"], "dy": size["Ly"], "dz": size["Lz"]},
    }]


def room_bounds(scene: dict) -> dict[str, float]:
    blocks = room_blocks(scene)
    return {
        "Lx": max(b["origin"]["x"] + b["size"]["dx"] for b in blocks),
        "Ly": max(b["origin"]["y"] + b["size"]["dy"] for b in blocks),
        "Lz": max(b["origin"]["z"] + b["size"]["dz"] for b in blocks),
    }


def cuboid_faces(x0, y0, z0, dx, dy, dz):
    x1, y1, z1 = x0 + dx, y0 + dy, z0 + dz
    return [
        [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0]],
        [[x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]],
        [[x0, y0, z0], [x1, y0, z0], [x1, y0, z1], [x0, y0, z1]],
        [[x0, y1, z0], [x1, y1, z0], [x1, y1, z1], [x0, y1, z1]],
        [[x0, y0, z0], [x0, y1, z0], [x0, y1, z1], [x0, y0, z1]],
        [[x1, y0, z0], [x1, y1, z0], [x1, y1, z1], [x1, y0, z1]],
    ]


def set_axes_equal(ax, bounds, scale=0.58):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    xmid, ymid, zmid = (xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2
    radius = scale * max(xmax - xmin, ymax - ymin, zmax - zmin)
    ax.set_xlim(xmid - radius, xmid + radius)
    ax.set_ylim(ymid - radius, ymid + radius)
    ax.set_zlim(zmid - radius, zmid + radius)
    ax.set_box_aspect((max(xmax - xmin, 1e-6), max(ymax - ymin, 1e-6), max(zmax - zmin, 1e-6)))


def style_3d_axes(ax):
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.line.set_color((1, 1, 1, 0))


def draw_room_3d(ax, scene: dict, *, face_alpha: float, face_color: str, edge_color: str, linewidth: float):
    for block in room_blocks(scene):
        x0, y0, z0 = block["origin"]["x"], block["origin"]["y"], block["origin"]["z"]
        dx, dy, dz = block["size"]["dx"], block["size"]["dy"], block["size"]["dz"]
        faces = cuboid_faces(x0, y0, z0, dx, dy, dz)
        if face_alpha > 0:
            ax.add_collection3d(
                Poly3DCollection(
                    faces,
                    facecolor=face_color,
                    edgecolor="none",
                    alpha=face_alpha,
                )
            )
        for face in faces:
            closed = face + [face[0]]
            ax.plot(*zip(*closed), color=edge_color, lw=linewidth, alpha=0.7)


def draw_obstacles_3d(ax, scene: dict, *, face_color: str, face_alpha: float, edge_color: str, linewidth: float):
    for obs in scene.get("obstacles", []):
        x0, y0, z0 = obs["min"]["x"], obs["min"]["y"], obs["min"]["z"]
        dx, dy, dz = obs["size"]["dx"], obs["size"]["dy"], obs["size"]["dz"]
        faces = cuboid_faces(x0, y0, z0, dx, dy, dz)
        if face_alpha > 0:
            ax.add_collection3d(
                Poly3DCollection(
                    faces,
                    facecolor=face_color,
                    edgecolor=edge_color,
                    linewidths=linewidth,
                    alpha=face_alpha,
                )
            )
        else:
            for face in faces:
                closed = face + [face[0]]
                ax.plot(*zip(*closed), color=edge_color, lw=linewidth, alpha=0.9)


def draw_openings_3d(ax, scene: dict):
    room = room_bounds(scene)
    for op in scene.get("openings", []):
        color = "#377eb8" if op["type"] == "inlet" else "#e41a1c"
        u, v = op["center"]["u"], op["center"]["v"]
        du, dv = op["size"]["du"], op["size"]["dv"]
        if op["wall"] == "west":
            poly = [[[0, u - du / 2, v - dv / 2], [0, u + du / 2, v - dv / 2], [0, u + du / 2, v + dv / 2], [0, u - du / 2, v + dv / 2]]]
        elif op["wall"] == "east":
            x = room["Lx"]
            poly = [[[x, u - du / 2, v - dv / 2], [x, u + du / 2, v - dv / 2], [x, u + du / 2, v + dv / 2], [x, u - du / 2, v + dv / 2]]]
        elif op["wall"] == "south":
            poly = [[[u - du / 2, 0, v - dv / 2], [u + du / 2, 0, v - dv / 2], [u + du / 2, 0, v + dv / 2], [u - du / 2, 0, v + dv / 2]]]
        elif op["wall"] == "north":
            y = room["Ly"]
            poly = [[[u - du / 2, y, v - dv / 2], [u + du / 2, y, v - dv / 2], [u + du / 2, y, v + dv / 2], [u - du / 2, y, v + dv / 2]]]
        else:
            continue
        ax.add_collection3d(Poly3DCollection(poly, facecolor=color, edgecolor="none", alpha=0.96))


def draw_floorplan(ax, scene: dict):
    room = room_bounds(scene)
    for block in room_blocks(scene):
        ox, oy = block["origin"]["x"], block["origin"]["y"]
        dx, dy = block["size"]["dx"], block["size"]["dy"]
        ax.add_patch(Rectangle((ox, oy), dx, dy, facecolor="#f7f8fb", edgecolor="#111111", linewidth=2.0))

    for obs in scene.get("obstacles", []):
        x, y = obs["min"]["x"], obs["min"]["y"]
        dx, dy = obs["size"]["dx"], obs["size"]["dy"]
        ax.add_patch(Rectangle((x, y), dx, dy, facecolor="#8f98a3", edgecolor="#30343a", alpha=0.95))

    for op in scene.get("openings", []):
        color = "#377eb8" if op["type"] == "inlet" else "#e41a1c"
        u, du = op["center"]["u"], op["size"]["du"]
        if op["wall"] == "west":
            ax.plot([0, 0], [u - du / 2, u + du / 2], color=color, linewidth=5)
        elif op["wall"] == "east":
            ax.plot([room["Lx"], room["Lx"]], [u - du / 2, u + du / 2], color=color, linewidth=5)
        elif op["wall"] == "south":
            ax.plot([u - du / 2, u + du / 2], [0, 0], color=color, linewidth=5)
        elif op["wall"] == "north":
            ax.plot([u - du / 2, u + du / 2], [room["Ly"], room["Ly"]], color=color, linewidth=5)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    ax.set_xlim(-0.04 * room["Lx"], 1.04 * room["Lx"])
    ax.set_ylim(-0.04 * room["Ly"], 1.04 * room["Ly"])


def infer_section_plane(scene: dict) -> tuple[str, float]:
    room = room_bounds(scene)
    openings = scene.get("openings", [])
    walls = {op.get("wall") for op in openings}
    if "west" in walls and "east" in walls:
        coords = [op["center"]["u"] for op in openings if op.get("wall") in {"west", "east"}]
        plane = sum(coords) / len(coords) if coords else room["Ly"] * 0.5
        return "y", plane
    if "south" in walls and "north" in walls:
        coords = [op["center"]["u"] for op in openings if op.get("wall") in {"south", "north"}]
        plane = sum(coords) / len(coords) if coords else room["Lx"] * 0.5
        return "x", plane
    room_kind = "composite" if "blocks" in scene.get("room", {}) else "rectangular"
    return ("y", room["Ly"] * 0.5) if room_kind == "composite" else ("x", room["Lx"] * 0.5)


def draw_section(ax, scene: dict):
    room = room_bounds(scene)
    axis, plane = infer_section_plane(scene)
    horizontal_extent = room["Lx"] if axis == "y" else room["Ly"]

    for block in room_blocks(scene):
        bx0, by0, bz0 = block["origin"]["x"], block["origin"]["y"], block["origin"]["z"]
        bx1 = bx0 + block["size"]["dx"]
        by1 = by0 + block["size"]["dy"]
        bz1 = bz0 + block["size"]["dz"]
        intersects = (by0 <= plane <= by1) if axis == "y" else (bx0 <= plane <= bx1)
        if not intersects:
            continue
        if axis == "y":
            ax.add_patch(Rectangle((bx0, bz0), bx1 - bx0, bz1 - bz0, facecolor="#f7f8fb", edgecolor="#111111", linewidth=2.0))
        else:
            ax.add_patch(Rectangle((by0, bz0), by1 - by0, bz1 - bz0, facecolor="#f7f8fb", edgecolor="#111111", linewidth=2.0))

    for obs in scene.get("obstacles", []):
        x0, y0, z0 = obs["min"]["x"], obs["min"]["y"], obs["min"]["z"]
        x1 = x0 + obs["size"]["dx"]
        y1 = y0 + obs["size"]["dy"]
        z1 = z0 + obs["size"]["dz"]
        intersects = (y0 <= plane <= y1) if axis == "y" else (x0 <= plane <= x1)
        if not intersects:
            continue
        if axis == "y":
            ax.add_patch(Rectangle((x0, z0), x1 - x0, z1 - z0, facecolor="#8f98a3", edgecolor="#30343a", alpha=0.95))
        else:
            ax.add_patch(Rectangle((y0, z0), y1 - y0, z1 - z0, facecolor="#8f98a3", edgecolor="#30343a", alpha=0.95))

    for op in scene.get("openings", []):
        color = "#377eb8" if op["type"] == "inlet" else "#e41a1c"
        u, v = op["center"]["u"], op["center"]["v"]
        du, dv = op["size"]["du"], op["size"]["dv"]
        if axis == "y" and op["wall"] in {"west", "east"}:
            if not (u - du / 2 <= plane <= u + du / 2):
                continue
            x = 0 if op["wall"] == "west" else room["Lx"]
            ax.plot([x, x], [v - dv / 2, v + dv / 2], color=color, linewidth=5)
        elif axis == "x" and op["wall"] in {"south", "north"}:
            if not (u - du / 2 <= plane <= u + du / 2):
                continue
            y = 0 if op["wall"] == "south" else room["Ly"]
            ax.plot([y, y], [v - dv / 2, v + dv / 2], color=color, linewidth=5)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    ax.set_xlim(-0.04 * horizontal_extent, 1.04 * horizontal_extent)
    ax.set_ylim(-0.05 * room["Lz"], 1.05 * room["Lz"])


def render_3d_scene(scene: dict, output: Path, *, elev: float, azim: float, room_face_alpha: float, obstacle_face_alpha: float, wireframe: bool = False):
    room = room_bounds(scene)
    bounds = (0, room["Lx"], 0, room["Ly"], 0, room["Lz"])
    fig = plt.figure(figsize=(6.4, 5.2), dpi=220)
    ax = fig.add_subplot(111, projection="3d")

    draw_room_3d(
        ax,
        scene,
        face_alpha=0.0 if wireframe else room_face_alpha,
        face_color="#dfe9f3",
        edge_color="#4d5661",
        linewidth=1.0 if wireframe else 1.2,
    )
    draw_obstacles_3d(
        ax,
        scene,
        face_color="#7e8791",
        face_alpha=0.0 if wireframe else obstacle_face_alpha,
        edge_color="#4b5158",
        linewidth=0.8,
    )
    draw_openings_3d(ax, scene)

    ax.view_init(elev=elev, azim=azim)
    set_axes_equal(ax, bounds, scale=0.55)
    style_3d_axes(ax)
    fig.tight_layout(pad=0.0)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_floorplan(scene: dict, output: Path):
    fig, ax = plt.subplots(figsize=(6.2, 5.0), dpi=220)
    draw_floorplan(ax, scene)
    fig.tight_layout(pad=0.0)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_section(scene: dict, output: Path):
    fig, ax = plt.subplots(figsize=(6.0, 4.8), dpi=220)
    draw_section(ax, scene)
    fig.tight_layout(pad=0.0)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_case(scene_path: Path, render_root: Path, views: list[str]) -> dict:
    scene = load_scene(scene_path)
    case_name = scene_path.stem if scene_path.stem.startswith("bench_") else f"bench_{scene_path.stem}"
    case_dir = render_root / case_name
    outputs: dict[str, str] = {}

    view_specs = {
        "perspective": lambda out: render_3d_scene(scene, out, elev=24, azim=-60, room_face_alpha=0.12, obstacle_face_alpha=0.98, wireframe=False),
        "wireframe": lambda out: render_3d_scene(scene, out, elev=24, azim=-60, room_face_alpha=0.0, obstacle_face_alpha=0.0, wireframe=True),
        "birdseye": lambda out: render_3d_scene(scene, out, elev=68, azim=-58, room_face_alpha=0.10, obstacle_face_alpha=0.95, wireframe=False),
        "floorplan": lambda out: render_floorplan(scene, out),
        "section": lambda out: render_section(scene, out),
    }

    for view in views:
        out = case_dir / view / f"{case_name}_{view}.png"
        view_specs[view](out)
        outputs[view] = str(out)

    section_axis, section_plane = infer_section_plane(scene)
    manifest = {
        "case_name": case_name,
        "source_scene": str(scene_path),
        "views": outputs,
        "room_kind": "composite" if "blocks" in scene.get("room", {}) else "rectangular",
        "obstacle_count": len(scene.get("obstacles", [])),
        "category": scene.get("meta", {}).get("benchmark", {}).get("category"),
        "section_view": {
            "axis": section_axis,
            "plane_coordinate": round(section_plane, 3),
        },
        "render_contract": {
            "layout": "benchmark/renderings/<case_name>/<view>/<case_name>_<view>.png",
            "note": "geometry-only benchmark input assets; distinct from reference CFD result visualizations",
        },
    }
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def collect_render_manifests(render_root: Path) -> list[dict]:
    manifests: list[dict] = []
    for manifest_path in sorted(render_root.glob("*/manifest.json")):
        manifests.append(json.loads(manifest_path.read_text(encoding="utf-8")))
    return manifests


def main() -> int:
    parser = argparse.ArgumentParser(description="Render benchmark 2D input views from scene JSON")
    parser.add_argument("--scenes", nargs="*", default=None, help="Scene JSON paths; default benchmark/scenes/*.json")
    parser.add_argument("--render-root", type=Path, default=DEFAULT_RENDER_ROOT)
    parser.add_argument("--views", nargs="+", default=["perspective", "birdseye", "floorplan", "wireframe", "section"], choices=["perspective", "birdseye", "floorplan", "wireframe", "section"])
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    scene_paths = [Path(p) for p in args.scenes] if args.scenes else sorted(DEFAULT_SCENES_DIR.glob("*.json"))
    if args.limit is not None:
        scene_paths = scene_paths[: args.limit]

    for scene_path in scene_paths:
        render_case(scene_path, args.render_root, args.views)

    manifests = collect_render_manifests(args.render_root)
    summary = {
        "ok": True,
        "count": len(manifests),
        "views": args.views,
        "render_root": str(args.render_root),
        "cases": manifests,
    }
    args.render_root.mkdir(parents=True, exist_ok=True)
    (args.render_root / "renderings_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
