#!/usr/bin/env python3
"""Visualize indoor scene geometry and OpenFOAM flow field side-by-side.

Outputs a 1x2 figure:
- left: geometry / obstacles / openings overview
- right: flow-field slice colored by velocity magnitude
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True


def load_scene(scene_json: Path) -> dict:
    with scene_json.open("r", encoding="utf-8") as f:
        return json.load(f)


def find_latest_vtk(case_dir: Path) -> Path:
    vtk_dir = case_dir / "VTK"
    if not vtk_dir.exists():
        raise FileNotFoundError(f"VTK directory not found: {vtk_dir}")

    candidates: list[tuple[int, Path]] = []
    for child in vtk_dir.iterdir():
        if child.is_dir():
            try:
                timestep = int(child.name.split("_")[-1])
            except ValueError:
                continue
            internal = child / "internal.vtu"
            if internal.exists():
                candidates.append((timestep, internal))

    if not candidates:
        raise FileNotFoundError(f"No internal.vtu found under {vtk_dir}")

    candidates.sort(key=lambda x: x[0])
    return candidates[-1][1]


def compute_velocity_point_data(mesh: pv.DataSet) -> pv.DataSet:
    if "U" in mesh.point_data:
        u = np.asarray(mesh.point_data["U"])
    elif "U" in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()
        u = np.asarray(mesh.point_data["U"])
    else:
        raise KeyError("Velocity field 'U' not found in mesh point_data or cell_data")

    mesh["Umag"] = np.linalg.norm(u, axis=1)
    return mesh


def room_blocks(scene: dict) -> list[dict]:
    room = scene["room"]
    if "blocks" in room:
        return room["blocks"]
    size = room["size"]
    return [{"origin": {"x": 0.0, "y": 0.0, "z": 0.0}, "size": {"dx": size["Lx"], "dy": size["Ly"], "dz": size["Lz"]}}]


def room_bounds(scene: dict) -> dict[str, float]:
    blocks = room_blocks(scene)
    return {
        "Lx": max(b["origin"]["x"] + b["size"]["dx"] for b in blocks),
        "Ly": max(b["origin"]["y"] + b["size"]["dy"] for b in blocks),
        "Lz": max(b["origin"]["z"] + b["size"]["dz"] for b in blocks),
    }


def draw_geometry_panel(ax, scene: dict) -> None:
    room = room_bounds(scene)
    Lx, Ly = room["Lx"], room["Ly"]

    for block in room_blocks(scene):
        ox = block["origin"]["x"]
        oy = block["origin"]["y"]
        dx = block["size"]["dx"]
        dy = block["size"]["dy"]
        ax.add_patch(Rectangle((ox, oy), dx, dy, fill=False, edgecolor="black", linewidth=2))

    for obs in scene["obstacles"]:
        x = obs["min"]["x"]
        y = obs["min"]["y"]
        dx = obs["size"]["dx"]
        dy = obs["size"]["dy"]
        ax.add_patch(Rectangle((x, y), dx, dy, facecolor="#9aa0a6", edgecolor="black", alpha=0.9))
        ax.text(x + dx / 2, y + dy / 2, obs["id"], ha="center", va="center", fontsize=7)

    for opening in scene["openings"]:
        wall = opening["wall"]
        u = opening["center"]["u"]
        du = opening["size"]["du"]
        color = "tab:blue" if opening["type"] == "inlet" else "tab:red"

        if wall == "west":
            ax.plot([0, 0], [u - du / 2, u + du / 2], color=color, linewidth=4)
        elif wall == "east":
            ax.plot([Lx, Lx], [u - du / 2, u + du / 2], color=color, linewidth=4)
        elif wall == "south":
            ax.plot([u - du / 2, u + du / 2], [0, 0], color=color, linewidth=4)
        elif wall == "north":
            ax.plot([u - du / 2, u + du / 2], [Ly, Ly], color=color, linewidth=4)

    ax.set_title("Gemini scene geometry (top view)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)
    ax.set_xlim(-0.05 * Lx, 1.05 * Lx)
    ax.set_ylim(-0.05 * Ly, 1.05 * Ly)


def draw_flow_panel(ax, vtk_path: Path, scene: dict) -> None:
    mesh = pv.read(vtk_path)
    mesh = compute_velocity_point_data(mesh)

    room = room_bounds(scene)
    center_z = 0.5 * room["Lz"]
    sl = mesh.slice(normal="z", origin=(0, 0, center_z))
    if sl.n_points == 0:
        raise RuntimeError("Mid-plane slice produced no points")

    pts = sl.points
    x = pts[:, 0]
    y = pts[:, 1]
    umag = np.asarray(sl["Umag"])

    tric = ax.tricontourf(x, y, umag, levels=30, cmap="turbo")
    ax.tricontour(x, y, umag, levels=10, colors="k", linewidths=0.3, alpha=0.4)

    for obs in scene["obstacles"]:
        ox = obs["min"]["x"]
        oy = obs["min"]["y"]
        odx = obs["size"]["dx"]
        ody = obs["size"]["dy"]
        ax.add_patch(Rectangle((ox, oy), odx, ody, fill=False, edgecolor="white", linewidth=1.0, linestyle="--"))

    cbar = plt.colorbar(tric, ax=ax)
    cbar.set_label("|U| [m/s]")
    ax.set_title(f"OpenFOAM velocity magnitude slice (z={center_z:.2f} m)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.15)
    ax.set_xlim(0, room["Lx"])
    ax.set_ylim(0, room["Ly"])


def main() -> int:
    parser = argparse.ArgumentParser(description="Create side-by-side indoor geometry/flow visualization")
    parser.add_argument("scene_json", type=Path)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("results/indoor_case_comparison.png"))
    parser.add_argument("--vtk", type=Path, default=None, help="Optional explicit VTK internal.vtu path")
    args = parser.parse_args()

    scene = load_scene(args.scene_json)
    vtk_path = args.vtk if args.vtk else find_latest_vtk(args.case_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    draw_geometry_panel(axes[0], scene)
    draw_flow_panel(axes[1], vtk_path, scene)
    fig.suptitle("Indoor CFD pipeline output", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
