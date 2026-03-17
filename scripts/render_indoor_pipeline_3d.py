#!/usr/bin/env python3
"""Render a 3D two-panel figure for the Gemini -> OpenFOAM indoor pipeline.

Panel (a): Gemini-generated 3D room geometry (cutaway)
Panel (b): OpenFOAM flow field in same room with multiple contour planes + streamlines
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True


def load_scene(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def load_flow(vtk_path: Path) -> pv.DataSet:
    mesh = pv.read(vtk_path)
    if "U" not in mesh.point_data and "U" in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()
    if "U" not in mesh.point_data:
        raise KeyError("Velocity field U not found")
    u = np.asarray(mesh.point_data["U"])
    mesh.point_data["Umag"] = np.linalg.norm(u, axis=1)
    return mesh


def find_latest_vtk(case_dir: Path) -> Path:
    vtk_dir = case_dir / "VTK"
    candidates: list[tuple[int, Path]] = []
    for child in vtk_dir.iterdir():
        if child.is_dir():
            try:
                t = int(child.name.split("_")[-1])
            except ValueError:
                continue
            vtu = child / "internal.vtu"
            if vtu.exists():
                candidates.append((t, vtu))
    if not candidates:
        raise FileNotFoundError(f"No internal.vtu under {vtk_dir}")
    candidates.sort(key=lambda x: x[0])
    return candidates[-1][1]


def set_axes_equal(ax, bounds, scale=0.58):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    xmid, ymid, zmid = (xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2
    radius = scale * max(xmax - xmin, ymax - ymin, zmax - zmin)
    ax.set_xlim(xmid - radius, xmid + radius)
    ax.set_ylim(ymid - radius, ymid + radius)
    ax.set_zlim(zmid - radius, zmid + radius)
    ax.set_box_aspect((xmax - xmin, ymax - ymin, zmax - zmin))


def style_axes(ax):
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.line.set_color((1, 1, 1, 0))


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


def draw_cutaway_room(ax, room, face_alpha=0.08, edge_color="0.62", lw=1.25):
    Lx, Ly, Lz = room["Lx"], room["Ly"], room["Lz"]
    # keep floor + back + left + right walls; omit front wall and ceiling for cutaway
    faces = [
        [[0, 0, 0], [Lx, 0, 0], [Lx, Ly, 0], [0, Ly, 0]],  # floor
        [[0, Ly, 0], [Lx, Ly, 0], [Lx, Ly, Lz], [0, Ly, Lz]],  # back wall
        [[0, 0, 0], [0, Ly, 0], [0, Ly, Lz], [0, 0, Lz]],  # west wall
        [[Lx, 0, 0], [Lx, Ly, 0], [Lx, Ly, Lz], [Lx, 0, Lz]],  # east wall
    ]
    coll = Poly3DCollection(faces, facecolor=(0.84, 0.90, 0.97, face_alpha), edgecolor="none")
    ax.add_collection3d(coll)

    # wire edges for cutaway frame
    pts = np.array([
        [0,0,0],[Lx,0,0],[Lx,Ly,0],[0,Ly,0],
        [0,0,Lz],[Lx,0,Lz],[Lx,Ly,Lz],[0,Ly,Lz]
    ])
    edges = [(0,1),(1,2),(2,3),(3,0),(3,7),(7,6),(6,2),(1,5),(5,6),(0,4),(4,7)]
    for i, j in edges:
        ax.plot(*zip(pts[i], pts[j]), color=edge_color, lw=lw, alpha=0.95)


def draw_obstacles(ax, obstacles, alpha=1.0, color="#7a7a7a"):
    for obs in obstacles:
        x0, y0, z0 = obs["min"]["x"], obs["min"]["y"], obs["min"]["z"]
        dx, dy, dz = obs["size"]["dx"], obs["size"]["dy"], obs["size"]["dz"]
        coll = Poly3DCollection(cuboid_faces(x0, y0, z0, dx, dy, dz), facecolor=color, edgecolor="0.35", linewidths=0.35, alpha=alpha)
        ax.add_collection3d(coll)


def draw_openings(ax, openings, room):
    Lx, Ly, Lz = room["Lx"], room["Ly"], room["Lz"]
    for op in openings:
        color = "#377eb8" if op["type"] == "inlet" else "#e41a1c"
        u, v = op["center"]["u"], op["center"]["v"]
        du, dv = op["size"]["du"], op["size"]["dv"]
        if op["wall"] == "west":
            face = [[[0, u-du/2, v-dv/2], [0, u+du/2, v-dv/2], [0, u+du/2, v+dv/2], [0, u-du/2, v+dv/2]]]
            ax.add_collection3d(Poly3DCollection(face, facecolor=color, edgecolor="none", alpha=0.95))
            ax.quiver(-0.7, u, v, 0.5, 0, 0, color=color, linewidth=2.4, arrow_length_ratio=0.25)
        elif op["wall"] == "east":
            face = [[[Lx, u-du/2, v-dv/2], [Lx, u+du/2, v-dv/2], [Lx, u+du/2, v+dv/2], [Lx, u-du/2, v+dv/2]]]
            ax.add_collection3d(Poly3DCollection(face, facecolor=color, edgecolor="none", alpha=0.95))
            ax.quiver(Lx+0.15, u, v, 0.5, 0, 0, color=color, linewidth=2.4, arrow_length_ratio=0.25)


def slice_to_tris(poly: pv.PolyData):
    if poly.n_cells == 0:
        return np.empty((0, 3, 3)), np.empty((0, 3), dtype=int)
    faces = poly.faces.reshape(-1, 4)[:, 1:4]
    return np.asarray(poly.points)[faces], faces


def lines_to_segments(poly: pv.PolyData):
    segs = []
    if poly.n_cells == 0 or poly.lines.size == 0:
        return segs
    pts = np.asarray(poly.points)
    arr = poly.lines
    i = 0
    while i < len(arr):
        n = int(arr[i])
        ids = arr[i + 1:i + 1 + n]
        xyz = pts[ids]
        if len(xyz) >= 2:
            segs.append(xyz)
        i += n + 1
    return segs


def make_inlet_seed(openings, nz=5, ny=5):
    inlet = next(op for op in openings if op["type"] == "inlet")
    u, v = inlet["center"]["u"], inlet["center"]["v"]
    du, dv = inlet["size"]["du"], inlet["size"]["dv"]
    ys = np.linspace(u - 0.40 * du, u + 0.40 * du, ny)
    zs = np.linspace(v - 0.40 * dv, v + 0.40 * dv, nz)
    yy, zz = np.meshgrid(ys, zs, indexing="xy")
    xx = np.full_like(yy, 0.08)
    return pv.PolyData(np.c_[xx.ravel(), yy.ravel(), zz.ravel()])


def render(scene: dict, mesh: pv.DataSet, output: Path, panel_a: Path | None = None, panel_b: Path | None = None):
    room = scene["room"]["size"]
    bounds = (0, room["Lx"], 0, room["Ly"], 0, room["Lz"])

    fig = plt.figure(figsize=(14.2, 5.7), dpi=220)
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")

    # panel a
    draw_cutaway_room(ax1, room, face_alpha=0.12, edge_color="0.58", lw=1.35)
    draw_obstacles(ax1, scene["obstacles"], alpha=1.0, color="#767676")
    draw_openings(ax1, scene["openings"], room)
    ax1.set_title("(a) Gemini-generated 3D scene", pad=10)

    # panel b
    draw_cutaway_room(ax2, room, face_alpha=0.06, edge_color="0.60", lw=1.25)
    draw_obstacles(ax2, scene["obstacles"], alpha=0.34, color="#7b7b7b")

    cmap = colormaps["turbo"]
    all_plane_vals = []
    plane_payloads = []

    # Plane 1: vertical longitudinal slice (x-z at y=center)
    ymid = 0.5 * room["Ly"]
    sl_y = mesh.slice(normal=(0, 1, 0), origin=(0, ymid, 0)).triangulate()
    tris_y, faces_y = slice_to_tris(sl_y)
    if len(faces_y):
        vals_y = np.asarray(sl_y["Umag"])[faces_y].mean(axis=1)
        all_plane_vals.extend(vals_y.tolist())
        plane_payloads.append((tris_y, vals_y, 0.72))

    # Plane 2: horizontal elevated slice to show room interior more clearly
    zcut = 0.55 * room["Lz"]
    sl_z = mesh.slice(normal=(0, 0, 1), origin=(0, 0, zcut)).triangulate()
    tris_z, faces_z = slice_to_tris(sl_z)
    if len(faces_z):
        vals_z = np.asarray(sl_z["Umag"])[faces_z].mean(axis=1)
        all_plane_vals.extend(vals_z.tolist())
        plane_payloads.append((tris_z, vals_z, 0.45))

    if all_plane_vals:
        norm = Normalize(vmin=min(all_plane_vals), vmax=max(all_plane_vals))
        for tris_i, vals_i, alpha_i in plane_payloads:
            colors_i = cmap(norm(vals_i))
            colors_i[:, 3] = alpha_i
            ax2.add_collection3d(Poly3DCollection(tris_i, facecolors=colors_i, edgecolor="none"))
    else:
        norm = Normalize(vmin=0.0, vmax=1.0)

    seeds = make_inlet_seed(scene["openings"], nz=5, ny=5)
    stream = mesh.streamlines_from_source(
        seeds,
        vectors="U",
        integration_direction="forward",
        initial_step_length=0.065,
        max_length=16.0,
        terminal_speed=1e-6,
        compute_vorticity=False,
    )
    segments = lines_to_segments(stream)
    if segments:
        metric_all = []
        metrics = []
        for seg in segments:
            m = np.linalg.norm(np.diff(seg, axis=0), axis=1)
            if len(m):
                metrics.append(m)
                metric_all.extend(m.tolist())
        if metric_all:
            line_norm = Normalize(vmin=min(metric_all), vmax=max(metric_all))
            for seg, m in zip(segments, metrics):
                for i in range(len(seg) - 1):
                    xyz = seg[i:i+2]
                    ax2.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=cmap(line_norm(m[i])), lw=1.0, alpha=0.82)

    draw_openings(ax2, scene["openings"], room)
    ax2.set_title("(b) OpenFOAM flow field with contour planes", pad=10)

    ax1.view_init(elev=24, azim=-62)
    ax2.view_init(elev=28, azim=-58)
    set_axes_equal(ax1, bounds, scale=0.52)
    set_axes_equal(ax2, bounds, scale=0.56)
    for ax in (ax1, ax2):
        style_axes(ax)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax2], fraction=0.022, pad=0.018)
    cbar.set_label("Velocity magnitude |U| [m/s]")

    # Framework identity / process annotations in figure coordinates
    fig.text(
        0.125, 0.90,
        'Prompt → Gemini\nStructured scene JSON',
        ha='left', va='top', fontsize=10.5,
        bbox=dict(boxstyle='round,pad=0.35', facecolor='#f5f7fb', edgecolor='#6e98ff', linewidth=1.15)
    )
    fig.text(
        0.515, 0.52,
        'Automated meshing\n+ OpenFOAM CFD',
        ha='center', va='center', fontsize=10.5, color='#4f4f4f'
    )
    arrow = FancyArrowPatch(
        (0.43, 0.52), (0.60, 0.52),
        transform=fig.transFigure,
        arrowstyle='simple', mutation_scale=16,
        fc='#c7cdd8', ec='#c7cdd8', alpha=0.85
    )
    fig.add_artist(arrow)

    fig.subplots_adjust(left=0.03, right=0.93, top=0.92, bottom=0.06, wspace=0.02)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    if output.suffix.lower() != ".pdf":
        fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # optional single-panel exports
    if panel_a or panel_b:
        for which, path in [("a", panel_a), ("b", panel_b)]:
            if path is None:
                continue
            f = plt.figure(figsize=(6.4, 5.6), dpi=220)
            ax = f.add_subplot(111, projection="3d")
            if which == "a":
                draw_cutaway_room(ax, room, face_alpha=0.12, edge_color="0.58", lw=1.35)
                draw_obstacles(ax, scene["obstacles"], alpha=1.0, color="#767676")
                draw_openings(ax, scene["openings"], room)
                ax.set_title("(a) Gemini-generated 3D scene", pad=10)
            else:
                draw_cutaway_room(ax, room, face_alpha=0.06, edge_color="0.60", lw=1.25)
                draw_obstacles(ax, scene["obstacles"], alpha=0.34, color="#7b7b7b")
                if all_plane_vals:
                    for tris_i, vals_i, alpha_i in plane_payloads:
                        colors_i = cmap(norm(vals_i))
                        colors_i[:, 3] = alpha_i
                        ax.add_collection3d(Poly3DCollection(tris_i, facecolors=colors_i, edgecolor="none"))
                if segments and metric_all:
                    for seg, m in zip(segments, metrics):
                        for i in range(len(seg) - 1):
                            xyz = seg[i:i+2]
                            ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=cmap(line_norm(m[i])), lw=1.0, alpha=0.82)
                draw_openings(ax, scene["openings"], room)
                ax.set_title("(b) OpenFOAM flow field with contour planes", pad=10)
            if which == "a":
                ax.view_init(elev=24, azim=-62)
                set_axes_equal(ax, bounds, scale=0.52)
            else:
                ax.view_init(elev=28, azim=-58)
                set_axes_equal(ax, bounds, scale=0.56)
            style_axes(ax)
            f.tight_layout()
            path.parent.mkdir(parents=True, exist_ok=True)
            f.savefig(path, bbox_inches="tight", facecolor="white")
            plt.close(f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render 3D indoor pipeline figure")
    parser.add_argument("--scene-json", type=Path, default=Path("/home/yhjoo/projects/OpenFOAM/generated/indoor_pipeline_g31.json"))
    parser.add_argument("--case-dir", type=Path, default=Path("/home/yhjoo/projects/OpenFOAM/cases/indoor_pipeline_g31"))
    parser.add_argument("--vtk", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=Path("/home/yhjoo/projects/OpenFOAM/results/indoor_pipeline_g31/indoor_pipeline_3d_comparison.png"))
    args = parser.parse_args()

    scene = load_scene(args.scene_json)
    vtk_path = args.vtk if args.vtk else find_latest_vtk(args.case_dir)
    mesh = load_flow(vtk_path)
    render(
        scene,
        mesh,
        args.output,
        panel_a=args.output.with_name("panel_geometry_3d.png"),
        panel_b=args.output.with_name("panel_flow_3d.png"),
    )
    print(args.output)
