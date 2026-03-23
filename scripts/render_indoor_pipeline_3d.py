#!/usr/bin/env python3
"""Render a 3D two-panel figure for the Gemini -> OpenFOAM indoor pipeline.

Panel (a): Gemini-generated 3D room geometry (walls, obstacles, openings)
Panel (b): OpenFOAM flow field with contour planes + streamlines

Uses PyVista off-screen rendering for publication-quality output.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

# Force software rendering for headless environments (WSL2, CI, etc.)
os.environ.setdefault("LIBGL_ALWAYS_SOFTWARE", "1")
os.environ.setdefault("MESA_GL_VERSION_OVERRIDE", "3.3")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np

import pyvista as pv

pv.OFF_SCREEN = True
pv.global_theme.anti_aliasing = "msaa"


# ---------------------------------------------------------------------------
# Scene loading helpers
# ---------------------------------------------------------------------------

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


def room_blocks(scene: dict) -> list[dict]:
    room = scene["room"]
    if "blocks" in room:
        return room["blocks"]
    size = room["size"]
    return [{"origin": {"x": 0.0, "y": 0.0, "z": 0.0},
             "size": {"dx": size["Lx"], "dy": size["Ly"], "dz": size["Lz"]}}]


def room_bounds(scene: dict) -> dict[str, float]:
    blocks = room_blocks(scene)
    return {
        "Lx": max(b["origin"]["x"] + b["size"]["dx"] for b in blocks),
        "Ly": max(b["origin"]["y"] + b["size"]["dy"] for b in blocks),
        "Lz": max(b["origin"]["z"] + b["size"]["dz"] for b in blocks),
    }


# ---------------------------------------------------------------------------
# PyVista geometry builders
# ---------------------------------------------------------------------------

def make_box(x0: float, y0: float, z0: float,
             dx: float, dy: float, dz: float) -> pv.PolyData:
    """Create a box as a PolyData surface mesh."""
    return pv.Box(bounds=(x0, x0 + dx, y0, y0 + dy, z0, z0 + dz))


def add_room_walls(plotter: pv.Plotter, scene: dict, *,
                   opacity: float = 0.15, color: str = "#d6e4f0",
                   edge_color: str = "#606060", edge_width: float = 2.0):
    """Add room walls as semi-transparent surfaces with solid edges."""
    blocks = room_blocks(scene)
    for block in blocks:
        x0 = block["origin"]["x"]
        y0 = block["origin"]["y"]
        z0 = block["origin"]["z"]
        dx, dy, dz = block["size"]["dx"], block["size"]["dy"], block["size"]["dz"]
        box = make_box(x0, y0, z0, dx, dy, dz)
        plotter.add_mesh(box, color=color, opacity=opacity,
                         show_edges=True, edge_color=edge_color,
                         line_width=edge_width)


def add_obstacles(plotter: pv.Plotter, obstacles: list[dict], *,
                  opacity: float = 1.0, color: str = "#7a7a7a",
                  edge_color: str = "#353535"):
    """Add obstacles as solid gray boxes."""
    for obs in obstacles:
        x0, y0, z0 = obs["min"]["x"], obs["min"]["y"], obs["min"]["z"]
        dx, dy, dz = obs["size"]["dx"], obs["size"]["dy"], obs["size"]["dz"]
        box = make_box(x0, y0, z0, dx, dy, dz)
        plotter.add_mesh(box, color=color, opacity=opacity,
                         show_edges=True, edge_color=edge_color,
                         line_width=1.5)


def _opening_patch(op: dict, room: dict) -> tuple[pv.PolyData, np.ndarray]:
    """Return a planar quad for the opening and the outward normal direction."""
    u, v = op["center"]["u"], op["center"]["v"]
    du, dv = op["size"]["du"], op["size"]["dv"]
    Lx, Ly, Lz = room["Lx"], room["Ly"], room["Lz"]
    wall = op["wall"]

    if wall == "west":
        pts = np.array([
            [0, u - du/2, v - dv/2],
            [0, u + du/2, v - dv/2],
            [0, u + du/2, v + dv/2],
            [0, u - du/2, v + dv/2],
        ])
        normal = np.array([-1, 0, 0], dtype=float)
    elif wall == "east":
        pts = np.array([
            [Lx, u - du/2, v - dv/2],
            [Lx, u + du/2, v - dv/2],
            [Lx, u + du/2, v + dv/2],
            [Lx, u - du/2, v + dv/2],
        ])
        normal = np.array([1, 0, 0], dtype=float)
    elif wall == "south":
        pts = np.array([
            [u - du/2, 0, v - dv/2],
            [u + du/2, 0, v - dv/2],
            [u + du/2, 0, v + dv/2],
            [u - du/2, 0, v + dv/2],
        ])
        normal = np.array([0, -1, 0], dtype=float)
    elif wall == "north":
        pts = np.array([
            [u - du/2, Ly, v - dv/2],
            [u + du/2, Ly, v - dv/2],
            [u + du/2, Ly, v + dv/2],
            [u - du/2, Ly, v + dv/2],
        ])
        normal = np.array([0, 1, 0], dtype=float)
    elif wall == "floor":
        pts = np.array([
            [u - du/2, v - dv/2, 0],
            [u + du/2, v - dv/2, 0],
            [u + du/2, v + dv/2, 0],
            [u - du/2, v + dv/2, 0],
        ])
        normal = np.array([0, 0, -1], dtype=float)
    elif wall == "ceiling":
        pts = np.array([
            [u - du/2, v - dv/2, Lz],
            [u + du/2, v - dv/2, Lz],
            [u + du/2, v + dv/2, Lz],
            [u - du/2, v + dv/2, Lz],
        ])
        normal = np.array([0, 0, 1], dtype=float)
    else:
        raise ValueError(f"Unknown wall: {wall}")

    faces = np.array([4, 0, 1, 2, 3])
    quad = pv.PolyData(pts, faces=faces)
    return quad, normal


def add_openings(plotter: pv.Plotter, openings: list[dict], room: dict, *,
                 arrow_scale: float = 0.8):
    """Add opening patches (blue=inlet, red=outlet) with flow-direction arrows."""
    for op in openings:
        color = "#377eb8" if op["type"] == "inlet" else "#e41a1c"
        quad, normal = _opening_patch(op, room)
        plotter.add_mesh(quad, color=color, opacity=0.92,
                         show_edges=True, edge_color=color, line_width=2.0)

        # Arrow showing flow direction
        center = quad.center
        # For inlet: arrow points inward (negate normal); for outlet: outward
        if op["type"] == "inlet":
            direction = -normal
            start = np.array(center) + normal * arrow_scale * 0.8
        else:
            direction = normal
            start = np.array(center)

        arrow = pv.Arrow(start=start, direction=direction,
                         scale=arrow_scale, tip_length=0.3,
                         tip_radius=0.10, shaft_radius=0.04)
        plotter.add_mesh(arrow, color=color, opacity=0.9)


def make_inlet_seed(openings: list[dict], room: dict,
                    nz: int = 5, ny: int = 5) -> pv.PolyData:
    """Create a grid of seed points just inside the inlet opening."""
    inlet = next(op for op in openings if op["type"] == "inlet")
    u, v = inlet["center"]["u"], inlet["center"]["v"]
    du, dv = inlet["size"]["du"], inlet["size"]["dv"]
    wall = inlet["wall"]
    Lx, Ly = room["Lx"], room["Ly"]

    us = np.linspace(u - 0.40 * du, u + 0.40 * du, ny)
    vs = np.linspace(v - 0.40 * dv, v + 0.40 * dv, nz)
    uu, vv = np.meshgrid(us, vs, indexing="xy")

    offset = 0.08  # slightly inside the room
    if wall == "west":
        xx = np.full_like(uu, offset)
        yy, zz = uu, vv
    elif wall == "east":
        xx = np.full_like(uu, Lx - offset)
        yy, zz = uu, vv
    elif wall == "south":
        yy = np.full_like(uu, offset)
        xx, zz = uu, vv
    elif wall == "north":
        yy = np.full_like(uu, Ly - offset)
        xx, zz = uu, vv
    else:
        # floor/ceiling: fall back to horizontal
        xx, yy = uu, vv
        zz = np.full_like(uu, offset)

    return pv.PolyData(np.c_[xx.ravel(), yy.ravel(), zz.ravel()])


# ---------------------------------------------------------------------------
# Camera helper
# ---------------------------------------------------------------------------

def setup_camera(plotter: pv.Plotter, room: dict, *,
                 elev: float = 30.0, azim: float = -55.0,
                 zoom: float = 1.0):
    """Set up a consistent isometric-like camera view."""
    Lx, Ly, Lz = room["Lx"], room["Ly"], room["Lz"]
    cx, cy, cz = Lx / 2, Ly / 2, Lz / 2
    plotter.camera.focal_point = (cx, cy, cz)

    r = max(Lx, Ly, Lz) * 2.2
    elev_rad = np.radians(elev)
    azim_rad = np.radians(azim)
    eye_x = cx + r * np.cos(elev_rad) * np.cos(azim_rad)
    eye_y = cy + r * np.cos(elev_rad) * np.sin(azim_rad)
    eye_z = cz + r * np.sin(elev_rad)
    plotter.camera.position = (eye_x, eye_y, eye_z)
    plotter.camera.up = (0, 0, 1)
    plotter.camera.zoom(zoom)


# ---------------------------------------------------------------------------
# Panel rendering
# ---------------------------------------------------------------------------

def render_panel_a(scene: dict, room: dict, img_path: Path,
                   width: int = 1800, height: int = 1500):
    """Render Panel (a): Gemini-generated 3D scene."""
    pl = pv.Plotter(off_screen=True, window_size=(width, height))
    pl.set_background("white")

    add_room_walls(pl, scene, opacity=0.15, color="#d6e4f0",
                   edge_color="#505050", edge_width=2.5)
    add_obstacles(pl, scene["obstacles"], opacity=1.0, color="#7a7a7a",
                  edge_color="#353535")
    add_openings(pl, scene["openings"], room, arrow_scale=0.7)

    setup_camera(pl, room, elev=28, azim=-58, zoom=1.05)
    pl.enable_anti_aliasing("msaa")

    img_path.parent.mkdir(parents=True, exist_ok=True)
    pl.screenshot(str(img_path), transparent_background=False, scale=2)
    pl.close()
    print(f"  Panel (a) saved: {img_path}")


def render_panel_b(scene: dict, room: dict, mesh: pv.DataSet,
                   img_path: Path, width: int = 1800, height: int = 1500):
    """Render Panel (b): OpenFOAM flow field with contour planes + streamlines."""
    pl = pv.Plotter(off_screen=True, window_size=(width, height))
    pl.set_background("white")

    # Room walls - more transparent to see internal flow
    add_room_walls(pl, scene, opacity=0.08, color="#dce8f4",
                   edge_color="#707070", edge_width=1.5)
    # Obstacles - semi-transparent
    add_obstacles(pl, scene["obstacles"], opacity=0.3, color="#8a8a8a",
                  edge_color="#555555")

    Lx, Ly, Lz = room["Lx"], room["Ly"], room["Lz"]

    # Compute global velocity range for consistent colormap
    umag = np.asarray(mesh.point_data["Umag"])
    vmin, vmax = float(umag.min()), float(umag.max())
    clim = [vmin, vmax]

    # Contour plane 1: horizontal at z = 0.5*Lz
    hz_origin = (Lx / 2, Ly / 2, 0.50 * Lz)
    sl_hz = mesh.slice(normal=(0, 0, 1), origin=hz_origin)
    if sl_hz.n_cells > 0:
        pl.add_mesh(sl_hz, scalars="Umag", cmap="turbo", clim=clim,
                    opacity=0.70, show_edges=False,
                    scalar_bar_args={"title": "|U| [m/s]", "n_labels": 3,
                                     "position_x": 0.83, "position_y": 0.22,
                                     "width": 0.10, "height": 0.40,
                                     "title_font_size": 34, "label_font_size": 26,
                                     "color": "black", "bold": True,
                                     "vertical": True,
                                     "fmt": "%.3f"})
    else:
        print("  Warning: horizontal slice has 0 cells")

    # Contour plane 2: vertical at y = 0.5*Ly
    vt_origin = (Lx / 2, 0.50 * Ly, Lz / 2)
    sl_vt = mesh.slice(normal=(0, 1, 0), origin=vt_origin)
    if sl_vt.n_cells > 0:
        pl.add_mesh(sl_vt, scalars="Umag", cmap="turbo", clim=clim,
                    opacity=0.65, show_edges=False, show_scalar_bar=False)
    else:
        print("  Warning: vertical slice has 0 cells")

    # Streamlines from inlet
    seeds = make_inlet_seed(scene["openings"], room, nz=6, ny=6)
    try:
        stream = mesh.streamlines_from_source(
            seeds,
            vectors="U",
            integration_direction="forward",
            initial_step_length=0.05,
            max_length=20.0,
            terminal_speed=1e-6,
            compute_vorticity=False,
        )
        if stream.n_cells > 0:
            # Color streamlines by velocity magnitude
            if "Umag" not in stream.point_data and "U" in stream.point_data:
                u_st = np.asarray(stream.point_data["U"])
                stream.point_data["Umag"] = np.linalg.norm(u_st, axis=1)
            pl.add_mesh(stream, scalars="Umag", cmap="turbo", clim=clim,
                        line_width=2.5, opacity=0.9, show_scalar_bar=False)
        else:
            print("  Warning: streamlines produced 0 cells")
    except Exception as e:
        print(f"  Warning: streamline generation failed: {e}")

    # Openings with arrows
    add_openings(pl, scene["openings"], room, arrow_scale=0.6)

    setup_camera(pl, room, elev=28, azim=-58, zoom=1.05)
    pl.enable_anti_aliasing("msaa")

    img_path.parent.mkdir(parents=True, exist_ok=True)
    pl.screenshot(str(img_path), transparent_background=False, scale=2)
    pl.close()
    print(f"  Panel (b) saved: {img_path}")


# ---------------------------------------------------------------------------
# Combined figure assembly
# ---------------------------------------------------------------------------

def combine_panels(panel_a_path: Path, panel_b_path: Path,
                   output_path: Path, dpi: int = 600):
    """Combine two panel images into a side-by-side publication figure."""
    img_a = plt.imread(str(panel_a_path))
    img_b = plt.imread(str(panel_b_path))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.2, 5.7))

    ax1.imshow(img_a)
    ax1.set_title("(a) Gemini-generated 3D scene", fontsize=13, pad=10)
    ax1.axis("off")

    ax2.imshow(img_b)
    ax2.set_title("(b) OpenFOAM flow field with contour planes", fontsize=13, pad=10)
    ax2.axis("off")

    # Process annotation between panels
    fig.text(
        0.125, 0.91,
        "Prompt \u2192 Gemini\nStructured scene JSON",
        ha="left", va="top", fontsize=10.5,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#f5f7fb",
                  edgecolor="#6e98ff", linewidth=1.15),
    )
    fig.text(
        0.515, 0.52,
        "Automated meshing\n+ OpenFOAM CFD",
        ha="center", va="center", fontsize=10.5, color="#4f4f4f",
    )
    arrow = FancyArrowPatch(
        (0.43, 0.52), (0.60, 0.52),
        transform=fig.transFigure,
        arrowstyle="simple", mutation_scale=16,
        fc="#c7cdd8", ec="#c7cdd8", alpha=0.85,
    )
    fig.add_artist(arrow)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.03, wspace=0.04)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save PNG
    png_path = output_path.with_suffix(".png")
    fig.savefig(str(png_path), dpi=dpi, bbox_inches="tight", facecolor="white")
    print(f"  Combined PNG: {png_path}")

    # Save PDF
    pdf_path = output_path.with_suffix(".pdf")
    fig.savefig(str(pdf_path), bbox_inches="tight", facecolor="white")
    print(f"  Combined PDF: {pdf_path}")

    plt.close(fig)


# ---------------------------------------------------------------------------
# Main render orchestrator
# ---------------------------------------------------------------------------

def render(scene: dict, mesh: pv.DataSet, output: Path,
           panel_a: Path | None = None, panel_b: Path | None = None):
    """Full rendering pipeline: individual panels then combined figure."""
    room = room_bounds(scene)

    # Determine panel output paths
    pa = panel_a or output.with_name("panel_geometry_3d.png")
    pb = panel_b or output.with_name("panel_flow_3d.png")

    print("Rendering Panel (a): 3D geometry ...")
    render_panel_a(scene, room, pa)

    print("Rendering Panel (b): flow field ...")
    render_panel_b(scene, room, mesh, pb)

    print("Combining panels ...")
    combine_panels(pa, pb, output)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Render 3D indoor pipeline figure (PyVista)")
    parser.add_argument(
        "--scene-json", "--scene", type=Path,
        default=Path("/home/yhjoo/projects/OpenFOAM/generated/"
                     "indoor_pipeline_g31.json"),
        help="Path to Gemini-generated scene JSON")
    parser.add_argument(
        "--case-dir", "--case", type=Path,
        default=Path("/home/yhjoo/projects/OpenFOAM/cases/"
                     "indoor_pipeline_g31"),
        help="Path to OpenFOAM case directory")
    parser.add_argument(
        "--vtk", type=Path, default=None,
        help="Direct path to VTK file (overrides auto-detection)")
    parser.add_argument(
        "--output", type=Path,
        default=Path("/home/yhjoo/projects/OpenFOAM/results/"
                     "indoor_pipeline_g31/indoor_pipeline_3d_comparison.png"),
        help="Output path for combined figure")
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
    print(f"Done: {args.output}")
