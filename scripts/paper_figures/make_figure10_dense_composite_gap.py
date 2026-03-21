#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import matplotlib.patheffects as pe

pv.OFF_SCREEN = True

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results/paper_figures"
PDF_OUT = OUT_DIR / "figure10_dense_composite_structure_vs_cfd_gap.pdf"
PNG_OUT = OUT_DIR / "figure10_dense_composite_structure_vs_cfd_gap.png"
PDF_RENDER_OUT = OUT_DIR / "figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png"

SELECTED = [
    {"case": "bench_a4_02", "view": "floorplan", "row_label": "(a) A4-02"},
    {"case": "bench_a4_04", "view": "floorplan", "row_label": "(b) A4-04"},
]

COLORS = {
    "title": "#111827",
    "subtitle": "#334155",
    "room_fill": "#EEF3F8",
    "room_edge": "#1F2937",
    "obstacle_fill": "#E5A657",
    "obstacle_edge": "#8A4B08",
    "inlet": "#2563EB",
    "outlet": "#C83E3A",
    "frame": "#94A3B8",
    "metric_box": "#F8FAFC",
    "metric_edge": "#CBD5E1",
}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def room_blocks(scene: dict[str, Any]) -> list[dict[str, float]]:
    room = scene["room"]
    if "blocks" in room:
        return [
            {
                "x": float(block["origin"]["x"]),
                "y": float(block["origin"]["y"]),
                "z": float(block["origin"].get("z", 0.0)),
                "dx": float(block["size"]["dx"]),
                "dy": float(block["size"]["dy"]),
                "dz": float(block["size"].get("dz", scene["room"].get("size", {}).get("Lz", 0.0))),
            }
            for block in room["blocks"]
        ]
    size = room["size"]
    return [{"x": 0.0, "y": 0.0, "z": 0.0, "dx": float(size["Lx"]), "dy": float(size["Ly"]), "dz": float(size["Lz"])}]


def room_extent(blocks: list[dict[str, float]]) -> tuple[float, float]:
    max_x = max(block["x"] + block["dx"] for block in blocks)
    max_y = max(block["y"] + block["dy"] for block in blocks)
    return max_x, max_y


def room_extent_3d(blocks: list[dict[str, float]]) -> tuple[float, float, float]:
    max_x = max(block["x"] + block["dx"] for block in blocks)
    max_y = max(block["y"] + block["dy"] for block in blocks)
    max_z = max(float(block.get("z", 0.0)) + float(block.get("dz", 0.0)) for block in blocks)
    return max_x, max_y, max_z


def prepare_opening_geometry(opening: dict[str, Any], x_max: float, y_max: float) -> tuple[str, dict[str, float], dict[str, float]]:
    wall = opening["wall"]
    center = {
        "u": float(opening["center"]["u"]),
        "wall_x_max": x_max,
        "wall_y_max": y_max,
    }
    size = {"du": float(opening["size"]["du"]), "dv": float(opening["size"]["dv"])}
    return wall, center, size


def draw_opening(ax: plt.Axes, wall: str, center: dict[str, float], size: dict[str, float], color: str, lw: float = 3.1) -> None:
    if wall in {"north", "south"}:
        x0 = center["u"] - size["du"] / 2.0
        x1 = center["u"] + size["du"] / 2.0
        y = 0.0 if wall == "south" else center["wall_y_max"]
        ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle="round", zorder=8)
    elif wall in {"west", "east"}:
        y0 = center["u"] - size["du"] / 2.0
        y1 = center["u"] + size["du"] / 2.0
        x = 0.0 if wall == "west" else center["wall_x_max"]
        ax.plot([x, x], [y0, y1], color=color, lw=lw, solid_capstyle="round", zorder=8)


def find_latest_vtk(case_dir: Path) -> Path:
    vtk_dir = case_dir / "VTK"
    candidates: list[tuple[int, Path]] = []
    for child in vtk_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            timestep = int(child.name.split("_")[-1])
        except ValueError:
            continue
        internal = child / "internal.vtu"
        if internal.exists():
            candidates.append((timestep, internal))
    if not candidates:
        raise FileNotFoundError(f"No internal.vtu found under {vtk_dir}")
    candidates.sort(key=lambda item: item[0])
    return candidates[-1][1]


def compute_velocity_point_data(mesh: pv.DataSet) -> pv.DataSet:
    if "U" in mesh.point_data:
        u = np.asarray(mesh.point_data["U"])
    elif "U" in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()
        u = np.asarray(mesh.point_data["U"])
    else:
        raise KeyError("Velocity field 'U' not found in mesh")
    mesh["Umag"] = np.linalg.norm(u, axis=1)
    return mesh


def extract_midplane_slice(scene: dict[str, Any], case_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vtk_path = find_latest_vtk(case_dir)
    mesh = compute_velocity_point_data(pv.read(vtk_path))
    blocks = room_blocks(scene)
    _, _, z_max = room_extent_3d(blocks)
    sl = mesh.slice(normal="z", origin=(0.0, 0.0, 0.5 * z_max)).triangulate()
    if sl.n_points == 0:
        raise RuntimeError(f"Mid-plane slice is empty for {case_dir}")
    points = sl.points
    return points[:, 0], points[:, 1], np.asarray(sl["Umag"])


def draw_geometry_panel(ax: plt.Axes, predicted_scene: dict[str, Any]) -> None:
    blocks = room_blocks(predicted_scene)
    x_max, y_max = room_extent(blocks)

    for block in blocks:
        ax.add_patch(
            Rectangle(
                (block["x"], block["y"]),
                block["dx"],
                block["dy"],
                facecolor=COLORS["room_fill"],
                edgecolor=COLORS["room_edge"],
                linewidth=1.7,
                zorder=2,
            )
        )

    for obstacle in predicted_scene.get("obstacles", []):
        ax.add_patch(
            Rectangle(
                (float(obstacle["min"]["x"]), float(obstacle["min"]["y"])),
                float(obstacle["size"]["dx"]),
                float(obstacle["size"]["dy"]),
                facecolor=COLORS["obstacle_fill"],
                edgecolor=COLORS["obstacle_edge"],
                linewidth=1.1,
                zorder=4,
            )
        )

    for opening in predicted_scene.get("openings", []):
        wall, center, size = prepare_opening_geometry(opening, x_max, y_max)
        color = COLORS["inlet"] if opening["type"] == "inlet" else COLORS["outlet"]
        draw_opening(ax, wall, center, size, color=color)

    ax.set_xlim(-0.03 * x_max, x_max * 1.02)
    ax.set_ylim(-0.03 * y_max, y_max * 1.03)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
        spine.set_edgecolor(COLORS["frame"])


def draw_cfd_panel(
    ax: plt.Axes,
    scene: dict[str, Any],
    case_dir: Path,
    norm: Normalize,
) -> None:
    x, y, umag = extract_midplane_slice(scene, case_dir)
    tric = ax.tricontourf(x, y, umag, levels=np.linspace(norm.vmin, norm.vmax, 18), cmap="viridis", norm=norm)
    ax.tricontour(x, y, umag, levels=np.linspace(norm.vmin, norm.vmax, 6), colors="white", linewidths=0.34, alpha=0.22)

    for obstacle in scene.get("obstacles", []):
        patch = Rectangle(
            (float(obstacle["min"]["x"]), float(obstacle["min"]["y"])),
            float(obstacle["size"]["dx"]),
            float(obstacle["size"]["dy"]),
            fill=False,
            edgecolor="white",
            linewidth=1.25,
            linestyle=(0, (3.0, 1.8)),
            alpha=0.98,
            zorder=4,
        )
        patch.set_path_effects([pe.Stroke(linewidth=2.2, foreground="#0F172A", alpha=0.55), pe.Normal()])
        ax.add_patch(patch)

    blocks = room_blocks(scene)
    x_max, y_max = room_extent(blocks)
    ax.set_xlim(0.0, x_max)
    ax.set_ylim(0.0, y_max)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.85)
        spine.set_edgecolor(COLORS["frame"])

    return tric


def metric_line(summary: dict[str, Any]) -> str:
    pred = summary["prediction_summary"]
    struct = float(pred["structural_score"])
    cfd = float(summary["cfd_summary"]["aggregate_score"]["cfd_score"])
    obstacle_f1 = float(pred["obstacle_match"]["f1"])
    return f"struct. {struct:.3f} · obs. {obstacle_f1:.2f} · CFD {cfd:.3f}"


def opening_state(summary: dict[str, Any]) -> str:
    pred = summary["prediction_summary"]
    opening_f1 = float(pred["opening_metrics"]["type_f1"])
    wall_ratio = float(pred["opening_metrics"].get("wall_match_ratio", 0.0))
    if opening_f1 >= 0.99 and wall_ratio >= 0.99:
        return "openings kept"
    if opening_f1 >= 0.66:
        return "openings shifted"
    return "openings degraded"


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"],
            "font.size": 10.0,
            "axes.titlesize": 10.2,
            "axes.labelsize": 9.6,
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    evaluation_root = PROJECT_ROOT / "benchmark/evaluations_posthoc_scaled_longest_span"
    rows: list[dict[str, Any]] = []
    all_umag: list[np.ndarray] = []
    for cfg in SELECTED:
        eval_dir = evaluation_root / cfg["case"] / cfg["view"]
        summary = load_json(eval_dir / "evaluation_summary.json")
        predicted_scene = load_json(eval_dir / "predicted_scene.json")
        reference_scene = load_json(eval_dir / "reference_scene.json")
        ref_case_dir = Path(summary["cfd_summary"]["reference"]["case"])
        pred_case_dir = Path(summary["cfd_summary"]["predicted"]["case"])
        ref_slice = extract_midplane_slice(reference_scene, ref_case_dir)
        pred_slice = extract_midplane_slice(predicted_scene, pred_case_dir)
        all_umag.extend([ref_slice[2], pred_slice[2]])
        rows.append(
            {
                "cfg": cfg,
                "summary": summary,
                "predicted_scene": predicted_scene,
                "reference_scene": reference_scene,
                "ref_case_dir": ref_case_dir,
                "pred_case_dir": pred_case_dir,
            }
        )

    vmax = max(float(np.percentile(values, 99.0)) for values in all_umag)
    norm = Normalize(vmin=0.0, vmax=vmax)

    fig = plt.figure(figsize=(7.35, 5.6), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        3,
        left=0.06,
        right=0.925,
        bottom=0.105,
        top=0.865,
        hspace=0.18,
        wspace=0.05,
        width_ratios=[1.18, 1.52, 1.52],
    )

    fig.text(0.205, 0.905, "predicted geometry", ha="center", va="center", fontsize=11.8, fontweight="bold", color=COLORS["title"])
    fig.text(0.535, 0.905, "reference CFD", ha="center", va="center", fontsize=11.8, fontweight="bold", color=COLORS["title"])
    fig.text(0.80, 0.905, "predicted CFD", ha="center", va="center", fontsize=11.8, fontweight="bold", color=COLORS["title"])
    fig.text(0.5, 0.948, "dense composite: strong structure, weaker flow fidelity", ha="center", va="center", fontsize=11.0, fontweight="semibold", color=COLORS["subtitle"])

    legend_handles = [
        Patch(facecolor=COLORS["obstacle_fill"], edgecolor=COLORS["obstacle_edge"], label="obstacle"),
        Line2D([0], [0], color=COLORS["inlet"], lw=3.0, label="inlet"),
        Line2D([0], [0], color=COLORS["outlet"], lw=3.0, label="outlet"),
        Line2D([0], [0], color="white", lw=1.3, linestyle=(0, (3.0, 1.8)), label="pred. obstacle"),
    ]

    first_geo_ax = None
    for row_idx, row in enumerate(rows):
        summary = row["summary"]
        predicted_scene = row["predicted_scene"]
        reference_scene = row["reference_scene"]

        ax_geo = fig.add_subplot(gs[row_idx, 0])
        if first_geo_ax is None:
            first_geo_ax = ax_geo
        draw_geometry_panel(ax_geo, predicted_scene)

        ax_ref = fig.add_subplot(gs[row_idx, 1])
        draw_cfd_panel(ax_ref, reference_scene, row["ref_case_dir"], norm)

        ax_pred = fig.add_subplot(gs[row_idx, 2])
        draw_cfd_panel(ax_pred, predicted_scene, row["pred_case_dir"], norm)

        geo_pos = ax_geo.get_position()
        row_y = max(ax_geo.get_position().y1, ax_ref.get_position().y1) + 0.006
        fig.text(
            geo_pos.x0,
            row_y,
            row["cfg"]["row_label"],
            ha="left",
            va="bottom",
            fontsize=11.7,
            fontweight="bold",
            color=COLORS["title"],
        )
        fig.text(
            geo_pos.x0 + 0.16,
            row_y,
            f"{metric_line(summary)} · {opening_state(summary)}",
            ha="left",
            va="bottom",
            fontsize=9.9,
            fontweight="semibold",
            color=COLORS["subtitle"],
        )

    if first_geo_ax is not None:
        first_geo_ax.legend(
            handles=legend_handles,
            loc="upper left",
            bbox_to_anchor=(0.01, 0.995),
            ncol=2,
            frameon=True,
            framealpha=0.94,
            facecolor="white",
            edgecolor=COLORS["metric_edge"],
            handlelength=1.8,
            columnspacing=0.9,
            fontsize=8.2,
            borderpad=0.35,
            labelspacing=0.4,
        )

    cax = fig.add_axes([0.935, 0.18, 0.018, 0.59])
    sm = ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
    cbar.set_label("shared |U| [m/s]", fontsize=9.8)
    ticks = np.linspace(norm.vmin, norm.vmax, 4)
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels([f"{tick:.3f}" for tick in ticks])
    cbar.ax.tick_params(labelsize=8.6, length=2.5)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT)
    fig.savefig(PNG_OUT, dpi=600)
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
