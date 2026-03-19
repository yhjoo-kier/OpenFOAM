#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results/paper_figures"
PDF_OUT = OUT_DIR / "figure9_obstacle_hallucination_limited_cfd_penalty.pdf"
PNG_OUT = OUT_DIR / "figure9_obstacle_hallucination_limited_cfd_penalty.png"

CASES = [
    {
        "case": "bench_a3_01",
        "view": "floorplan",
        "row_title": "A3-01 empty composite",
        "tagline": "",
    },
    {
        "case": "bench_a3_05",
        "view": "floorplan",
        "row_title": "A3-05 one-obstacle composite",
        "tagline": "",
    },
]

COLORS = {
    "room_fill": "#EEF3F8",
    "room_edge": "#1F2937",
    "obstacle_fill": "#E5A657",
    "obstacle_edge": "#8A4B08",
    "reference_obstacle_outline": "#2563EB",
    "inlet": "#2563EB",
    "outlet": "#C83E3A",
    "title": "#111827",
    "subtitle": "#334155",
    "metric_box": "#F8FAFC",
    "grid": "#E5E7EB",
    "callout": "#FFF7ED",
    "callout_edge": "#D97706",
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
                "dx": float(block["size"]["dx"]),
                "dy": float(block["size"]["dy"]),
            }
            for block in room["blocks"]
        ]
    size = room["size"]
    return [{"x": 0.0, "y": 0.0, "dx": float(size["Lx"]), "dy": float(size["Ly"])}]


def room_extent(blocks: list[dict[str, float]]) -> tuple[float, float]:
    max_x = max(block["x"] + block["dx"] for block in blocks)
    max_y = max(block["y"] + block["dy"] for block in blocks)
    return max_x, max_y


def prepare_opening_geometry(opening: dict[str, Any], x_max: float, y_max: float) -> tuple[str, dict[str, float], dict[str, float]]:
    wall = opening["wall"]
    center = {
        "u": float(opening["center"]["u"]),
        "wall_x_max": x_max,
        "wall_y_max": y_max,
    }
    size = {"du": float(opening["size"]["du"]), "dv": float(opening["size"]["dv"])}
    return wall, center, size


def draw_opening(ax: plt.Axes, wall: str, center: dict[str, float], size: dict[str, float], color: str, lw: float = 2.6) -> None:
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


def draw_scene(
    ax: plt.Axes,
    scene: dict[str, Any],
    *,
    title: str,
    panel_extent: tuple[float, float],
    reference_scene: dict[str, Any] | None = None,
    metric_text: str | None = None,
    callout_text: str | None = None,
) -> None:
    x_max, y_max = panel_extent
    for block in room_blocks(scene):
        ax.add_patch(
            Rectangle(
                (block["x"], block["y"]),
                block["dx"],
                block["dy"],
                facecolor=COLORS["room_fill"],
                edgecolor=COLORS["room_edge"],
                linewidth=1.3,
                zorder=2,
            )
        )

    if reference_scene is not None:
        for obstacle in reference_scene.get("obstacles", []):
            ax.add_patch(
                Rectangle(
                    (float(obstacle["min"]["x"]), float(obstacle["min"]["y"])),
                    float(obstacle["size"]["dx"]),
                    float(obstacle["size"]["dy"]),
                    facecolor="none",
                    edgecolor=COLORS["reference_obstacle_outline"],
                    linewidth=1.4,
                    linestyle=(0, (4, 2)),
                    zorder=5,
                )
            )

    for obstacle in scene.get("obstacles", []):
        ax.add_patch(
            Rectangle(
                (float(obstacle["min"]["x"]), float(obstacle["min"]["y"])),
                float(obstacle["size"]["dx"]),
                float(obstacle["size"]["dy"]),
                facecolor=COLORS["obstacle_fill"],
                edgecolor=COLORS["obstacle_edge"],
                linewidth=0.95,
                zorder=4,
            )
        )

    for opening in scene.get("openings", []):
        wall, center, size = prepare_opening_geometry(opening, x_max, y_max)
        color = COLORS["inlet"] if opening["type"] == "inlet" else COLORS["outlet"]
        draw_opening(ax, wall, center, size, color=color)

    ax.set_xlim(-0.08 * x_max, x_max * 1.04)
    ax.set_ylim(-0.10 * y_max, y_max * 1.08)
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.grid(True, color=COLORS["grid"], linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(title, loc="left", fontsize=9.5, fontweight="bold", color=COLORS["subtitle"], pad=6)
    ax.tick_params(labelsize=8.0)

    if metric_text:
        ax.text(
            0.98,
            0.06,
            metric_text,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.8,
            color=COLORS["subtitle"],
            bbox={"boxstyle": "round,pad=0.30", "facecolor": COLORS["metric_box"], "edgecolor": "#CBD5E1", "linewidth": 0.9},
            zorder=10,
        )
    if callout_text:
        ax.text(
            0.02,
            0.98,
            callout_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.15,
            color=COLORS["callout_edge"],
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.28", "facecolor": COLORS["callout"], "edgecolor": COLORS["callout_edge"], "linewidth": 0.9},
            zorder=10,
        )


def metric_text(summary: dict[str, Any]) -> str:
    pred = summary["prediction_summary"]
    cfd = summary["cfd_summary"]["aggregate_score"]["cfd_score"]
    openings_ref = "/".join(w[0].upper() for w in pred["reference_opening_walls"])
    openings_pred = "/".join(w[0].upper() for w in pred["predicted_opening_walls"])
    return (
        f"S = {pred['structural_score']:.3f}\n"
        f"CFD = {cfd:.3f}\n"
        f"obstacles {pred['reference_obstacle_count']} → {pred['predicted_obstacle_count']}\n"
        f"openings {openings_ref} → {openings_pred}"
    )


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.4,
            "axes.titlesize": 9.7,
            "axes.labelsize": 9.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(7.40, 5.80), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        3,
        left=0.06,
        right=0.988,
        bottom=0.14,
        top=0.84,
        wspace=0.22,
        hspace=0.44,
        width_ratios=[0.94, 1.0, 1.08],
    )

    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

    for row_idx, case_cfg in enumerate(CASES):
        case = case_cfg["case"]
        view = case_cfg["view"]
        eval_dir = PROJECT_ROOT / "benchmark" / "evaluations" / case / view
        task = load_json(eval_dir / "task.json")
        summary = load_json(eval_dir / "evaluation_summary.json")
        reference_scene = load_json(Path(task["reference_scene"]))
        predicted_scene = load_json(eval_dir / "predicted_scene.json")
        panel_extent = tuple(
            max(a, b) for a, b in zip(room_extent(room_blocks(reference_scene)), room_extent(room_blocks(predicted_scene)))
        )

        ax_img = fig.add_subplot(gs[row_idx, 0])
        image = mpimg.imread(task["input_image"])
        ax_img.imshow(image)
        ax_img.set_xticks([])
        ax_img.set_yticks([])
        for spine in ax_img.spines.values():
            spine.set_visible(False)
        ax_img.set_title(f"{panel_labels[row_idx * 3]} {case_cfg['row_title']}\nInput floor-plan image", loc="left", fontsize=9.3, fontweight="bold", color=COLORS["subtitle"], pad=5)

        ax_ref = fig.add_subplot(gs[row_idx, 1])
        draw_scene(ax_ref, reference_scene, title=f"{panel_labels[row_idx * 3 + 1]} Reference geometry", panel_extent=panel_extent)

        ax_pred = fig.add_subplot(gs[row_idx, 2])
        callout = "extra obstacles hallucinated\nwhile room topology stays composite"
        if row_idx == 1:
            callout = "obstacle count inflates to three\nbut opening topology still matches"
        draw_scene(
            ax_pred,
            predicted_scene,
            title=f"{panel_labels[row_idx * 3 + 2]} Predicted geometry",
            panel_extent=panel_extent,
            reference_scene=reference_scene,
            metric_text=metric_text(summary),
            callout_text=callout,
        )

    fig.suptitle(
        "Figure 9. A3 composite cases can hallucinate obstacles while keeping strong CFD fidelity",
        x=0.5,
        y=0.965,
        fontsize=10.9,
        fontweight="bold",
        color=COLORS["title"],
    )
    fig.text(
        0.5,
        0.924,
        "Both examples preserve the opening-wall topology (blue/red) and keep CFD score near 0.69 even though obstacle count is inflated.",
        ha="center",
        va="center",
        fontsize=8.6,
        color=COLORS["subtitle"],
    )
    fig.text(
        0.5,
        0.080,
        "Predicted panels use orange filled boxes for inferred obstacles; dashed blue outlines indicate reference obstacles when present.",
        ha="center",
        va="center",
        fontsize=8.1,
        color=COLORS["subtitle"],
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
