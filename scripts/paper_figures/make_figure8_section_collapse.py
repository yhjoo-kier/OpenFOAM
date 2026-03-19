#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results/paper_figures"
PDF_OUT = OUT_DIR / "figure8_section_view_composite_collapse.pdf"
PNG_OUT = OUT_DIR / "figure8_section_view_composite_collapse.png"

CASES = [
    {
        "case": "bench_a3_04",
        "row_title": "A3-04 sparse composite",
        "summary": PROJECT_ROOT / "benchmark/evaluations/bench_a3_04/section/evaluation_summary.json",
        "reference": PROJECT_ROOT / "benchmark/evaluations/bench_a3_04/section/reference_scene.json",
        "predicted": PROJECT_ROOT / "benchmark/evaluations/bench_a3_04/section/predicted_scene.json",
    },
    {
        "case": "bench_a4_05",
        "row_title": "A4-05 dense composite",
        "summary": PROJECT_ROOT / "benchmark/evaluations/bench_a4_05/section/evaluation_summary.json",
        "reference": PROJECT_ROOT / "benchmark/evaluations/bench_a4_05/section/reference_scene.json",
        "predicted": PROJECT_ROOT / "benchmark/evaluations/bench_a4_05/section/predicted_scene.json",
    },
]

COLORS = {
    "room_fill": "#EEF3F8",
    "room_edge": "#1F2937",
    "obstacle_fill": "#9AA5B1",
    "obstacle_edge": "#4B5563",
    "inlet": "#2C7FB8",
    "outlet": "#CB3A31",
    "reference_overlay": "#C0392B",
    "missing_fill": "#FDECEC",
    "title": "#111827",
    "subtitle": "#334155",
    "metric_box": "#F8FAFC",
    "grid": "#E5E7EB",
}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def room_blocks(scene: dict[str, Any]) -> list[dict[str, float]]:
    room = scene["room"]
    if "blocks" in room:
        blocks = room["blocks"]
        return [
            {
                "x": float(block["origin"]["x"]),
                "y": float(block["origin"]["y"]),
                "dx": float(block["size"]["dx"]),
                "dy": float(block["size"]["dy"]),
            }
            for block in blocks
        ]
    size = room["size"]
    return [{"x": 0.0, "y": 0.0, "dx": float(size["Lx"]), "dy": float(size["Ly"])}]


def room_extent(blocks: list[dict[str, float]]) -> tuple[float, float]:
    max_x = max(block["x"] + block["dx"] for block in blocks)
    max_y = max(block["y"] + block["dy"] for block in blocks)
    return max_x, max_y


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


def prepare_opening_geometry(opening: dict[str, Any], x_max: float, y_max: float) -> tuple[str, dict[str, float], dict[str, float]]:
    wall = opening["wall"]
    center = {"u": float(opening["center"]["u"]), "v": float(opening["center"]["v"]), "wall_x_max": x_max, "wall_y_max": y_max}
    size = {"du": float(opening["size"]["du"]), "dv": float(opening["size"]["dv"])}
    return wall, center, size


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=9.8, fontweight="bold", color=COLORS["title"], ha="left", va="bottom")


def draw_scene(ax: plt.Axes, scene: dict[str, Any], *, title: str, overlay_scene: dict[str, Any] | None = None, annotate_collapse: bool = False, metric_text: str | None = None, panel_extent: tuple[float, float] | None = None) -> None:
    blocks = room_blocks(scene)
    x_max, y_max = room_extent(blocks)
    if panel_extent is not None:
        x_max, y_max = panel_extent

    for block in blocks:
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

    for obstacle in scene.get("obstacles", []):
        x = float(obstacle["min"]["x"])
        y = float(obstacle["min"]["y"])
        dx = float(obstacle["size"]["dx"])
        dy = float(obstacle["size"]["dy"])
        ax.add_patch(
            Rectangle(
                (x, y),
                dx,
                dy,
                facecolor=COLORS["obstacle_fill"],
                edgecolor=COLORS["obstacle_edge"],
                linewidth=0.9,
                zorder=4,
            )
        )

    if overlay_scene is not None:
        overlay_blocks = room_blocks(overlay_scene)
        for block in overlay_blocks:
            ax.add_patch(
                Rectangle(
                    (block["x"], block["y"]),
                    block["dx"],
                    block["dy"],
                    facecolor="none",
                    edgecolor=COLORS["reference_overlay"],
                    linewidth=1.4,
                    linestyle=(0, (4, 2)),
                    alpha=0.95,
                    zorder=6,
                )
            )

    for opening in scene.get("openings", []):
        wall, center, size = prepare_opening_geometry(opening, x_max, y_max)
        color = COLORS["inlet"] if opening["type"] == "inlet" else COLORS["outlet"]
        draw_opening(ax, wall, center, size, color=color)

    ax.set_xlim(-0.18 * x_max, x_max * 1.05)
    ax.set_ylim(-0.10 * y_max, y_max * 1.10)
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

    if annotate_collapse:
        ax.text(
            0.02,
            0.98,
            "section-view prediction\ncomposite → rectangular",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.3,
            color=COLORS["reference_overlay"],
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": COLORS["reference_overlay"], "linewidth": 0.9},
            zorder=10,
        )
    if metric_text:
        ax.text(
            0.98,
            0.08,
            metric_text,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.7,
            color=COLORS["subtitle"],
            bbox={"boxstyle": "round,pad=0.28", "facecolor": COLORS["metric_box"], "edgecolor": "#CBD5E1", "linewidth": 0.8},
            zorder=10,
        )


def collapse_metric_text(summary: dict[str, Any]) -> str:
    pred = summary["prediction_summary"]
    cfd = summary["cfd_summary"]["aggregate_score"]["cfd_score"]
    ref_open = "/".join(w[0].upper() for w in pred["reference_opening_walls"])
    pred_open = "/".join(w[0].upper() for w in pred["predicted_opening_walls"])
    return (
        f"S = {pred['structural_score']:.3f}\n"
        f"CFD = {cfd:.3f}\n"
        f"blocks {pred['reference_room_block_count']} → {pred['predicted_room_block_count']}\n"
        f"openings {ref_open} → {pred_open}"
    )


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.6,
            "axes.titlesize": 9.8,
            "axes.labelsize": 9.1,
            "xtick.labelsize": 8.4,
            "ytick.labelsize": 8.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(7.30, 5.65), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.985, bottom=0.10, top=0.88, wspace=0.18, hspace=0.32)

    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for row_idx, case in enumerate(CASES):
        summary = load_json(case["summary"])
        reference_scene = load_json(case["reference"])
        predicted_scene = load_json(case["predicted"])
        ref_extent = room_extent(room_blocks(reference_scene))
        pred_extent = room_extent(room_blocks(predicted_scene))
        panel_extent = (max(ref_extent[0], pred_extent[0]), max(ref_extent[1], pred_extent[1]))

        ax_ref = fig.add_subplot(gs[row_idx, 0])
        draw_scene(
            ax_ref,
            reference_scene,
            title="Reference composite geometry",
            panel_extent=panel_extent,
        )
        ax_ref.text(0.0, 1.11, case["row_title"], transform=ax_ref.transAxes, ha="left", va="bottom", fontsize=9.8, fontweight="bold", color=COLORS["title"])
        add_panel_label(ax_ref, panel_labels[row_idx * 2])

        ax_pred = fig.add_subplot(gs[row_idx, 1])
        draw_scene(
            ax_pred,
            predicted_scene,
            title="Prediction from section view",
            overlay_scene=reference_scene,
            annotate_collapse=True,
            metric_text=collapse_metric_text(summary),
            panel_extent=panel_extent,
        )
        add_panel_label(ax_pred, panel_labels[row_idx * 2 + 1])

    fig.suptitle(
        "Figure 8. Section-view input caused the only composite→rectangular room collapses",
        x=0.5,
        y=0.965,
        fontsize=10.9,
        fontweight="bold",
        color=COLORS["title"],
    )
    fig.text(
        0.5,
        0.925,
        "Dashed red outlines mark the missing composite arm in the predicted top-view abstraction.",
        ha="center",
        va="center",
        fontsize=8.6,
        color=COLORS["subtitle"],
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
