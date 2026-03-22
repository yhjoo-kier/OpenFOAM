#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "results/paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_discuss_section_collapse.pdf"
PNG_OUT = OUT_DIR / "fig_discuss_section_collapse.png"
META_OUT = OUT_DIR / "fig_discuss_section_collapse_meta.json"
AGG_PATH = PROJECT_ROOT / "benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

CASES = [
    {
        "case": "bench_a3_04",
        "view": "section",
        "row_title": "A3-04 sparse composite",
        "summary": PROJECT_ROOT / "benchmark/evaluations/bench_a3_04/section/evaluation_summary.json",
        "reference": PROJECT_ROOT / "benchmark/evaluations/bench_a3_04/section/reference_scene.json",
        "predicted": PROJECT_ROOT / "benchmark/evaluations/bench_a3_04/section/predicted_scene.json",
    },
    {
        "case": "bench_a4_05",
        "view": "section",
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
    "title": "#111827",
    "subtitle": "#334155",
    "metric_box": "#F8FAFC",
    "grid": "#E5E7EB",
    "panel": "#34495E",
}


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return "sans-serif"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def build_scaled_metric_lookup() -> dict[tuple[str, str], dict[str, float]]:
    payload = load_json(AGG_PATH)
    lookup: dict[tuple[str, str], dict[str, float]] = {}
    for entry in payload["derived_tags"]["task_examples"]["section_room_kind_collapse"]:
        lookup[(entry["case"], entry["view"])] = {
            "structural_score": float(entry["structural_score"]),
            "cfd_score": float(entry["cfd_score"]),
        }
    return lookup


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


def draw_opening(ax: plt.Axes, wall: str, center: dict[str, float], size: dict[str, float], color: str, lw: float = 2.5) -> None:
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
    center = {
        "u": float(opening["center"]["u"]),
        "v": float(opening["center"]["v"]),
        "wall_x_max": x_max,
        "wall_y_max": y_max,
    }
    size = {"du": float(opening["size"]["du"]), "dv": float(opening["size"]["dv"])}
    return wall, center, size


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=10.4,
        fontweight="bold",
        color=COLORS["panel"],
        ha="left",
        va="bottom",
    )


def draw_scene(
    ax: plt.Axes,
    scene: dict[str, Any],
    *,
    title: str,
    overlay_scene: dict[str, Any] | None = None,
    collapse_note: str | None = None,
    metric_text: str | None = None,
    panel_extent: tuple[float, float] | None = None,
) -> None:
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
        for block in room_blocks(overlay_scene):
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

    ax.set_xlim(-0.16 * x_max, x_max * 1.04)
    ax.set_ylim(-0.09 * y_max, y_max * 1.09)
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.grid(True, color=COLORS["grid"], linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(title, loc="left", fontsize=9.8, fontweight="bold", color=COLORS["subtitle"], pad=6)
    ax.tick_params(labelsize=8.4)

    if collapse_note:
        ax.text(
            0.03,
            0.97,
            collapse_note,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.4,
            color=COLORS["reference_overlay"],
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.26", "facecolor": "white", "edgecolor": COLORS["reference_overlay"], "linewidth": 0.9},
            zorder=10,
        )

    if metric_text:
        ax.text(
            0.98,
            0.09,
            metric_text,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.7,
            color=COLORS["subtitle"],
            bbox={"boxstyle": "round,pad=0.28", "facecolor": COLORS["metric_box"], "edgecolor": "#CBD5E1", "linewidth": 0.8},
            zorder=10,
        )


def collapse_metric_text(summary: dict[str, Any], scaled_metrics: dict[str, float]) -> str:
    pred = summary["prediction_summary"]
    ref_open = "/".join(w[0].upper() for w in pred["reference_opening_walls"])
    pred_open = "/".join(w[0].upper() for w in pred["predicted_opening_walls"])
    return (
        f"S = {scaled_metrics['structural_score']:.3f}\n"
        f"CFD = {scaled_metrics['cfd_score']:.3f}\n"
        f"blocks {pred['reference_room_block_count']} → {pred['predicted_room_block_count']}\n"
        f"openings {ref_open} → {pred_open}"
    )


def main() -> None:
    selected_font = pick_font()
    plt.rcParams.update(
        {
            "font.family": selected_font,
            "font.size": 9.8,
            "axes.titlesize": 9.8,
            "axes.labelsize": 9.2,
            "xtick.labelsize": 8.4,
            "ytick.labelsize": 8.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    scaled_lookup = build_scaled_metric_lookup()
    fig = plt.figure(figsize=(7.25, 5.90), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.985, bottom=0.14, top=0.94, wspace=0.19, hspace=0.40)
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for row_idx, case in enumerate(CASES):
        summary = load_json(case["summary"])
        reference_scene = load_json(case["reference"])
        predicted_scene = load_json(case["predicted"])
        scaled_metrics = scaled_lookup[(case["case"], case["view"])]

        ref_extent = room_extent(room_blocks(reference_scene))
        pred_extent = room_extent(room_blocks(predicted_scene))
        panel_extent = (max(ref_extent[0], pred_extent[0]), max(ref_extent[1], pred_extent[1]))

        ax_ref = fig.add_subplot(gs[row_idx, 0])
        draw_scene(ax_ref, reference_scene, title="Reference", panel_extent=panel_extent)
        ax_ref.text(
            0.0,
            1.11,
            case["row_title"],
            transform=ax_ref.transAxes,
            ha="left",
            va="bottom",
            fontsize=10.0,
            fontweight="bold",
            color=COLORS["title"],
        )
        add_panel_label(ax_ref, panel_labels[row_idx * 2])

        ax_pred = fig.add_subplot(gs[row_idx, 1])
        draw_scene(
            ax_pred,
            predicted_scene,
            title="Prediction + GT outline",
            overlay_scene=reference_scene,
            collapse_note="Composite → rectangular",
            metric_text=collapse_metric_text(summary, scaled_metrics),
            panel_extent=panel_extent,
        )
        add_panel_label(ax_pred, panel_labels[row_idx * 2 + 1])

    # --- Legend -----------------------------------------------------------
    legend_handles = [
        Patch(facecolor=COLORS["room_fill"], edgecolor=COLORS["room_edge"],
              linewidth=1.0, label="Predicted room"),
        Patch(facecolor="none", edgecolor=COLORS["reference_overlay"],
              linewidth=1.4, linestyle=(0, (4, 2)), label="Reference (GT) room"),
        Patch(facecolor=COLORS["obstacle_fill"], edgecolor=COLORS["obstacle_edge"],
              linewidth=0.9, hatch="///", label="Obstacle"),
        Line2D([], [], color=COLORS["inlet"], linewidth=2.5,
               solid_capstyle="round", label="Inlet"),
        Line2D([], [], color=COLORS["outlet"], linewidth=2.5,
               solid_capstyle="round", label="Outlet"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=5,
        fontsize=11,
        frameon=True,
        edgecolor="#CBD5E1",
        fancybox=True,
        borderpad=0.5,
        columnspacing=1.2,
        handletextpad=0.5,
        bbox_to_anchor=(0.5, 0.005),
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)

    meta = {
        "source_artifacts": {
            "scaled_manifest": str(AGG_PATH),
            "task_summaries": [str(case["summary"]) for case in CASES],
            "reference_scenes": [str(case["reference"]) for case in CASES],
            "predicted_scenes": [str(case["predicted"]) for case in CASES],
        },
        "font_family_requested": FONT_CANDIDATES,
        "font_family_selected": selected_font,
        "intended_width": "double-column",
        "panel_layout": "2x2",
        "subfigure_labels": panel_labels,
        "png_dpi": 600,
        "pdf_vector": True,
        "internal_caption_text": False,
        "notes": {
            "metric_setting": "posthoc_scaled_longest_span",
            "cases_locked": [f"{case['case']}/{case['view']}" for case in CASES],
        },
    }
    META_OUT.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"Selected font: {selected_font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {META_OUT}")


if __name__ == "__main__":
    main()
