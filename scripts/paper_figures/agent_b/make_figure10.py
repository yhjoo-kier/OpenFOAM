#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "results/paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_discuss_structure_cfd_gap.pdf"
PNG_OUT = OUT_DIR / "fig_discuss_structure_cfd_gap.png"

SELECTED = [
    {
        "case": "bench_a4_02",
        "view": "floorplan",
        "geometry_title": "(a) A4-02",
        "flow_title": "(b) A4-02",
        "geometry_callout": "openings kept",
        "flow_callout": "moderate CFD",
    },
    {
        "case": "bench_a4_04",
        "view": "floorplan",
        "geometry_title": "(c) A4-04",
        "flow_title": "(d) A4-04",
        "geometry_callout": "strong structure",
        "flow_callout": "moderate CFD",
    },
]

COLORS = {
    "title": "#111827",
    "subtitle": "#334155",
    "grid": "#E5E7EB",
    "room_fill": "#EEF3F8",
    "room_edge": "#1F2937",
    "obstacle_fill": "#E5A657",
    "obstacle_edge": "#8A4B08",
    "reference_outline": "#2563EB",
    "inlet": "#2563EB",
    "outlet": "#C83E3A",
    "metric_box": "#F8FAFC",
    "callout_bg": "#FFF7ED",
    "callout_edge": "#D97706",
    "bar_struct": "#475569",
    "bar_cfd": "#0F766E",
    "divider": "#CBD5E1",
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


def draw_opening(ax: plt.Axes, wall: str, center: dict[str, float], size: dict[str, float], color: str, lw: float = 2.8) -> None:
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


def crop_nonwhite(img: np.ndarray, threshold: float = 0.985, pad: int = 10) -> np.ndarray:
    if img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)
    if img.shape[2] == 4:
        img = img[:, :, :3]
    mask = np.any(img < threshold, axis=2)
    if not mask.any():
        return img
    ys, xs = np.where(mask)
    y0 = max(int(ys.min()) - pad, 0)
    y1 = min(int(ys.max()) + pad + 1, img.shape[0])
    x0 = max(int(xs.min()) - pad, 0)
    x1 = min(int(xs.max()) + pad + 1, img.shape[1])
    return img[y0:y1, x0:x1]


def pad_to_height(img: np.ndarray, target_h: int) -> np.ndarray:
    if img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)
    if img.shape[2] == 4:
        img = img[:, :, :3]
    h, _, _ = img.shape
    if h == target_h:
        return img
    pad_top = (target_h - h) // 2
    pad_bottom = target_h - h - pad_top
    return np.pad(img, ((pad_top, pad_bottom), (0, 0), (0, 0)), mode="constant", constant_values=1.0)


def score_bar_inset(ax: plt.Axes, structural: float, cfd: float) -> None:
    inset = ax.inset_axes([0.03, 0.03, 0.44, 0.19])
    inset.barh([1, 0], [structural, cfd], color=[COLORS["bar_struct"], COLORS["bar_cfd"]], height=0.42)
    inset.set_xlim(0, 1.0)
    inset.set_yticks([1, 0])
    inset.set_yticklabels(["Structural", "CFD"], fontsize=9.4)
    inset.set_xticks([])
    inset.grid(False)
    inset.set_axisbelow(True)
    inset.spines["top"].set_visible(False)
    inset.spines["right"].set_visible(False)
    inset.text(structural + 0.02, 1, f"{structural:.3f}", va="center", fontsize=8.2, color=COLORS["title"])
    inset.text(cfd + 0.02, 0, f"{cfd:.3f}", va="center", fontsize=8.2, color=COLORS["title"])


def draw_geometry_panel(
    ax: plt.Axes,
    reference_scene: dict[str, Any],
    predicted_scene: dict[str, Any],
    input_img: np.ndarray,
    *,
    title: str,
    callout_text: str,
    summary: dict[str, Any],
) -> None:
    panel_extent = tuple(
        max(a, b) for a, b in zip(room_extent(room_blocks(reference_scene)), room_extent(room_blocks(predicted_scene)))
    )
    x_max, y_max = panel_extent

    for block in room_blocks(predicted_scene):
        ax.add_patch(
            Rectangle(
                (block["x"], block["y"]),
                block["dx"],
                block["dy"],
                facecolor=COLORS["room_fill"],
                edgecolor=COLORS["room_edge"],
                linewidth=1.5,
                zorder=2,
            )
        )

    for block in room_blocks(reference_scene):
        ax.add_patch(
            Rectangle(
                (block["x"], block["y"]),
                block["dx"],
                block["dy"],
                facecolor="none",
                edgecolor=COLORS["reference_outline"],
                linewidth=1.45,
                linestyle=(0, (4, 2)),
                zorder=6,
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
                linewidth=1.0,
                zorder=4,
            )
        )

    for obstacle in reference_scene.get("obstacles", []):
        ax.add_patch(
            Rectangle(
                (float(obstacle["min"]["x"]), float(obstacle["min"]["y"])),
                float(obstacle["size"]["dx"]),
                float(obstacle["size"]["dy"]),
                facecolor="none",
                edgecolor=COLORS["reference_outline"],
                linewidth=1.15,
                linestyle=(0, (4, 2)),
                zorder=5,
            )
        )

    for opening in predicted_scene.get("openings", []):
        wall, center, size = prepare_opening_geometry(opening, x_max, y_max)
        color = COLORS["inlet"] if opening["type"] == "inlet" else COLORS["outlet"]
        draw_opening(ax, wall, center, size, color=color)

    ax.set_xlim(-0.06 * x_max, x_max * 1.03)
    ax.set_ylim(-0.10 * y_max, y_max * 1.10)
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.tick_params(labelsize=9.0)
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", color=COLORS["subtitle"], pad=10)

    pred = summary["prediction_summary"]
    structural = float(pred["structural_score"])
    metric_text = (
        f"struct. {structural:.3f}\n"
        f"obs F1 {pred['obstacle_match']['f1']:.2f}\n"
        f"openings {'yes' if pred['opening_wall_match'] else 'no'}"
    )
    ax.text(
        0.66,
        0.10,
        metric_text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        color=COLORS["subtitle"],
        bbox={"boxstyle": "round,pad=0.34", "facecolor": COLORS["metric_box"], "edgecolor": COLORS["divider"], "linewidth": 0.95, "alpha": 0.98},
        zorder=10,
    )
    ax.text(
        0.02,
        0.92,
        callout_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        color=COLORS["callout_edge"],
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.28", "facecolor": COLORS["callout_bg"], "edgecolor": COLORS["callout_edge"], "linewidth": 0.9},
        zorder=10,
    )

    inset = ax.inset_axes([0.62, 0.60, 0.30, 0.26])
    inset.imshow(crop_nonwhite(input_img, threshold=0.995, pad=4))
    inset.set_xticks([])
    inset.set_yticks([])
    for spine in inset.spines.values():
        spine.set_edgecolor(COLORS["divider"])
        spine.set_linewidth(0.95)


def draw_flow_panel(ax: plt.Axes, ref_img: np.ndarray, pred_img: np.ndarray, *, title: str, cfd_score: float, callout_text: str) -> None:
    ref_img = crop_nonwhite(ref_img, threshold=0.992, pad=8)
    pred_img = crop_nonwhite(pred_img, threshold=0.992, pad=8)
    top_trim_ref = int(ref_img.shape[0] * 0.04)
    top_trim_pred = int(pred_img.shape[0] * 0.04)
    ref_img = ref_img[top_trim_ref:, :, :]
    pred_img = pred_img[top_trim_pred:, :, :]
    target_h = max(ref_img.shape[0], pred_img.shape[0])
    ref_pad = pad_to_height(ref_img, target_h)
    pred_pad = pad_to_height(pred_img, target_h)
    gap = np.ones((target_h, 16, 3), dtype=float)
    combo = np.concatenate([ref_pad, gap, pred_pad], axis=1)
    ax.imshow(combo)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", color=COLORS["subtitle"], pad=10)
    ax.text(0.03, 0.97, callout_text, transform=ax.transAxes, ha="left", va="top", fontsize=9.4, fontweight="bold", color=COLORS["title"], bbox={"boxstyle": "round,pad=0.22", "facecolor": "#F1F5F9", "edgecolor": COLORS["subtitle"], "linewidth": 0.85})
    ax.text(0.24, 0.88, "Reference flow", transform=ax.transAxes, ha="center", va="top", fontsize=9.1, fontweight="bold", color=COLORS["subtitle"], bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": COLORS["divider"], "linewidth": 0.8})
    ax.text(0.77, 0.88, "Predicted flow", transform=ax.transAxes, ha="center", va="top", fontsize=9.1, fontweight="bold", color=COLORS["subtitle"], bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": COLORS["divider"], "linewidth": 0.8})
    ax.plot([0.5, 0.5], [0.04, 0.94], transform=ax.transAxes, color=COLORS["divider"], lw=1.0)
    ax.text(
        0.5,
        0.06,
        f"CFD {cfd_score:.3f}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10.0,
        fontweight="bold",
        color=COLORS["subtitle"],
        bbox={"boxstyle": "round,pad=0.26", "facecolor": "white", "edgecolor": COLORS["divider"], "linewidth": 0.9},
    )


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"],
            "font.size": 9.7,
            "axes.titlesize": 9.8,
            "axes.labelsize": 9.3,
            "xtick.labelsize": 8.3,
            "ytick.labelsize": 8.3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    aggregate = load_json(PROJECT_ROOT / "benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json")
    a4 = aggregate["by_category"]["A4"]

    fig = plt.figure(figsize=(7.35, 7.85), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, left=0.065, right=0.988, bottom=0.135, top=0.905, hspace=0.30, wspace=0.24, width_ratios=[1.00, 1.34])

    evaluation_root = PROJECT_ROOT / "benchmark/evaluations_posthoc_scaled_longest_span"

    for row, cfg in enumerate(SELECTED):
        eval_dir = evaluation_root / cfg["case"] / cfg["view"]
        task = load_json(eval_dir / "task.json")
        summary = load_json(eval_dir / "evaluation_summary.json")
        reference_scene = load_json(Path(task["reference_scene"]))
        predicted_scene = load_json(eval_dir / "predicted_scene.json")
        input_img = mpimg.imread(task["input_image"])
        ref_flow = mpimg.imread(Path(task["reference_results"]) / "panel_flow_3d.png")
        pred_flow = mpimg.imread(Path(task["actual_outputs"]["predicted_results_dir"]) / "panel_flow_3d.png")

        ax_geo = fig.add_subplot(gs[row, 0])
        draw_geometry_panel(
            ax_geo,
            reference_scene,
            predicted_scene,
            input_img,
            title=cfg["geometry_title"],
            callout_text=cfg["geometry_callout"],
            summary=summary,
        )

        ax_flow = fig.add_subplot(gs[row, 1])
        draw_flow_panel(
            ax_flow,
            ref_flow,
            pred_flow,
            title=cfg["flow_title"],
            cfd_score=float(summary["cfd_summary"]["aggregate_score"]["cfd_score"]),
            callout_text=cfg["flow_callout"],
        )

    fig.text(
        0.285,
        0.955,
        "geometry overlap + input",
        ha="center",
        va="center",
        fontsize=10.0,
        fontweight="bold",
        color=COLORS["title"],
    )
    fig.text(
        0.748,
        0.955,
        "reference vs predicted flow",
        ha="center",
        va="center",
        fontsize=10.0,
        fontweight="bold",
        color=COLORS["title"],
    )
    fig.text(
        0.5,
        0.928,
        f"A4 mean  struct. {a4['mean_structural_score']:.3f}  CFD {a4['mean_cfd_score']:.3f}",
        ha="center",
        va="center",
        fontsize=9.4,
        color=COLORS["subtitle"],
    )

    # --- Legend below the figure ---
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_handles = [
        Line2D([0], [0], color=COLORS["reference_outline"], lw=1.5, linestyle=(0, (4, 2)), label="Reference room / obstacle"),
        Patch(facecolor=COLORS["obstacle_fill"], edgecolor=COLORS["obstacle_edge"], linewidth=1.0, label="Predicted obstacle"),
        Line2D([0], [0], color=COLORS["inlet"], lw=2.8, solid_capstyle="round", label="Inlet"),
        Line2D([0], [0], color=COLORS["outlet"], lw=2.8, solid_capstyle="round", label="Outlet"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        frameon=True,
        fontsize=9.0,
        edgecolor=COLORS["divider"],
        facecolor="white",
        bbox_to_anchor=(0.5, 0.005),
        handlelength=2.0,
        columnspacing=1.5,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT)
    fig.savefig(PNG_OUT, dpi=600)
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
