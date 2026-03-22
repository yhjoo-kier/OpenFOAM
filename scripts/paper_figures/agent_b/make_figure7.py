#!/usr/bin/env python3
"""Figure 7 – Representative cross-view outcome panel.

Compares prediction quality across 3 key view types (floorplan, wireframe,
section) for two representative cases (A2-03 and A3-03).
Layout: 2 rows (cases) x 3 columns (views).

Each cell shows reference geometry (dashed purple outlines) overlaid with
prediction (coloured fill: matched tan, hallucinated orange-hatched).
Structural and CFD scores annotated per panel.

Output: figure7_cross_view_outcome.pdf / .png (600 dpi)
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_result_crossview_outcome.pdf"
PNG_OUT = OUT_DIR / "fig_result_crossview_outcome.png"

# ---------------------------------------------------------------------------
# Configuration – reduced to 3 most informative views
# ---------------------------------------------------------------------------
CASES = [
    {"case": "bench_a2_03", "label": "A2-03"},
    {"case": "bench_a3_03", "label": "A3-03"},
]
VIEWS = ["floorplan", "birdseye", "section"]
VIEW_LABELS = {
    "floorplan": "Floorplan",
    "birdseye": "Bird's-eye",
    "section": "Section",
}

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

COLORS = {
    # Room
    "room_fill": "#EEF3F8",
    "room_edge": "#334155",
    # Reference geometry (drawn as outlines / hatched)
    "ref_obstacle_fill": "none",
    "ref_obstacle_edge": "#7C3AED",
    # Predicted geometry
    "pred_obstacle_fill": "#F5E6CC",
    "pred_obstacle_edge": "#8B6A2F",
    "extra_obstacle_fill": "#E8863C",
    "extra_obstacle_edge": "#8B3A12",
    "missing_obstacle_fill": "none",
    "missing_obstacle_edge": "#7C3AED",
    # Openings
    "ref_inlet": "#93C5FD",
    "ref_outlet": "#FCA5A5",
    "pred_inlet": "#2563EB",
    "pred_outlet": "#DC2626",
    # Text / chrome
    "subtitle": "#1E293B",
    "panel_frame": "#94A3B8",
    "score_good": "#047857",
    "score_mid": "#B45309",
    "score_bad": "#B91C1C",
    "row_label_bg": "#F1F5F9",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "DejaVu Sans"


def load_json(p: Path) -> Any:
    return json.loads(p.read_text(encoding="utf-8"))


def room_blocks(scene: dict[str, Any]) -> list[dict[str, float]]:
    room = scene["room"]
    if "blocks" in room:
        return [
            {
                "x": float(b["origin"]["x"]),
                "y": float(b["origin"]["y"]),
                "dx": float(b["size"]["dx"]),
                "dy": float(b["size"]["dy"]),
            }
            for b in room["blocks"]
        ]
    s = room["size"]
    return [{"x": 0.0, "y": 0.0, "dx": float(s["Lx"]), "dy": float(s["Ly"])}]


def room_extent(blocks: list[dict[str, float]]) -> tuple[float, float]:
    return (
        max(b["x"] + b["dx"] for b in blocks),
        max(b["y"] + b["dy"] for b in blocks),
    )


def obstacle_boxes(scene: dict[str, Any]) -> list[tuple[float, float, float, float]]:
    result = []
    for o in scene.get("obstacles", []):
        x0 = float(o["min"]["x"])
        y0 = float(o["min"]["y"])
        dx = float(o["size"]["dx"])
        dy = float(o["size"]["dy"])
        result.append((x0, y0, dx, dy))
    return result


def classify_obstacles(
    pred_scene: dict[str, Any], ref_scene: dict[str, Any],
) -> tuple[
    list[tuple[float, float, float, float]],  # matched
    list[tuple[float, float, float, float]],  # extra (hallucinated)
    list[tuple[float, float, float, float]],  # missing (unmatched ref)
]:
    ref_boxes = obstacle_boxes(ref_scene)
    pred_boxes = obstacle_boxes(pred_scene)

    matched_pred: list[tuple[float, float, float, float]] = []
    extra_pred: list[tuple[float, float, float, float]] = []
    ref_matched_idx: set[int] = set()

    for pb in pred_boxes:
        px0, py0, pdx, pdy = pb
        best_iou = 0.0
        best_idx = -1
        for ri, rb in enumerate(ref_boxes):
            rx0, ry0, rdx, rdy = rb
            ix = max(0.0, min(px0 + pdx, rx0 + rdx) - max(px0, rx0))
            iy = max(0.0, min(py0 + pdy, ry0 + rdy) - max(py0, ry0))
            inter = ix * iy
            union = pdx * pdy + rdx * rdy - inter
            iou = inter / union if union > 0 else 0.0
            if iou > best_iou:
                best_iou = iou
                best_idx = ri
        if best_iou >= 0.15 and best_idx not in ref_matched_idx:
            matched_pred.append(pb)
            ref_matched_idx.add(best_idx)
        else:
            extra_pred.append(pb)

    missing_ref = [rb for ri, rb in enumerate(ref_boxes) if ri not in ref_matched_idx]
    return matched_pred, extra_pred, missing_ref


def draw_opening(
    ax: plt.Axes,
    wall: str,
    u: float,
    du: float,
    x_max: float,
    y_max: float,
    color: str,
    lw: float = 7.0,
) -> None:
    if wall in ("north", "south"):
        x0, x1 = u - du / 2, u + du / 2
        y = 0.0 if wall == "south" else y_max
        line = ax.plot(
            [x0, x1], [y, y], color=color, lw=lw,
            solid_capstyle="butt", zorder=12,
        )[0]
    elif wall in ("west", "east"):
        y0, y1 = u - du / 2, u + du / 2
        x = 0.0 if wall == "west" else x_max
        line = ax.plot(
            [x, x], [y0, y1], color=color, lw=lw,
            solid_capstyle="butt", zorder=12,
        )[0]
    else:
        return
    line.set_path_effects(
        [pe.Stroke(linewidth=lw + 3.0, foreground="white"), pe.Normal()]
    )


def score_color(val: float) -> str:
    if val >= 0.8:
        return COLORS["score_good"]
    elif val >= 0.5:
        return COLORS["score_mid"]
    return COLORS["score_bad"]


# ---------------------------------------------------------------------------
# Panel drawing – overlay style (ref outlines + pred fills)
# ---------------------------------------------------------------------------

def draw_overlay_panel(
    ax: plt.Axes,
    ref_scene: dict[str, Any],
    pred_scene: dict[str, Any],
    panel_extent: tuple[float, float],
    pad_factor: float = 0.10,
) -> None:
    x_max, y_max = panel_extent

    # Room blocks (use predicted room shape)
    for b in room_blocks(pred_scene):
        ax.add_patch(
            Rectangle(
                (b["x"], b["y"]), b["dx"], b["dy"],
                facecolor=COLORS["room_fill"],
                edgecolor=COLORS["room_edge"],
                linewidth=1.8, zorder=2,
            )
        )

    # Classify obstacles
    matched, extra, missing = classify_obstacles(pred_scene, ref_scene)

    # Draw matched predicted obstacles
    for x0, y0, dx, dy in matched:
        ax.add_patch(
            Rectangle(
                (x0, y0), dx, dy,
                facecolor=COLORS["pred_obstacle_fill"],
                edgecolor=COLORS["pred_obstacle_edge"],
                linewidth=1.8, alpha=0.80, zorder=4,
            )
        )

    # Draw extra (hallucinated) obstacles
    for x0, y0, dx, dy in extra:
        ax.add_patch(
            Rectangle(
                (x0, y0), dx, dy,
                facecolor=COLORS["extra_obstacle_fill"],
                edgecolor=COLORS["extra_obstacle_edge"],
                linewidth=2.2, alpha=0.85, hatch="xx", zorder=5,
            )
        )

    # Draw missing reference obstacles as dashed outlines
    for x0, y0, dx, dy in missing:
        ax.add_patch(
            Rectangle(
                (x0, y0), dx, dy,
                facecolor="none",
                edgecolor=COLORS["missing_obstacle_edge"],
                linewidth=3.5, linestyle="--", zorder=6,
            )
        )

    # Draw reference obstacle outlines (dashed purple) – thicker for visibility
    for o in ref_scene.get("obstacles", []):
        ax.add_patch(
            Rectangle(
                (float(o["min"]["x"]), float(o["min"]["y"])),
                float(o["size"]["dx"]), float(o["size"]["dy"]),
                facecolor="none",
                edgecolor=COLORS["ref_obstacle_edge"],
                linewidth=3.5, linestyle="--", zorder=7,
            )
        )

    # Openings – reference (lighter, thinner) then predicted (solid, thick)
    for op in ref_scene.get("openings", []):
        wall = op["wall"]
        u = float(op["center"]["u"])
        du = float(op["size"]["du"])
        ref_x, ref_y = room_extent(room_blocks(ref_scene))
        color = COLORS["ref_inlet"] if op["type"] == "inlet" else COLORS["ref_outlet"]
        draw_opening(ax, wall, u, du, ref_x, ref_y, color, lw=4.0)

    for op in pred_scene.get("openings", []):
        wall = op["wall"]
        u = float(op["center"]["u"])
        du = float(op["size"]["du"])
        color = COLORS["pred_inlet"] if op["type"] == "inlet" else COLORS["pred_outlet"]
        draw_opening(ax, wall, u, du, x_max, y_max, color, lw=7.0)

    # Axis styling
    ax.set_xlim(-pad_factor * x_max, x_max * (1 + pad_factor))
    ax.set_ylim(-pad_factor * y_max, y_max * (1 + pad_factor))
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(1.0)
        sp.set_color(COLORS["panel_frame"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font],
        "font.size": 11.0,
        "axes.titlesize": 12.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    n_rows = len(CASES)
    n_cols = len(VIEWS)

    # Larger figure for 2x3 layout
    fig_w = 7.2
    fig_h = 6.0
    fig = plt.figure(figsize=(fig_w, fig_h))

    # Grid: 2 rows x 3 cols, room for row labels on left and legend below
    gs = fig.add_gridspec(
        n_rows, n_cols,
        left=0.08, right=0.995, top=0.90, bottom=0.13,
        wspace=0.10, hspace=0.20,
    )

    # Column headers (view type names)
    for ci, view in enumerate(VIEWS):
        pos = gs[0, ci].get_position(fig)
        x_center = pos.x0 + pos.width / 2
        fig.text(
            x_center, 0.94, VIEW_LABELS[view],
            ha="center", va="bottom", fontsize=13, fontweight="bold",
            color=COLORS["subtitle"],
        )

    panel_idx = 0
    panel_letters = "abcdef"

    for ri, case_cfg in enumerate(CASES):
        case_name = case_cfg["case"]
        case_label = case_cfg["label"]

        # Row label on the left
        pos0 = gs[ri, 0].get_position(fig)
        y_center = pos0.y0 + pos0.height / 2
        fig.text(
            0.025, y_center, case_label,
            ha="center", va="center", fontsize=13, fontweight="bold",
            color=COLORS["subtitle"], rotation=90,
        )

        for ci, view in enumerate(VIEWS):
            eval_dir = EVAL_ROOT / case_name / view
            ref_scene = load_json(eval_dir / "reference_scene.json")
            pred_scene = load_json(eval_dir / "predicted_scene.json")
            summary = load_json(eval_dir / "evaluation_summary.json")

            # Compute shared extent
            pe_ref = room_extent(room_blocks(ref_scene))
            pe_pred = room_extent(room_blocks(pred_scene))
            ext = (max(pe_ref[0], pe_pred[0]), max(pe_ref[1], pe_pred[1]))

            # Scores
            cfd_score = float(
                summary["cfd_summary"]["aggregate_score"]["cfd_score"]
            )
            struct_score = float(
                summary["prediction_summary"]["structural_score"]
            )

            ax = fig.add_subplot(gs[ri, ci])
            draw_overlay_panel(ax, ref_scene, pred_scene, ext, pad_factor=0.10)

            # Panel label – larger font
            letter = panel_letters[panel_idx]
            panel_idx += 1
            ax.text(
                -0.02, 1.06, f"({letter})",
                transform=ax.transAxes, fontsize=14, fontweight="bold",
                color=COLORS["subtitle"], ha="left", va="bottom",
            )

            # Score badge – larger, with clear white background
            s_col = score_color(struct_score)
            c_col = score_color(cfd_score)
            score_txt = f"S={struct_score:.2f}  C={cfd_score:.2f}"
            ax.text(
                0.03, 0.04, score_txt,
                transform=ax.transAxes, fontsize=11, fontweight="bold",
                color="#334155", ha="left", va="bottom", zorder=20,
                bbox=dict(
                    facecolor="white", edgecolor="#94A3B8",
                    alpha=0.92, pad=2.5, boxstyle="round,pad=0.35",
                ),
            )

    # --- Legend ---
    legend_handles = [
        Patch(
            facecolor=COLORS["room_fill"], edgecolor=COLORS["room_edge"],
            linewidth=1.4, label="Room",
        ),
        Patch(
            facecolor="none", edgecolor=COLORS["ref_obstacle_edge"],
            linewidth=3.5, linestyle="--", label="Ref. obstacle",
        ),
        Patch(
            facecolor=COLORS["pred_obstacle_fill"],
            edgecolor=COLORS["pred_obstacle_edge"],
            linewidth=1.4, label="Matched pred.",
        ),
        Patch(
            facecolor=COLORS["extra_obstacle_fill"],
            edgecolor=COLORS["extra_obstacle_edge"],
            linewidth=1.8, hatch="xx", label="Hallucinated",
        ),
        Line2D([0], [0], color=COLORS["pred_inlet"], lw=5.0, label="Inlet"),
        Line2D([0], [0], color=COLORS["pred_outlet"], lw=5.0, label="Outlet"),
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=6,
        bbox_to_anchor=(0.53, 0.005),
        frameon=True,
        fancybox=True,
        shadow=False,
        edgecolor="#CBD5E1",
        facecolor="white",
        fontsize=11,
        handlelength=2.8,
        handletextpad=0.7,
        columnspacing=1.8,
        borderaxespad=0.8,
    )

    # Save
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
