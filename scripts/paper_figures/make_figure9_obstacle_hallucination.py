#!/usr/bin/env python3
"""Figure 9: Obstacle hallucination with limited CFD penalty.

Split-asset rebuild.  Layout: 2x2 geometry panels (reference / prediction)
+ bottom summary strip.  Wireframe inputs are omitted from panels
(referenced in caption) to eliminate the persistent wireframe-vs-plan-view
balance instability and give each prediction panel ~2x the area of the
previous 4-column layout.

Representative cases (locked):
  - bench_a3_01 / wireframe  (0 -> 3 obstacles, CFD 0.604, N/W preserved)
  - bench_a3_03 / wireframe  (1 -> 3 obstacles, CFD 0.600, E/S preserved)
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

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures"
PDF_OUT = OUT_DIR / "figure9_obstacle_hallucination_limited_cfd_penalty.pdf"
PNG_OUT = OUT_DIR / "figure9_obstacle_hallucination_limited_cfd_penalty.png"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

CASES = [
    {
        "case": "bench_a3_01",
        "view": "wireframe",
        "label": "A3-01",
        "obstacle_gt": 0,
        "obstacle_pred": 3,
        "topology_note": "N / W preserved",
    },
    {
        "case": "bench_a3_03",
        "view": "wireframe",
        "label": "A3-03",
        "obstacle_gt": 1,
        "obstacle_pred": 3,
        "topology_note": "E / S preserved",
    },
]

COLORS = {
    "room_fill": "#EEF3F8",
    "room_edge": "#334155",
    "ref_obstacle_fill": "#D2D8E1",
    "ref_obstacle_edge": "#667085",
    "matched_fill": "#F5E6CC",
    "matched_edge": "#8B5A1F",
    "extra_fill": "#C45A20",
    "extra_edge": "#7C2D12",
    "gt_outline": "#7C3AED",
    "inlet": "#2563EB",
    "outlet": "#C83E3A",
    "subtitle": "#334155",
    "panel_frame": "#CBD5E1",
    "summary_bg": "#F8FAFC",
    "summary_edge": "#94A3B8",
    "ok_green": "#047857",
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


def find_extra_obstacles(
    pred_scene: dict[str, Any], ref_scene: dict[str, Any]
) -> set[tuple[float, float, float, float]]:
    ref_boxes: list[tuple[float, float, float, float]] = []
    for o in ref_scene.get("obstacles", []):
        x0 = float(o["min"]["x"])
        y0 = float(o["min"]["y"])
        dx = float(o["size"]["dx"])
        dy = float(o["size"]["dy"])
        ref_boxes.append((x0, y0, x0 + dx, y0 + dy))

    extras: set[tuple[float, float, float, float]] = set()
    for o in pred_scene.get("obstacles", []):
        x0 = float(o["min"]["x"])
        y0 = float(o["min"]["y"])
        dx = float(o["size"]["dx"])
        dy = float(o["size"]["dy"])
        matched = False
        for rx0, ry0, rx1, ry1 in ref_boxes:
            ix = max(0.0, min(x0 + dx, rx1) - max(x0, rx0))
            iy = max(0.0, min(y0 + dy, ry1) - max(y0, ry0))
            inter = ix * iy
            union = dx * dy + (rx1 - rx0) * (ry1 - ry0) - inter
            if union > 0 and inter / union >= 0.18:
                matched = True
                break
        if not matched:
            extras.add((x0, y0, dx, dy))
    return extras


def draw_opening(
    ax: plt.Axes,
    wall: str,
    u: float,
    du: float,
    x_max: float,
    y_max: float,
    color: str,
    lw: float = 6.0,
) -> None:
    line = None
    if wall in ("north", "south"):
        x0, x1 = u - du / 2, u + du / 2
        y = 0.0 if wall == "south" else y_max
        line = ax.plot(
            [x0, x1], [y, y], color=color, lw=lw,
            solid_capstyle="round", zorder=12,
        )[0]
    elif wall in ("west", "east"):
        y0, y1 = u - du / 2, u + du / 2
        x = 0.0 if wall == "west" else x_max
        line = ax.plot(
            [x, x], [y0, y1], color=color, lw=lw,
            solid_capstyle="round", zorder=12,
        )[0]
    if line is not None:
        line.set_path_effects(
            [pe.Stroke(linewidth=lw + 3.0, foreground="white"), pe.Normal()]
        )


# ---------------------------------------------------------------------------
# Panel drawing
# ---------------------------------------------------------------------------

def draw_panel(
    ax: plt.Axes,
    scene: dict[str, Any],
    *,
    title: str,
    panel_extent: tuple[float, float],
    ref_scene: dict[str, Any] | None = None,
    style: str = "reference",
) -> None:
    x_max, y_max = panel_extent

    # Room blocks
    for b in room_blocks(scene):
        ax.add_patch(
            Rectangle(
                (b["x"], b["y"]), b["dx"], b["dy"],
                facecolor=COLORS["room_fill"],
                edgecolor=COLORS["room_edge"],
                linewidth=1.8, zorder=2,
            )
        )

    if style == "reference":
        for o in scene.get("obstacles", []):
            ax.add_patch(
                Rectangle(
                    (float(o["min"]["x"]), float(o["min"]["y"])),
                    float(o["size"]["dx"]), float(o["size"]["dy"]),
                    facecolor=COLORS["ref_obstacle_fill"],
                    edgecolor=COLORS["ref_obstacle_edge"],
                    linewidth=1.8, hatch="///", zorder=4,
                )
            )
    else:
        extras = find_extra_obstacles(scene, ref_scene) if ref_scene else set()

        # GT obstacle outlines (from reference)
        if ref_scene:
            for o in ref_scene.get("obstacles", []):
                x0 = float(o["min"]["x"])
                y0 = float(o["min"]["y"])
                dx = float(o["size"]["dx"])
                dy = float(o["size"]["dy"])
                p = Rectangle(
                    (x0, y0), dx, dy,
                    facecolor=(0.93, 0.88, 1.0, 0.12),
                    edgecolor=COLORS["gt_outline"],
                    linewidth=3.0,
                    linestyle=(0, (5.5, 2.0)),
                    zorder=7,
                )
                p.set_path_effects(
                    [pe.Stroke(linewidth=4.5, foreground="white"), pe.Normal()]
                )
                ax.add_patch(p)

        # Predicted obstacles
        for o in scene.get("obstacles", []):
            x0 = float(o["min"]["x"])
            y0 = float(o["min"]["y"])
            dx = float(o["size"]["dx"])
            dy = float(o["size"]["dy"])
            is_extra = (x0, y0, dx, dy) in extras
            ax.add_patch(
                Rectangle(
                    (x0, y0), dx, dy,
                    facecolor=COLORS["extra_fill"] if is_extra else COLORS["matched_fill"],
                    edgecolor=COLORS["extra_edge"] if is_extra else COLORS["matched_edge"],
                    linewidth=2.4 if is_extra else 1.4,
                    alpha=0.92 if is_extra else 0.45,
                    hatch="xx" if is_extra else None,
                    zorder=5 if is_extra else 4,
                )
            )

    # Openings
    for op in scene.get("openings", []):
        wall = op["wall"]
        u = float(op["center"]["u"])
        du = float(op["size"]["du"])
        color = COLORS["inlet"] if op["type"] == "inlet" else COLORS["outlet"]
        draw_opening(ax, wall, u, du, x_max, y_max, color, lw=6.5)

    # Panel styling — generous padding
    pad = 0.07
    ax.set_xlim(-pad * x_max, x_max * (1 + pad))
    ax.set_ylim(-pad * y_max, y_max * (1 + pad))
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        title, loc="left", fontsize=13.5, fontweight="bold",
        color=COLORS["subtitle"], pad=6,
    )
    for sp in ax.spines.values():
        sp.set_linewidth(1.0)
        sp.set_color(COLORS["panel_frame"])


# ---------------------------------------------------------------------------
# Summary strip
# ---------------------------------------------------------------------------

def draw_summary_box(
    ax: plt.Axes,
    label: str,
    obstacle_gt: int,
    obstacle_pred: int,
    topology_note: str,
    cfd_score: float,
) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_facecolor(COLORS["summary_bg"])
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(1.1)
        sp.set_color(COLORS["summary_edge"])

    # Case label (left edge)
    ax.text(
        0.03, 0.50, label, fontsize=16.5, fontweight="bold",
        color=COLORS["subtitle"], va="center", rotation=90,
    )

    # Three metric blocks evenly spaced
    col_xs = [0.18, 0.50, 0.82]

    # Obstacle count
    ax.text(
        col_xs[0], 0.76, "obstacles", fontsize=12.0, fontweight="bold",
        color=COLORS["subtitle"], va="center", ha="center",
    )
    ax.text(
        col_xs[0], 0.34, f"{obstacle_gt} \u2192 {obstacle_pred}",
        fontsize=18.5, fontweight="bold", color=COLORS["extra_edge"],
        va="center", ha="center",
    )

    # Opening topology
    topo_value = topology_note.replace(" / ", "/").replace(" preserved", "\npreserved")
    ax.text(
        col_xs[1], 0.76, "openings", fontsize=12.0, fontweight="bold",
        color=COLORS["subtitle"], va="center", ha="center",
    )
    ax.text(
        col_xs[1], 0.34, topo_value,
        fontsize=13.0, fontweight="bold", color=COLORS["ok_green"],
        va="center", ha="center", linespacing=1.0,
    )

    # CFD cue with explicit interpretation
    ax.text(
        col_xs[2], 0.76, "CFD response", fontsize=12.0, fontweight="bold",
        color=COLORS["subtitle"], va="center", ha="center",
    )
    ax.text(
        col_xs[2], 0.38, "moderate", fontsize=16.5, fontweight="bold",
        color=COLORS["subtitle"], va="center", ha="center",
    )
    ax.text(
        col_xs[2], 0.15, f"score {cfd_score:.2f}",
        fontsize=11.2, color=COLORS["subtitle"], va="center", ha="center",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font],
        "font.size": 12.0,
        "axes.titlesize": 13.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(11.6, 9.3), constrained_layout=False)

    # --- Main 2x2 geometry grid ---
    gs_main = fig.add_gridspec(
        2, 2,
        left=0.05, right=0.95, top=0.90, bottom=0.34,
        wspace=0.14, hspace=0.20,
        width_ratios=[1.0, 1.15],
    )

    # Column headers
    fig.text(
        0.28, 0.935, "reference", ha="center", fontsize=14.5,
        fontweight="bold", color=COLORS["subtitle"],
    )
    fig.text(
        0.72, 0.935, "prediction", ha="center",
        fontsize=15.0, fontweight="bold", color=COLORS["subtitle"],
    )

    panel_labels = ["(a)", "(b)", "(c)", "(d)"]
    cfd_scores: list[float] = []

    for row, cfg in enumerate(CASES):
        eval_dir = EVAL_ROOT / cfg["case"] / cfg["view"]
        task = load_json(eval_dir / "task.json")
        summary = load_json(eval_dir / "evaluation_summary.json")
        ref_scene = load_json(Path(task["reference_scene"]))
        pred_scene = load_json(eval_dir / "predicted_scene.json")
        pe_val = tuple(
            max(a, b)
            for a, b in zip(
                room_extent(room_blocks(ref_scene)),
                room_extent(room_blocks(pred_scene)),
            )
        )
        cfd_score = float(
            summary["cfd_summary"]["aggregate_score"]["cfd_score"]
        )
        cfd_scores.append(cfd_score)

        # Reference panel
        ax_ref = fig.add_subplot(gs_main[row, 0])
        draw_panel(
            ax_ref, ref_scene,
            title=f"{panel_labels[row * 2]}  {cfg['label']}",
            panel_extent=pe_val,
            style="reference",
        )

        # Prediction panel
        ax_pred = fig.add_subplot(gs_main[row, 1])
        draw_panel(
            ax_pred, pred_scene,
            title=f"{panel_labels[row * 2 + 1]}  {cfg['label']}",
            panel_extent=pe_val,
            ref_scene=ref_scene,
            style="prediction",
        )

    # --- Bottom summary strip ---
    gs_summary = fig.add_gridspec(
        1, 2,
        left=0.05, right=0.95, top=0.275, bottom=0.125,
        wspace=0.14,
    )

    for i, cfg in enumerate(CASES):
        ax_sum = fig.add_subplot(gs_summary[0, i])
        draw_summary_box(
            ax_sum,
            label=cfg["label"],
            obstacle_gt=cfg["obstacle_gt"],
            obstacle_pred=cfg["obstacle_pred"],
            topology_note=cfg["topology_note"],
            cfd_score=cfd_scores[i],
        )

    # --- Legend ---
    legend_handles = [
        Patch(
            facecolor=COLORS["ref_obstacle_fill"],
            edgecolor=COLORS["ref_obstacle_edge"],
            hatch="///", label="ref. obstacle",
        ),
        Patch(
            facecolor=COLORS["matched_fill"],
            edgecolor=COLORS["matched_edge"],
            label="matched pred.",
        ),
        Patch(
            facecolor=COLORS["extra_fill"],
            edgecolor=COLORS["extra_edge"],
            hatch="xx", label="extra pred. (hallucinated)",
        ),
        Rectangle(
            (0, 0), 1, 1,
            facecolor="none", edgecolor=COLORS["gt_outline"],
            linewidth=3.0, linestyle=(0, (8, 3)),
            label="GT obstacle outline",
        ),
        Line2D([0], [0], color=COLORS["inlet"], lw=5.5, label="inlet"),
        Line2D([0], [0], color=COLORS["outlet"], lw=5.5, label="outlet"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.50, 0.020),
        frameon=False,
        fontsize=13.8,
        handlelength=3.2,
        handletextpad=1.0,
        columnspacing=2.1,
        borderaxespad=0.8,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
