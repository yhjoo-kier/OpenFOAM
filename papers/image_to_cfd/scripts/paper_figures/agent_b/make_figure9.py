#!/usr/bin/env python3
"""Figure 9 – Obstacle hallucination with limited CFD penalty.

Fresh rebuild: maximally clean 2-row x 2-col (reference / prediction) layout.
No cue cards, no evidence lane, no summary strips.  The visual comparison
carries the story; the caption supplies the interpretation.

Each prediction panel encodes:
  - matched obstacle  -> light tan, thin edge
  - extra (hallucinated) obstacle -> strong orange, XX hatch, thick edge
  - GT obstacle outline -> dashed purple (only in prediction)
  - inlet / outlet -> blue / red thick lines

A slim metric annotation in each prediction panel shows the CFD score,
and a compact shared legend sits below.

Representative cases (locked):
  A3-01 / wireframe  -- 0 -> 3 obstacles, N/W preserved, CFD 0.604
  A3-03 / wireframe  -- 1 -> 3 obstacles, E/S preserved, CFD 0.600

Agent-B copy: outputs to results/paper_figures_agent_b/
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[3]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_agent_b"
PDF_OUT = OUT_DIR / "fig_discuss_obstacle_hallucination.pdf"
PNG_OUT = OUT_DIR / "fig_discuss_obstacle_hallucination.png"

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

CASES = [
    {
        "case": "bench_a3_01",
        "view": "wireframe",
        "label": "A3-01",
        "obstacle_gt": 0,
        "obstacle_pred": 3,
        "topology_note": "N/W",
    },
    {
        "case": "bench_a3_03",
        "view": "wireframe",
        "label": "A3-03",
        "obstacle_gt": 1,
        "obstacle_pred": 3,
        "topology_note": "E/S",
    },
]

COLORS = {
    "room_fill": "#EEF3F8",
    "room_edge": "#334155",
    "ref_obstacle_fill": "#C8CED8",
    "ref_obstacle_edge": "#5B6573",
    "matched_fill": "#F5E6CC",
    "matched_edge": "#8B6A2F",
    "extra_fill": "#E8863C",
    "extra_edge": "#8B3A12",
    "gt_outline": "#7C3AED",
    "inlet": "#2563EB",
    "outlet": "#DC2626",
    "subtitle": "#1E293B",
    "panel_frame": "#94A3B8",
    "metric_bg": "#F1F5F9",
    "metric_edge": "#CBD5E1",
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
    pred_scene: dict[str, Any], ref_scene: dict[str, Any],
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
    lw: float = 7.0,
) -> None:
    """Draw a thick colored line segment for an opening."""
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
    pad_factor: float = 0.10,
) -> None:
    x_max, y_max = panel_extent

    # Room blocks
    for b in room_blocks(scene):
        ax.add_patch(
            Rectangle(
                (b["x"], b["y"]), b["dx"], b["dy"],
                facecolor=COLORS["room_fill"],
                edgecolor=COLORS["room_edge"],
                linewidth=2.0, zorder=2,
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
                    linewidth=2.0, hatch="///", zorder=4,
                )
            )
    else:
        extras = find_extra_obstacles(scene, ref_scene) if ref_scene else set()

        # Predicted obstacles (GT outline removed -- reference panel provides
        # the ground-truth comparison directly via side-by-side layout)
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
                    linewidth=2.8 if is_extra else 1.6,
                    alpha=0.90 if is_extra else 0.65,
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
        draw_opening(ax, wall, u, du, x_max, y_max, color, lw=7.0)

    # Panel styling -- generous padding
    ax.set_xlim(-pad_factor * x_max, x_max * (1 + pad_factor))
    ax.set_ylim(-pad_factor * y_max, y_max * (1 + pad_factor))
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        title, loc="left", fontsize=14.0, fontweight="bold",
        color=COLORS["subtitle"], pad=8,
    )
    for sp in ax.spines.values():
        sp.set_linewidth(1.0)
        sp.set_color(COLORS["panel_frame"])


def add_metric_badge(
    ax: plt.Axes,
    cfd_score: float,
    n_extra: int,
    topology: str,
) -> None:
    """Add a compact metric badge inside the prediction panel (top-right)."""
    txt_obs = f"+{n_extra} hallucinated"
    txt_open = f"openings preserved ({topology})"
    txt_cfd = f"CFD score {cfd_score:.2f}"

    lines = [txt_obs, txt_open, txt_cfd]
    colors = [COLORS["extra_edge"], COLORS["ok_green"], COLORS["subtitle"]]
    sizes = [12.0, 11.0, 11.5]

    y_start = 0.97
    for i, (line, color, sz) in enumerate(zip(lines, colors, sizes)):
        y = y_start - i * 0.09
        t = ax.text(
            0.97, y, line,
            transform=ax.transAxes, fontsize=sz, fontweight="bold",
            color=color, ha="right", va="top", zorder=20,
        )
        t.set_path_effects([
            pe.Stroke(linewidth=4.5, foreground="white"),
            pe.Normal(),
        ])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font],
        "font.size": 12.5,
        "axes.titlesize": 14.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    # Figure: double-column, 2 rows x 2 cols + legend
    fig = plt.figure(figsize=(12.0, 10.2), constrained_layout=False)

    gs = fig.add_gridspec(
        2, 2,
        left=0.04, right=0.96, top=0.90, bottom=0.13,
        wspace=0.10, hspace=0.14,
        width_ratios=[1.0, 1.15],
        height_ratios=[0.85, 1.0],
    )

    # Column headers
    fig.text(
        0.26, 0.935, "reference", ha="center", fontsize=15.5,
        fontweight="bold", color=COLORS["subtitle"],
    )
    fig.text(
        0.73, 0.935, "prediction", ha="center",
        fontsize=15.5, fontweight="bold", color=COLORS["subtitle"],
    )

    panel_labels = [
        ["(a)", "(b)"],
        ["(c)", "(d)"],
    ]

    for row, cfg in enumerate(CASES):
        eval_dir = EVAL_ROOT / cfg["case"] / cfg["view"]
        task = load_json(eval_dir / "task.json")
        summary = load_json(eval_dir / "evaluation_summary.json")
        ref_scene = load_json(eval_dir / "reference_scene.json")
        pred_scene = load_json(eval_dir / "predicted_scene.json")

        # Shared extent for consistent scaling
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

        # Extra padding for row 2 (A3-03) to avoid E/S crowding
        pad = 0.18 if row == 1 else 0.09

        # Reference panel
        ax_ref = fig.add_subplot(gs[row, 0])
        draw_panel(
            ax_ref, ref_scene,
            title=f"{panel_labels[row][0]}  {cfg['label']}",
            panel_extent=pe_val,
            style="reference",
            pad_factor=pad,
        )

        # Prediction panel
        ax_pred = fig.add_subplot(gs[row, 1])
        draw_panel(
            ax_pred, pred_scene,
            title=f"{panel_labels[row][1]}  {cfg['label']}",
            panel_extent=pe_val,
            ref_scene=ref_scene,
            style="prediction",
            pad_factor=pad,
        )

        # Compact metric badge in prediction panel
        n_extra = cfg["obstacle_pred"] - cfg["obstacle_gt"]
        add_metric_badge(ax_pred, cfd_score, n_extra, cfg["topology_note"])

    # --- Legend ---
    legend_handles = [
        Patch(
            facecolor=COLORS["ref_obstacle_fill"],
            edgecolor=COLORS["ref_obstacle_edge"],
            linewidth=1.5, hatch="///",
            label="reference obstacle",
        ),
        Patch(
            facecolor=COLORS["matched_fill"],
            edgecolor=COLORS["matched_edge"],
            linewidth=1.2,
            label="matched prediction",
        ),
        Patch(
            facecolor=COLORS["extra_fill"],
            edgecolor=COLORS["extra_edge"],
            linewidth=2.0, hatch="xx",
            label="hallucinated (extra)",
        ),
        Line2D([0], [0], color=COLORS["inlet"], lw=6.0, label="inlet"),
        Line2D([0], [0], color=COLORS["outlet"], lw=6.0, label="outlet"),
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=5,
        bbox_to_anchor=(0.50, 0.015),
        frameon=True,
        fancybox=True,
        shadow=False,
        edgecolor=COLORS["metric_edge"],
        facecolor="white",
        fontsize=14.0,
        handlelength=3.5,
        handletextpad=0.9,
        columnspacing=2.5,
        borderaxespad=1.0,
    )

    # Save
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
