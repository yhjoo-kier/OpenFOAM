#!/usr/bin/env python3
"""Fig 14: Application to real architectural floor plans — 2 demo cases."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.font_manager as fm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_demo_floorplan_application.pdf"
PNG_OUT = OUT_DIR / "fig_demo_floorplan_application.png"

DEMO_DIR = PROJECT_ROOT / "benchmark" / "real_image_demo"
FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]

CASES = [
    {
        "input_img": DEMO_DIR / "case1_simple_rectangular.png",
        "scene_json": PROJECT_ROOT / "generated" / "demo_case1_rectangular.json",
        "case_dir": PROJECT_ROOT / "cases" / "demo_case1_rectangular",
        "label_row": "Case 1: Rectangular office",
    },
    {
        "input_img": DEMO_DIR / "case2_lshaped_composite.png",
        "scene_json": PROJECT_ROOT / "generated" / "demo_case2_lshaped.json",
        "case_dir": PROJECT_ROOT / "cases" / "demo_case2_lshaped",
        "label_row": "Case 2: L-shaped composite",
    },
]

COLORS = {
    "room_fill": "#EEF3F8",
    "room_edge": "#334155",
    "obstacle_fill": "#F5E6CC",
    "obstacle_edge": "#8B6A2F",
    "inlet": "#2563EB",
    "outlet": "#DC2626",
}


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def load_json(p: Path) -> Any:
    return json.loads(p.read_text(encoding="utf-8"))


def room_blocks(scene: dict) -> list[dict]:
    room = scene["room"]
    if "blocks" in room:
        return [{"x": float(b["origin"]["x"]), "y": float(b["origin"]["y"]),
                 "dx": float(b["size"]["dx"]), "dy": float(b["size"]["dy"])} for b in room["blocks"]]
    s = room["size"]
    return [{"x": 0.0, "y": 0.0, "dx": float(s["Lx"]), "dy": float(s["Ly"])}]


def obstacles(scene: dict) -> list[dict]:
    result = []
    for o in scene.get("obstacles", []):
        orig = o.get("origin") or o.get("min") or {"x": 0, "y": 0}
        size = o.get("size", {})
        result.append({"x": float(orig["x"]), "y": float(orig["y"]),
                       "dx": float(size.get("dx", 0)), "dy": float(size.get("dy", 0))})
    return result


def openings(scene: dict) -> list[dict]:
    result = []
    for o in scene.get("openings", []):
        result.append({
            "wall": o.get("wall", ""),
            "type": o.get("type", ""),
            "center": o.get("center", {}),
            "width": float(o.get("width", 0.5)),
        })
    return result


def draw_scene(ax, scene: dict, title: str):
    blocks = room_blocks(scene)
    obs = obstacles(scene)
    opens = openings(scene)

    for b in blocks:
        ax.add_patch(Rectangle((b["x"], b["y"]), b["dx"], b["dy"],
                                facecolor=COLORS["room_fill"], edgecolor=COLORS["room_edge"],
                                linewidth=1.5, zorder=1))

    for o in obs:
        ax.add_patch(Rectangle((o["x"], o["y"]), o["dx"], o["dy"],
                                facecolor=COLORS["obstacle_fill"], edgecolor=COLORS["obstacle_edge"],
                                linewidth=0.8, zorder=2))

    # Draw openings as colored lines on walls
    for op in opens:
        color = COLORS["inlet"] if op["type"] == "inlet" else COLORS["outlet"]
        wall = op["wall"]
        center = op.get("center", {})
        w = op["width"]
        cx = float(center.get("x", 0))
        cy = float(center.get("y", 0))
        cz = float(center.get("z", 0))

        if wall in ("north", "south"):
            ax.plot([cx - w / 2, cx + w / 2], [cy, cy], color=color, linewidth=3, zorder=5)
        elif wall in ("east", "west"):
            ax.plot([cx, cx], [cy - w / 2, cy + w / 2], color=color, linewidth=3, zorder=5)

    all_blocks = blocks
    x_max = max(b["x"] + b["dx"] for b in all_blocks) if all_blocks else 10
    y_max = max(b["y"] + b["dy"] for b in all_blocks) if all_blocks else 10
    margin = max(x_max, y_max) * 0.1
    ax.set_xlim(-margin, x_max + margin)
    ax.set_ylim(-margin, y_max + margin)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]", fontsize=8)
    ax.set_ylabel("y [m]", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=9, pad=4)


def main() -> None:
    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))

    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for row, case_cfg in enumerate(CASES):
        # Left: input image
        ax_img = axes[row, 0]
        if case_cfg["input_img"].exists():
            img = mpimg.imread(str(case_cfg["input_img"]))
            ax_img.imshow(img)
        ax_img.set_title(f"Input floor plan", fontsize=9, pad=4)
        ax_img.axis("off")

        # Right: VLM-extracted geometry
        ax_geo = axes[row, 1]
        if case_cfg["scene_json"].exists():
            scene = load_json(case_cfg["scene_json"])
            draw_scene(ax_geo, scene, "VLM-extracted geometry")
        else:
            ax_geo.text(0.5, 0.5, "Scene JSON\nnot found", transform=ax_geo.transAxes,
                        ha="center", fontsize=10, color="red")
            ax_geo.axis("off")

        # Row label
        axes[row, 0].text(-0.15, 0.5, case_cfg["label_row"], transform=axes[row, 0].transAxes,
                          fontsize=9.5, fontweight="bold", color="#334155",
                          ha="center", va="center", rotation=90)

        # Panel labels
        axes[row, 0].text(-0.05, 1.05, panel_labels[row * 2], transform=axes[row, 0].transAxes,
                          fontsize=11, fontweight="bold", color="#34495E")
        axes[row, 1].text(-0.12, 1.05, panel_labels[row * 2 + 1], transform=axes[row, 1].transAxes,
                          fontsize=11, fontweight="bold", color="#34495E")

    fig.tight_layout(rect=[0.05, 0.0, 1.0, 1.0], h_pad=2.0, w_pad=1.5)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02, bbox_inches="tight")
    print(f"Font: {font}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
