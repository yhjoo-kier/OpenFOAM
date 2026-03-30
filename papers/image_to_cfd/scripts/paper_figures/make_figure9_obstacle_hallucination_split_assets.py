#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from PIL import Image, ImageEnhance, ImageFilter, ImageOps

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
ASSET_DIR = PROJECT_ROOT / "results" / "paper_figures" / "figure9_split_assets"
MANIFEST_OUT = ASSET_DIR / "figure9_split_manifest.json"

CASES = [
    {
        "case": "bench_a3_01",
        "view": "wireframe",
        "label": "A3-01",
        "panel_prefix": "(a)",
        "obstacle_note": "0 -> 3",
        "topology_note": "N/W preserved",
        "detail_crop": None,
        "detail_note": None,
    },
    {
        "case": "bench_a3_03",
        "view": "wireframe",
        "label": "A3-03",
        "panel_prefix": "(d)",
        "obstacle_note": "1 -> 3",
        "topology_note": "E/S preserved",
        "detail_crop": (5.95, 8.35, 0.0, 2.20),
        "detail_note": "detail near E/S openings",
    },
]

COLORS = {
    "room_fill": "#EEF3F8",
    "room_edge": "#334155",
    "reference_obstacle_fill": "#D2D8E1",
    "reference_obstacle_edge": "#667085",
    "matched_obstacle_fill": "#EBCB9A",
    "matched_obstacle_edge": "#8B5A1F",
    "extra_obstacle_fill": "#E6A04A",
    "extra_obstacle_edge": "#9A3412",
    "reference_outline": "#7C3AED",
    "inlet": "#2563EB",
    "outlet": "#C83E3A",
    "subtitle": "#334155",
    "panel_frame": "#CBD5E1",
    "summary_box_face": "#F8FAFC",
    "summary_box_edge": "#94A3B8",
    "summary_ok": "#047857",
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


def draw_opening(ax: plt.Axes, wall: str, center: dict[str, float], size: dict[str, float], color: str, lw: float) -> None:
    line = None
    if wall in {"north", "south"}:
        x0 = center["u"] - size["du"] / 2.0
        x1 = center["u"] + size["du"] / 2.0
        y = 0.0 if wall == "south" else center["wall_y_max"]
        line = ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle="round", zorder=12)[0]
    elif wall in {"west", "east"}:
        y0 = center["u"] - size["du"] / 2.0
        y1 = center["u"] + size["du"] / 2.0
        x = 0.0 if wall == "west" else center["wall_x_max"]
        line = ax.plot([x, x], [y0, y1], color=color, lw=lw, solid_capstyle="round", zorder=12)[0]
    if line is not None:
        line.set_path_effects([pe.Stroke(linewidth=lw + 2.4, foreground="white"), pe.Normal()])


def load_enhanced_input(image_path: Path) -> np.ndarray:
    image = Image.open(image_path).convert("RGBA")
    array = np.asarray(image)
    alpha = array[..., 3] > 0
    visible = np.any(array[..., :3] < 249, axis=-1) & alpha
    if visible.any():
        ys, xs = np.where(visible)
        pad_x = max(5, int(0.015 * (xs.max() - xs.min() + 1)))
        pad_y = max(5, int(0.015 * (ys.max() - ys.min() + 1)))
        left = max(0, xs.min() - pad_x)
        right = min(image.width, xs.max() + pad_x + 1)
        top = max(0, ys.min() - pad_y)
        bottom = min(image.height, ys.max() + pad_y + 1)
        image = image.crop((left, top, right, bottom))
    image = ImageOps.expand(image, border=8, fill=(255, 255, 255, 255)).convert("L")
    image = ImageOps.autocontrast(image, cutoff=0.5)
    image = ImageEnhance.Contrast(image).enhance(2.1)
    image = image.resize((image.width * 2, image.height * 2), resample=Image.Resampling.LANCZOS)
    image = image.filter(ImageFilter.UnsharpMask(radius=1.2, percent=115, threshold=2))
    image = image.filter(ImageFilter.MinFilter(3))
    arr = np.asarray(image).astype(np.uint8)
    arr = np.where(arr < 232, 26, 255).astype(np.uint8)
    image = Image.fromarray(arr, mode="L").resize((image.width // 2, image.height // 2), resample=Image.Resampling.NEAREST)
    arr = np.asarray(image).astype(np.uint8)
    return np.stack([arr, arr, arr], axis=-1)


def find_extra_obstacles(predicted_scene: dict[str, Any], reference_scene: dict[str, Any]) -> list[dict[str, float]]:
    ref_boxes = []
    for obstacle in reference_scene.get("obstacles", []):
        x0 = float(obstacle["min"]["x"])
        y0 = float(obstacle["min"]["y"])
        dx = float(obstacle["size"]["dx"])
        dy = float(obstacle["size"]["dy"])
        ref_boxes.append((x0, y0, x0 + dx, y0 + dy))

    extras: list[dict[str, float]] = []
    for obstacle in predicted_scene.get("obstacles", []):
        x0 = float(obstacle["min"]["x"])
        y0 = float(obstacle["min"]["y"])
        dx = float(obstacle["size"]["dx"])
        dy = float(obstacle["size"]["dy"])
        box = (x0, y0, x0 + dx, y0 + dy)
        matched = False
        for rx0, ry0, rx1, ry1 in ref_boxes:
            inter_x = max(0.0, min(box[2], rx1) - max(box[0], rx0))
            inter_y = max(0.0, min(box[3], ry1) - max(box[1], ry0))
            inter = inter_x * inter_y
            union = dx * dy + (rx1 - rx0) * (ry1 - ry0) - inter
            if union > 0 and inter / union >= 0.18:
                matched = True
                break
        if not matched:
            extras.append({"x": x0, "y": y0, "dx": dx, "dy": dy})
    return extras


def draw_scene(
    ax: plt.Axes,
    scene: dict[str, Any],
    *,
    title: str,
    panel_extent: tuple[float, float],
    reference_scene: dict[str, Any] | None = None,
    style: str = "prediction",
    crop_variant: str = "default",
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
                linewidth=1.45,
                zorder=2,
            )
        )

    if style == "reference":
        for obstacle in scene.get("obstacles", []):
            ax.add_patch(
                Rectangle(
                    (float(obstacle["min"]["x"]), float(obstacle["min"]["y"])),
                    float(obstacle["size"]["dx"]),
                    float(obstacle["size"]["dy"]),
                    facecolor=COLORS["reference_obstacle_fill"],
                    edgecolor=COLORS["reference_obstacle_edge"],
                    linewidth=1.5,
                    hatch="///",
                    zorder=4,
                )
            )
    else:
        extra_boxes = {
            (extra["x"], extra["y"], extra["dx"], extra["dy"])
            for extra in find_extra_obstacles(scene, reference_scene or {"obstacles": []})
        }
        for obstacle in scene.get("obstacles", []):
            x0 = float(obstacle["min"]["x"])
            y0 = float(obstacle["min"]["y"])
            dx = float(obstacle["size"]["dx"])
            dy = float(obstacle["size"]["dy"])
            key = (x0, y0, dx, dy)
            is_extra = key in extra_boxes
            ax.add_patch(
                Rectangle(
                    (x0, y0),
                    dx,
                    dy,
                    facecolor=COLORS["extra_obstacle_fill"] if is_extra else COLORS["matched_obstacle_fill"],
                    edgecolor=COLORS["extra_obstacle_edge"] if is_extra else COLORS["matched_obstacle_edge"],
                    linewidth=1.95 if is_extra else 1.45,
                    alpha=0.88 if is_extra else 0.65,
                    hatch="xx" if is_extra else None,
                    zorder=5 if is_extra else 4,
                )
            )

    for opening in scene.get("openings", []):
        wall, center, size = prepare_opening_geometry(opening, x_max, y_max)
        color = COLORS["inlet"] if opening["type"] == "inlet" else COLORS["outlet"]
        draw_opening(ax, wall, center, size, color=color, lw=6.3 if style == "prediction" else 5.5)

    if style == "prediction":
        if crop_variant == "focus_right":
            ax.set_xlim(-0.03 * x_max, x_max * 1.16)
            ax.set_ylim(-0.12 * y_max, y_max * 1.02)
        else:
            ax.set_xlim(-0.03 * x_max, x_max * 1.08)
            ax.set_ylim(-0.10 * y_max, y_max * 1.04)
    else:
        ax.set_xlim(-0.04 * x_max, x_max * 1.05)
        ax.set_ylim(-0.06 * y_max, y_max * 1.03)
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, loc="left", fontsize=12.0, fontweight="bold", color=COLORS["subtitle"], pad=4)
    for spine in ax.spines.values():
        spine.set_linewidth(0.85)
        spine.set_color(COLORS["panel_frame"])


def apply_panel_crop(ax: plt.Axes, crop: tuple[float, float, float, float]) -> None:
    x0, x1, y0, y1 = crop
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)


def add_evidence_asset(ax: plt.Axes, *, label: str, obstacle_note: str, topology_note: str, cfd_score: float) -> None:
    ax.set_axis_off()
    ax.set_facecolor("white")
    ax.text(0.0, 0.965, f"{label} evidence", ha="left", va="bottom", fontsize=13.2, fontweight="bold", color=COLORS["subtitle"], transform=ax.transAxes)
    cards = [
        (0.69, "obstacles", obstacle_note, COLORS["extra_obstacle_edge"], 20.0),
        (0.39, "openings", topology_note, COLORS["summary_ok"], 18.2),
        (0.09, "CFD", f"{cfd_score:.2f}  moderate", COLORS["subtitle"], 17.2),
    ]
    for y, title, value, color, value_size in cards:
        box = Rectangle(
            (0.0, y),
            1.0,
            0.22,
            transform=ax.transAxes,
            facecolor=COLORS["summary_box_face"],
            edgecolor=COLORS["summary_box_edge"],
            linewidth=1.0,
        )
        ax.add_patch(box)
        ax.text(0.05, y + 0.16, title, ha="left", va="center", fontsize=11.4, fontweight="bold", color=COLORS["subtitle"], transform=ax.transAxes)
        ax.text(0.05, y + 0.075, value, ha="left", va="center", fontsize=value_size, fontweight="bold", color=color, transform=ax.transAxes)


def add_compact_preview_evidence(ax: plt.Axes, *, label: str, obstacle_note: str, topology_note: str, cfd_score: float) -> None:
    ax.set_axis_off()
    ax.set_facecolor("white")
    ax.text(0.0, 0.965, f"{label} evidence", ha="left", va="bottom", fontsize=12.6, fontweight="bold", color=COLORS["subtitle"], transform=ax.transAxes)
    box = Rectangle(
        (0.0, 0.06),
        1.0,
        0.80,
        transform=ax.transAxes,
        facecolor=COLORS["summary_box_face"],
        edgecolor=COLORS["summary_box_edge"],
        linewidth=1.0,
    )
    ax.add_patch(box)
    row_ys = [0.72, 0.49, 0.26]
    labels = ["obs.", "open.", "CFD"]
    values = [obstacle_note, topology_note, f"{cfd_score:.2f}  moderate"]
    colors = [COLORS["extra_obstacle_edge"], COLORS["summary_ok"], COLORS["subtitle"]]
    value_sizes = [17.5, 16.5, 15.6]
    for y, row_label, value, color, value_size in zip(row_ys, labels, values, colors, value_sizes):
        ax.text(0.05, y, row_label, ha="left", va="center", fontsize=11.0, fontweight="bold", color=COLORS["subtitle"], transform=ax.transAxes)
        ax.text(0.28, y, value, ha="left", va="center", fontsize=value_size, fontweight="bold", color=color, transform=ax.transAxes)


def save_dual(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def setup_rcparams() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 11.0,
            "axes.titlesize": 11.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def build_case_bundles() -> list[dict[str, Any]]:
    bundles: list[dict[str, Any]] = []
    for case_cfg in CASES:
        eval_dir = EVAL_ROOT / case_cfg["case"] / case_cfg["view"]
        task = load_json(eval_dir / "task.json")
        summary = load_json(eval_dir / "evaluation_summary.json")
        reference_scene = load_json(Path(task["reference_scene"]))
        predicted_scene = load_json(eval_dir / "predicted_scene.json")
        bundles.append(
            {
                "cfg": case_cfg,
                "task": task,
                "summary": summary,
                "reference_scene": reference_scene,
                "predicted_scene": predicted_scene,
                "panel_extent": tuple(
                    max(a, b)
                    for a, b in zip(room_extent(room_blocks(reference_scene)), room_extent(room_blocks(predicted_scene)))
                ),
                "input_image": Path(task["input_image"]),
                "cfd_score": float(summary["cfd_summary"]["aggregate_score"]["cfd_score"]),
            }
        )
    return bundles


def render_input_asset(bundle: dict[str, Any], stem_name: str) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 2.6), constrained_layout=False)
    ax.imshow(load_enhanced_input(bundle["input_image"]), interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        f"{bundle['cfg']['panel_prefix']} {bundle['cfg']['label']}",
        loc="left",
        fontsize=12.3,
        fontweight="bold",
        color=COLORS["subtitle"],
        pad=4,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.85)
        spine.set_color(COLORS["panel_frame"])
    save_dual(fig, ASSET_DIR / stem_name)


def render_reference_asset(bundle: dict[str, Any], stem_name: str, panel_label: str) -> None:
    fig, ax = plt.subplots(figsize=(4.5, 2.6), constrained_layout=False)
    draw_scene(
        ax,
        bundle["reference_scene"],
        title=panel_label,
        panel_extent=bundle["panel_extent"],
        style="reference",
    )
    save_dual(fig, ASSET_DIR / stem_name)


def render_prediction_asset(bundle: dict[str, Any], stem_name: str, panel_label: str) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 2.9), constrained_layout=False)
    draw_scene(
        ax,
        bundle["predicted_scene"],
        title=panel_label,
        panel_extent=bundle["panel_extent"],
        reference_scene=bundle["reference_scene"],
        style="prediction",
        crop_variant="focus_right" if bundle["cfg"]["detail_crop"] else "default",
    )
    save_dual(fig, ASSET_DIR / stem_name)


def render_evidence_asset(bundle: dict[str, Any], stem_name: str) -> None:
    fig, ax = plt.subplots(figsize=(4.7, 3.3), constrained_layout=False)
    add_evidence_asset(
        ax,
        label=bundle["cfg"]["label"],
        obstacle_note=bundle["cfg"]["obstacle_note"],
        topology_note=bundle["cfg"]["topology_note"],
        cfd_score=bundle["cfd_score"],
    )
    save_dual(fig, ASSET_DIR / stem_name)


def render_detail_pair_asset(bundle: dict[str, Any], stem_name: str) -> None:
    crop = bundle["cfg"]["detail_crop"]
    if crop is None:
        return
    fig = plt.figure(figsize=(5.9, 3.2), constrained_layout=False)
    gs = fig.add_gridspec(1, 2, left=0.06, right=0.98, bottom=0.12, top=0.82, wspace=0.10)
    ax_pred = fig.add_subplot(gs[0, 0])
    ax_gt = fig.add_subplot(gs[0, 1])

    draw_scene(
        ax_pred,
        bundle["predicted_scene"],
        title="(g) pred.",
        panel_extent=bundle["panel_extent"],
        reference_scene=bundle["reference_scene"],
        style="prediction",
    )
    apply_panel_crop(ax_pred, crop)
    for spine in ax_pred.spines.values():
        spine.set_color(COLORS["extra_obstacle_edge"])
        spine.set_linewidth(1.35)
    ax_pred.text(
        0.04,
        0.08,
        "extra pred.",
        ha="left",
        va="bottom",
        transform=ax_pred.transAxes,
        fontsize=9.6,
        fontweight="bold",
        color=COLORS["extra_obstacle_edge"],
        bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=COLORS["extra_obstacle_edge"], lw=0.8),
        zorder=20,
    )

    draw_scene(
        ax_gt,
        bundle["reference_scene"],
        title="GT",
        panel_extent=bundle["panel_extent"],
        style="reference",
    )
    apply_panel_crop(ax_gt, crop)
    for spine in ax_gt.spines.values():
        spine.set_color(COLORS["reference_outline"])
        spine.set_linewidth(1.35)
    ax_gt.text(
        0.04,
        0.08,
        "true obstacle",
        ha="left",
        va="bottom",
        transform=ax_gt.transAxes,
        fontsize=9.6,
        fontweight="bold",
        color=COLORS["reference_outline"],
        bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=COLORS["reference_outline"], lw=0.8),
        zorder=20,
    )
    fig.suptitle(bundle["cfg"]["detail_note"], y=0.97, fontsize=12.0, fontweight="bold", color=COLORS["subtitle"])
    save_dual(fig, ASSET_DIR / stem_name)


def render_legend_asset(stem_name: str) -> None:
    fig = plt.figure(figsize=(6.8, 1.2), constrained_layout=False)
    handles = [
        Patch(facecolor=COLORS["reference_obstacle_fill"], edgecolor=COLORS["reference_obstacle_edge"], hatch="///", label="ref. obstacle"),
        Patch(facecolor=COLORS["matched_obstacle_fill"], edgecolor=COLORS["matched_obstacle_edge"], label="matched pred."),
        Patch(facecolor=COLORS["extra_obstacle_fill"], edgecolor=COLORS["extra_obstacle_edge"], hatch="xx", label="extra pred."),
        Line2D([0], [0], color=COLORS["inlet"], lw=5.2, label="inlet"),
        Line2D([0], [0], color=COLORS["outlet"], lw=5.2, label="outlet"),
    ]
    fig.legend(
        handles=handles,
        loc="center",
        ncol=3,
        frameon=False,
        fontsize=12.0,
        handlelength=2.5,
        handletextpad=0.8,
        columnspacing=1.5,
    )
    save_dual(fig, ASSET_DIR / stem_name)


def render_geometry_grid(bundles: list[dict[str, Any]], stem_name: str) -> None:
    fig = plt.figure(figsize=(11.5, 5.7), constrained_layout=False)
    gs = fig.add_gridspec(2, 3, left=0.04, right=0.99, bottom=0.10, top=0.90, wspace=0.12, hspace=0.18, width_ratios=[1.0, 0.92, 1.48])
    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

    for row_idx, bundle in enumerate(bundles):
        ax_input = fig.add_subplot(gs[row_idx, 0])
        ax_input.imshow(load_enhanced_input(bundle["input_image"]), interpolation="nearest")
        ax_input.set_xticks([])
        ax_input.set_yticks([])
        ax_input.set_title(
            f"{panel_labels[row_idx * 3]} {bundle['cfg']['label']}",
            loc="left",
            fontsize=12.1,
            fontweight="bold",
            color=COLORS["subtitle"],
            pad=4,
        )
        for spine in ax_input.spines.values():
            spine.set_linewidth(0.85)
            spine.set_color(COLORS["panel_frame"])

        ax_ref = fig.add_subplot(gs[row_idx, 1])
        draw_scene(
            ax_ref,
            bundle["reference_scene"],
            title=panel_labels[row_idx * 3 + 1],
            panel_extent=bundle["panel_extent"],
            style="reference",
        )

        ax_pred = fig.add_subplot(gs[row_idx, 2])
        draw_scene(
            ax_pred,
            bundle["predicted_scene"],
            title=panel_labels[row_idx * 3 + 2],
            panel_extent=bundle["panel_extent"],
            reference_scene=bundle["reference_scene"],
            style="prediction",
            crop_variant="focus_right" if bundle["cfg"]["detail_crop"] else "default",
        )

    fig.text(0.133, 0.93, "input wireframe", ha="center", va="bottom", fontsize=12.8, fontweight="bold", color=COLORS["subtitle"])
    fig.text(0.432, 0.93, "reference", ha="center", va="bottom", fontsize=12.8, fontweight="bold", color=COLORS["subtitle"])
    fig.text(0.790, 0.93, "prediction", ha="center", va="bottom", fontsize=12.8, fontweight="bold", color=COLORS["subtitle"])
    save_dual(fig, ASSET_DIR / stem_name)


def render_split_preview(bundles: list[dict[str, Any]], stem_name: str) -> None:
    fig = plt.figure(figsize=(15.4, 8.8), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        4,
        left=0.040,
        right=0.970,
        bottom=0.158,
        top=0.930,
        wspace=0.11,
        hspace=0.18,
        width_ratios=[1.00, 0.90, 1.55, 1.88],
    )
    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

    for row_idx, bundle in enumerate(bundles):
        ax_input = fig.add_subplot(gs[row_idx, 0])
        ax_input.imshow(load_enhanced_input(bundle["input_image"]), interpolation="nearest")
        ax_input.set_xticks([])
        ax_input.set_yticks([])
        ax_input.set_title(
            f"{panel_labels[row_idx * 3]} {bundle['cfg']['label']}",
            loc="left",
            fontsize=12.2,
            fontweight="bold",
            color=COLORS["subtitle"],
            pad=4,
        )
        for spine in ax_input.spines.values():
            spine.set_linewidth(0.85)
            spine.set_color(COLORS["panel_frame"])

        ax_ref = fig.add_subplot(gs[row_idx, 1])
        draw_scene(
            ax_ref,
            bundle["reference_scene"],
            title=panel_labels[row_idx * 3 + 1],
            panel_extent=bundle["panel_extent"],
            style="reference",
        )

        ax_pred = fig.add_subplot(gs[row_idx, 2])
        draw_scene(
            ax_pred,
            bundle["predicted_scene"],
            title=panel_labels[row_idx * 3 + 2],
            panel_extent=bundle["panel_extent"],
            reference_scene=bundle["reference_scene"],
            style="prediction",
            crop_variant="focus_right" if bundle["cfg"]["detail_crop"] else "default",
        )

        if row_idx == 0:
            ax_evidence = fig.add_subplot(gs[row_idx, 3])
            add_evidence_asset(
                ax_evidence,
                label=bundle["cfg"]["label"],
                obstacle_note=bundle["cfg"]["obstacle_note"],
                topology_note=bundle["cfg"]["topology_note"],
                cfd_score=bundle["cfd_score"],
            )
        else:
            side_gs = gs[row_idx, 3].subgridspec(2, 1, height_ratios=[0.46, 0.54], hspace=0.12)
            ax_evidence = fig.add_subplot(side_gs[0, 0])
            add_compact_preview_evidence(
                ax_evidence,
                label=bundle["cfg"]["label"],
                obstacle_note=bundle["cfg"]["obstacle_note"],
                topology_note=bundle["cfg"]["topology_note"],
                cfd_score=bundle["cfd_score"],
            )
            detail_gs = side_gs[1, 0].subgridspec(1, 2, wspace=0.09)
            ax_pred_detail = fig.add_subplot(detail_gs[0, 0])
            ax_gt_detail = fig.add_subplot(detail_gs[0, 1])
            crop = bundle["cfg"]["detail_crop"]
            draw_scene(
                ax_pred_detail,
                bundle["predicted_scene"],
                title="(g) pred.",
                panel_extent=bundle["panel_extent"],
                reference_scene=bundle["reference_scene"],
                style="prediction",
            )
            apply_panel_crop(ax_pred_detail, crop)
            for spine in ax_pred_detail.spines.values():
                spine.set_color(COLORS["extra_obstacle_edge"])
                spine.set_linewidth(1.35)
            ax_pred_detail.text(
                0.04,
                0.08,
                "extra pred.",
                ha="left",
                va="bottom",
                transform=ax_pred_detail.transAxes,
                fontsize=9.4,
                fontweight="bold",
                color=COLORS["extra_obstacle_edge"],
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=COLORS["extra_obstacle_edge"], lw=0.8),
                zorder=20,
            )

            draw_scene(
                ax_gt_detail,
                bundle["reference_scene"],
                title="GT",
                panel_extent=bundle["panel_extent"],
                style="reference",
            )
            apply_panel_crop(ax_gt_detail, crop)
            for spine in ax_gt_detail.spines.values():
                spine.set_color(COLORS["reference_outline"])
                spine.set_linewidth(1.35)
            ax_gt_detail.text(
                0.04,
                0.08,
                "true obstacle",
                ha="left",
                va="bottom",
                transform=ax_gt_detail.transAxes,
                fontsize=9.4,
                fontweight="bold",
                color=COLORS["reference_outline"],
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=COLORS["reference_outline"], lw=0.8),
                zorder=20,
            )

    fig.text(0.125, 0.944, "input wireframe", ha="center", va="bottom", fontsize=12.8, fontweight="bold", color=COLORS["subtitle"])
    fig.text(0.370, 0.944, "reference", ha="center", va="bottom", fontsize=12.8, fontweight="bold", color=COLORS["subtitle"])
    fig.text(0.660, 0.944, "prediction", ha="center", va="bottom", fontsize=12.8, fontweight="bold", color=COLORS["subtitle"])
    fig.text(0.882, 0.944, "evidence / detail", ha="center", va="bottom", fontsize=12.8, fontweight="bold", color=COLORS["subtitle"])

    handles = [
        Patch(facecolor=COLORS["reference_obstacle_fill"], edgecolor=COLORS["reference_obstacle_edge"], hatch="///", label="ref. obstacle"),
        Patch(facecolor=COLORS["matched_obstacle_fill"], edgecolor=COLORS["matched_obstacle_edge"], label="matched pred."),
        Patch(facecolor=COLORS["extra_obstacle_fill"], edgecolor=COLORS["extra_obstacle_edge"], hatch="xx", label="extra pred."),
        Line2D([0], [0], color=COLORS["inlet"], lw=5.2, label="inlet"),
        Line2D([0], [0], color=COLORS["outlet"], lw=5.2, label="outlet"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.50, 0.040),
        frameon=False,
        fontsize=12.0,
        handlelength=2.5,
        handletextpad=0.8,
        columnspacing=1.5,
    )
    save_dual(fig, ASSET_DIR / stem_name)


def write_manifest(bundles: list[dict[str, Any]]) -> None:
    asset_specs = [
        {"name": "figure9_row1_input", "role": "panel asset", "message": "wireframe input for A3-01"},
        {"name": "figure9_row1_reference", "role": "panel asset", "message": "reference geometry for A3-01"},
        {"name": "figure9_row1_prediction", "role": "panel asset", "message": "prediction geometry for A3-01 without GT overlay"},
        {"name": "figure9_row1_evidence", "role": "summary asset", "message": "row-1 obstacle/opening/CFD summary"},
        {"name": "figure9_row2_input", "role": "panel asset", "message": "wireframe input for A3-03"},
        {"name": "figure9_row2_reference", "role": "panel asset", "message": "reference geometry for A3-03"},
        {"name": "figure9_row2_prediction", "role": "panel asset", "message": "prediction geometry for A3-03 without GT overlay"},
        {"name": "figure9_row2_evidence", "role": "summary asset", "message": "row-2 obstacle/opening/CFD summary"},
        {"name": "figure9_row2_detail_pair", "role": "detail asset", "message": "dominant row-2 blocker resolver: pred-vs-GT crop"},
        {"name": "figure9_legend", "role": "shared asset", "message": "shared legend for split-asset assembly"},
        {"name": "figure9_geometry_grid", "role": "assembly support", "message": "2x3 geometry-only board for panel-level QC"},
        {"name": "figure9_split_asset_preview", "role": "assembly preview", "message": "candidate rebuilt composition from split assets"},
    ]
    manifest = {
        "source_manifest": str(PROJECT_ROOT / "benchmark" / "manifests" / "evaluation_aggregate_summary_posthoc_scaled_longest_span.json"),
        "representative_cases": [
            {
                "case": bundle["cfg"]["case"],
                "view": bundle["cfg"]["view"],
                "label": bundle["cfg"]["label"],
                "obstacles": bundle["cfg"]["obstacle_note"],
                "openings": bundle["cfg"]["topology_note"],
                "cfd_score": bundle["cfd_score"],
            }
            for bundle in bundles
        ],
        "message_allocation": {
            "geometry_panels": "show obstacle hallucination and opening preservation without GT overlays in the main prediction panels",
            "summary_assets": "carry row-level obstacle/opening/CFD interpretation in standalone cards",
            "row2_detail_pair": "carry the crowded panel-(f) comparison as a direct pred-vs-GT crop instead of a stacked overlay",
            "legend": "only explain ref obstacle, matched pred., extra pred., inlet, outlet",
        },
        "qc_gate": [
            "local panel-level QC on row2 prediction and row2 detail pair",
            "local panel-level QC on row1/row2 evidence assets",
            "local assembled QC on geometry grid and split preview",
            "external subagent QC on one exact rebuilt preview revision",
            "external Gemini QC on the same exact rebuilt preview revision",
        ],
        "assets": [
            {
                **spec,
                "pdf": str((ASSET_DIR / spec["name"]).with_suffix(".pdf")),
                "png": str((ASSET_DIR / spec["name"]).with_suffix(".png")),
            }
            for spec in asset_specs
        ],
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    setup_rcparams()
    bundles = build_case_bundles()
    ASSET_DIR.mkdir(parents=True, exist_ok=True)

    render_input_asset(bundles[0], "figure9_row1_input")
    render_reference_asset(bundles[0], "figure9_row1_reference", "(b)")
    render_prediction_asset(bundles[0], "figure9_row1_prediction", "(c)")
    render_evidence_asset(bundles[0], "figure9_row1_evidence")

    render_input_asset(bundles[1], "figure9_row2_input")
    render_reference_asset(bundles[1], "figure9_row2_reference", "(e)")
    render_prediction_asset(bundles[1], "figure9_row2_prediction", "(f)")
    render_evidence_asset(bundles[1], "figure9_row2_evidence")
    render_detail_pair_asset(bundles[1], "figure9_row2_detail_pair")

    render_legend_asset("figure9_legend")
    render_geometry_grid(bundles, "figure9_geometry_grid")
    render_split_preview(bundles, "figure9_split_asset_preview")
    write_manifest(bundles)

    print(f"Wrote split assets to {ASSET_DIR}")
    print(f"Wrote manifest to {MANIFEST_OUT}")


if __name__ == "__main__":
    main()
