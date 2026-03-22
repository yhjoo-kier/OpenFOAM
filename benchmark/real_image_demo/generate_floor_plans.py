#!/usr/bin/env python3
"""Generate realistic synthetic floor plan images for the real-image demo section.

Produces two floor plans that mimic architectural drawings:
  1. case1_simple_rectangular.png  -- 6m x 4m rectangular room (A1/A2 complexity)
  2. case2_lshaped_composite.png   -- L-shaped room (A3/A4 complexity)

Features: hatched walls, door arcs, window symbols, furniture outlines,
dimension lines with arrows, and a title block.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Arc, FancyArrowPatch, Rectangle
import matplotlib.lines as mlines
import numpy as np
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent

# ── drawing helpers ──────────────────────────────────────────────────────────

WALL_THICKNESS = 0.18          # metres (drawn to scale)
WALL_COLOR = "#222222"
WALL_HATCH = "///"
FURNITURE_COLOR = "#888888"
DIM_COLOR = "#0055AA"
FONT = {"family": "monospace"}


def _draw_wall_segment(ax, x0, y0, x1, y1, t=WALL_THICKNESS):
    """Draw a thick hatched wall between two points (axis-aligned only)."""
    dx, dy = x1 - x0, y1 - y0
    if abs(dy) < 1e-9:          # horizontal wall
        rect = Rectangle((min(x0, x1), y0 - t / 2), abs(dx), t,
                          linewidth=1.2, edgecolor=WALL_COLOR,
                          facecolor="#CCCCCC", hatch=WALL_HATCH, zorder=2)
    else:                        # vertical wall
        rect = Rectangle((x0 - t / 2, min(y0, y1)), t, abs(dy),
                          linewidth=1.2, edgecolor=WALL_COLOR,
                          facecolor="#CCCCCC", hatch=WALL_HATCH, zorder=2)
    ax.add_patch(rect)


def _draw_walls_from_polygon(ax, pts, t=WALL_THICKNESS, gaps=None):
    """Draw walls along a closed polygon, optionally leaving gaps.

    *gaps* is a list of (seg_index, gap_start_frac, gap_end_frac) that
    suppress part of a wall segment (used for doors/windows).
    """
    gaps = gaps or []
    gap_map = {}
    for seg_idx, gs, ge in gaps:
        gap_map.setdefault(seg_idx, []).append((gs, ge))

    n = len(pts)
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        segs_to_draw = [(0.0, 1.0)]
        if i in gap_map:
            for gs, ge in gap_map[i]:
                new_segs = []
                for s, e in segs_to_draw:
                    if ge <= s or gs >= e:
                        new_segs.append((s, e))
                    else:
                        if s < gs:
                            new_segs.append((s, gs))
                        if ge < e:
                            new_segs.append((ge, e))
                segs_to_draw = new_segs
        for s, e in segs_to_draw:
            sx = x0 + (x1 - x0) * s
            sy = y0 + (y1 - y0) * s
            ex = x0 + (x1 - x0) * e
            ey = y0 + (y1 - y0) * e
            _draw_wall_segment(ax, sx, sy, ex, ey, t)


def _draw_door(ax, cx, cy, width=0.9, orientation="up"):
    """Draw a door symbol (arc + line)."""
    angles = {"up": (0, 90), "down": (180, 270), "left": (90, 180), "right": (270, 360)}
    a1, a2 = angles[orientation]
    arc = Arc((cx, cy), width * 2, width * 2, angle=0, theta1=a1, theta2=a2,
              linewidth=1.0, edgecolor=WALL_COLOR, linestyle="--", zorder=3)
    ax.add_patch(arc)
    # door leaf line
    rad = np.radians(a1)
    ax.plot([cx, cx + width * np.cos(rad)], [cy, cy + width * np.sin(rad)],
            color=WALL_COLOR, linewidth=1.5, zorder=3)
    rad2 = np.radians(a2)
    ax.plot([cx, cx + width * np.cos(rad2)], [cy, cy + width * np.sin(rad2)],
            color=WALL_COLOR, linewidth=1.5, zorder=3)


def _draw_window(ax, x0, y0, x1, y1):
    """Draw a window symbol (double-line with glass panes)."""
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    dx, dy = x1 - x0, y1 - y0
    length = np.hypot(dx, dy)
    if abs(dy) < 1e-9:  # horizontal
        # two parallel lines
        for off in [-0.06, 0.06]:
            ax.plot([x0, x1], [y0 + off, y1 + off],
                    color="#0077CC", linewidth=2.0, zorder=4)
        # glass pane ticks
        n_panes = max(2, int(length / 0.4))
        for frac in np.linspace(0.05, 0.95, n_panes):
            px = x0 + dx * frac
            ax.plot([px, px], [y0 - 0.08, y0 + 0.08],
                    color="#0077CC", linewidth=0.8, zorder=4)
    else:  # vertical
        for off in [-0.06, 0.06]:
            ax.plot([x0 + off, x1 + off], [y0, y1],
                    color="#0077CC", linewidth=2.0, zorder=4)
        n_panes = max(2, int(length / 0.4))
        for frac in np.linspace(0.05, 0.95, n_panes):
            py = y0 + dy * frac
            ax.plot([x0 - 0.08, x0 + 0.08], [py, py],
                    color="#0077CC", linewidth=0.8, zorder=4)


def _dim_line(ax, x0, y0, x1, y1, label, offset=0.35, side="outside"):
    """Draw an architectural dimension line with arrows and label."""
    dx, dy = x1 - x0, y1 - y0
    length = np.hypot(dx, dy)
    if abs(dy) < 1e-9:  # horizontal
        yoff = -offset if side == "outside" else offset
        # extension lines
        ax.plot([x0, x0], [y0, y0 + yoff], color=DIM_COLOR, linewidth=0.5, zorder=1)
        ax.plot([x1, x1], [y1, y1 + yoff], color=DIM_COLOR, linewidth=0.5, zorder=1)
        # dimension line
        ax.annotate("", xy=(x1, y0 + yoff), xytext=(x0, y0 + yoff),
                     arrowprops=dict(arrowstyle="<->", color=DIM_COLOR, lw=0.8))
        ax.text((x0 + x1) / 2, y0 + yoff + 0.08, label,
                ha="center", va="bottom", fontsize=8, color=DIM_COLOR, **FONT)
    else:  # vertical
        xoff = -offset if side == "outside" else offset
        ax.plot([x0, x0 + xoff], [y0, y0], color=DIM_COLOR, linewidth=0.5, zorder=1)
        ax.plot([x1, x1 + xoff], [y1, y1], color=DIM_COLOR, linewidth=0.5, zorder=1)
        ax.annotate("", xy=(x0 + xoff, y1), xytext=(x0 + xoff, y0),
                     arrowprops=dict(arrowstyle="<->", color=DIM_COLOR, lw=0.8))
        ax.text(x0 + xoff - 0.08, (y0 + y1) / 2, label,
                ha="right", va="center", fontsize=8, color=DIM_COLOR, rotation=90, **FONT)


def _furniture_rect(ax, x, y, w, h, label=None):
    """Simple rectangular furniture outline."""
    rect = Rectangle((x, y), w, h, linewidth=0.8, edgecolor=FURNITURE_COLOR,
                      facecolor="#F0F0F0", linestyle="-", zorder=1)
    ax.add_patch(rect)
    if label:
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=6, color=FURNITURE_COLOR, **FONT)


def _title_block(ax, title, subtitle=""):
    """Add a title block in the lower-right."""
    ax.text(0.98, 0.02, title, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=10, fontweight="bold", **FONT)
    if subtitle:
        ax.text(0.98, 0.06, subtitle, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7, color="#666666", **FONT)


def _north_arrow(ax, x, y, size=0.3):
    """Draw a north arrow."""
    ax.annotate("", xy=(x, y + size), xytext=(x, y),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5))
    ax.text(x, y + size + 0.08, "N", ha="center", va="bottom",
            fontsize=9, fontweight="bold", **FONT)


# ── Case 1: Simple Rectangular Room ─────────────────────────────────────────

def make_case1():
    W, H = 6.0, 4.0
    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5))
    ax.set_aspect("equal")
    ax.set_xlim(-1.0, W + 1.2)
    ax.set_ylim(-1.0, H + 1.0)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # floor
    floor = Rectangle((0, 0), W, H, linewidth=0, facecolor="#FAFAF5", zorder=0)
    ax.add_patch(floor)

    # walls -- polygon with gaps for door and window
    pts = [(0, 0), (W, 0), (W, H), (0, H)]
    # seg 0: south wall (0,0)->(6,0) -- door gap at ~40-55%
    # seg 2: north wall (6,4)->(0,4) -- window gap at ~25-65%
    door_start, door_end = 0.38, 0.53   # ~0.9m door centred at x=2.7
    # Window on north wall: 1.5m wide, centred at x=3.0 (mid-wall)
    # Seg 2 goes (6,4)->(0,4), so frac f maps to x = 6 - 6*f.
    # We want window from x=2.25 to x=3.75 => frac_start=(6-3.75)/6=0.375, frac_end=(6-2.25)/6=0.625
    win_frac_start, win_frac_end = 0.375, 0.625  # 1.5m window centred on north wall
    gaps = [
        (0, door_start, door_end),
        (2, win_frac_start, win_frac_end),  # seg 2 goes right-to-left
    ]
    _draw_walls_from_polygon(ax, pts, gaps=gaps)

    # door on south wall
    door_cx = W * door_start
    _draw_door(ax, door_cx, 0.0, width=0.9, orientation="up")

    # window on north wall — centred, 1.5m wide
    wx0 = W * (1 - win_frac_end)   # = 6 * 0.375 = 2.25
    wx1 = W * (1 - win_frac_start) # = 6 * 0.625 = 3.75
    _draw_window(ax, wx0, H, wx1, H)

    # furniture
    # desk (1.2m x 0.6m) upper-left area
    _furniture_rect(ax, 0.4, 2.8, 1.2, 0.6, "Desk")
    # chair
    _furniture_rect(ax, 0.7, 2.2, 0.5, 0.5, "Chair")
    # bookshelf along east wall
    _furniture_rect(ax, 5.1, 0.5, 0.7, 2.0, "Shelf")

    # dimension lines
    _dim_line(ax, 0, 0, W, 0, f"{W:.1f} m", offset=0.55, side="outside")
    _dim_line(ax, W, 0, W, H, f"{H:.1f} m", offset=0.55, side="outside")

    # room label
    ax.text(W / 2, H / 2, "OFFICE\n24.0 m²", ha="center", va="center",
            fontsize=12, color="#444444", **FONT)

    _north_arrow(ax, -0.5, H - 0.3)
    _title_block(ax, "Case 1 — Simple Rectangular Room",
                 "6.0 m \u00d7 4.0 m  |  A1/A2 complexity")

    # scale bar
    ax.plot([0.3, 1.3], [-0.75, -0.75], color="black", linewidth=2)
    ax.text(0.8, -0.85, "1 m", ha="center", va="top", fontsize=7, **FONT)

    out = OUTPUT_DIR / "case1_simple_rectangular.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {out}")


# ── Case 2: L-shaped Composite Room ─────────────────────────────────────────

def make_case2():
    # Main block: 5m x 4m, extension: 3m x 2.5m (attached on the east side, lower portion)
    fig, ax = plt.subplots(1, 1, figsize=(11, 8))
    ax.set_aspect("equal")
    ax.set_xlim(-1.2, 8.5)
    ax.set_ylim(-1.2, 5.0)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # L-shape vertices (counter-clockwise)
    #   D─────────────C
    #   |             |
    #   |     main    F───────G
    #   |             |  ext  |
    #   A─────────────B───────H  (note: H is at bottom-right of extension)
    # wait, let me define clearly:
    # Main: (0,0) to (5,4)
    # Extension: (5,0) to (8, 2.5)
    A = (0, 0)
    B = (5, 0)
    C = (5, 4)
    D = (0, 4)
    # extension
    E = (5, 2.5)  # internal corner (already on main east wall)
    F = (8, 2.5)
    G = (8, 0)
    # Full polygon: A -> B -> G -> F -> E -> C -> D
    # But E is at (5,2.5), C is at (5,4) -- the east wall of main from 2.5 to 4 is exposed
    pts = [A, B, G, F, E, C, D]
    # segments: 0:A-B, 1:B-G, 2:G-F, 3:F-E, 4:E-C, 5:C-D, 6:D-A

    # floor
    from matplotlib.patches import Polygon
    floor = Polygon(pts, closed=True, facecolor="#FAFAF5", edgecolor="none", zorder=0)
    ax.add_patch(floor)

    # gaps: door on west wall (seg 6: D(0,4)->A(0,0)), window on seg 4: E(5,2.5)->C(5,4)
    # seg 6 goes from (0,4) to (0,0) -- 4m long. Door at ~50-72.5% => y from 2.0 to 1.1
    door_s, door_e = 0.50, 0.725  # ~0.9m door
    # seg 1: B(5,0)->G(8,0) 3m long. Window at 20-70% => x from 5.6 to 7.1 (~1.5m window)
    win_s, win_e = 0.20, 0.70

    gaps = [
        (6, door_s, door_e),
        (1, win_s, win_e),
    ]
    _draw_walls_from_polygon(ax, pts, gaps=gaps)

    # door on west wall -- seg 6 goes (0,4) to (0,0)
    # frac 0.50 => y = 4 - 4*0.50 = 2.0, frac 0.725 => y = 4 - 4*0.725 = 1.1
    door_y_top = 4.0 - 4.0 * door_s  # 2.0
    _draw_door(ax, 0.0, door_y_top, width=0.9, orientation="right")

    # window on south wall of extension -- seg 1: (5,0)->(8,0)
    wx0 = 5.0 + 3.0 * win_s  # 5.6
    wx1 = 5.0 + 3.0 * win_e  # 7.1
    _draw_window(ax, wx0, 0.0, wx1, 0.0)

    # furniture -- sofa in main room
    _furniture_rect(ax, 1.0, 0.4, 2.0, 0.8, "Sofa")

    # dimension lines -- main block
    _dim_line(ax, 0, 0, 5, 0, "5.0 m", offset=0.55, side="outside")
    _dim_line(ax, 0, 4, 0, 0, "4.0 m", offset=0.55, side="outside")
    # extension
    _dim_line(ax, 5, 0, 8, 0, "3.0 m", offset=0.90, side="outside")
    _dim_line(ax, 8, 0, 8, 2.5, "2.5 m", offset=0.55, side="outside")

    # room labels
    ax.text(2.5, 2.5, "LIVING ROOM\n20.0 m²", ha="center", va="center",
            fontsize=11, color="#444444", **FONT)
    ax.text(6.5, 1.25, "ALCOVE\n7.5 m²", ha="center", va="center",
            fontsize=9, color="#444444", **FONT)

    # dashed line at internal boundary
    ax.plot([5, 5], [0, 2.5], color="#999999", linewidth=0.6, linestyle=":", zorder=1)

    _north_arrow(ax, -0.7, 3.5)
    _title_block(ax, "Case 2 — L-shaped Composite Room",
                 "5.0\u00d74.0 m + 3.0\u00d72.5 m  |  A3/A4 complexity")

    # scale bar
    ax.plot([0.3, 1.3], [-0.90, -0.90], color="black", linewidth=2)
    ax.text(0.8, -1.0, "1 m", ha="center", va="top", fontsize=7, **FONT)

    out = OUTPUT_DIR / "case2_lshaped_composite.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    make_case1()
    make_case2()
    print("Done.")
