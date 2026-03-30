#!/usr/bin/env python3
"""Scatter plot: structural score vs CFD fidelity score for all evaluation cases."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
STATS_FILE = Path(__file__).resolve().parents[3] / "benchmark" / "manifests" / "evaluation_statistics.json"
OUTPUT_DIR = Path(__file__).resolve().parents[3] / "results" / "paper_figures_agent_b"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

data = json.load(open(STATS_FILE))
records = data["records"]

# ---------------------------------------------------------------------------
# Visual encoding
# ---------------------------------------------------------------------------
CATEGORY_COLORS = {
    "A1": "#1f77b4",  # blue
    "A2": "#ff7f0e",  # orange
    "A3": "#2ca02c",  # green
    "A4": "#d62728",  # red
}

VIEW_MARKERS = {
    "perspective":  "o",  # circle
    "birdseye":     "s",  # square
    "floorplan":    "D",  # diamond
    "wireframe":    "^",  # triangle-up
    "section":      "v",  # triangle-down
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 5.0))

for rec in records:
    ax.scatter(
        rec["structural_score"],
        rec["cfd_score"],
        c=CATEGORY_COLORS[rec["category"]],
        marker=VIEW_MARKERS[rec["view"]],
        s=48,
        alpha=0.80,
        edgecolors="white",
        linewidths=0.4,
        zorder=3,
    )

# ---------------------------------------------------------------------------
# Legend — category colors + view markers
# ---------------------------------------------------------------------------
cat_handles = [
    mlines.Line2D([], [], color=CATEGORY_COLORS[c], marker="o", linestyle="None",
                  markersize=7, label=c)
    for c in ["A1", "A2", "A3", "A4"]
]
view_handles = [
    mlines.Line2D([], [], color="gray", marker=VIEW_MARKERS[v], linestyle="None",
                  markersize=7, label=v)
    for v in ["perspective", "birdseye", "floorplan", "wireframe", "section"]
]

leg1 = ax.legend(handles=cat_handles, title="Category", loc="upper left",
                 fontsize=10, title_fontsize=11, framealpha=0.9)
ax.add_artist(leg1)
ax.legend(handles=view_handles, title="View", loc="lower right",
          fontsize=10, title_fontsize=11, framealpha=0.9)

# ---------------------------------------------------------------------------
# Axes
# ---------------------------------------------------------------------------
ax.set_xlabel("Structural score", fontsize=12)
ax.set_ylabel("CFD fidelity score", fontsize=12)
ax.tick_params(labelsize=11)
ax.set_xlim(-0.02, 1.05)
ax.set_ylim(-0.02, 1.05)
ax.grid(True, alpha=0.25, linewidth=0.5)

fig.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
fig.savefig(OUTPUT_DIR / "fig_discuss_scatter_structural_cfd.pdf", dpi=600, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig_discuss_scatter_structural_cfd.png", dpi=600, bbox_inches="tight")
plt.close(fig)

print(f"Saved to {OUTPUT_DIR / 'fig_discuss_scatter_structural_cfd.pdf'}")
print(f"Saved to {OUTPUT_DIR / 'fig_discuss_scatter_structural_cfd.png'}")
