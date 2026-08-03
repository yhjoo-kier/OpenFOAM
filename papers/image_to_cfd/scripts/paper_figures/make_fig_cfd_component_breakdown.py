#!/usr/bin/env python3
"""Fig: CFD agreement score component breakdown (Phase 2, preset-matched).

Revision round 1 (2026-07): rebuilt to read the preset-matched aggregate
manifest (evaluation_statistics_phase2.json) instead of the stale per-case
posthoc evaluation tree, which reflected unmatched inlet-velocity data and
produced an incorrect overall mean (0.453 instead of 0.477).

The figure now shows the overall component decomposition (n = 97 converged
cases) with SD error bars, plus the overall composite CFD agreement score
as a reference line. Note: component means are computed over cases with
available field data, while the composite score assigns 0.0 to unavailable
components, so the composite mean (0.477) sits slightly below the average
of the four component means.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

PAPER_ROOT = Path(__file__).resolve().parents[2]   # papers/image_to_cfd
REPO_ROOT = Path(__file__).resolve().parents[4]    # repository root
STATS_PATH = REPO_ROOT / "benchmark" / "manifests" / "evaluation_statistics_phase2.json"
OUT_DIR = PAPER_ROOT / "results" / "paper_figures_phase2"
PDF_OUT = OUT_DIR / "fig_result_cfd_component_breakdown.pdf"
PNG_OUT = OUT_DIR / "fig_result_cfd_component_breakdown.png"
META_OUT = OUT_DIR / "fig_result_cfd_component_breakdown_meta.json"

COMPONENTS = ["overlap", "vel_mag", "vel_dir", "pressure"]
COMP_LABELS = ["Overlap\nratio", "Velocity\nmagnitude", "Velocity\ndirection", "Pressure"]
COMP_COLORS = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759"]

FONT_CANDIDATES = ["Arial", "Liberation Sans", "DejaVu Sans"]


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    for c in FONT_CANDIDATES:
        if c in available:
            return c
    return "sans-serif"


def main() -> None:
    stats = json.loads(STATS_PATH.read_text(encoding="utf-8"))
    overall = stats["overall"]
    comps = overall["components"]
    cfd = overall["cfd"]

    means = np.array([comps[c]["mean"] for c in COMPONENTS])
    sds = np.array([comps[c]["sd"] for c in COMPONENTS])
    n = cfd["n"]
    overall_mean = cfd["mean"]

    font = pick_font()
    plt.rcParams.update({
        "font.family": font,
        "font.size": 10,
        "axes.labelsize": 10.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(5.6, 3.6))

    x = np.arange(len(COMPONENTS))
    bars = ax.bar(x, means, width=0.6, color=COMP_COLORS, zorder=3,
                  yerr=sds, error_kw=dict(ecolor="#4B5563", elinewidth=0.8,
                                          capsize=3, capthick=0.8))

    for xi, m in zip(x, means):
        ax.text(xi, 0.02, f"{m:.3f}", ha="center", va="bottom",
                fontsize=9, color="#1F2937", zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(COMP_LABELS, fontsize=9.5)
    ax.set_ylabel("Component score")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.3, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Overall composite score reference line
    ax.axhline(overall_mean, color="#666666", linestyle="--", linewidth=1.0,
               zorder=2, alpha=0.8)
    ax.text(len(COMPONENTS) - 0.55, overall_mean + 0.02,
            f"Overall CFD agreement score: {overall_mean:.3f} (n = {n})",
            fontsize=8.5, color="#444444", ha="right", va="bottom",
            fontstyle="italic")

    fig.tight_layout(pad=1.0)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT, pad_inches=0.02)
    fig.savefig(PNG_OUT, dpi=600, pad_inches=0.02)

    meta = {
        "source": str(STATS_PATH.relative_to(REPO_ROOT)),
        "n": n,
        "overall_cfd_mean": overall_mean,
        "components": {c: {"mean": comps[c]["mean"], "sd": comps[c]["sd"]} for c in COMPONENTS},
        "font": font,
        "note": "Component means over cases with available field data; composite assigns 0.0 to unavailable components.",
    }
    META_OUT.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    print(f"Font: {font}")
    print(f"Overall CFD mean: {overall_mean} (n={n})")
    print(f"Component means: {dict(zip(COMPONENTS, means.round(4).tolist()))}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
