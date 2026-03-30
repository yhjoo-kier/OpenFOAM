# Figure 9 QC Log — obstacle hallucination with limited CFD penalty

> Date: 2026-03-19
> Figure: 9
> Status: PASS (self QC + subagent QC + Gemini QC)

## Figure intent

A3 composite cases can hallucinate extra obstacles while preserving opening/topology cues strongly enough that CFD fidelity remains high.

## Source artifact / case lock

Selected representative tasks:

1. `bench_a3_01 / floorplan`
   - reference obstacle count `0`
   - predicted obstacle count `3`
   - opening walls `N/W -> N/W`
   - structural score `0.750`
   - CFD score `0.697`
2. `bench_a3_05 / floorplan`
   - reference obstacle count `1`
   - predicted obstacle count `3`
   - opening walls `E/W -> E/W`
   - structural score `0.750`
   - CFD score `0.693`

Rationale:
- both are clean A3 examples of the `opening_topology_preserved_with_obstacle_hallucination` pattern
- both keep the room composite
- both preserve opening-wall topology
- both keep headline CFD fidelity near `0.69`, which makes the message visually legible without overloading the panel set

## Layout / manuscript decision

- intended width: **double-column**
- reason: a `2 x 3` comparative layout was needed to show `(input image / reference geometry / predicted geometry)` for two representative A3 cases without collapsing font size or annotation readability
- panel layout: **2 x 3**
- panel labels: integrated directly into panel titles as `(a)`–`(f)`
- subfigure labels: yes, embedded in title strings rather than detached floating labels
- subcaption handling: no separate subcaptions; each panel title carries the local role, while the main caption carries the manuscript claim
- export size: approximately `7.40 x 5.80 in`
- PNG dpi: **600**
- PDF backend: matplotlib PDF (vector text / vector primitives)

## Output artifacts

- script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure9_obstacle_hallucination.py`
- PDF: `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.pdf`
- PNG: `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.png`
- Gemini QC artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure9-visual-qc-20260319-114625.md`
- QC image copy: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/figure9_obstacle_hallucination_qc_input.png`

## Self visual inspection

### First-pass issues found

Initial draft had two blocking layout problems:
- row-title / panel-title / panel-label collisions in the upper title band
- left-panel explanatory text overlapped adjacent panels and was unsafe at journal reduction

### Revision applied

- removed long row-level explanation lines from the panel belt
- shortened center/right panel titles to `Reference geometry` / `Predicted geometry`
- merged panel labels directly into the panel titles
- retained the manuscript message in the figure-level subtitle, predicted-panel callouts, and metric boxes
- increased bottom margin slightly so the footer note no longer sat on the canvas edge

### Final self-QC verdict

**PASS**

Checked items:
- no visible clipping on PDF/PNG export
- title belt no longer collides with panel callouts
- metric boxes remain legible in double-column layout
- callouts communicate the exact claim without requiring caption rescue
- row-wise comparison is easy to scan: input -> reference -> prediction
- PDF/PNG layout was generated from the same script and appears consistent

## Independent subagent visual QC

Session result (`agent:main:subagent:e81f5667-d47e-41cf-9abe-52f02adf7495`):

- Verdict: **PASS**
- Blocking issues: none
- Non-blocking polish notes:
  - metric-box text and footer note are near the lower safe bound for aggressive reduction
  - callout boxes could breathe slightly more
  - left image panels are visually lighter than geometry panels
- Manuscript-readiness summary:
  - figure is manuscript-ready and communicates the target claim effectively

Interpretation:
- independent QC confirmed that the earlier title-belt issue was resolved in the revised version
- remaining comments are polish-level only

## Gemini CLI visual QC

Artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure9-visual-qc-20260319-114625.md`

Gemini verdict: **PASS**

Key notes:
- no blocking issues
- suggested mild symmetry improvement for left image panels
- warned that metric-box contrast should remain print-safe
- noted interior `y [m]` labels are slightly crowded between columns
- panel `(f)` callout text is a little tight but still acceptable

Interpretation:
- Gemini agrees the figure is publication-ready after only minor optional polish
- no mandatory regeneration was requested

## Final release decision

**Approved as Figure 9 production candidate.**

This figure passes the required three-stage QC gate:
1. self visual inspection — PASS
2. separate subagent visual QC — PASS
3. Gemini CLI visual QC — PASS

Therefore Figure 9 can be marked complete in the production plan.
