# Figure 8 QC Log — Section-view composite-room collapse

- Date: 2026-03-19
- Figure: Figure 8
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure8_section_collapse.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure8_section_view_composite_collapse.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure8_section_view_composite_collapse.png`
- Source artifacts:
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a3_04/section/reference_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a3_04/section/predicted_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a3_04/section/evaluation_summary.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a4_05/section/reference_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a4_05/section/predicted_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a4_05/section/evaluation_summary.json`
  - supporting interpretation note: `/home/yhjoo/projects/OpenFOAM/docs/26-03-19_cli_eval_aggregate_results.md`

## Figure intent
Show the paper’s clearest topology-failure mode: under section-view input, the only two room-kind collapses in the frozen-20 benchmark both converted composite rooms into rectangular predictions. The figure emphasizes the lost composite arm directly in top-view geometry rather than relying on the blank section-render PNGs.

## Source artifact / case lock
- Cases fixed to the only two `section_room_kind_collapse` tasks in the frozen-20 aggregate manifest:
  - `bench_a3_04 / section`
  - `bench_a4_05 / section`
- Each row uses:
  - left: reference composite room geometry
  - right: predicted top-view abstraction from section-view input
- Metrics shown in prediction panels:
  - structural score
  - CFD score
  - room-block count collapse (`2 -> 1`)
  - opening-wall collapse (`N/S -> E/W`)

## Layout / typography decision
- Initial width candidate: `single-column`
- Final intended width: `double-column`
- Decision reason: single-column reduction made the L-shape vs rectangular collapse, dashed reference overlay, and metric box simultaneously too weak. The figure needs a 2x2 layout with readable axes and annotations.
- Panel layout: `2x2`
- Subfigure labels: yes — `(a)` to `(d)` inside the figure
- Subcaption handling: no separate subcaption blocks inside the image; main manuscript caption should explain that left panels are reference composite rooms and right panels are section-view predictions
- Row semantics:
  - top row: `A3-04 sparse composite`
  - bottom row: `A4-05 dense composite`
- Export size: `7.30 x 5.65 in`
- PNG dpi: `600`
- PDF backend: matplotlib PDF (vector)

## QC protocol completion
1. Source artifact/case fixed — PASS
2. Single/double-column decision + panel design fixed — PASS
3. Script written — PASS
4. PDF + PNG generated — PASS
5. Self visual inspection — PASS
6. Independent subagent visual QC — PASS
7. Gemini CLI visual QC — PASS

## Iteration notes
- First draft used longer titles and footer text; QC flagged top/bottom crowding and borderline small text.
- Revised version shortened the figure headline, removed nonessential footer text, increased figure height, synchronized row-wise axis extents between reference and prediction panels, and enlarged the metric-box text.
- The prediction panels now show the reference composite extent with dashed red outlines so the missing arm remains visually explicit.

## Self visual QC
- Verdict: PASS
- Method: direct inspection of the exported PNG plus rasterized PDF rendering (`pdftoppm`) for layout consistency
- Summary:
  - no visible clipping in either PDF-rendered or PNG output
  - 2x2 panel alignment is stable
  - dashed red overlay clearly shows the missing composite arm in both prediction panels
  - metric boxes remain readable at the intended double-column layout
  - figure message is understandable without extra prose in the plot area

## Subagent visual QC
- Verdict: PASS
- Session: `agent:main:subagent:50ad1e1c-b5e4-464a-aa09-6470cd51c02b`
- Summary:
  - no obvious clipping; only minor proximity of right-side callouts to the edge
  - small text is borderline but still readable at double-column width
  - panel spacing is clean and consistent
  - red dashed overlays explain the composite-to-rectangular collapse clearly
  - judged publication-ready

## Gemini CLI visual QC
- Verdict: PASS
- Artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure8-visual-qc-20260319-104508.md`
- Summary:
  - Gemini reported no blocking issues
  - headline claim is clearly supported by the dashed red overlays
  - only non-blocking polish notes were slight row-spacing tightness, metric font size caution under aggressive downscaling, and left y-axis gutter awareness

## Final disposition
- Figure 8 is approved for current manuscript production use.
- Manuscript insertion recommendation: use the PDF version in LaTeX as a `double-column` discussion/result figure.
- If the final figure set becomes more tightly packed later, only minor row-spacing polish may be considered; this is non-blocking.
