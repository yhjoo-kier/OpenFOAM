# Figure 5 QC Log — Aggregate performance across input views (scale-calibrated main-setting refresh)

- Date: 2026-03-20
- Figure: 5
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure5_view_aggregate.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure5_view_aggregate_performance.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure5_view_aggregate_performance.png`
- PDF render check:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/figure5_view_aggregate_pdf_render.png`
- Export meta:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure5_view_aggregate_performance_meta.json`

## 1. Source artifact / metric selection

Current manuscript target was refreshed from the older no-scaling aggregate to the scale-calibrated full aggregate.

- main source artifact:
  - `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
- supporting comparison note:
  - `/home/yhjoo/projects/OpenFOAM/docs/26-03-20_posthoc_scaled_full_comparison.md`
- selected view-level metrics:
  - mean structural score
  - mean CFD score
  - room-kind match rate
  - opening-wall match rate

Rationale:
- user/topic rule says the paper main pipeline is the single-factor scale-calibrated setting
- therefore Figure 5 should no longer point to the legacy no-scaling aggregate as its source artifact

## 2. Layout / export decision

- intended width: `double-column`
- panel layout: `1 x 2`
  - `(a)` score metrics
  - `(b)` agreement rates
- export size: `7.20 in x 3.95 in`
- PNG dpi: `600`
- PDF backend: matplotlib vector PDF
- font family request: `Arial -> Liberation Sans -> DejaVu Sans`
- selected font in this environment: `DejaVu Sans`
- internal caption / figure-number rule:
  - confirmed no `Figure 5`, `Fig. 5`, or caption-style sentence exists inside the artwork

## 3. Current design state

Latest refresh focuses on manuscript-safe readability rather than dense numeric annotation.

- removed the earlier per-bar numeric labels to reduce reduction-time crowding
- retained short panel titles, axis labels, and panel labels only
- latest variant moves each legend above the plotting region to eliminate compact-render crowding
- PNG/PDF export now avoids tight-bbox compression so both formats stay visually aligned
- source artifact now reflects the scale-calibrated full aggregate

## 4. QC status in this turn

### Self visual QC
- verdict: `PASS`
- summary:
  - revised legend placement removes overlap risk from the data area
  - PNG and PDF renders are now visually aligned with no obvious tight-bbox drift
  - top spacing, panel labels, and axis text remain readable at the intended double-column size

### Subagent visual QC
- verdict: `PASS`
- note file:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/26-03-20_figure5_subagent_visual_qc.md`
- subagent summary:
  - PNG/PDF consistency is acceptable after the legend/top-margin revision
  - legends no longer crowd the bars
  - no clipping or forbidden in-figure prose detected

### Gemini visual QC
- verdict: `PASS`
- artifact:
  - `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure5-visual-qc-20260320-rerun.md`
- Gemini summary:
  - judged the rerun manuscript-ready with no clipping, acceptable readability, and no remaining blockers

## 5. Final disposition for this turn

- overall figure completion verdict: `COMPLETE`
- reason:
  - 3-stage QC unanimity was achieved on the rerun (`self PASS`, `subagent PASS`, `Gemini PASS`)
- consequence:
  - Figure 5 can be marked complete in the v2 checklist
  - Figure 5 can be stored in Google Drive under `Fig05_view_aggregate/`

## 6. Delivery status

- Google Drive folder:
  - `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig05_view_aggregate/`
- delivered files:
  - `figure5_view_aggregate_performance.pdf`
  - `figure5_view_aggregate_performance.png`
  - `qc_summary.md`
