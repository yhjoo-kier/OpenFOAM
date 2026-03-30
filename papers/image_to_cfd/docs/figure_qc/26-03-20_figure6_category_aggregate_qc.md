# Figure 6 QC Log — Aggregate performance across benchmark categories (scale-calibrated main-setting refresh)

- Date: 2026-03-20
- Figure: 6
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure6_category_aggregate.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure6_category_aggregate_performance.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure6_category_aggregate_performance.png`
- PDF render check:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/figure6_category_aggregate_pdf_render.png`
- Export meta:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure6_category_aggregate_performance_meta.json`

## 1. Source artifact / metric selection

Current manuscript target was refreshed from the older no-scaling aggregate to the scale-calibrated full aggregate.

- main source artifact:
  - `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
- selected category order:
  - `A1`, `A2`, `A3`, `A4`
- selected metrics:
  - panel (a): mean structural score, mean CFD score
  - panel (b): room-kind match rate, opening-wall match rate

Rationale:
- topic rule says the paper main pipeline is the single-factor scale-calibrated setting
- Figure 6 therefore had to be regenerated from the scale-calibrated aggregate rather than the legacy no-scaling manifest

## 2. Layout / export decision

- intended width: `double-column`
- panel layout: `1 x 2`
  - `(a)` score metrics
  - `(b)` agreement rates
- export size: `7.30 in x 4.15 in`
- PNG dpi: `600`
- PDF backend: matplotlib vector PDF
- font family request: `Arial -> Liberation Sans -> DejaVu Sans`
- selected font in this environment: `DejaVu Sans`
- internal caption / figure-number rule:
  - confirmed no `Figure 6`, `Fig. 6`, or caption-style sentence exists inside the artwork

## 3. Iteration notes

- initial refresh still inherited the older layout margin assumptions and duplicated y labels on the right panel; these were removed
- a first rerun fixed the y-label duplication but left the left margin too tight for subagent QC
- final rerun widened the left margin and slightly increased overall figure/font sizing to improve double-column readability while keeping the panel structure unchanged

## 4. QC status in this turn

### Self visual QC
- verdict: `PASS`
- method:
  - multimodal inspection of the exported PNG and a PDF-rendered PNG copy
- summary:
  - no clipping or truncation remained after the final left-margin fix
  - legends stay outside the plotting area and do not crowd the bars
  - PNG and PDF renders are materially consistent
  - typography is readable at intended double-column size

### Subagent visual QC
- verdict: `PASS`
- session:
  - `agent:main:subagent:27bf20df-1c53-43fa-8d55-9969438525cc`
- summary:
  - subplot alignment and export consistency judged clean
  - no forbidden in-figure prose detected
  - final caution was only that text sizing is somewhat tight, but still acceptable/manuscript-ready

### Gemini visual QC
- verdict: `PASS`
- artifact:
  - `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure6-visual-qc-20260320-231055.md`
- summary:
  - Gemini judged the figure professionally formatted and ready for double-column publication
  - no clipping, no material PNG/PDF mismatch, and no forbidden in-figure prose were reported

## 5. Final disposition for this turn

- overall figure completion verdict: `COMPLETE`
- reason:
  - 3-stage QC unanimity achieved (`self PASS`, `subagent PASS`, `Gemini PASS`)
- consequence:
  - Figure 6 can be marked complete in the v2 checklist
  - Figure 6 can be stored in Google Drive under `Fig06_category_aggregate/`
