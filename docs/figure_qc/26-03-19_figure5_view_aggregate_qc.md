# Figure 5 QC Log — Aggregate performance across input views

- Date: 2026-03-19
- Figure: 5
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure5_view_aggregate.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure5_view_aggregate_performance.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure5_view_aggregate_performance.png`
- PDF render check:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/figure5_view_aggregate_pdf_render.png`
- Export meta:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure5_view_aggregate_performance_meta.json`
- Source artifacts:
  - `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary.json`
  - `/home/yhjoo/projects/OpenFOAM/docs/26-03-19_cli_eval_aggregate_results.md`

## 1. Source artifact / metric selection

Selected aggregate by-view metrics from the completed frozen-20 benchmark (`100` tasks):

- mean structural score
- mean CFD score
- room-kind match rate
- opening-wall match rate

Headline message preserved for manuscript use:
- `floorplan` is strongest overall
- `section` is the clearest stress view
- opening-wall agreement drops most sharply for `section`

## 2. Layout / panel / caption design

- intended width: `double-column`
- rationale:
  - two metric families must remain readable together without shrinking labels/annotations too aggressively
  - single-column stacking was rejected up front because bar labels and category labels would become too small for stable print use
- panel layout: `1 x 2`
  - `(a)` score metrics
  - `(b)` agreement rates
- subfigure labels: `(a)`, `(b)` only
- subcaption handling:
  - no internal subcaption sentences
  - manuscript caption will explain the metric-color mapping and interpretation outside the figure
- export size: `7.10 in x 3.20 in`
- PNG dpi: `600`
- PDF backend / vector status:
  - matplotlib PDF backend
  - axes/text/vector bars preserved as vector in PDF
- font family request: `Arial -> Liberation Sans -> DejaVu Sans`
- selected font in this environment: `DejaVu Sans`
- Arial availability note:
  - Arial was not available in the current environment
  - a sans-serif fallback was used and is recorded explicitly here for reproducibility
- internal caption / figure-number rule:
  - confirmed no `Figure 5`, `Fig. 5`, or caption-style sentence exists inside the artwork

## 3. Self visual QC

### Checks
- [x] no in-figure figure number or caption sentence
- [x] panel titles remain short and non-caption-like
- [x] labels/ticks/values readable at intended double-column width
- [x] no clipping or overlap in final export
- [x] no panel misalignment
- [x] PNG and PDF-rendered PNG are visually consistent
- [x] sans-serif typography is consistent

### Notes
- Earlier regenerated draft failed Gemini QC because the right-edge `100%` labels in panel `(b)` were too close to the boundary.
- Final script revision expanded panel `(b)` x-limit and shifted the in-bar text left slightly.
- After that correction, the final PNG/PDF pair no longer shows right-edge clipping.

## 4. Self verdict

- self QC verdict: `PASS`
- summary:
  - final regenerated version satisfies the new manuscript rule set: no internal figure number/caption text, Arial-first fallback recorded, and 600 dpi PNG + vector PDF both produced.

## 5. Independent QC

### Subagent QC
- status: `PASS`
- session: `agent:main:subagent:c937f603-9f5e-4341-a3f6-a8ee9f059396`
- note file: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/26-03-19_figure5_subagent_visual_qc.md`
- summary:
  - corrected version shows no remaining clipping/overlap/alignment issue
  - PNG/PDF consistency acceptable
  - DejaVu Sans fallback explicitly noted

### Gemini CLI QC
- status: `PASS`
- artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure5-visual-qc-20260319-123950.md`
- summary:
  - final version judged aligned and readable
  - no forbidden internal figure number/caption text detected
  - no remaining clipping in annotations after the right-edge fix

## 6. Final disposition

- final verdict: `PASS`
- production status: approved for manuscript production use under the tightened figure rules
- delivery status:
  - production plan may remain checked complete
  - Google Drive QC-pass folder may contain this regenerated final version
