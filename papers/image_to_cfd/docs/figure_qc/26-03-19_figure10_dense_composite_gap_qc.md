# Figure 10 QC Log — Dense composite structure-vs-CFD gap

- Date: 2026-03-19
- Figure: 10
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure10_dense_composite_gap.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.png`
- PDF render check:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png`

## 1. Source artifact / case selection

Selected dense-composite floor-plan cases with correct room kind + opening walls but only moderate CFD fidelity:

1. `bench_a4_02 / floorplan`
   - structural = `0.8125`
   - CFD = `0.4023`
   - room kind match = `True`
   - opening-wall match = `True`
2. `bench_a4_04 / floorplan`
   - structural = `1.0000`
   - CFD = `0.4596`
   - room kind match = `True`
   - opening-wall match = `True`

Aggregate context used in header:
- category `A4` structural = `0.7124`
- category `A4` CFD = `0.4266`
- room-kind match = `96%`
- opening-wall match = `80%`

## 2. Layout / panel / caption design

- intended width: `double-column`
- rationale:
  - Figure needs two representative cases, each with geometry evidence and enlarged flow comparison.
  - Initial denser 3-column layout failed QC because lower-panel readability degraded at two-column scale.
  - Current 2x2 layout gives each flow comparison substantially more space.
- panel layout: `2 x 2`
  - `(a)` A4-02 geometry + input
  - `(b)` A4-02 flow comparison
  - `(c)` A4-04 geometry + input
  - `(d)` A4-04 flow comparison
- subfigure labels: embedded directly into each panel title
- subcaption handling:
  - no separate external subcaptions
  - panel titles carry panel identity and role
  - figure-level title + one short aggregate header line provide the general claim and category context
- main caption role (for manuscript later): explain that dense composite layouts can preserve coarse structure/openings while still degrading downstream CFD fidelity
- export size: `7.35 in x 6.65 in`
- PNG dpi: `600`
- PDF backend / vector status:
  - matplotlib PDF backend
  - text/axes/vector primitives preserved
  - embedded flow/input imagery remain raster inside vector PDF

## 3. Self visual QC

### Pass criteria checked

- [x] panel alignment is stable
- [x] no annotation clipping observed
- [x] geometry overlays are readable
- [x] flow panels enlarged enough to read as qualitative comparison
- [x] structure-vs-CFD gap is explicit via score bars + CFD badges + header context
- [x] PDF render is visually consistent with direct PNG export
- [x] footer/title/header text are readable at figure scale

### Notes

- Earlier denser layout was rejected internally after subagent QC because flow panels and metric text were too small.
- Revised layout removes the crowded top scatter and moves score-gap evidence into per-case inset bars, which materially improves readability.
- Geometry panels use dashed blue overlays for reference footprints and orange filled obstacles for predicted blockage placement.

## 4. Self verdict

- self QC verdict: `PASS`
- remaining caution:
  - manuscript caption should explicitly state that the figure shows representative floor-plan exemplars within the broader A4 regime, while the aggregate category statistics are summarized in the header.

## 5. Independent QC

### Subagent QC
- status: `PASS`
- session: `agent:main:subagent:0217a215-18b4-44a8-a401-1c7d353bde55`
- summary:
  - clean 2x2 layout and clear panel hierarchy
  - no clipping observed
  - PNG/PDF visually consistent
  - claim reads clearly; only minor residual crowding in geometry panels

### Gemini CLI QC
- status: `PASS`
- artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-20260319-122200.md`
- summary:
  - panel alignment stable
  - fonts readable at two-column scale
  - no meaningful PNG-vs-PDF mismatch
  - structural-vs-CFD contrast reads clearly from the figure alone
