# Figure 8 QC Log — Section-view composite-room collapse (v2 rerun, scale-calibrated main-setting)

- Date: 2026-03-20
- Figure: 8
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure8_section_collapse.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure8_section_view_composite_collapse.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure8_section_view_composite_collapse.png`
- PDF render check:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/figure8_section_view_composite_collapse_pdf_render.png`
- Gemini input image:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/figure8_section_view_composite_collapse_gemini_input.png`
- Export meta:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure8_section_view_composite_collapse_meta.json`

## 1. Source artifact / case lock

Figure 8 was refreshed under the v2 rule that the manuscript main setting must use the single-factor scale-calibrated pipeline.

- scale-calibrated metric source:
  - `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
- geometry / scene source artifacts:
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a3_04/section/reference_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a3_04/section/predicted_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a3_04/section/evaluation_summary.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a4_05/section/reference_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a4_05/section/predicted_scene.json`
  - `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/bench_a4_05/section/evaluation_summary.json`
- representative cases were re-locked to the only two `section_room_kind_collapse` examples in the scale-calibrated aggregate manifest:
  - `bench_a3_04 / section`
  - `bench_a4_05 / section`

Metric box values now use the scale-calibrated aggregate entries (`S=0.417`, `CFD=0.334` and `0.327`) instead of the older no-scaling CFD values.

## 2. Layout / export decision

- intended width: `double-column`
- panel layout: `2 x 2`
  - `(a)` A3-04 reference
  - `(b)` A3-04 prediction + GT outline
  - `(c)` A4-05 reference
  - `(d)` A4-05 prediction + GT outline
- export size: `7.25 in x 5.20 in`
- PNG dpi: `600`
- PDF backend: matplotlib vector PDF
- font family request: `Arial -> Liberation Sans -> DejaVu Sans`
- selected font in this environment: `DejaVu Sans`
- internal caption / figure-number rule:
  - confirmed no `Figure 8`, `Fig. 8`, or caption-style sentence exists inside the artwork

## 3. Rerun changes vs old 2026-03-19 version

- removed the forbidden in-figure headline/subtitle that explicitly said `Figure 8 ...`
- refreshed metric sourcing so the figure aligns with the scale-calibrated main-setting rule
- kept the same two representative section-collapse cases because they remain the only room-kind collapse tasks in the scale-calibrated aggregate
- added Arial-first fallback logging via the meta file
- regenerated PDF/PNG and rechecked PNG/PDF consistency

## 4. QC status in this turn

### Self visual QC
- verdict: `PASS`
- method:
  - direct inspection of regenerated PNG
  - separate PDF raster render inspection via `pdftoppm`
- summary:
  - no forbidden in-figure number/caption prose remains
  - panel ordering and spacing are stable at double-column width
  - dashed GT outlines clearly expose the missing composite arm in both prediction panels
  - PNG and PDF render keep the same layout without reflow/clipping
  - only minor non-blocking polish remains (tight left gutter and metric-box crowding in dense panel)

### Subagent visual QC
- verdict: `PASS`
- session:
  - `agent:main:subagent:8632df27-f0b9-4399-ad41-53296439e768`
- summary:
  - no forbidden in-figure prose detected
  - PNG/PDF layout consistency acceptable; no actual clipping
  - message is clear, with dashed GT outlines and collapse badge readable at double-column scale
  - non-blocking polish: slightly more top margin, slightly stronger metric-box / dashed-outline print robustness, and less overlap from the metric box if space allows

### Gemini visual QC
- verdict: `PASS`
- artifact:
  - `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure8-visual-qc-20260320-235409.md`
- summary:
  - Gemini reported no blocking issues
  - collapse message was judged visually clear and manuscript-ready for double-column use
  - non-blocking polish: left y-axis gutter is a bit tight, panel (d) stats box could be cleaner, and the x-axis `0` tick overlap is cosmetic only

## 5. Final disposition for this turn

- overall figure completion verdict: `COMPLETE`
- reason:
  - all three QC gates passed on the v2 rerun (`self PASS`, `subagent PASS`, `Gemini PASS`)
- consequence:
  - Figure 8 can be marked complete in `/home/yhjoo/projects/OpenFOAM/docs/26-03-20_paper_figure_production_checklist_v2.md`
  - Figure 8 can be stored in Google Drive under `Fig08_section_collapse/`

## 6. Delivery status

- Google Drive target folder:
  - `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig08_section_collapse/`
- delivered files:
  - `figure8_section_view_composite_collapse.pdf`
  - `figure8_section_view_composite_collapse.png`
  - `qc_summary.md`
