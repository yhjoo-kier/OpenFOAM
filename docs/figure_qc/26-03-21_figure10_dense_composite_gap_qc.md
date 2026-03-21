# Figure 10 QC Log — Dense composite structure-vs-CFD gap (2026-03-21 refresh)

- Date: 2026-03-21
- Figure: 10
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure10_dense_composite_gap.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.png`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png`

## 1. Source artifact / case lock

Main setting was refreshed to the scale-calibrated manifest:
- aggregate source: `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
- evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span/`

Representative dense-composite cases remained:
1. `bench_a4_02 / floorplan`
   - structural = `0.9375`
   - CFD = `0.4086`
   - derived tag = `dense_composite_structure_physics_gap`
2. `bench_a4_04 / floorplan`
   - structural = `1.0000`
   - CFD = `0.4584`
   - derived tag = `dense_composite_structure_physics_gap`

Aggregate context used in the figure header:
- `A4` mean structural = `0.791`
- `A4` mean CFD = `0.458`
- room-kind match = `88%`
- opening-wall match = `72%`

## 2. Layout / export decisions

- intended width: `double-column`
- panel layout: `2 x 2`
  - `(a)` A4-02 geometry
  - `(b)` A4-02 flow
  - `(c)` A4-04 geometry
  - `(d)` A4-04 flow
- top column headers only: `geometry + input` / `reference vs predicted flow`
- forbidden in-figure `Figure 10 ...` headline and footer sentence were removed on this refresh
- export size: `7.35 in x 6.65 in`
- PNG dpi: `600`
- font family request: `Arial` first
- actual fallback on this host: `DejaVu Sans` (Matplotlib reported `Arial` and `Liberation Sans` unavailable)
- save path was changed away from `bbox_inches="tight"` so PNG and PDF keep a more stable outer frame

## 3. 2026-03-21 refresh summary

This refresh updated the figure from the older production candidate to the current scale-calibrated main-setting rule and tightened rule compliance:
- switched aggregate source from legacy manifest to `evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
- switched per-case evaluation root to `benchmark/evaluations_posthoc_scaled_longest_span/`
- removed the explicit `Figure 10` headline and the sentence-style footer
- shortened panel titles and callouts to rule-safe short annotations
- regenerated PDF / PNG / PDF-render check on the exact new revision

## 4. 2026-03-21 noon declutter rerender

A focused readability pass was applied to the same representative pair and same scale-calibrated source.

Changes on this exact revision:
- removed the extra structural-vs-CFD score-bar inset from the geometry panels
- shortened panel titles to bare case IDs (`(a) A4-02`, `(b) A4-02`, ...)
- simplified the geometry-side callouts to single short cues (`openings kept`, `strong structure`)
- moved the main CFD emphasis into larger flow-side badges (`CFD 0.409`, `CFD 0.458`) plus a short `moderate CFD` cue
- removed the geometry-panel background grid to reduce clutter
- simplified the top aggregate line to `A4 mean  struct. 0.791  CFD 0.458`
- kept the input thumbnail but slightly reduced its footprint relative to the main geometry overlay

## 5. Exact current QC state

### Self visual QC
- verdict: `PASS`
- basis:
  - left geometry panels are materially cleaner after removing the bar inset and background grid
  - the structural-vs-CFD contrast is more direct because structure stays on the left metric box while CFD is emphasized on the right flow badge
  - panel spacing / trim safety look acceptable on the exact regenerated PNG/PDF pair

### Independent subagent multimodal QC
- verdict: `FAIL`
- session: `agent:main:subagent:de501d4f-aa4a-4c07-9128-316c87707b68`
- main blockers recovered from the exact current revision:
  - small metric boxes, axis/tick labels, and several panel annotations are still hard to read at conservative double-column print scale
  - the reference-vs-predicted CFD panels still do not show a strong enough visual contrast on first glance, so `moderate CFD` still leans on the score badge
  - callouts, dashed overlays, inset thumbnails, axes, and metric text still compete with each other instead of producing a fully dominant message hierarchy
  - the caption-free takeaway is improved but still not immediate enough for the strict external gate
  - the CFD score labels are clearer than before but still not visually dominant enough to fully carry the figure on their own

### Gemini multimodal QC
- verdict: `PASS`
- artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-rerun-20260321-121015.md`
- accepted points:
  - orange callouts now communicate the main idea more directly
  - the 2x2 grid feels clean and trim-safe
  - metric boxes and CFD badges were considered readable on Gemini's gate
  - the structure-vs-flow gap is judged identifiable without caption help

## 6. Current decision

- Figure 10 is **not complete yet**.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.

## 7. Immediate next fixes

1. Push one more print-safety pass on the smallest geometry-panel text (`struct./obs/openings`, ticks, labels).
2. Make the flow-side `moderate CFD` contrast more visually obvious without depending mainly on the numeric badge.
3. Reduce one more annotation layer on the left if needed (likely the inset thumbnail or some axis/detail load).
4. Re-run the same three-way QC on the exact next rerender only.

## 8. 2026-03-21 afternoon simplification rerender

A further simplification pass was applied on the same representative pair and same scale-calibrated source.

Changes on this exact revision:
- geometry-panel axes/ticks and input thumbnails were removed to free print-scale space for the room/obstacle overlays
- geometry metric boxes were reduced to two-line badges (`struct.` + `openings kept/shifted`)
- top headers were rewritten to the direct contrast pair `strong structure` / `moderate CFD`
- flow-column width was increased (`width_ratios=[0.92, 1.48]`, `wspace=0.16`) and the inter-image gap was reduced to enlarge the CFD renders
- micro labels (`Reference`, `Predicted`, CFD badge, top aggregate line) were enlarged relative to the noon rerender

Exact current QC split on this latest rerender:
- self visual QC: `FAIL`
  - smallest labels and badges are still borderline for conservative double-column print
  - CFD-side emphasis remains weaker than desired, and the badge/label styling still does too much work relative to the render itself
  - PNG/PDF are close, but text/shadow rendering still differs enough to keep print conservatism high
- independent subagent multimodal QC: `FAIL`
  - `Reference` / `Predicted` microlabels and top summary remain too small/low-contrast
  - right-column CFD panels still read slightly underscaled versus the available whitespace, so the `moderate CFD fidelity` takeaway is not immediate enough without reading badges
  - PNG/PDF are materially consistent, but the consistency preserves the same readability weakness
- Gemini multimodal QC: `NOT YET RECOVERED` on this exact rerender
  - a local Gemini CLI rerun was attempted, but it did not return a usable verdict in this turn; previous Gemini PASS applies only to the prior noon rerender and must not be reused for this exact asset

Decision on this exact rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.

## 9. 2026-03-21 late rerender — 2D flow-slice pivot

A more structural simplification pass was applied on the same representative pair and the same scale-calibrated source, but the **flow-side evidence asset** was changed from the previous 3D contour render crop to the `comparison_1x2.png` CFD-slice crop so the right-column message relies less on tiny 3D details.

Exact changes on this latest asset:
- retained the same representative pair: `bench_a4_02/floorplan` and `bench_a4_04/floorplan`
- retained the same scale-calibrated aggregate source: `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
- kept the left geometry panels minimal (no axes/ticks/input thumbnail) and preserved only the short structural badge
- replaced the right flow panels with tighter crops from each case's `comparison_1x2.png` so the 2D velocity slice pattern is visually larger and more direct than the prior 3D view
- kept only short in-panel labels (`Reference`, `Predicted`, `CFD 0.409`, `CFD 0.458`) and preserved the caption-safe top headers `strong structure` / `moderate CFD`
- slightly reduced overall figure height and increased the right-column width ratio to enlarge the flow evidence inside the same double-column footprint

Exact current QC split on this latest rerender:
- self visual QC: `FAIL`
  - the right-column evidence is materially clearer than the earlier 3D version, but there is still too much unused whitespace between quadrants relative to actual panel area
  - top summary text remains conservative-print fragile, and the whole figure still shrinks the evidence more than necessary for a double-column slot
  - the CFD badges are cleaner now, but they still occupy the lower flow field rather than reading as fully non-invasive support
- independent subagent multimodal QC: `FAIL`
  - exact verdict: enlarge all evidence panels and text by roughly 20–30% while tightening whitespace; the structure-vs-CFD gap is only partly self-evident without caption help
  - main blockers: smallest text still fragile, message still partly badge-dependent, and contour panels will lose legibility margin after journal downscaling
- Gemini multimodal QC: `FAIL`
  - exact artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-rerun-20260321-124412.md`
  - main blockers: excessive negative space, top summary/tick labels too small, CFD badges obscure important flow-field area, and subfigure labels feel too detached from the panel bodies

Decision on this latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Next likely fix: aggressively tighten the 2x2 spacing/margins first, then retest before changing content again.

## 10. 2026-03-21 late-afternoon rerender — tighter spacing + off-data badges

A follow-up rerender targeted the exact blockers from the previous FAIL split, still on the same representative pair and same scale-calibrated source.

Exact changes on this latest asset:
- reduced outer margins and inter-panel gutters to scale up the four evidence panels within the same double-column footprint
- moved the structural / CFD numeric badges out of the flow-field interiors and into below-panel figure text so the CFD maps are less occluded
- moved the `(a)–(d)` case labels out of the panel interiors into figure-level labels above each panel
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same manifest root (`evaluation_aggregate_summary_posthoc_scaled_longest_span.json`)
- kept Arial-first font request, but the host still fell back to `DejaVu Sans`

Exact current self-QC state on this latest rerender:
- verdict: `FAIL`
- main blockers:
  - moving panel labels outside removed direct data occlusion, but the outside label stack now collides locally with `Reference` / `Predicted` headings near the flow panels
  - the right-column CFD evidence still does not have enough visual weight relative to the geometry column, so the structure-vs-CFD gap is still not immediate enough without caption support
  - a shared color scale / colorbar is still missing, which weakens first-glance interpretability of the CFD comparison
  - typography is cleaner than before, but label hierarchy is still inconsistent and some text sits too close to panel borders for conservative print safety

Decision on this latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix now is a layout rebalance that gives the CFD panels clearly more area and introduces a shared color scale, rather than further micro-tweaking badge placement.

## 11. 2026-03-21 evening rerender — 3-column simplification (geometry / reference CFD / predicted CFD)

A more aggressive simplification pass was applied on the same representative pair and the same scale-calibrated source.

Exact changes on this latest asset:
- rewrote the layout from the previous 2x2-with-combined-flow structure into a `2 x 3` grid: predicted geometry / reference CFD / predicted CFD
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same aggregate source (`evaluation_aggregate_summary_posthoc_scaled_longest_span.json`)
- removed the separate outside panel-label stack that had collided with `Reference` / `Predicted` in the prior revision
- replaced per-panel CFD badges with one per-row metric footer (`struct. ... · openings kept · CFD ...`) so the flow fields themselves stay unobscured
- retained Arial-first font request, but the host again fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` remain unavailable
- regenerated both PNG and PDF, then rendered the PDF back to PNG for exact layout consistency checking

Exact current QC split on this latest rerender:
- self visual QC: `FAIL`
  - the core message is clearer than the prior combined-flow revision, but the exact crop still leaves stray CFD-source axis/title fragments in the center/right panels, which reads as unpolished at manuscript quality
  - the structure-vs-CFD contrast is visible, yet the figure still does not explain the blue dashed reference geometry / opening cues quickly enough without caption support
  - the bold per-row metric footer is readable but now competes too strongly with the CFD panels
- independent subagent multimodal QC: `FAIL`
  - exact verdict: the story is only partly inferable without caption; bottom-row center CFD source title remains clipped/crowded; smallest text remains marginal at double-column downscaling; PNG/PDF layout is consistent but PDF-rendered CFD panels look slightly softer
  - session: `agent:main:subagent:f7963117-c6a6-4355-99ae-4f4a223b9583`
- Gemini multimodal QC: `FAIL`
  - exact artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-rerun-20260321-130913.md`
  - main blockers: broken / inconsistent CFD axis labeling on the cropped panels, clipped internal CFD-source title on the lower reference panel, excessive vertical whitespace, and unexplained geometry-side dashed/colored cues that still lean on caption knowledge

Decision on this latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Next likely fix: stop reusing the pre-rendered `comparison_1x2.png` crops and instead regenerate clean CFD slice panels directly from the underlying data so titles/axes/color scale can be controlled consistently.

## 12. 2026-03-21 evening rerender — cleaner slice crop + simplified geometry cue

A new rerender tested the most conservative simplification path on the same representative pair and the same scale-calibrated source.

Exact changes on this latest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same aggregate source (`evaluation_aggregate_summary_posthoc_scaled_longest_span.json`)
- removed the dashed reference overlays from the geometry panels and showed only the predicted geometry / obstacles / inlet-outlet markers to avoid unexplained blue cue lines
- replaced the previous dirty CFD crop with a wider crop from each `comparison_1x2.png` right-half so axis labels and the per-panel colorbar survive cleanly while the long internal title is cropped away
- kept the `2 x 3` structure (`predicted geometry / reference CFD / predicted CFD`) but widened the two CFD columns and tightened inter-panel spacing
- moved the row metrics into a lighter rounded badge above each row (`struct. ... · openings kept · CFD ...`) instead of the previous heavy footer
- retained Arial-first font request, but the host still fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` remain unavailable

Exact current self-QC state on this latest rerender:
- verdict: `FAIL`
- main blockers:
  - conservative double-column print readability is still not closed because the CFD-axis ticks / colorbar ticks inherited from the source images remain too small after journal downscaling
  - the row-metric badge still competes with the panel stack and locally crowds the top border area instead of disappearing into the layout
  - the CFD evidence is cleaner than the prior clipped crop, but there is still no trustworthy shared normalization visible across all four CFD panels, so quantitative comparison remains visually ambiguous
  - geometry-side cues are simpler than before, but the obstacle/opening encoding is still not fully self-explanatory without either a legend or a stronger in-figure cue system
  - panel-area usage improved only partially; the layout is cleaner, yet still does not feel decisively publication-ready

Decision on this latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now to regenerate the CFD slice panels from underlying simulation data (or another title-free source asset) so label size, shared color normalization, and colorbar design can be controlled directly instead of inherited from `comparison_1x2.png`.

## 13. 2026-03-21 night rerender — direct VTK slice + shared colorbar

A fresh rerender was built by abandoning `comparison_1x2.png` reuse and drawing the CFD panels directly from each case's VTK field.

Exact changes on this latest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and same scale-calibrated source / evaluation root
- rewrote `scripts/paper_figures/make_figure10_dense_composite_gap.py` so the reference/predicted CFD panels are rendered from `internal.vtu` mid-plane slices via PyVista rather than cropped raster screenshots
- introduced one shared `turbo` normalization / colorbar across all four CFD maps and removed the inherited per-panel OpenFOAM titles, axis ticks, and per-panel colorbars
- kept the `2 x 3` layout (`predicted geometry / reference CFD / predicted CFD`), retained only short row labels plus one compact metric line per row, and added a tiny geometry legend (`obstacle`, `inlet`, `outlet`)
- generated a fresh PDF and a PDF-rendered PNG check (`figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png`) for exact cross-format QC
- typography still requested `Arial` first, but both host and container again fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` were unavailable

Exact current QC split on this latest rerender:
- self visual QC: `FAIL`
  - the core story is much clearer than the prior crop-based variants, but the smallest remaining text still does not clear conservative journal print margins
  - shared-colorbar tick labels are too tight/colliding, and the top summary / legend cluster still feels cramped relative to the available header space
  - geometry-side context is still undersized versus the CFD evidence, so the structure side of the contrast is weaker than it should be
- independent subagent multimodal QC: `FAIL`
  - exact verdict: top summary labels, per-row metric text, legend, and colorbar text remain too small/light; geometry panels are undersized; dashed white overlays may lose contrast in print; PDF render is slightly thinner/softer than PNG
  - session: `agent:main:subagent:e6cd1ab1-258f-41e7-9499-409f7af5d198`
- Gemini multimodal QC: `FAIL`
  - exact artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-current-20260321-134113.md`
  - main blockers: colorbar tick-label collision, cramped top legend/header alignment, and row-metric text that remains risky at double-column downscaling

Decision on this latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now a typography/layout pass: enlarge or simplify the top legend/header cluster, enlarge geometry panels (or shrink whitespace), thicken low-contrast dashed overlays, and redesign the shared colorbar ticks so PDF/PNG remain equally readable.

## 14. 2026-03-21 night rerender — larger panels + in-panel legend + print-safe colormap attempt

A fresh typography/layout pass was applied on the same representative pair and same scale-calibrated source, still targeting the locked Figure 10 brief only.

Exact changes on this latest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root
- widened the geometry column (`width_ratios=[1.18, 1.52, 1.52]`) and tightened row/column gaps to enlarge all six evidence panels inside the same double-column footprint
- replaced the previous `turbo` shared CFD colormap with `viridis` and simplified the shared colorbar to four ticks in a right-side vertical bar
- moved the global legend into the first geometry panel and added an explicit `pred. obstacle` key for the dashed CFD overlay cue
- enlarged row labels / per-row metric text and thickened the CFD obstacle overlays using a dark halo path effect for better print contrast
- retained Arial-first font request, but this host still fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` remain unavailable

Exact current self-QC state on this latest rerender:
- verdict: `FAIL`
- main blockers:
  - the in-panel legend solved the floating-header issue but now directly occludes the top-left geometry panel and makes the left column feel visually unbalanced
  - the per-row metric line remains too long for the available width and still rides too close to the panels for conservative double-column print safety
  - the shared colorbar is cleaner than before, but its label/ticks are still too small relative to journal reduction and the CFD-vs-geometry message still depends on reading supporting text
  - PDF-rendered text/panel sharpness remains slightly softer than the PNG version, so the exact asset is not yet print-safe enough

Decision on this latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- External QC was intentionally not advanced on this exact asset because the new rerender did not clear the self-gate.
- Highest-leverage next fix is likely to move the legend completely outside the data region, compress the per-row metrics into a smaller dedicated strip or footer, and give the colorbar/headers one cleaner typography lane instead of stacking them around the panel grid.
