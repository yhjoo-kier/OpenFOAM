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

## 15. 2026-03-21 late-night rerender — outside legend + bottom horizontal colorbar

A new declutter pass was applied on the same representative pair and the same scale-calibrated source, still only for Figure 10.

Exact changes on this latest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and same scale-calibrated evaluation root
- moved the legend out of the data region into a dedicated bottom-left lane and moved the shared colorbar to a dedicated bottom horizontal lane under the CFD columns
- changed the per-row metric line from one long strip to a shorter multi-line block in the left-panel title area to reduce horizontal crowding
- slightly increased row-title / metric / colorbar tick typography and increased bottom-lane height for print safety
- retained direct VTK slice rendering with shared `viridis` normalization; Arial-first was requested again, but the host still fell back to `DejaVu Sans`

Exact current self-QC state on this latest rerender:
- verdict: `FAIL`
- main blockers on the exact current asset:
  - top-left header stack still collides with the upper row title area, and the collision shifts slightly between PNG and PDF renderings
  - bottom legend + horizontal colorbar no longer occlude panel interiors, but the combined bottom lane is still too tight for conservative print and remains visually busy
  - legend text (`dashed obstacle mask`) and per-row metric text are still borderline small at strict double-column print scale
  - the structure-vs-CFD message is cleaner than before, but the weaker CFD fidelity is still not sufficiently self-evident without reading the metric text and colorbar carefully
  - PDF vs PNG still show small typography placement/weight differences, so exact cross-format stability is not fully closed yet

Decision on this latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- External multimodal QC was not advanced on this exact asset because it still failed the self-gate.
- Highest-leverage next fix is now to simplify the top lane (remove or relocate the top-left header entirely), shorten the legend vocabulary, and either demote the metric text into a cleaner side/bottom strip or make the visual CFD mismatch itself more explicit with one stronger cue.

## 16. 2026-03-21 latest rerender — legend-free CFD overlay simplification

A new simplification pass targeted only the exact blockers from section 15 while staying on the same representative pair and same scale-calibrated source.

Exact changes on this latest asset:
- removed the separate bottom legend lane entirely, so the bottom row is now only a shared horizontal colorbar
- removed the top-left aggregate header (`A4 dense composite`) to free the top lane for panel content only
- changed the CFD obstacle cue from dashed white outlines to semi-transparent obstacle-colored overlays with white edges, so the geometry-to-flow correspondence is more intuitive and less legend-dependent
- kept the direct VTK mid-plane slice rendering and shared `viridis` normalization
- rewrote the row metadata into a lighter two-line header (`(a) A4-02`, metrics line below) placed above each geometry panel instead of using the old crowded title stack
- Arial was still requested first, but the host again fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` remain unavailable

Exact current self-QC state on this latest rerender:
- verdict: `FAIL`
- main blockers on the exact current asset:
  - the top lane is materially cleaner, but the per-row metric line is still too small / too long for strict conservative double-column print
  - the figure is less legend-dependent than before, yet the structure-vs-CFD gap is still not fully self-evident without reading the `CFD` numbers in the row header
  - semi-transparent obstacle overlays are more intuitive than dashed masks, but they still partially veil the flow field in the predicted panels
  - overall panel sizing is improved, but publication-grade message hierarchy is not fully closed yet, so this exact asset should not advance to external QC

Decision on this exact latest rerender:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- External multimodal QC was intentionally not advanced because the exact current asset still fails the self-gate.
- Highest-leverage next fix is now to shorten or relocate the row metrics and to make the CFD mismatch read more directly from the panel bodies rather than from header text.

## 17. 2026-03-21 final afternoon pass — message-first header simplification

A new self-only rerender tested a more message-first hierarchy on the same locked representative pair and same scale-calibrated source.

Exact changes on this newest asset:
- rewrote the top lane from neutral column names into a direct contrast pair: `high structural match` vs `weaker CFD fidelity`
- shortened the old long row-metric sentence into two compact badges only (`struct. ...`, `CFD ...`)
- removed the separate `openings kept` line from the latest exported asset after it proved non-essential and clutter-prone during self review
- slightly widened the geometry column and reduced CFD obstacle overlay opacity so the predicted flow is less veiled
- regenerated PNG, PDF, and PDF-render check on the exact same script/output paths

Exact current self-QC state on this newest asset:
- verdict: `FAIL`
- main blockers on the exact current asset:
  - the figure-level message is clearer than the prior legend-free revision, but the row label / badge lane above the geometry panels is still cramped and currently collides with the case-label area
  - the structure-vs-CFD gap is more legible at first glance, yet it still depends too much on the explicit `CFD` score badges rather than reading decisively enough from the panel bodies alone
  - PDF and PNG remain broadly consistent, but the remaining top-lane crowding is still too risky for conservative manuscript print

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- External multimodal QC was intentionally not advanced because the exact newest rerender still fails the self-gate.
- Highest-leverage next fix is now to redesign the row-identification lane itself (for example, move case identity fully into caption/panel labels and reserve the top-of-row strip only for short score badges), or to switch to a layout where the CFD mismatch reads without any score badges.

## 18. 2026-03-21 late-afternoon simplification pass — neutral headers + left-margin row IDs

A cleaner exact-asset rerender was built by removing the entire message/badge lane and letting the panels carry the contrast more directly.

Exact changes on this newest asset:
- removed the figure-level `high structural match` / `weaker CFD fidelity` headline and the two row score badges entirely
- switched to three plain column headers only: `predicted geometry`, `reference |U|`, `predicted |U|`
- moved row identifiers (`(a) A4-02`, `(b) A4-04`) into the left margin so the top of each row is no longer stacked with badges
- kept the shared bottom colorbar and verified PNG/PDF consistency by rendering the PDF back to `figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png`
- reduced predicted-panel obstacle overlay opacity further (`alpha=0.14`) so the CFD field stays more legible
- exact review copies for Gemini were staged at:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/fig10_exact_asset_20260321.png`
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/fig10_exact_asset_20260321_pdf_render.png`

Exact current QC split on this newest asset:
- self visual QC: `PASS`
  - the top lane is materially calmer and no longer competes with the panel bodies
  - the structure-vs-CFD contrast now reads mainly from the panel stack rather than from score badges
  - PNG/PDF layout consistency is acceptable on the exact rendered pair, with no new clipping or bottom-lane collisions
- independent subagent multimodal QC: `NOT RECOVERED`
  - an in-turn subagent QC attempt was made, but the verdict could not be recovered under the current cron-session visibility limits, so this gate remains open
- Gemini multimodal QC: `FAIL`
  - exact artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-150355.md`
  - Gemini accepted that the core message is visible and PNG/PDF are consistent, but still flagged two blockers on the exact asset: perceived row-label overlap near the left geometry panels and a row-alignment / panel-size mismatch between `reference |U|` and `predicted |U|`

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Next pass should target only the two surviving external concerns: make the left-margin row IDs even less invasive, and verify/normalize the apparent row-wise panel height alignment between the middle and right CFD columns before re-running external QC.

## 19. 2026-03-21 evening exact-asset rerender — row labels moved above panels + shared row extents

A narrow follow-up pass targeted only the two surviving blockers from section 18 on the same locked representative pair and same scale-calibrated source.

Exact changes on this newest asset:
- moved the row identifiers from the left margin into a small above-panel position over each geometry panel (`(a) A4-02`, `(b) A4-04`) so they no longer visually intrude into the data region
- unified the row-wise plotting extents across predicted geometry / reference `|U|` / predicted `|U|` by using the same max `(x, y)` room bounds within each row for all three panels
- reduced the left outer margin and regenerated PNG / PDF / PDF-render check on the exact same output paths
- font request remained Arial-first, but the host still fell back to `DejaVu Sans`

Exact current QC split on this newest asset:
- self visual QC: `PASS`
  - row identifiers are now cleanly separated from the geometry panels
  - PNG and PDF-render remain visually consistent on the exact exported pair
  - the previous apparent middle/right panel-size mismatch is materially reduced by forcing shared row extents
- independent subagent multimodal QC: `FAIL`
  - exact verdict: row labels no longer interfere and PNG/PDF consistency is good, but the figure still does not make the `structure vs CFD gap` sufficiently self-evident without caption help; the left geometry column still feels underweighted relative to the two CFD panels
  - session: `agent:main:subagent:30224e04-f13a-44c6-bb30-22b12f4b295f`
- Gemini multimodal QC: `NO VALID EXACT-ASSET VERDICT YET`
  - a local Gemini CLI rerun was attempted and artifact was captured at `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-152259.md`
  - however, the returned critique referenced stale/removed elements (`global legend`, `aggregate metrics`) that are not present in the exact current asset, so this output is not reliable enough to count as the required exact-asset Gemini gate

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is no longer label placement; it is now message strength: make the structure-vs-CFD contrast read more directly from the panel bodies (or rebalance geometry vs CFD visual weight) before re-running both external QC gates.

## 20. 2026-03-21 late-evening exact-asset rerender — geometry reweight + lighter CFD overlays

A narrow, exact-asset rerender targeted the remaining message-strength blockers while preserving the same overall concept and the same locked representative pair.

Exact changes on this newest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root / manifest
- increased the geometry column share (`width_ratios=[1.24, 1.34, 1.34]`) so the left structure evidence no longer reads as a visibly minor side lane
- thickened/darkened the geometry room + obstacle strokes slightly for better print-scale salience
- reduced predicted-CFD obstacle overlay opacity further (`alpha=0.04`) and softened the halo so flow peaks remain more readable under the obstacle cue
- enlarged column-header / row-label / colorbar typography slightly and thinned the colorbar outline for cleaner PNG/PDF export
- regenerated the exact current outputs at the same paths:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.png`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png`

Exact current QC split on this newest asset:
- self visual QC: `PASS`
  - after the geometry-column reweight and lighter obstacle overlays, the three-column stack reads more directly as `structure evidence` vs `reference/predicted CFD evidence`
  - no new label collisions or trim issues were observed on the exact PNG/PDF-render pair
  - remaining self note is minor only: the left geometry column could still carry slightly more salience, but it is no longer the main blocking issue
- independent subagent multimodal QC: `FAIL`
  - exact session: `agent:main:subagent:51ba40cb-a052-48c8-a06d-7d2e553bf3f7`
  - main blockers on the exact current asset:
    - the intended causal read (`dense composite with high structural match yet weaker CFD fidelity`) is still not fully self-evident without caption help; first glance still reads more generically as `prediction differs`
    - geometry column remains somewhat underweighted relative to the two CFD panels, even after the latest rebalance
    - colorbar max tick looked too close to the right edge / borderline clipped in the external reading
    - colorbar typography remains near the lower bound for conservative double-column print
- Gemini multimodal QC: `PENDING / NO VALID EXACT-ASSET VERDICT YET`
  - one Gemini rerun inside an external workspace produced a stale-path critique referencing removed elements and therefore does **not** count as a valid gate (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-153818.md`)
  - a corrected repo-local exact-asset Gemini rerun was launched against the current outputs, but no usable verdict was recovered within this turn

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now message explicitness rather than generic cleanup: either strengthen the geometry-side cue that the structure is close to reference, or add one very small caption-safe in-panel cue that makes the `structure ok / CFD still off` contrast more self-evident without reintroducing clutter.

## 21. 2026-03-21 final exact-asset rerender — reference-outline revival pass

A narrow, self-only rerender was applied on the same locked representative pair and same scale-calibrated source to make the structural-match side more self-evident without reintroducing the old score-badge lane.

Exact changes on this newest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root / manifest
- revived the geometry-side reference cue by adding short `gray dashed` reference room/obstacle outlines directly into the predicted-geometry panels
- rebalanced the three columns slightly toward parity (`width_ratios=[1.30, 1.29, 1.29]`) so the left geometry lane no longer grows relative to the CFD evidence
- reduced geometry stroke heaviness modestly and simplified the shared colorbar ticks to the two endpoints only to lower bottom-lane crowding
- regenerated the exact current outputs at the same paths:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap.png`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png`

Exact current self-QC state on this newest asset:
- verdict: `FAIL`
- basis on the exact PNG/PDF-render pair:
  - the revived dashed reference outline is directionally helpful, but it is still not consistently strong enough to make `structure is close to reference` read immediately at print scale
  - a new cross-format issue appeared in the top lane: PNG and PDF-render do not match cleanly enough around the left column header / top spacing, so exact export stability is not closed
  - the bottom colorbar is cleaner than before, but the endpoint tick labels still sit too close to the bar edges in the PDF-rendered check
  - the overall message remains only partly self-evident; without caption help the figure still reads too much like a generic prediction comparison rather than a decisive `structure plausible / CFD weaker` contrast

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- External QC was intentionally **not** advanced on this exact asset because the rerender still fails the self-gate.
- Highest-leverage next fix is now a more structural redesign rather than another micro-pass: likely either (a) add a caption-safe direct CFD mismatch cue (difference map or tightly targeted mismatch annotation), or (b) switch the left column to a more explicit predicted-vs-reference structural comparison panel so the `structure close` claim stops relying on faint dashed outlines.

## 22. 2026-03-21 late-night exact-asset rerender — direct message headers + exact external QC refresh

A narrow rerender targeted the two remaining themes from the prior pass: message explicitness and colorbar / header print safety, while staying on the same locked representative pair and the same scale-calibrated source.

Exact changes on this newest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root / manifest
- increased the geometry-column share again (`width_ratios=[1.38, 1.24, 1.24]`) so the left structural evidence reads less like a side lane
- replaced the neutral top labels with a more direct, caption-safe contrast pair: `close geometry` / `reference |U|` / `shifted |U|`
- rewrote the left subtitle from `gray dashed = ref. outline` to the shorter legend-safe cue `solid pred. · dashed ref.`
- strengthened the dashed reference outlines slightly and pushed the colorbar end tick alignment inward by adjusting tick-label anchoring
- regenerated PNG / PDF / PDF-render check on the exact same output paths; host font fallback remains `DejaVu Sans` because `Arial` / `Liberation Sans` were unavailable

Exact current QC split on this newest asset:
- self visual QC: `PASS`
  - the intended read is more direct than the prior neutral-header version, and the endpoint colorbar labels no longer crowd the bar edges
  - PNG and PDF-render remain materially consistent on the exact exported pair
  - no new clipping was observed on the exact current asset
- independent subagent multimodal QC: `FAIL`
  - exact latest verdict: the figure still does not make `close geometry but shifted CFD` sufficiently self-evident without caption help; the smallest header / legend / colorbar text remains weak for conservative double-column print; and the top header/subtitle lane still feels cramped against the geometry panels
- Gemini multimodal QC: `FAIL`
  - valid exact-current artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-161354.md`
  - exact latest blockers recovered from the usable terminal verdict on the exact-root review copies: left header / subtitle spacing is still too tight, the first-column header block is vertically misaligned against the other two column titles, and Gemini still does not accept the current geometry-to-reference comparison as publication-safe

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now no longer micro-typography; it is a column-1 redesign that makes the structural comparison explicit without a cramped subtitle (for example, a more direct predicted-vs-reference structural comparison panel or a different header architecture with the first-column cue moved out of the title lane).

## 23. 2026-03-21 16:23 exact-current rerender — explicit reference/predicted geometry pair

A new exact-current rerender replaced the old single geometry-overlay column with an explicit side-by-side geometry comparison, while keeping the same locked representative pair and the same scale-calibrated source.

Exact changes on this newest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root / manifest
- rewrote the layout from `2 x 3` into `2 x 4`: `ref. geometry` / `pred. geometry` / `reference |U|` / `predicted |U|`
- removed the cramped first-column subtitle lane entirely; structural similarity is now shown directly by the two left panels rather than by dashed overlays
- kept the shared bottom colorbar only under the two CFD columns and rechecked exact PNG/PDF consistency with `pdftoppm`
- Arial was still requested first, but the host again fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` were unavailable
- exact review copies used for external QC:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/fig10_exact_asset_20260321_v2.png`
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/fig10_exact_asset_20260321_v2_pdf_render.png`

Exact current QC split on this newest asset:
- self visual QC: `FAIL`
  - the figure message is more direct than the prior overlay-based version because geometry closeness is now read from a literal `ref.` vs `pred.` pair
  - PNG/PDF layout consistency is good and no clipping/collision was observed on the exact exported pair
  - however, for conservative double-column print the top column headers and colorbar endpoint ticks remain too small, and the geometry panels are still visually underweighted relative to the CFD panels; the message is clearer but not yet manuscript-safe enough
- independent subagent multimodal QC: `FAIL`
  - exact session: `agent:main:subagent:ab76e513-1b83-4949-8179-bac472c212e1`
  - exact verdict: the message is only partly immediate at print scale; the geometry panels and headers are too small relative to the CFD panels; and the colorbar ticks / fine contour detail remain below a conservative publication-safe threshold
- Gemini multimodal QC: `PASS`
  - exact artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-162350.md`
  - Gemini accepted the new 2x4 structure as caption-free readable, well-aligned, and PNG/PDF-consistent

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now a print-scale rebalance rather than concept selection: enlarge the two geometry panels and the top/colorbar typography without sacrificing the explicit four-column comparison.

## 24. 2026-03-21 16:34 exact-current rerender — print-scale rebalance pass

A narrow follow-up pass stayed on the same explicit `2 x 4` concept and the same locked representative pair, but targeted only conservative print readability.

Exact changes on this newest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root / manifest
- increased overall figure height slightly and rebalanced the four columns closer to parity (`width_ratios=[1.16, 1.16, 1.18, 1.18]`) so the two geometry panels no longer lag as much behind the CFD pair
- enlarged the top column headers / row labels / shared colorbar typography and thickened geometry room + obstacle strokes for better print survival
- lowered the top header lane and added a bit more bottom space for the shared colorbar while preserving the same output paths
- Arial was still requested first, but the host again fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` remain unavailable

Exact current self-QC state on this newest asset:
- verdict: `FAIL`
- basis on the exact PNG/PDF-render pair:
  - the concept remains the right one, but the geometry pair is still visually weaker than the CFD pair, so `structure close / CFD shifted` is not yet self-evident enough at first glance
  - top/header and colorbar readability improved, but PNG vs PDF-render still show small vertical placement drift, especially in the header-to-panel gap and bottom colorbar position
  - overall vertical whitespace remains slightly uneven, with too much empty space above the bottom colorbar lane relative to the panel stack
  - a temporary colorbar-end spacing experiment produced an incorrect endpoint-tick read during self review and was reverted; the figure remains in FAIL state and external QC was intentionally not advanced on this exact asset

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now still concept-preserving but more structural: either further strengthen the geometry-side contrast (stronger edge/fill hierarchy or less dominant CFD rendering) or redesign the bottom colorbar lane so PDF/PNG export stability closes before external QC rerun.

## 25. 2026-03-21 16:56 exact-current rerender — geometry-up / colorbar-fix print-safety pass

A narrow print-safety rerender stayed on the same explicit `2 x 4` concept and the same locked representative pair, and targeted the two most obvious remaining issues from the prior exact asset: the underweighted geometry pair and the misleading/too-coarse shared colorbar ticks.

Exact changes on this newest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root / manifest
- increased overall figure height slightly and shifted column balance toward the two geometry panels (`width_ratios=[1.26, 1.26, 1.10, 1.10]`) so the left structural comparison reads less like a side lane
- tightened outer margins / gutters slightly while giving the bottom shared colorbar lane a bit more height for print safety
- enlarged top headers / row labels / colorbar label modestly and replaced the old endpoint-only `%.1f` colorbar ticks with small-range-aware ticks (`0.000`, midpoint, `vmax`) so the scale no longer collapses to visually duplicated `0.0` endpoints on low-velocity slices
- regenerated the exact current outputs at the same paths and restaged exact-review copies at:
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/fig10_exact_asset_20260321_v3.png`
  - `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/fig10_exact_asset_20260321_v3_pdf_render.png`
- Arial was still requested first, but the host again fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` were unavailable

Exact current QC split on this newest asset:
- self visual QC: `FAIL`
  - the colorbar-value bug is fixed and the geometry pair is stronger than in the prior exact asset
  - however, on strict double-column review the shared colorbar numbers / label still feel borderline small, and the PNG vs PDF-render pair still shows slight top-lane / colorbar-lane spacing drift that keeps the exact self-gate conservative
  - the figure message is now clearer, but I still do not want to over-call it manuscript-safe on self review
- independent subagent multimodal QC: `PASS`
  - exact session: `agent:main:subagent:cf70134f-8bdb-4c63-9280-12092455abac`
  - returned verdict: `PASS`; row/column message was judged immediately readable and PNG/PDF consistency acceptable on the exact asset
- Gemini multimodal QC: `PASS`
  - exact artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-165608.md`
  - returned verdict: `PASS`; exact-current prompt referenced only the repo-local v3 review copies and Gemini accepted both layout consistency and caption-free message readability

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now very narrow: either make the shared colorbar typography decisively print-safe or simplify/remove it if the panel-to-panel contrast already carries the message.

## 26. 2026-03-21 17:xx exact-current rerender — spacing/colorbar print-safety micro-pass

A final narrow micro-pass stayed on the same explicit `2 x 4` concept and the same locked representative pair, and touched only print-safety details that were still bothering the self gate.

Exact changes on this newest asset:
- kept the same representative pair (`bench_a4_02/floorplan`, `bench_a4_04/floorplan`) and the same scale-calibrated evaluation root / manifest
- slightly increased the geometry-column share again (`width_ratios=[1.28, 1.28, 1.08, 1.08]`) and strengthened the geometry room/obstacle strokes plus inlet/outlet markers
- tightened the top lane by lowering the column headers and row labels, while slightly increasing the shared colorbar label/tick typography
- increased the bottom colorbar lane height a bit and regenerated the exact PNG / PDF / PDF-render check directly from the script (the script now also refreshes `figure10_dense_composite_structure_vs_cfd_gap_pdf_render.png` automatically via `pdftoppm`)
- Arial was still requested first, but the host again fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` were unavailable

Exact current QC state on this newest asset:
- self visual QC: `FAIL`
  - the PDF-render exact check looks materially stable and the shared colorbar is safer than before, but I am keeping the self gate conservative because the geometry pair still reads weaker than the CFD pair at strict journal print review
  - top/header spacing no longer looks broken, yet the overall top lane still feels a bit heavy for the amount of evidence it carries
  - the figure is close, but not clearly enough above the manuscript-safe line to over-call it complete
- external multimodal QC on this exact asset: `NOT ADVANCED`
  - I intentionally did **not** advance subagent/Gemini exact-current QC on this rerender because the self gate did not close cleanly
  - quick local visual review of the PDF-rendered exact asset suggests the remaining blocker is no longer PNG/PDF drift, but rather geometry-vs-CFD visual balance under conservative print assumptions

Decision on this newest asset:
- Figure 10 remains incomplete.
- Do **not** mark the checklist complete.
- Do **not** copy to Google Drive.
- Highest-leverage next fix is now to rebalance evidence rather than typography: either enlarge/salience-boost the geometry pair one more step, or demote/simplify the shared colorbar so the left-side structural comparison carries more weight without needing header support.
