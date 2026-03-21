# Figure 9 QC Log — obstacle hallucination with limited CFD penalty

> Date: 2026-03-21
> Figure: 9
> Status: IN PROGRESS (self QC PASS, subagent QC FAIL, Gemini QC FAIL)

## Source artifact / case lock (current revision)

Scale-calibrated manifest:
- `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`

Locked representative tasks for the current revision:
1. `bench_a3_01 / wireframe`
   - structural score: `0.750`
   - CFD score: `0.604`
   - obstacle count: `0 -> 3`
   - opening walls: `N/W -> N/W`
2. `bench_a3_03 / wireframe`
   - structural score: `0.875`
   - CFD score: `0.600`
   - obstacle count: `1 -> 3`
   - opening walls: `E/S -> E/S`

Rationale:
- both belong to the scale-calibrated `opening_topology_preserved_with_obstacle_hallucination` pattern,
- both preserve opening-wall topology despite obstacle inflation,
- wireframe inputs were selected after birdseye inputs failed independent QC for being too faint at print scale.

## Layout decision

- intended width: **double-column**
- panel layout: **2 x 3**
- panel sequence per row: `input / reference / prediction`
- export: `PDF + PNG(600 dpi)`
- font family requested: `Arial`
- actual fallback available on this host: **DejaVu Sans** (`Arial` unavailable)

## Current outputs

- script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure9_obstacle_hallucination.py`
- PDF: `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.pdf`
- PNG: `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.png`
- staged PNG for external QC: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/assets/figure9_obstacle_hallucination_limited_cfd_penalty.png`
- Gemini artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r3-20260321-0112.md`

## Revision summary for this turn

The figure was revised again after the previous external FAILs as follows:
- prediction badges were rewritten from shorthand to explicit compact wording (`# obstacles ...`, `openings match`, `CFD score ...`),
- GT outline color was moved away from inlet blue to a purple dashed style so crowded prediction panels separate more cleanly,
- figure-level typography, legend size, room/opening strokes, and dashed-outline thickness were increased,
- wireframe inputs were recropped more aggressively and passed through stronger autocontrast / contrast / sharpness enhancement,
- panel width ratios and bottom legend spacing were adjusted to give the input column slightly more visual weight.

The rerender materially improved readability and GT-vs-inlet separation, but external QC still did not clear the figure for final manuscript delivery.

## Self visual QC

Verdict: **PASS**

Self-check notes on the current revision:
- no forbidden figure-number / caption-style prose remains inside the figure,
- prediction badges are now less shorthand-heavy and the GT outline no longer collides chromatically with inlet blue,
- PDF/PNG export succeeded with matching layout and the bottom legend is less trim-risky than before,
- however, the wireframe-vs-plan-view stroke balance and the standalone interpretation of `CFD score 0.60` still look borderline enough that external QC remains the decisive gate.

## Independent subagent visual QC (current revision)

Verdict: **FAIL**

Blocking issues reported:
- `limited CFD penalty` is still not self-explanatory from the figure alone; `CFD score 0.60` lacks enough on-figure interpretation,
- print-scale readability is still risky for the prediction callout text, legend items, and smaller markers,
- GT outline vs opening markers is improved but still not robust enough in crowded panel `(f)`,
- the prediction text boxes continue to compete with geometry, and the left wireframe column remains visually weaker than the plan-view columns.

## Gemini multimodal QC (current revision)

Verdict: **FAIL**

Blocking issues reported:
- the main remaining blocking issue is **stroke-weight imbalance**: 2D plan-view boundaries now feel heavy while the 3D wireframe inputs still look too thin/light for conservative print reduction,
- panel `(c)` still has slight top-left badge crowding near the inlet marker,
- the `CFD score 0.60` cue is clearer than before but still depends on reader intuition to understand why it represents a limited penalty.

## 2026-03-21 late revision update (same turn)

A further repair pass was attempted after the earlier FAILs:
- metric badges were moved out of the prediction geometry region and re-anchored in figure-space above the prediction panels,
- wireframe preprocessing was strengthened again (grayscale autocontrast + stronger contrast + unsharp mask + min-filter dark-line expansion),
- prediction-column width was trimmed slightly while the input column was widened,
- GT outline color/linewidth/dash pattern were strengthened again,
- figure height / inter-row spacing were increased to reclaim non-panel whitespace for annotations.

### Self visual re-check on the late revision

Verdict: **FAIL**

Blocking findings from the conservative print-scale recheck:
- the metric badges still read as too close to panel frames rather than clearly separated margin annotations,
- the wireframe inputs remain inconsistent in stroke visibility after reduction; some room/inner edges still feel too light relative to the plan-view panels,
- the purple GT outline is stronger than before but still not yet robust enough across all panels.

### External QC rerun status on the late revision

- Subagent multimodal QC: **NOT RUN** on this revision (held back because self re-check still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (held back because self re-check still failed)

## 2026-03-21 latest revision update (annotation-lane rerender)

A further repair pass was completed on the same day:
- the metric callouts were moved into a dedicated right-side annotation lane instead of hovering near the prediction-panel frame,
- the layout was rebalanced so the reference/prediction plan panels gained width and the large inter-row whitespace was reduced,
- wireframe preprocessing was tightened again to force near-binary dark linework,
- predicted/reference obstacle edges were thickened,
- opening markers and GT outline were rendered with stronger white-halo separation and slightly heavier strokes,
- legend font size / handle spacing / bottom margin were increased for safer trim.

### Self visual re-check on the annotation-lane rerender

Verdict: **PASS**

Self-check notes:
- the badge/frame collision from the previous late revision is gone,
- wireframe contrast is now more uniform across the two rows,
- plan-view comparison panels have more visual weight,
- PDF/PNG remain layout-matched and trim-safe on the current local check.

### Gemini multimodal QC on the annotation-lane rerender

Verdict: **FAIL (conditional / provisional pass wording, but not a clean pass)**

Blocking findings carried by the latest Gemini artifact (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r6-20260321-021130.md`):
- panel `(f)` still leaves the purple GT outline semantically ambiguous when it clusters near the openings,
- inlet/outlet markers remain smaller than ideal for conservative print reduction,
- callout boxes still sit visually close to the room boundary region,
- wireframe/floorplan line-weight mismatch improved but was still flagged for additional polish.

### Independent subagent QC on the annotation-lane rerender

Verdict: **PENDING at end of this turn**

- a new rerun (`fig9-subagent-qc-rerun2`) was launched on the latest revision,
- the run had not returned a final verdict before this cron turn ended,
- therefore Figure 9 remains incomplete regardless of the improved layout.

## 2026-03-21 latest repair pass after annotation-lane rerender

A further same-night repair pass was applied to the annotation-lane revision:
- legend text was clarified from `GT outline` to `GT obstacle outline`,
- reference obstacle fill/edge contrast was darkened slightly for grayscale safety,
- inlet/outlet strokes were enlarged again with a stronger white halo,
- panel `(f)` gained a short purple `GT obst.` arrow label so the dashed reference obstacle is not mistaken for a wall mismatch,
- the right annotation lane width / inter-panel spacing were increased slightly.

### Self visual re-check on this repair pass

Verdict: **PASS**

Self-check notes:
- panel `(f)` now explains the purple dashed geometry more explicitly than the previous revision,
- the dedicated right-side callout lane remains outside the prediction geometry and trim-safe,
- opening markers and legend semantics are clearer than the previous iteration,
- PDF/PNG layout match remains intact on local inspection.

### Independent subagent QC on this repair pass

Verdict: **FAIL / REVISE**

Blocking findings from the latest independent rerun:
- panel `(f)` GT obstacle outline is still borderline at conservative print scale because the purple outline/label competes with nearby obstacles and opening markers,
- annotation text size is still a little tight for print reduction, especially the purple GT label and smaller legend/callout text,
- the intended message (obstacle hallucination while wall topology is preserved) is more visible now but still not fully self-evident without caption support.

### Gemini multimodal QC on this repair pass

Verdict: **FAIL**

Blocking findings from `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r7-20260321-023923.md`:
- GT obstacle outline in panel `(f)` is still visually obscured by the heavy wall boundary / nearby predicted obstacle,
- overall line-weight balance is still imperfect because the plan-view walls remain heavier than the wireframe inputs,
- legend text remains slightly smaller than ideal for conservative print-scale readability.

## 2026-03-21 latest rerender (external-callout revision)

A further same-turn rerender was produced after the annotation-lane version:
- wireframe crops were tightened again and linework was thickened slightly to raise print-scale visibility,
- plan-view wall / obstacle strokes were reduced modestly to soften the wall-dominance imbalance,
- legend and callout typography were enlarged,
- the panel `(f)` `GT obst.` cue was moved out of the prediction panel into the right annotation lane with an external purple leader arrow.

### Self visual re-check on the external-callout revision

Verdict: **PASS**

Self-check notes:
- wireframe panels are larger and darker than the previous revision,
- panel `(f)` no longer contains an in-panel purple label competing directly with the predicted geometry,
- legend readability improved and the callout lane remains trim-safe on the local inspection.

### Independent subagent QC on the external-callout revision

Verdict: **FAIL / REVISE**

Blocking findings:
- panel `(f)` GT obstacle cue is still borderline at conservative print scale because the purple outline remains small and close to the inlet/outlet cluster,
- inlet/outlet markers in panel `(f)` are still partially occluded / too small for comfortable print reduction,
- wireframe-vs-plan-view line-weight balance improved but is still not fully closed,
- right-side callout overlap / edge safety improved, but the figure is not yet clean enough for final manuscript delivery.

### Gemini multimodal QC on the external-callout revision

Verdict: **FAIL**

Blocking findings from `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r8-20260321-030943.md`:
- the external purple `GT obst.` label + leader still collides visually with the gray metric bubble in row 2,
- callout text appears vertically too low inside the rounded boxes,
- panel `(f)` still has insufficient separation among GT cue / opening markers / nearby predicted geometry.

## 2026-03-21 latest rerender (GT inset revision)

A further repair pass replaced the external boxed `GT obst.` tag with a dedicated inset:
- the row-2 GT cue was converted from a boxed text tag into a small right-lane inset that shows the GT obstacle alone,
- metric box text was raised and rewritten from `CFD 0.60 (moderate)` to `CFD stays 0.60` to make the limited-penalty message slightly more self-explanatory,
- wireframe column width was increased while plan-view room / obstacle / GT-outline strokes were softened modestly to reduce wall-dominance,
- legend typography and inlet/outlet legend strokes were enlarged,
- overall figure height increased slightly and row spacing was tightened.

### Self visual re-check on the GT inset revision

Verdict: **PASS**

Self-check notes:
- the previous direct collision between the row-2 purple cue and the gray metric bubble is resolved,
- panel `(f)` now explains the GT obstacle through a dedicated inset instead of an in-panel boxed label,
- PDF/PNG layout match remains intact on the local inspection,
- however, conservative external print-scale QC still remains the gate.

### Independent subagent QC on the GT inset revision

Verdict: **FAIL / REVISE**

Blocking findings:
- the intended message is still not fully self-evident without caption support because the obstacle discrepancy depends on a small GT cue and compact plan views,
- panel `(f)` remains crowded around openings / nearby geometry even after moving the GT cue into an inset,
- wireframe lines are still too thin relative to the plan-view walls for conservative 2-column reduction,
- right-side callouts and legend text are improved but still marginal at print scale.

### Gemini multimodal QC on the GT inset revision

Verdict: **FAIL**

Blocking findings from `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r9-20260321-034323.md`:
- wireframe line weights in `(a)` and `(d)` are still too thin and risk fading under 2-column reduction,
- the new `GT obstacle` inset/connector still sits too far outside the natural panel cluster and the purple connector remains too thin,
- vertical whitespace between rows is still larger than ideal, so the figure reads as slightly disconnected.

## Carry-over repair targets

1. materially thicken the wireframe stroke visibility in panels `(a)` and `(d)` rather than relying only on contrast preprocessing,
2. pull the GT obstacle inset/cue closer to panel `(f)` and thicken the purple connector / cue geometry,
3. reduce the inter-row whitespace further so the 2-row composition feels denser and more unified,
4. keep Figure 9 **incomplete** until self + subagent + Gemini all pass on the same revision.

## 2026-03-21 dawn revision update (count-badge / de-crowding pass)

A further same-night repair pass was applied after the GT-inset branch stalled:
- the prediction panels were rerendered with stronger room / obstacle / GT-outline / opening strokes,
- the wireframe preprocessing was thickened again and the panel grid was rebalanced toward more equal visual scale across input / reference / prediction,
- the previous GT inset and leader-arrow treatment was removed entirely,
- the ambiguous `extra obs.` arrow callout was replaced by a compact top-right count badge inside each prediction panel (`extra +N`, `pred / GT`),
- the right annotation lane wording was expanded from compass abbreviations to full `North/West` and `East/South` wording,
- the reference obstacle fill gained hatching for grayscale robustness,
- the bottom legend was split into two rows to reduce horizontal compression risk.

### Self visual re-check on the count-badge revision

Verdict: **FAIL**

Blocking findings from the latest conservative self/vision check:
- wireframe strokes in `(a)` and `(d)` are still too light relative to the plan-view panels for safe 2-column reduction,
- prediction panel `(f)` remains crowded because the GT dashed obstacle, predicted obstacle, and inlet/outlet markers are still too compressed in the same local region,
- the new in-panel count badge is clearer than the previous arrow callout, but the key message still relies on small text and color cues rather than instantly legible geometry,
- legend and callout typography improved, but not enough yet for a safe manuscript final.

### External QC rerun status on the count-badge revision

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)

## Updated carry-over repair targets

1. either materially enlarge / zoom the crowded prediction region in `(f)` or introduce a local inset so GT-vs-predicted obstacle mismatch and opening markers separate cleanly,
2. further thicken wireframe and GT-outline strokes for conservative print reduction,
3. enlarge the count/callout typography again and make the limited-CFD cue more immediately interpretable at a glance,
4. keep Figure 9 **incomplete** and do **not** copy it to Drive until a single revision passes self + subagent + Gemini QC together.

## 2026-03-21 pre-dawn revision update (row-summary + zoom-inset pass)

A further repair pass was produced after the count-badge revision:
- wireframe preprocessing was strengthened again with heavier dilation/contrast,
- the layout was rebalanced toward a wider input column and a wider right annotation lane,
- in-panel `extra +N / pred / GT` count badges were removed to reduce prediction-panel clutter,
- the right annotation lane was rewritten into explicit per-row summary cards (`A3-01 summary`, `A3-03 summary`) with `False positives`, `Openings preserved`, and `CFD score ... (higher is better)` wording,
- panel `(f)` gained a local zoom inset tied to the crowded right-side opening/obstacle region.

### Self visual re-check on the row-summary + zoom-inset revision

Verdict: **FAIL**

Blocking findings from the latest strict self/vision check:
- panel `(f)` is still too crowded at conservative manuscript scale even with the inset; the dense overlap problem was reduced but not actually closed,
- the inset itself is not yet clean enough to serve as a decisive de-crowding device at print scale,
- bottom legend readability is still borderline-small,
- overall 2-column reduction safety remains insufficient, so external QC was withheld on this revision.

### External QC rerun status on the row-summary + zoom-inset revision

- Subagent multimodal QC: **NOT RUN** on this revision (held back because latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (held back because latest self gate still failed)

## 2026-03-21 dawn-late revision update (matched-vs-extra split pass)

A further repair pass was produced after the older GT-inset / badge variants stalled:
- prediction obstacles were split into two visual roles so the figure now distinguishes `matched pred.` from `extra pred.` in both the panels and legend,
- the older row-2 zoom/inset device was removed to declutter panel `(f)` and reclaim panel area,
- the overall layout was widened/rebalanced toward a larger wireframe column and a larger prediction column with tighter row spacing,
- wireframe preprocessing was thickened again to improve conservative print-scale visibility,
- summary boxes were repadded and shortened (`N/W` / `E/S` wording) after Gemini flagged border-overlap / cramped text issues,
- the latest staged asset is `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/assets/figure9_obstacle_hallucination_limited_cfd_penalty.png`.

### Self visual re-check on the matched-vs-extra split revision

Verdict: **FAIL**

Blocking findings from the latest conservative self gate:
- the matched-vs-extra split materially clarifies the false-positive story, but the summary cards and bottom legend are still too dense for comfortable 2-column reduction,
- panel `(f)` is cleaner than the previous inset branch but still visually crowded near the east/south opening cluster,
- the wireframe column is darker now, yet it remains raster-like / alias-prone compared with the cleaner plan-view panels,
- the `CFD stays 0.60` claim is better framed than before, but it is still primarily text-driven rather than visually supported.

### Independent subagent QC on the matched-vs-extra split revision

Verdict: **FAIL**

Blocking findings:
- summary panels and the bottom legend remain too dense for conservative print-scale readability,
- the semantic story is still not fully self-evident without caption help because the key claims are carried mainly by the summary cards,
- panel `(f)` is still crowded around the outlet / GT / matched-prediction region,
- the limited-CFD message is still numerical-textual rather than visually demonstrated inside the panel set.

### Gemini multimodal QC on the matched-vs-extra split revision

Verdict: **FAIL**

Blocking findings from the latest Gemini CLI rerun:
- line-weight / style consistency is still off: thick purple GT segments and rasterized wireframe strokes compete awkwardly with the cleaner 2D wall lines,
- panel `(f)` remains cluttered because the GT obstacle cue still overlaps visually with nearby predicted geometry,
- row alignment / whitespace still feels somewhat disjointed,
- summary-box typography is still cramped, and the legend icon for `GT obstacle` does not perfectly match how the GT cue reads inside panel `(f)`.

## 2026-03-21 latest revision update (compact cue-card + column-header pass)

A new cleanup pass was produced on the same morning:
- panel titles were shortened to bare panel labels while the semantic headers moved to figure-level column headers (`input wireframe`, `reference`, `prediction`, `key cues`),
- the right summary lane was compressed into lighter cue cards with fewer lines and larger numeric emphasis (`Obstacles ...`, opening-wall match, `CFD 0.60`, `limited penalty`),
- panel borders were softened slightly and the width balance was shifted toward the reference/prediction columns,
- the figure was rerendered to reduce text density and make the key-cue lane less caption-like.

### Self visual re-check on the compact cue-card revision

Verdict: **FAIL**

Blocking findings from the latest strict self/vision gate:
- panel `(f)` is still too crowded around the east/south opening cluster; the purple GT outline, extra-prediction patch, and inlet/outlet markers still merge under conservative 2-column reduction,
- the core message (`extra obstacles but limited CFD penalty because openings are preserved`) is still not fully self-evident without caption support,
- line-weight balance is still off because the wireframe column remains visually stronger than the plan-view panels,
- the cue cards and bottom legend are lighter than before, but they still remain borderline-dense for conservative print reduction.

### External QC rerun status on the compact cue-card revision

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)

## 2026-03-21 sunrise revision update (wider cue-card + retained-CFD bar pass)

A new cleanup pass was produced on the same morning:
- the overall double-column canvas was widened slightly and the row spacing was tightened,
- the right-side cue-card lane was widened and rewritten from sentence-heavy notes into more compact `0 → 3 obstacles`, `N/W openings match`, `CFD retained` style cards,
- the prediction panels gained slightly stronger opening markers and a semi-transparent purple GT-obstacle fill to separate GT from the tan predicted obstacle,
- room-wall strokes were softened slightly while keeping the plan-view panels larger than the earlier versions,
- the latest staged asset was refreshed at `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/assets/figure9_obstacle_hallucination_limited_cfd_penalty.png`.

### Self visual re-check on the sunrise revision

Verdict: **PASS (local only)**

Self-check notes:
- the older key-card text-overlap failure mode is materially reduced,
- the cue cards are now shorter and the `CFD retained` bar reads more directly than the earlier sentence-style cards,
- PDF/PNG remain layout-matched on the local inspection,
- however, the panel `(f)` east/south cluster and conservative trim margin remain close enough that external QC still decides the gate.

### Independent subagent QC on the sunrise revision

Verdict: **FAIL**

Blocking findings from the latest independent rerun:
- the core message is still not fully self-evident without caption help because small legend / cue semantics carry too much of the explanation,
- panel `(f)` remains crowded around the east/south opening cluster,
- wireframe inputs `(a)` and `(d)` are still visually overweight relative to the plan-view panels,
- legend and key-cue symbols remain borderline-small for conservative print reduction,
- bottom/right trim safety is improved but still not comfortably closed.

### Gemini multimodal QC on the sunrise revision

Verdict: **FAIL**

Blocking findings from `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r10c-20260321-061125.md`:
- wireframe line weight in `(a)` and `(d)` is still too heavy and risks collapsing into dark blobs after two-column reduction,
- panel `(f)` still has a blocking local crowding issue: the south inlet touches the room-corner region and the tan / purple / east-wall cluster lacks enough separation,
- the right-side cue cards remain too close to the right trim boundary for a comfortable journal margin.

### Carry-over repair targets after the sunrise revision

1. materially reduce wireframe stroke weight in `(a)` and `(d)` rather than only widening the canvas,
2. de-crowd panel `(f)` by moving the south inlet away from the corner and increasing separation between the tan predicted obstacle, purple GT obstacle cue, and east wall,
3. add a slightly larger right outer margin / breathing room for the cue-card lane,
4. keep Figure 9 **incomplete** and do **not** copy it to Drive until a single revision passes self + subagent + Gemini QC together.

## 2026-03-21 06:31 KST latest revision update (margin + balance rebake)

A further self-only repair pass was produced after the sunrise revision:
- wireframe preprocessing was intentionally softened (lower contrast / less dilation) to reduce the overweight left-column appearance flagged by the latest external checks,
- the prediction column was widened while the right cue-card lane also gained width, and the overall canvas/right margin were adjusted to create more breathing room near trim,
- prediction/reference panels now use larger outer padding so the east/south edge cluster has more room inside panel `(f)`,
- the purple GT-obstacle cue was softened (lighter fill + thinner dashed outline) to reduce direct competition with the tan predicted obstacles,
- cue-card typography and the retained-CFD bar/value were enlarged so the right lane reads faster at a glance,
- the bottom legend was enlarged slightly and respaced.

### Self visual re-check on the margin + balance rebake

Verdict: **FAIL**

Blocking findings from the latest conservative self/vision gate:
- the core message is still not visually self-sufficient enough; readers still have to decode the legend/cue cards before the `extra obstacles but CFD retained` claim lands,
- panel `(f)` remains the primary blocker because the purple GT cue, tan extra obstacle, and opening markers are still too compressed in the east/south region,
- right cue cards are improved but still sit closer to the trim boundary than is comfortable for final journal submission,
- legend/cue text density remains borderline-small for conservative two-column print reduction, especially if grayscale reproduction is considered.

### External QC rerun status on the margin + balance rebake

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)

### Updated carry-over repair targets after the margin + balance rebake

1. make the limited-penalty message more on-panel and less legend-dependent (likely by simplifying the cue lane and/or adding a stronger per-row retained-CFD cue),
2. de-crowd panel `(f)` more aggressively, potentially with a cleaner local zoom or stronger geometry separation rather than additional small symbols,
3. push the cue-card lane further inward / increase the right outer margin again,
4. enlarge or simplify legend/text encodings further and reduce color-only dependence for grayscale safety,
5. keep Figure 9 **incomplete** and do **not** copy it to Drive until a single revision passes self + subagent + Gemini QC together.

## 2026-03-21 07:0x KST latest revision update (cue-card de-overlap rebalance)

A further local repair pass was completed on the same morning:
- the right cue cards were rebuilt to eliminate the visible `CFD retained` / `0.60` overlap from the previous rebake branch,
- the overall canvas was widened slightly again and the prediction / cue-card columns gained width,
- prediction-panel padding was increased so the east/south cluster in panel `(f)` sits farther from the frame,
- the bottom legend spacing was enlarged modestly.

### Self visual re-check on the cue-card de-overlap rebalance

Verdict: **FAIL**

Blocking findings from the latest conservative self/image check:
- the cue-card overlap bug is largely fixed, but the key claim is still not self-evident enough without caption help,
- panel `(f)` remains the dominant blocker: the purple GT cue, orange obstacle patch, and east/south opening markers are still too tightly clustered,
- the right cue-card lane is cleaner than before but still borderline-tight for print reduction,
- wireframe panels remain visually lighter/smaller than the plan-view panels, so input-vs-output balance is not yet convincingly closed.

### External QC rerun status on the cue-card de-overlap rebalance

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)

### Updated carry-over repair targets after the cue-card de-overlap rebalance

1. de-crowd panel `(f)` more aggressively, likely with either a cleaner local zoom / split cue or stronger geometry separation around the east/south opening cluster,
2. make the `limited CFD penalty` message land more directly on-figure instead of through card/legend decoding,
3. enlarge or rebalance the input wireframe panels so they no longer look underweighted against the plan-view panels,
4. keep Figure 9 **incomplete** and do **not** copy it to Drive until a single revision passes self + subagent + Gemini QC together.

## 2026-03-21 07:3x KST latest local rebalance (inward cue-lane + softer wireframe + obstacle→CFD link)

A new self-only repair pass was produced this turn:
- wireframe preprocessing was softened again to avoid the left-column overweight look from earlier sunrise revisions,
- the whole canvas was widened slightly while the right cue-card lane was pulled inward (`right=0.914`) to create more trim safety,
- the prediction column was widened and prediction-panel padding increased again so the east/south cluster in panel `(f)` sits farther from the frame,
- cue cards were simplified further (`3 extra obs.` / opening-wall line / big `CFD 0.60 retained`) and a thin orange connector was added to make the obstacle→CFD relation slightly more explicit,
- the latest staged asset was refreshed at `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/assets/figure9_obstacle_hallucination_limited_cfd_penalty.png`.

### Self visual re-check on the inward cue-lane rebalance

Verdict: **FAIL**

Blocking findings from the latest strict self/vision gate:
- panel `(f)` is still the dominant blocker: the purple GT cue, orange extra obstacle, and blue/red opening markers remain too compressed around the east/south cluster for conservative two-column reduction,
- right cue-card trim safety improved, but the `0.60` value / bar / new connector still sit too close to the card edge for a safe final export,
- the new obstacle→CFD connector helps slightly, but the `extra obstacles but CFD retained` message is still not visually self-evident enough to justify external QC yet.

### External QC rerun status on the inward cue-lane rebalance

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)

## 2026-03-21 08:0x KST latest revision update (prediction-width + simplified cue-card pass)

A further self-only repair pass was produced this turn:
- the prediction column was widened again and the right outer margin increased slightly to buy more trim safety,
- prediction-panel x/y padding was increased further so the east/south cluster in panel `(f)` sits farther from the frame,
- the right lane was simplified from the overloaded `retained` bar layout into a lighter card with three cues only: `extra obstacles`, `openings match`, and `CFD 0.60 stays moderate`,
- the figure-level right header was shortened from `obstacles / openings / CFD` to `key cues` to reduce title pressure,
- the latest staged output remains `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.png`.

### Self visual re-check on the prediction-width + simplified cue-card pass

Verdict: **FAIL**

Blocking findings from the latest strict self/vision gate:
- the cue-card overlap bug is materially reduced, but the right lane is still text-led and the `limited CFD penalty` message is not yet self-evident enough without caption support,
- panel `(f)` remains the dominant blocker: the purple GT cue, tan obstacle, and east/south opening markers are still too compressed for a conservative two-column print reduction,
- the right cue card is cleaner than before but still feels tight around the `CFD 0.60` cluster, so external QC was withheld again on this revision.

### External QC rerun status on the prediction-width + simplified cue-card pass

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because latest self gate still failed)

## 2026-03-21 08:0x KST strategy pivot — split-asset rebuild for Figure 9

Following repeated local-polish loops, the working strategy for Figure 9 is now explicitly changed from `single monolithic rerender with denser cue patches` to a **split-asset rebuild**.

### Diagnosis behind the pivot

The current blocker is no longer treated as a minor typography / spacing issue. The root cause is the combination of:
- a **crowded `(f)` prediction panel** where GT cue, extra predicted obstacle, and inlet/outlet markers collide in the east/south region,
- a **message-architecture mismatch** where the intended claim (`extra obstacles but limited CFD penalty because openings are preserved`) depends too heavily on cue cards / legend decoding rather than panel evidence,
- persistent **wireframe vs plan-view balance instability** under conservative 2-column print reduction.

### New production direction (approved pivot)

Adopt the following structure for the next Figure 9 rebuild attempt:

1. **Split-asset production first**
   - stop treating Figure 9 as a single all-in-one rendering target,
   - generate panel assets and cue assets as separate files first,
   - assemble only after panel-level readability is confirmed.

2. **Geometry evidence vs consequence evidence separation**
   - the main geometry panel set should focus on the visual fact that obstacle hallucination occurs while opening topology is preserved,
   - the `limited CFD penalty` message should be carried by a companion cue/summary asset rather than by dense in-panel overlays.

3. **Panel `(f)` gets special handling**
   - treat `(f)` as the dominant blocker and produce it as an independent asset,
   - permit larger internal padding, a local crop/zoom strategy, or stricter cue-role separation specifically for `(f)`,
   - do not reintroduce dense GT/pred/opening cue stacking inside `(f)` unless it clearly clears print-scale readability.

### Asset plan for the next rebuild

Recommended intermediate outputs:
- `fig9_row1_input`
- `fig9_row1_reference`
- `fig9_row1_prediction`
- `fig9_row2_input`
- `fig9_row2_reference`
- `fig9_row2_prediction`
- `fig9_row1_summary`
- `fig9_row2_summary`
- `fig9_legend`
- optional: `fig9_row2_prediction_zoom` if `(f)` still requires local de-crowding

### Information-allocation rules for the rebuild

- **Inside reference panels**: show GT/reference geometry only.
- **Inside prediction panels**: keep only the minimum geometry/cue set needed to identify predicted obstacles and openings.
- **Do not force the full `limited CFD penalty` explanation into `(f)` itself.**
- **Row-summary assets** should carry compact, row-level interpretation such as:
  - `obstacles: 1 -> 3`
  - `openings: preserved`
  - `CFD: 0.60 retained`
- Keep summary wording short and verdict-like; avoid sentence-style caption fragments.

### QC gate update after the pivot

For the next Figure 9 iteration:
- first clear a **panel-level local QC** for `(f)` and the row-summary assets,
- then clear an **assembled-figure local QC**,
- only then rerun subagent and Gemini external QC.

Until a single rebuilt revision passes self + subagent + Gemini together, Figure 9 remains **incomplete**.

## 2026-03-21 08:1x KST split-asset-informed rebalance (tighter prediction crop + shorter cue cards)

A further local repair pass was completed after the split-asset pivot was recorded:
- the assembled figure kept the same high-level 2x3+cue layout, but the rebuild adopted the split-asset intent by reallocating information away from sentence-like cue cards and back into larger geometry panels,
- the prediction column was widened again while inter-row spacing was tightened, and the right cue lane was narrowed slightly to buy more usable width for the plan-view evidence,
- prediction-panel crops were tightened materially (especially for the row-2 prediction panel) so empty whitespace was reduced and the crowded geometry cluster occupies a larger fraction of the panel,
- the cue cards were shortened from `extra obstacles / openings match / stays moderate` wording into more compact `extra obs. / preserved / moderate` verdict-style labels,
- the latest outputs were regenerated at `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.pdf` and `.png`.

### Self visual re-check on the tighter-crop revision

Verdict: **FAIL**

Blocking findings from the latest conservative self gate:
- panel `(f)` is larger and less whitespace-wasteful than before, but the east/south cluster is still the dominant blocker; the purple GT cue, tan obstacle patch, and opening markers remain too compressed for a conservative two-column reduction,
- the figure is cleaner than the earlier sentence-heavy cue-card branches, but the `limited CFD penalty` claim is still carried too much by the right cue cards rather than being visually self-evident from the panel set,
- legend readability improved only marginally and the bottom legend still feels borderline-small relative to the dense prediction evidence,
- the revision is useful progress, but not strong enough to justify a fresh subagent/Gemini rerun yet.

### External QC rerun status on the tighter-crop revision

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)

## 2026-03-21 08:34 KST redesign pass (2x3 + row-summary + panel-f zoom-inset)

A more structural rebuild was attempted this turn to break the repeated cue-lane polishing loop:
- the old right-side cue-card lane was removed entirely,
- the figure was rebuilt as a `2x3` main panel grid (`input / reference / prediction`) with a dedicated row-summary strip under each case,
- the summary strip now carries the row-level evidence in three separate boxes (`obstacles`, `openings`, `CFD`) instead of sentence-like cards,
- extra predicted obstacles were re-encoded with stronger orange fill plus `xx` hatching so hallucinated obstacles are less color-only,
- panel `(f)` received a dedicated local zoom inset focused on the crowded east/south opening + GT/prediction region,
- the shared legend was retained but simplified to stay compatible with the new panel grammar.

### Self visual re-check on the redesign pass

Verdict: **FAIL**

Blocking findings from the latest strict local self/image gate:
- the redesign direction is better architecturally, but the current execution still fails at conservative print scale because the summary-strip typography is too small,
- panel `(f)` remains the dominant blocker: the new inset is conceptually right but still too small / too weakly linked to the source region to truly resolve the east/south crowding,
- the bottom legend is still too dense and too small relative to the amount of symbol decoding required,
- the wireframe input panels are still low-utility at this scale and remain visually unbalanced relative to the plan-view evidence,
- line weight / hatch / dashed-outline visibility still needs another pass for grayscale-safe journal reduction.

### External QC rerun status on the redesign pass

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)

### Updated carry-over repair targets after the redesign pass

1. materially enlarge the row-summary typography and simplify the secondary text so the summary row survives conservative print reduction,
2. turn the current `(f)` inset into either a much larger zoom device or a dedicated subpanel so the GT/predicted/opening cluster separates cleanly,
3. either enlarge or demote/remove the small wireframe thumbnails so the plan-view evidence becomes the clear primary visual channel,
4. simplify/enlarge the bottom legend and reduce symbol decoding pressure,
5. keep Figure 9 **incomplete** and do **not** copy it to Drive until a single revision passes self + subagent + Gemini QC together.

## 2026-03-21 08:5x KST latest revision update (evidence-column + dedicated `(g)` detail panel)

A new structural rebuild was completed to push Figure 9 away from the failed `2x3 + thin summary strip` branch:
- the old per-row bottom summary strips were removed,
- the figure was rebuilt into a `2 x 4` composition: `input / reference / prediction / evidence`,
- row 1 now uses a full-height evidence column with larger `obstacles / openings / CFD` cards,
- row 2 uses a compact evidence block plus a dedicated `(g)` detail panel instead of an in-panel inset,
- the row-2 zoom crop was tightened to focus more directly on the east/south opening cluster,
- the prediction panel keeps the boxed zoom-source region, while the separate `(g)` panel reduces direct overlay pressure inside `(f)`.

### Self visual re-check on the evidence-column rebuild

Verdict: **PASS (local only)**

Self-check notes:
- the previous tiny row-summary strip failure mode is materially reduced,
- the figure grammar is clearer now because row-level evidence is spatially separated from the main geometry panels,
- the `(f)` blocker is improved: the dedicated `(g)` panel is cleaner than the prior inset branch and the zoom-source box is now explicit,
- PDF/PNG remain layout-matched on local inspection.

### Independent subagent QC on the evidence-column rebuild

Verdict: **FAIL**

Blocking findings:
- the core message is still not self-evident enough without caption support,
- panel `(f)` remains too crowded on the right side after downscaling, with matched/extra/GT/opening cues still merging,
- panel `(g)` is cleaner than the old inset but still too small to fully resolve the ambiguity,
- evidence boxes and the legend remain borderline-dense for conservative two-column print,
- line-weight balance still favors the plan/annotation side over the wireframe inputs.

### Gemini multimodal QC on the evidence-column rebuild

Verdict: **PASS**

Artifact:
- `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r11-20260321-085x.md`

Pass reasons highlighted by Gemini:
- hallucination evidence reads clearly through `extra pred.` hatching plus the explicit obstacle-count card,
- the `CFD 0.60 / moderate` evidence card is legible,
- panel `(g)` is considered helpful for alignment verification,
- legend spacing and plan-view readability were judged acceptable on this pass.

### Current gate state after the latest rebuild

- self QC: **PASS**
- subagent QC: **FAIL**
- Gemini QC: **PASS**

Figure 9 therefore remains **incomplete**, must **not** be marked done in the checklist, and must **not** be copied to Drive yet.

## 2026-03-21 09:0x KST follow-up revision update (larger `(g)` + direct zoom annotations)

A further focused repair pass was attempted on the same morning, targeting only the persistent `(f)/(g)` blocker:
- the rightmost column was widened and the row-2 split was rebalanced so panel `(g)` receives more vertical area,
- panel `(g)` was upgraded from a plain detail crop to a labeled detail panel with direct short callouts (`extra pred.`, `GT obst.`, `openings`) instead of relying only on legend decoding,
- the source zoom box in panel `(f)` was strengthened and given a short `see (g)` cue,
- compact evidence typography for row 2 was enlarged slightly, and the overall canvas / margins were relaxed for better trim safety,
- refreshed outputs were written again to:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.png`

### Self visual re-check on the larger-`(g)` annotated revision

Verdict: **FAIL**

Blocking findings from the latest strict self/image gate:
- panel `(g)` is larger and more interpretable than the prior branch, but the GT obstacle cue is still too thin / low-contrast relative to the tan predicted obstacle and does not survive conservative print reduction safely,
- panel `(f)` still leaves the row-2 east/south cluster too compressed; the new `see (g)` cue helps navigation but does not fully close the local ambiguity,
- the new zoom annotations and bottom legend remain too small/thin for a confident final manuscript send,
- the central message is improved but still not sufficiently self-evident without caption support because too much decoding still depends on small legend / color / dashed-outline semantics.

### External QC rerun status on the larger-`(g)` annotated revision

- Subagent multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)
- Gemini multimodal QC: **NOT RUN** on this revision (withheld because the latest self gate still failed)

### Updated carry-over repair targets after the annotated-zoom revision

1. materially thicken and simplify the GT obstacle encoding in both `(f)` and `(g)` so it survives conservative print reduction,
2. enlarge `(g)` again or promote it to a more dominant comparison device rather than a supporting detail,
3. simplify/enlarge legend symbols and annotation strokes so the figure stops depending on tiny color-coded cues,
4. keep Figure 9 **incomplete** and do **not** copy it to Drive until a single revision passes self + subagent + Gemini QC together.

## 2026-03-21 09:2x KST follow-up revision update (split `(g)` pred-vs-GT detail pair)

A further focused rebuild was completed on the same morning to reduce caption/legend dependence in the stubborn row-2 crowded region:
- the old single `(g)` annotated zoom panel was removed,
- the row-2 detail area was rebuilt into a side-by-side pair: `(g) pred.` and `GT`, both using the same crop so the hallucinated prediction can be compared against the true obstacle configuration directly,
- the fourth-column width was increased and the row-2 detail region was given more vertical area,
- the compact row-2 evidence card was simplified into a single larger box (`1 → 3`, `E/S preserved`, `CFD 0.60`) to avoid the prior overlap failure mode,
- the row-2 source box and `see (g)` cue in panel `(f)` were retained so the crowded region maps cleanly into the new direct-comparison detail pair.

### Self visual re-check on the split-detail revision

Verdict: **PASS (local only)**

Self-check notes:
- the previous row-2 evidence overlap bug is resolved,
- `(g)` is now structurally clearer because the reader can compare prediction and GT with the same crop instead of decoding multiple small arrows inside a single inset,
- the central row-2 message (`extra predicted obstacles while opening topology is preserved`) is more self-evident than the prior single-zoom branch,
- PDF/PNG remain layout-matched on local inspection.

### External QC rerun status on the split-detail revision

- Subagent multimodal QC: **NOT RUN** on this revision yet
- Gemini multimodal QC: **NOT RUN** on this revision yet

### Carry-over status after the split-detail revision

- Figure 9 remains **incomplete**,
- Drive copy remains **forbidden**,
- next needed step is a fresh external QC rerun on this exact revision before any completion decision.

## 2026-03-21 09:5x KST follow-up refinement (wider detail/evidence column + row-2 evidence reflow)

A further same-morning repair loop targeted only the remaining external blockers from the split-detail branch:
- the fourth column was widened again and the row-2 detail pair was given slightly more visual priority,
- panel `(g)` bottom note tags (`extra obstacles`, `true obstacle`) were lifted upward for better trim safety,
- the panel `(f)` source box was pushed behind the opening markers so the dashed source rectangle no longer dominates the red outlet marker,
- the row-2 evidence box was first rewritten into compact aligned `obs. / open. / CFD` rows, then rerendered again with a separated `0.60` + `moderate` stack after the first reflow still looked too dense.

### Self visual re-check on the latest follow-up refinement

Verdict: **FAIL**

Blocking findings from the latest strict self gate:
- the row-2 evidence box remains visually fragile; after the second reflow it is no longer a hard collision locally, but the value/qualifier stack is still cramped enough to look borderline at conservative double-column reduction,
- panel `(f)` is improved but still locally crowded because the source box, dashed GT cue, and opening markers all compete in the same east/south region,
- the bottom legend still sits close enough to the lower trim edge that it does not yet feel comfortably final.

### Independent subagent QC on the latest follow-up refinement

Verdict: **FAIL**

Blocking findings:
- row 2 is still not self-evident enough without caption support; the reviewer still feels the reader must decode legend/text to understand the full `extra predicted obstacles + preserved openings + limited CFD penalty` claim,
- panel `(f)` remains too crowded at print scale,
- panel `(g)` is improved but still judged too small to carry the argument cleanly on its own,
- the bottom legend / right-side detail area still feels somewhat margin-tight.

### Gemini multimodal QC on the latest follow-up refinement

Verdict: **FAIL**

Blocking findings from `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r14-20260321-0951.md`:
- the newest row-2 evidence-card reflow is still being read as overlapping / insufficiently separated at print scale,
- the bottom legend still wants a bit more vertical gutter,
- panel `(f)` crowding is reduced but remains visibly tight near the `see (g)` source region.

### Updated carry-over repair targets after the latest follow-up refinement

1. rebuild the row-2 evidence box with materially more vertical room (not another micro-adjustment) so external QC stops reading any overlap/alignment ambiguity,
2. give the bottom legend a larger explicit lower margin,
3. either simplify panel `(f)` source-box signaling further or promote the row-2 detail comparison more aggressively so the crowded region stops carrying so much decoding burden,
4. keep Figure 9 **incomplete** and do **not** copy it to Drive until one exact revision passes self + subagent + Gemini together.

## 2026-03-21 implementation note — split-asset builder committed to repo

A dedicated split-asset production script now exists:
- `scripts/paper_figures/make_figure9_obstacle_hallucination_split_assets.py`

Companion rebuild note:
- `docs/figure_qc/26-03-21_figure9_split_asset_rebuild_path.md`

Implemented outputs:
- `results/paper_figures/figure9_split_assets/figure9_row1_*`
- `results/paper_figures/figure9_split_assets/figure9_row2_*`
- `results/paper_figures/figure9_split_assets/figure9_row2_detail_pair.*`
- `results/paper_figures/figure9_split_assets/figure9_geometry_grid.*`
- `results/paper_figures/figure9_split_assets/figure9_split_asset_preview.*`
- `results/paper_figures/figure9_split_assets/figure9_split_manifest.json`

Design rule of the implemented path:
- main prediction panels drop the GT overlay,
- row-2 GT comparison is promoted into a separate `pred.` vs `GT` detail asset,
- row-level CFD interpretation is kept in standalone evidence assets rather than the overloaded cue lane.

QC status on this implementation note:
- assets and preview render successfully,
- self / subagent / Gemini were **not** rerun on this new exact preview revision yet,
- Figure 9 therefore remains **incomplete**.

## 2026-03-21 09:5x KST follow-up refinement 2 (row-2 declutter + staged Gemini rerun)

A further structural cleanup pass was completed on the same morning:
- row-2 prediction panel `(f)` no longer shows the full purple GT-obstacle outline overlay; instead the crowded region is delegated primarily to the dedicated `(g) pred.` vs `GT` detail pair,
- the old `see (g)` pill label was removed from `(f)` to reduce local annotation clutter, while the source crop box was retained in a lighter form,
- the fourth column was widened again and the row-2 evidence/detail stack was rebalanced toward a taller evidence box (`0.33 / 0.67` split),
- the row-2 evidence card was rewritten into three cleaner single-line rows (`obs.`, `open.`, `CFD`) to eliminate the previous value/qualifier stacking failure,
- the overall figure bottom margin and legend anchor were lifted slightly to create more gutter below the legend,
- the refreshed exact revision was restaged to `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/assets/figure9_obstacle_hallucination_limited_cfd_penalty.png` and rerun through both independent subagent QC and Gemini CLI QC.

### Self visual re-check on the row-2 declutter revision

Verdict: **PASS (local only)**

Self-check notes:
- the earlier row-2 evidence-card overlap is resolved on the exact staged revision,
- panel `(f)` is materially cleaner because the GT-vs-predicted comparison burden shifted into the dedicated `(g)` pair,
- the legend now has slightly more lower gutter than the previous revision,
- PDF/PNG remain layout-matched on local inspection.

### Independent subagent QC on the row-2 declutter revision

Verdict: **FAIL**

Blocking findings:
- the row-2 story is still not self-evident enough without caption support; the reader still has to infer too much from the A3-03 evidence/detail region,
- panel `(f)` is cleaner but remains crowded at conservative print scale,
- the row-2 evidence box is no longer overlapping, but it still reads as cramped / shorthand-heavy for a final two-column reduction,
- the `(g) pred.` / `GT` detail pair is only moderately clear and still a bit too small,
- the bottom legend gutter is improved but still judged marginal rather than comfortable.

### Gemini multimodal QC on the row-2 declutter revision

Verdict: **FAIL**

Artifact:
- `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-r15-20260321-095x.md`

Blocking findings:
- row-wise evidence formatting is inconsistent because row 2 uses a single compressed box while row 1 uses three spacious boxes,
- panel `(g)` internal labels (`extra obstacles`, `true obstacle`) are too small for conservative print reduction,
- the `(g)` title area is vertically cramped against the A3-03 evidence box above,
- panel `(f)` still has residual crowding where the dashed purple crop/GT cue and nearby geometry compete,
- external QC also notes that the reference obstacle role in panel `(e)` is not reinforced locally relative to the stronger `GT` labeling used in `(g)`.

### Updated carry-over repair targets after the row-2 declutter revision

1. make row-2 evidence formatting more consistent with row 1 (likely separate cards instead of one compressed box),
2. enlarge `(g)` labels and create more vertical gutter between the A3-03 evidence box and the detail-pair titles,
3. either further simplify the purple cueing in `(f)` or enlarge the row-2 detail comparison so the crowded region stops carrying so much semantic load,
4. keep Figure 9 **incomplete** and do **not** copy it to Drive until one exact revision passes self + subagent + Gemini together.

## 2026-03-21 10:0x KST repo implementation note (split-asset production path committed to repo)

A concrete repo-side split-asset production path is now implemented to preserve the locked representative pair while breaking the monolithic rerender loop:
- script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure9_obstacle_hallucination_split_assets.py`
- asset root: `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_split_assets/`
- manifest: `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure9_split_assets/figure9_split_manifest.json`

Generated asset set on this path:
- `figure9_row1_input`, `figure9_row1_reference`, `figure9_row1_prediction`
- `figure9_row2_input`, `figure9_row2_reference`, `figure9_row2_prediction`
- `figure9_row1_evidence`, `figure9_row2_evidence`
- `figure9_row2_detail_pair`
- `figure9_legend`
- `figure9_geometry_grid`
- `figure9_split_asset_preview`

Message allocation on the implemented path:
- main geometry panels carry only `obstacle hallucination + opening preservation`,
- evidence assets carry the row-level `obstacles / openings / CFD` interpretation,
- the row-2 `pred.` vs `GT` detail pair is the dominant blocker-resolver for crowded panel `(f)`,
- the main row-2 prediction panel is intentionally kept lighter than the older overlay-heavy revisions.

Expected QC implication of this implementation:
1. local panel-level QC should now inspect `row2_prediction`, `row2_detail_pair`, and both evidence assets before touching the assembled preview,
2. local assembled QC should then inspect `figure9_geometry_grid` and `figure9_split_asset_preview`,
3. only after those pass should subagent + Gemini reruns be launched on the exact same preview revision,
4. QC standards are unchanged: Figure 9 remains **incomplete** until one exact rebuilt revision passes self + subagent + Gemini together.

## 2026-03-21 radical layout rebuild — 2x2 + summary strip

### Diagnosis that motivated this rebuild

After 15+ revision cycles (all failing QC), the root causes were identified as architectural rather than cosmetic:

1. **Panel (f) crowding** — the 4-column layout gave the A3-03 prediction panel only ~2.5 effective inches. No amount of stroke/padding/cue tweaking resolved the east/south GT/extra/opening overlap.
2. **Message-architecture mismatch** — the `limited CFD penalty` claim was carried by cue cards and legend decoding rather than by panel evidence. Every cue-card variant failed QC.
3. **Wireframe vs plan-view format mismatch** — raster wireframe inputs vs vector plan-view drawings created inherent visual weight instability that oscillated between "too thin" and "too heavy" across 15+ attempts.

### Architectural change

The figure was **completely restructured** from `2x4 (input/reference/prediction/evidence)` to `2x2 (reference/prediction) + bottom summary strip`:

- **Wireframe input column: dropped entirely**. Figures 3/5/7 already show multi-view examples; the input view type is referenced in the caption. This eliminates root cause #3.
- **Right-side cue-card lane: dropped entirely**. Quantitative evidence moved to a compact bottom summary strip. This eliminates root cause #2.
- **Each prediction panel now gets ~2x the area** of the previous layout (from ~2.5 to ~5+ effective inches at double-column width). This directly addresses root cause #1.
- Panel labels simplified: `(a)` reference A3-01, `(b)` prediction A3-01, `(c)` reference A3-03, `(d)` prediction A3-03.
- Bottom summary strip: per-case obstacle count (`0 → 3`, `1 → 3`), opening match (`N/W ✓`, `E/S ✓`), CFD score (`0.60`).

### Representative cases

Unchanged from the locked set:
1. `bench_a3_01 / wireframe` — 0 → 3 obstacles, CFD 0.604, N/W preserved
2. `bench_a3_03 / wireframe` — 1 → 3 obstacles, CFD 0.600, E/S preserved

### Current outputs

- script: `scripts/paper_figures/make_figure9_obstacle_hallucination.py` (rewritten)
- PDF: `results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.pdf`
- PNG: `results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.png`
- staged QC asset: `docs/figure_qc/assets/figure9_obstacle_hallucination_limited_cfd_penalty.png`

### Self visual QC on the 2x2 rebuild

Verdict: **PASS**

Self-check findings:
- panel (d) (former blocker panel (f)) now has clear separation between GT outline, extra predicted obstacles, matched obstacles, and opening markers in the east/south region,
- the summary strip is clean with no text overlap and clear per-case metrics,
- legend is large, well-spaced, and all 6 items are clearly distinguishable,
- no figure-internal caption text, sentence-style prose, or figure numbers,
- PDF/PNG remain layout-matched,
- the message (obstacle hallucination + opening preservation + limited CFD penalty) is visually self-evident from the reference-vs-prediction comparison plus the summary strip,
- wireframe balance instability is eliminated entirely.

### Independent subagent QC on the 2x2 rebuild

Verdict: **PASS**

The first subagent run returned CONDITIONAL PASS with two blocking items:
1. (HIGH) grayscale robustness — matched pred vs extra pred too similar
2. (MEDIUM) CFD score lacked interpretive anchor

Both were fixed:
- matched pred lightened (α=0.45, pale cream `#F5E6CC`) vs extra pred darkened (α=0.92, burnt orange `#C45A20` + crosshatch)
- CFD anchor added: "score 0.60 / 1.00" → later refined to "moderate" + "score 0.60"

A second subagent run on the fixed version returned **PASS** on all 8 criteria with zero blocking issues. Minor suggestions (legend font borderline, panel (d) density, summary label size) were all non-blocking.

A subsequent polish pass (by user/linter) center-aligned all summary text and used a two-line topology format ("N/W\npreserved"). The final rendered figure retains all improvements.

### Gemini multimodal QC on the 2x2 rebuild

Verdict: **NOT YET RUN**

### Current gate state

- self QC: **PASS**
- subagent QC: **PASS**
- Gemini QC: **PENDING**

Figure 9 remains **incomplete** until Gemini QC passes on the same revision.

## 2026-03-21 10:1x KST follow-up regression check (explicit-summary-strip refinement)

A same-morning follow-up pass was attempted before launching Gemini on the 2x2 branch:
- the column header was simplified to `prediction`,
- the summary strip was rewritten to make the CFD interpretation more explicit (`moderate` + `score 0.60`),
- the canvas / summary-strip height / legend sizing were increased again.

### Self visual QC on the explicit-summary-strip refinement

Verdict: **FAIL**

Blocking findings from the newest exact local/image gate:
- the 2x2 architecture is still the strongest direction so far, but the main claim remains too **summary-strip dependent**; `limited CFD penalty` still does not land from the panels alone,
- summary-strip wording/abbreviations (`N/W`, `E/S`, `score 0.60`) remain borderline at conservative double-column reduction,
- the legend and symbol grammar are still too subtle for grayscale-safe print reduction (matched vs extra vs reference still depends on fine fill/hatch/color distinctions),
- the figure is cleaner than the older 2x4/cue-lane branches, but four panels + two summary strips + six legend items still feel somewhat dense for final manuscript delivery.

### External QC rerun status on the explicit-summary-strip refinement

- Independent subagent QC: **NOT RUN** on this exact revision (withheld because the newest self/image gate failed)
- Gemini multimodal QC: **NOT RUN** on this exact revision (withheld because the newest self/image gate failed)

### Updated carry-over repair targets after the explicit-summary-strip refinement

1. make the `limited CFD penalty` cue less dependent on reading small summary-strip text,
2. simplify/enlarge the legend and reduce dependence on subtle hatch/color distinctions for grayscale safety,
3. either enlarge or simplify the per-case summary grammar so opening preservation and CFD interpretation survive conservative print reduction,
4. keep Figure 9 **incomplete** and do **not** copy it to Drive until one exact revision passes self + subagent + Gemini together.
