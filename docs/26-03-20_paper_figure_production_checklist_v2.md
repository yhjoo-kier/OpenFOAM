# Paper Figure Production Checklist v2

> 작성일: 2026-03-20
> 목적: OpenFOAM Image-to-CFD 논문의 메인 figure 11종을 논문 투고 품질로 순차 생산하고, 각 figure마다 자체 QC + 서브에이전트 멀티모달 QC + Gemini 멀티모달 QC를 모두 통과한 경우에만 완료 처리 및 Google Drive 저장을 수행한다.

## 0. 기본 방침

- 현재 논문 메인 파이프라인은 **single-factor scale calibration이 포함된 image-to-CFD pipeline** 이다.
- 기존 no-scaling baseline은 논문 메인 결과로 사용하지 않는다.
- Figure는 "성능 개선 트릭" 서사가 아니라, **scale-calibrated pipeline의 benchmarked performance / failure structure / robustness** 를 전달해야 한다.
- figure 내부에는 caption성 문장, figure 번호(Figure 5, Fig. 5 등), 장문 설명을 넣지 않는다.
- figure 내부 허용 요소는 축 라벨, 범례, 짧은 annotation, panel label `(a), (b)` 및 매우 짧은 panel title 수준으로 제한한다.
- 모든 figure는 반드시 먼저 `single-column` 또는 `double-column` 여부를 결정하고 기록한다.
- 산출물은 항상 `PDF + PNG` 동시 생성, PNG는 **600 dpi 이상**이어야 한다.
- 가능하면 Arial을 우선 사용하고, 불가 시 가장 가까운 sans fallback을 QC log에 기록한다.
- PDF와 PNG는 동일 레이아웃/annotation/여백을 유지해야 한다.
- 완료 처리(`[x]`)는 **자체 QC + 서브에이전트 멀티모달 QC + Gemini 멀티모달 QC 모두 통과한 뒤에만** 허용한다.
- 완료본만 Google Drive에 복사한다.
- Google Drive 기존 폴더의 과거 산출물은 새 체계 시작 전에 정리하고, 이후에는 **이번 v2 체크리스트에서 완료 처리된 figure만** 저장한다.

## 1. 저장 경로 / 상태 원본

- 상태 원본 문서: `/home/yhjoo/projects/OpenFOAM/docs/26-03-20_paper_figure_production_checklist_v2.md`
- figure outputs: `/home/yhjoo/projects/OpenFOAM/results/paper_figures/`
- figure scripts: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/`
- QC logs: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/`
- Google Drive root: `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/`

## 2. 생산 우선순위 (효율 우선)

### Current autopilot-lane override (2026-03-21)

This checklist remains the state source, but the current figure-autopilot lane is operating under the following explicit override:

1. Figure 10 — Dense composite structure-vs-CFD gap
2. Figure 2 — Benchmark dataset design and difficulty matrix
3. Figure 3 — Multi-view rendering protocol
4. Figure 4 — Reference benchmark and evaluation pathway
5. Figure 7 — Representative cross-view outcome panel
6. Figure 11 — Robustness and convergence summary
7. Figure 1 — Overall framework with mandatory scale calibration (PaperBanana)

Additional lane rules now in force:
- Figure 9 is being handled in a separate lane and must be skipped here. Do not edit/re-render/re-QC Figure 9 from this autopilot lane.
- Figures 5, 6, and 8 stay closed as completed figures; do not reopen them.
- Each turn should advance only the single highest-priority unfinished figure in the list above.

## 3. 공통 작업 루프

각 cron 턴에서 다음을 수행한다.

1. 이 체크리스트와 기존 산출물/로그를 읽는다.
2. 아직 완료되지 않은 figure 중 우선순위가 가장 높은 1개를 선택한다.
3. source artifact / metrics / representative case를 확정한다.
4. single-column / double-column 여부를 먼저 결정하고 QC log 또는 체크리스트 메모에 남긴다.
5. figure script 작성 또는 갱신 후 PDF/PNG를 생성한다.
6. **자체 visual QC** 수행
7. **서브에이전트 멀티모달 QC** 수행
8. **Gemini 멀티모달 QC** 수행
9. 세 QC가 모두 통과하면 `[x]` 처리하고 Drive로 복사한다.
10. QC가 하나라도 실패하면 실패 원인과 수정 계획을 짧게 남기고, 다음 턴에서 이어서 보수한다.

## 4. Google Drive 정리 규칙

- v2 체계 시작 시, `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/` 아래의 **기존 figure 산출물은 삭제 후 재구성**한다.
- 단, 삭제는 이 figure 전용 폴더 범위 내에서만 수행한다.
- 완료한 figure는 아래와 같이 figure별 폴더로 저장한다.
  - `Fig05_view_aggregate/`
  - `Fig06_category_aggregate/`
  - ...
- 각 figure 폴더에는 최소한 아래를 포함한다.
  - 최종 PDF
  - 최종 PNG (>=600 dpi)
  - `qc_summary.md`

---

## Figure 5 — Aggregate performance across input views
- [x] single/double-column 결정
- [x] single-factor scaling 포함 메인 aggregate 결과 반영 확인
- [x] plotting script 작성/갱신
- [x] PDF/PNG 생성
- [x] 자체 visual QC
- [x] 서브에이전트 멀티모달 QC
- [x] Gemini 멀티모달 QC
- [x] Drive 저장
  - 2026-03-20 final note: source artifact was refreshed from legacy no-scaling aggregate to `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json` and Figure 5 was regenerated under the scale-calibrated main-setting rule.
  - Final layout decision: `double-column`, `1x2`, legends moved above each panel to eliminate compact-render crowding; PNG/PDF render consistency rechecked.
  - Final QC state: self `PASS`, subagent `PASS`, Gemini `PASS`.
  - Final QC log: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/26-03-20_figure5_view_aggregate_qc.md`
  - Drive folder: `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig05_view_aggregate/`

## Figure 6 — Aggregate performance across benchmark categories
- [x] single/double-column 결정
- [x] single-factor scaling 포함 category aggregate 반영 확인
- [x] plotting script 작성/갱신
- [x] PDF/PNG 생성
- [x] 자체 visual QC
- [x] 서브에이전트 멀티모달 QC
- [x] Gemini 멀티모달 QC
- [x] Drive 저장
  - 2026-03-20 final note: source artifact was refreshed from the legacy no-scaling aggregate to `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json` and Figure 6 was regenerated under the scale-calibrated main-setting rule.
  - Final layout decision: `double-column`, `1x2`; left margin and font sizing were increased during rerun so category labels stay unclipped and readable at print scale.
  - Final QC state: self `PASS`, subagent `PASS`, Gemini `PASS`.
  - Final QC log: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/26-03-20_figure6_category_aggregate_qc.md`
  - Drive folder: `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig06_category_aggregate/`

## Figure 8 — Composite-room collapse under section view
- [x] representative cases 재확인
- [x] single/double-column 결정
- [x] panel/annotation 정리
- [x] PDF/PNG 생성 또는 재생성
- [x] 자체 visual QC
- [x] 서브에이전트 멀티모달 QC
- [x] Gemini 멀티모달 QC
- [x] Drive 저장
  - 2026-03-20 final note: representative cases remain the only two `section_room_kind_collapse` tasks in the scale-calibrated aggregate (`bench_a3_04/section`, `bench_a4_05/section`), but the figure was rerun to remove forbidden in-figure caption text and to refresh the metric box from `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`.
  - Final layout decision: `double-column`, `2x2`; prediction panels use dashed GT outlines plus a short `Composite → rectangular` cue, with no embedded `Figure 8` / caption-style prose.
  - Final QC state: self `PASS`, subagent `PASS`, Gemini `PASS`.
  - Final QC log: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/26-03-20_figure8_section_collapse_qc.md`
  - Drive folder: `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig08_section_collapse/`

## Figure 9 — Obstacle hallucination with limited CFD penalty
- [x] representative cases 재확인
- [x] single/double-column 결정
- [x] annotation 최소화 정리
- [x] PDF/PNG 생성 또는 재생성
- [x] 자체 visual QC
- [x] 서브에이전트 멀티모달 QC
- [x] Gemini 멀티모달 QC
- [x] Drive 저장
  - 2026-03-21 final note: after 15+ revision rounds, the figure was rebuilt from scratch using `make_figure9_v2.py` with a radically simplified 2×2 (reference / prediction) layout. Key breakthroughs: (1) GT obstacle outlines removed from prediction panels — reference panels provide comparison via side-by-side layout, eliminating panel (d) SE corner crowding; (2) metric badges made self-explanatory ("+N hallucinated", "openings preserved (X/Y)", "CFD score 0.60"); (3) legend simplified to 5-item single row; (4) wireframe inputs omitted (referenced in caption).
  - Final layout decision: `double-column`, `2×2` grid with column headers, asymmetric height ratios (row 2 gets more space), row 2 pad_factor=0.18 for E/S breathing room.
  - Final QC state: self `PASS`, subagent (Opus) `PASS` (9/9 criteria), Gemini `PASS` (6/6 criteria).
  - Final script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure9_v2.py`
  - Final QC logs: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig9-visual-qc-v3-final.md`
  - Drive folder: `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig09_obstacle_hallucination/`

## Figure 10 — Dense composite structure-vs-CFD gap
- [x] representative cases 재확인
- [x] single/double-column 결정
- [x] aggregate/discussion 메시지 정합 확인
- [x] PDF/PNG 생성 또는 재생성
- [ ] 자체 visual QC
- [ ] 서브에이전트 멀티모달 QC
- [ ] Gemini 멀티모달 QC
- [ ] Drive 저장
  - 2026-03-21 progress note: Figure 10 was refreshed onto the scale-calibrated main setting. Aggregate source is now `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`, and per-case assets were rerouted to `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span/`.
  - Representative pair stayed locked to `bench_a4_02/floorplan` + `bench_a4_04/floorplan`, which remain the cleanest `dense_composite_structure_physics_gap` exemplars under the current manifest.
  - Current layout decision remains `double-column`, `2x2`, with short top headers only and no in-figure `Figure 10`/caption-style prose.
  - Typography note: Arial was requested first, but this host fell back to `DejaVu Sans` because `Arial` / `Liberation Sans` were unavailable in Matplotlib.
  - 2026-03-21 late follow-up: the right-column evidence pivoted from the previous 3D contour render crop to a tighter crop of each case's `comparison_1x2.png` CFD slice, while the left geometry panels stayed minimal. This made the flow mismatch more direct, but the exact current asset still failed external gates.
  - Exact current QC split on the latest rerender: self `FAIL`, independent subagent `FAIL`, Gemini `FAIL` (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-rerun-20260321-124412.md`).
  - Current blockers: too much negative space for the delivered evidence area, top summary / smallest text still conservative-print fragile, CFD badges still intrude on the flow field, and the structure-vs-CFD gap is clearer than before but still not fully self-evident without badge support.
  - Working QC log: `/home/yhjoo/projects/OpenFOAM/docs/figure_qc/26-03-21_figure10_dense_composite_gap_qc.md`
  - 2026-03-21 late-afternoon follow-up: tightened margins/gutters and moved structural/CFD badges plus `(a)–(d)` labels out of the panel interiors to reduce data occlusion. This improved rule compliance, but the exact new rerender still self-fails because the flow column remains visually undersized, the outside label stack now locally collides with `Reference` / `Predicted`, and the figure still lacks a shared color scale.
  - 2026-03-21 evening follow-up: tried a more radical `2 x 3` simplification (`predicted geometry / reference CFD / predicted CFD`) and replaced the intrusive CFD badges with one per-row metric footer. The story became cleaner than the prior combined-flow layout, but the exact current asset still failed all three gates because cropped `comparison_1x2.png` panels leave inconsistent/clipped source-axis text, the geometry-side dashed cues are still not fully self-explanatory without caption support, and the per-row footer now competes too strongly with the CFD evidence.
  - Exact current QC split on this latest rerender: self `FAIL`, subagent `FAIL`, Gemini `FAIL` (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-rerun-20260321-130913.md`).
  - 2026-03-21 late-evening follow-up: tested a cleaner-source variant that removed the dashed reference overlays from geometry, widened the CFD columns, and switched to a less aggressive crop of each `comparison_1x2.png` right-half so the slice axes and per-panel colorbar survive while the long internal title is cropped away. This fixed the worst clipping from the prior revision, but the exact new asset still self-fails: inherited axis/colorbar tick text remains too small for conservative print, the row-metric badge still crowds the top border area, there is still no trustworthy shared CFD normalization visible across panels, and the simplified geometry cues are still not fully self-explanatory without a legend.
  - Exact current QC split on this latest rerender: self `FAIL`; external QC was not advanced because the new asset did not clear the self-gate.
  - Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 night follow-up: replaced the crop-based CFD panels with direct mid-plane VTK slice renders plus a shared colorbar by rewriting `scripts/paper_figures/make_figure10_dense_composite_gap.py`. This removed the inherited OpenFOAM titles/per-panel colorbars and made the structure-vs-CFD gap materially clearer.
  - Exact current QC split on this newest rerender: self `FAIL`, subagent `FAIL`, Gemini `FAIL` (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-current-20260321-134113.md`).
  - Current blockers on the exact latest asset: shared colorbar tick labels still collide at print scale; the top legend/header cluster is cramped; row-metric text remains too small/light; geometry panels are undersized versus CFD panels; dashed white overlays may lose contrast in print; PDF-rendered text looks slightly thinner/softer than PNG.
  - Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 late-night follow-up: tried a typography/layout pass that widened the geometry column, switched the shared CFD map from `turbo` to `viridis`, thickened dashed obstacle overlays, and moved the legend into panel (a) to free the top lane. The exact new rerender still self-fails because the in-panel legend now occludes the geometry evidence, the per-row metric string is still too long/tight for conservative print, and the right-side colorbar typography remains small. External QC was not advanced on this exact asset because the new revision did not clear the self-gate.
  - 2026-03-21 latest follow-up: rerouted the legend and colorbar fully outside the panel interiors (bottom-left legend lane + bottom horizontal shared colorbar), and split the row metrics into shorter multi-line blocks. This removed the worst in-panel occlusion, but the exact current asset still self-fails because the top-left header stack still collides with the upper row title area, the bottom legend/colorbar lane remains too tight/busy at conservative print scale, and the figure still needs too much metric/colorbar reading before the CFD-fidelity gap becomes obvious. Figure 10 stays incomplete; external QC was intentionally not advanced on this exact rerender.
  - 2026-03-21 latest-latest follow-up: removed the bottom legend lane entirely, removed the top-left aggregate header, and replaced the dashed CFD obstacle mask with semi-transparent obstacle-colored overlays so the figure depends less on explicit legend reading. The exact newest asset is cleaner, but it still self-fails because the per-row metric line remains too small/long for conservative print, the structure-vs-CFD gap still depends too much on header metrics, and the CFD overlays still slightly veil the predicted flow field. Figure 10 remains incomplete and was not advanced to external QC.
  - 2026-03-21 final afternoon pass: rewrote the top lane into the more direct contrast pair `high structural match` / `weaker CFD fidelity`, collapsed the old long row-metric sentence into two short score badges, and reduced CFD overlay opacity. This improved first-glance message hierarchy, but the exact current asset still self-fails because the row label/badge strip above the geometry panels is still cramped/colliding and the figure still leans too much on explicit `CFD` badges instead of letting the mismatch read from the panel bodies alone. External QC was intentionally not advanced on this exact asset.
  - 2026-03-21 late-afternoon simplification pass: removed the figure-level message/badge lane entirely, switched to plain column headers (`predicted geometry`, `reference |U|`, `predicted |U|`), moved row IDs into the left margin, and kept only the shared bottom colorbar. The exact new PNG/PDF pair is materially cleaner and cross-format-consistent; self-QC is now `PASS` on this exact asset. Gemini exact-asset QC was rerun on copied, non-ignored review images (`docs/figure_qc/fig10_exact_asset_20260321*.png`) and returned `FAIL`, mainly flagging panel-label/geometry overlap and row-alignment concerns. A subagent QC was attempted in-turn but the verdict was not recoverable under current session-visibility limits, so Figure 10 remains incomplete and must not advance or be copied to Drive yet.
  - 2026-03-21 evening follow-up: addressed exactly those two external concerns by moving row IDs above the geometry panels and forcing shared row-wise plotting extents across geometry / reference `|U|` / predicted `|U|`. The exact new asset self-passes and the prior label-overlap / row-alignment complaints are materially reduced; independent subagent exact-asset QC now explicitly confirms the labels are clean and PNG/PDF consistency is good, but still returns `FAIL` because the `structure vs CFD gap` remains too implicit without caption support and the geometry column still feels underweighted relative to the CFD columns. A local Gemini CLI rerun was captured (`.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-152259.md`), but it referenced stale removed elements (`global legend`, `aggregate metrics`) and therefore does not count as a valid exact-current-asset Gemini gate. Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 late-evening follow-up: pushed one more exact-asset rerender aimed only at the surviving message-strength blockers. The geometry column was reweighted (`width_ratios=[1.24, 1.34, 1.34]`), room/obstacle strokes were slightly strengthened, predicted-CFD obstacle overlays were lightened again (`alpha=0.04`), and colorbar/header typography was nudged up for print safety. Self-QC on the exact new asset remains `PASS`, but independent subagent QC still returns `FAIL`: the figure is cleaner and less occluded, yet the intended read `structure is close / CFD still off` is still not fully self-evident without caption help, and the colorbar maximum tick remains borderline near the right edge. One Gemini rerun again produced a stale-path / non-exact critique and was discarded; a corrected repo-local Gemini exact-asset rerun was launched but no valid verdict was recovered within this turn. Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 final rerender in this turn: tested a reference-outline revival pass on the same exact asset. The geometry panels now include short gray dashed reference room/obstacle outlines, the three columns were rebalanced closer to parity (`width_ratios=[1.30, 1.29, 1.29]`), geometry strokes were softened slightly, and the shared colorbar was reduced to endpoint ticks only. This did **not** close the self-gate: the dashed reference cue is still too faint/inconsistent to carry the `structure close` claim at print scale, PNG vs PDF-render diverged again around the top-left header/spacing, and the endpoint tick labels still sit too close to the colorbar edges in the PDF-render check. External QC was intentionally not advanced on this exact asset. Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 late-night exact-current rerender: pushed one more message-first pass on the same locked pair by widening the geometry column (`width_ratios=[1.38, 1.24, 1.24]`), switching the column headers to `close geometry` / `reference |U|` / `shifted |U|`, shortening the left cue to `solid pred. · dashed ref.`, and re-anchoring the colorbar endpoint labels. Exact latest split is still open-fail: self `PASS`, subagent `FAIL`, Gemini `FAIL` (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-161354.md`). Newly stabilized point: PNG/PDF consistency and colorbar edge clipping are acceptable again. Remaining blockers: the top-left header/subtitle lane is still cramped, the first-column header block is vertically misaligned against the other columns, and the figure still does not make the `close geometry but shifted CFD` contrast sufficiently self-evident for strict external gates. Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 16:23 exact-current rerender: replaced the old single overlay-style geometry column with an explicit `ref. geometry` / `pred. geometry` pair, yielding a `2 x 4` comparison (`ref. geometry`, `pred. geometry`, `reference |U|`, `predicted |U|`) while keeping the same locked A4 pair and the same scale-calibrated source. This materially improved caption-free message clarity and removed the old cramped subtitle/header problem. Exact current split on the new asset: self `FAIL`, subagent `FAIL`, Gemini `PASS` (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-162350.md`). Remaining blockers are now print-scale rather than conceptual: the two geometry panels and the top/colorbar typography are still too small relative to the CFD panels for a conservative double-column gate. Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 16:34 exact-current rerender: stayed on the same `2 x 4` concept and pushed only a print-scale rebalance pass. The overall figure height increased slightly, the four columns were rebalanced closer to parity (`width_ratios=[1.16, 1.16, 1.18, 1.18]`), top headers / row labels / shared colorbar typography were enlarged, and the geometry room/obstacle strokes were strengthened. This improved local print readability, but the exact current asset still self-fails: the geometry pair remains visually weaker than the CFD pair, PNG vs PDF-render still show small vertical placement drift around the header-to-panel gap / bottom colorbar lane, and overall vertical whitespace is still slightly uneven. A temporary colorbar-end spacing experiment was reverted after it produced an incorrect endpoint-tick read during self review. External QC was intentionally not advanced on this exact asset. Figure 10 remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 16:56 exact-current rerender: pushed one narrow print-safety pass on the same locked `2 x 4` concept. The geometry pair was enlarged (`width_ratios=[1.26, 1.26, 1.10, 1.10]`), the figure got slightly taller/tighter, and the shared colorbar ticks were rewritten for the actual low-speed range so the old misleading `0.0 ... 0.0` endpoint read is gone. Exact latest QC split: self `FAIL`, subagent `PASS`, Gemini `PASS` (`/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-fig10-visual-qc-exact-20260321-165608.md`). I am keeping the self gate conservative because the shared colorbar label/ticks still feel borderline small at strict double-column print and the PNG/PDF pair still shows slight top-lane / colorbar-lane spacing drift. Figure 10 therefore remains incomplete and must **not** be copied to Drive yet.
  - 2026-03-21 17:xx exact-current micro-pass: kept the same locked `2 x 4` concept and only touched print-safety details. Geometry strokes/opening markers were strengthened, column widths shifted to `[1.28, 1.28, 1.08, 1.08]`, the top lane was compressed slightly, and shared colorbar typography was increased; the script now also refreshes the PDF-render check automatically via `pdftoppm`. Exact current self state is still `FAIL`: cross-format stability is materially better, but the geometry pair still reads a bit too weak relative to the CFD pair under conservative double-column print review. External QC was intentionally not advanced on this exact rerender because the self gate did not fully close. Figure 10 remains incomplete and must **not** be copied to Drive yet.

## Figure 2 — Benchmark dataset design and difficulty matrix
- [ ] representative A1/A2/A3/A4 cases 선정
- [ ] single/double-column 결정
- [ ] 2x2 matrix layout script 작성
- [ ] PDF/PNG 생성
- [ ] 자체 visual QC
- [ ] 서브에이전트 멀티모달 QC
- [ ] Gemini 멀티모달 QC
- [ ] Drive 저장

## Figure 3 — Multi-view rendering protocol
- [ ] representative case 선정
- [ ] rectangular/composite 2-row 필요 여부 결정
- [ ] single/double-column 결정
- [ ] 5-view panel 조립 script 작성
- [ ] PDF/PNG 생성
- [ ] 자체 visual QC
- [ ] 서브에이전트 멀티모달 QC
- [ ] Gemini 멀티모달 QC
- [ ] Drive 저장

## Figure 4 — Reference benchmark and evaluation pathway
- [ ] reference/evaluation flow 구성요소 확정
- [ ] single/double-column 결정
- [ ] pipeline diagram script 작성
- [ ] PDF/PNG 생성
- [ ] 자체 visual QC
- [ ] 서브에이전트 멀티모달 QC
- [ ] Gemini 멀티모달 QC
- [ ] Drive 저장

## Figure 7 — Representative cross-view outcome panel
- [ ] representative cases 선정
- [ ] good/weak/topology-sensitive/CFD-sensitive panel 구성 설계
- [ ] single/double-column 결정
- [ ] composition script 작성
- [ ] PDF/PNG 생성
- [ ] 자체 visual QC
- [ ] 서브에이전트 멀티모달 QC
- [ ] Gemini 멀티모달 QC
- [ ] Drive 저장

## Figure 11 — Robustness and convergence summary
- [ ] 포함 지표(성공/실패/stabilization) 확정
- [ ] single/double-column 결정
- [ ] plotting/diagram script 작성
- [ ] PDF/PNG 생성
- [ ] 자체 visual QC
- [ ] 서브에이전트 멀티모달 QC
- [ ] Gemini 멀티모달 QC
- [ ] Drive 저장

## Figure 1 — Overall framework with mandatory scale calibration (PaperBanana)
- [ ] PaperBanana용 spec 문서 작성
- [ ] mandatory scale calibration block 반영
- [ ] single/double-column 결정
- [ ] PaperBanana 생성 실행
- [ ] PDF/PNG 정리
- [ ] 자체 visual QC
- [ ] 서브에이전트 멀티모달 QC
- [ ] Gemini 멀티모달 QC
- [ ] Drive 저장
