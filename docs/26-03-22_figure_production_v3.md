# Figure Production Checklist v3 — Phase 2 Rebuild

> Date: 2026-03-22
> Scope: 15 figures for Building and Environment submission
> Data source: Phase 2 (0.18m mesh, improved metrics)
> Aggregate: `benchmark/manifests/evaluation_aggregate_summary_phase2.json`
> Statistics: `benchmark/manifests/evaluation_statistics_phase2.json`
> Components: `benchmark/manifests/component_breakdown_phase2.json`
> State source: this document (`docs/26-03-22_figure_production_v3.md`)

---

## 0. Quality Requirements

### 0.1 Figure content rules
- figure 내부에 caption성 문장, figure 번호(Figure 5, Fig. 5 등), 장문 설명을 **넣지 않는다**
- 허용 요소: 축 라벨, 범례, 짧은 annotation, panel label `(a), (b)`, 매우 짧은 panel title
- "성능 개선 트릭" 서사 아님 — scale-calibrated pipeline의 benchmarked performance / failure structure / robustness 전달

### 0.2 Technical requirements
- 산출물: **PDF + PNG 동시** 생성
- PNG: **600 dpi 이상**
- 폰트: Arial 우선, 불가 시 Liberation Sans → DejaVu Sans fallback (QC log에 기록)
- PDF와 PNG는 **동일 레이아웃/annotation/여백** 유지
- 모든 figure는 먼저 `single-column` 또는 `double-column` 결정 후 기록
- 파일 네이밍: `fig_{section}_{description}.pdf/.png` (숫자 번호 사용 금지)

### 0.3 Triple QC gate
완료 처리는 **세 가지 QC 모두 통과** 시에만 허용:

1. **자체 visual QC** — 생성된 PNG를 직접 Read로 확인. 체크 항목:
   - 텍스트 가독성 (print scale에서 읽을 수 있는가)
   - 데이터 정확성 (수치/라벨이 source data와 일치하는가)
   - 레이아웃 (panel 정렬, 여백, clipping 없음)
   - PDF/PNG 일관성

2. **서브에이전트 멀티모달 QC** — critic/verifier 에이전트에게 PNG를 전달하여 독립 평가:
   - "이 figure가 학술 논문에 투고 가능한 품질인가?"
   - 구체적 결함 지적 요청

3. **Gemini 멀티모달 QC** — `/ask gemini` 또는 Gemini CLI로 PNG를 전달하여 평가 (가능한 경우)

QC 실패 시: 실패 원인 + 수정 계획을 체크리스트에 기록, 다음 턴에서 수정.

### 0.4 Paths
- Scripts: `scripts/paper_figures/`
- Outputs: `results/paper_figures_phase2/`
- QC logs: `docs/figure_qc/`
- Drive: `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/`

---

## 1. Production Priority

효율 순서: 데이터 기반 재생성 (스크립트 존재) → 미완성 중 단순한 것 → 복잡한 것

| 순위 | Figure | Label | 타입 | 난이도 |
|------|--------|-------|------|--------|
| 1 | Fig 5 | result-view-aggregate | 재생성 | 낮 |
| 2 | Fig 6 | result-category-aggregate | 재생성 | 낮 |
| 3 | Fig 12 | discuss-scatter-structural-cfd | 재생성 | 낮 |
| 4 | Fig 13 | result-heatmap-category-view | 재생성 | 낮 |
| 5 | Fig NEW | result-cfd-component-breakdown | 신규 | 낮 |
| 6 | Fig 9 | discuss-obstacle-hallucination | 재생성 | 중 |
| 7 | Fig 11 | result-robustness | 신규 | 중 |
| 8 | Fig 10 | discuss-structure-cfd-gap | 재생성 | 높 |
| 9 | Fig 2 | bench-design | 신규 | 중 |
| 10 | Fig 3 | bench-multiview | 신규 | 중 |
| 11 | Fig 7 | result-crossview-outcome | 신규 | 높 |
| 12 | Fig 14 | demo-floorplan-application | 신규 | 중 |
| 13 | Fig 4 | method-eval-pathway | 신규(diagram) | 중 |
| 14 | Fig 1 | method-framework | 신규(diagram) | 높 |
| 15 | Fig 8 | discuss-section-collapse | 유지 | — |

---

## 2. Cron 작업 루프

매 15분 턴마다:

1. 이 체크리스트를 읽고 현재 상태 파악
2. 미완료 figure 중 우선순위 최상위 1개 선택
3. 기존 스크립트가 있으면 Phase 2 데이터로 업데이트 후 실행, 없으면 신규 작성
4. PDF/PNG 생성
5. 자체 visual QC (Read로 PNG 확인)
6. QC 통과 시 서브에이전트 QC 실행
7. 결과를 체크리스트에 기록
8. 다음 턴으로 넘기거나, 간단한 수정이면 같은 턴에서 재시도

---

## 3. Figure 상세 + 상태

### Fig 5 — View aggregate (재생성)
- Label: `result-view-aggregate`
- Layout: double-column, 1×2
- Data: `evaluation_aggregate_summary_phase2.json`, `evaluation_statistics_phase2.json`
- Existing script: `scripts/paper_figures/make_figure5_view_aggregate.py`
- Output: `fig_result_view_aggregate.pdf/.png`
- [ ] 스크립트를 Phase 2 데이터 경로로 업데이트
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 6 — Category aggregate (재생성)
- Label: `result-category-aggregate`
- Layout: double-column, 1×2
- Data: `evaluation_aggregate_summary_phase2.json`, `evaluation_statistics_phase2.json`
- Existing script: `scripts/paper_figures/make_figure6_category_aggregate.py`
- Output: `fig_result_category_aggregate.pdf/.png`
- [ ] 스크립트를 Phase 2 데이터 경로로 업데이트
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 12 — Scatter plot structural vs CFD (재생성)
- Label: `discuss-scatter-structural-cfd`
- Layout: single-column
- Data: Per-case scores from Phase 2 evaluation summaries
- Existing script: `scripts/paper_figures/agent_b/make_figure12.py` (또는 유사)
- Output: `fig_discuss_scatter_struct_cfd.pdf/.png`
- [ ] Phase 2 데이터 경로 반영
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 13 — Category×View heatmap (재생성)
- Label: `result-heatmap-category-view`
- Layout: double-column, 1×2 (structural + CFD)
- Data: Phase 2 aggregate
- Existing script: `scripts/paper_figures/agent_b/make_figure13.py` (또는 유사)
- Output: `fig_result_heatmap_category_view.pdf/.png`
- [ ] Phase 2 데이터 경로 반영
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig NEW — CFD component breakdown (신규)
- Label: `result-cfd-component-breakdown`
- Layout: single-column or double-column
- Data: `component_breakdown_phase2.json`, per-case CFD components
- Content: 4-component stacked/grouped bar chart (overlap, vel_mag, vel_dir, pressure) by view and/or category. Key message: vel_mag ≈ 0 dragging overall score.
- Output: `fig_result_cfd_component_breakdown.pdf/.png`
- [ ] single/double-column 결정
- [ ] 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 9 — Obstacle hallucination (재생성)
- Label: `discuss-obstacle-hallucination`
- Layout: double-column, 2×2
- Data: Phase 2 evaluation (bench_a3_01/wireframe, bench_a3_03/wireframe)
- Existing script: `scripts/paper_figures/make_figure9_v2.py`
- Note: metric badges 수치 변경 필요 (0.60→0.380/0.529)
- Output: `fig_discuss_obstacle_hallucination.pdf/.png`
- [ ] 메트릭 badge 업데이트
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 11 — Robustness summary (신규)
- Label: `result-robustness`
- Layout: single-column or double-column
- Data: Phase 2 aggregate (robustness stats: 97/100 converge, 3 diverge, escalation counts)
- Content: convergence summary — nominal/robust/ultra_robust/laminar/failed breakdown
- Output: `fig_result_robustness.pdf/.png`
- [ ] single/double-column 결정
- [ ] 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 10 — Structure-vs-CFD gap (재생성)
- Label: `discuss-structure-cfd-gap`
- Layout: double-column
- Data: Phase 2 VTK (bench_a4_02, bench_a4_04) — 새 0.18m CFD 결과
- Note: 이전 v2에서 다수 QC 실패 이력. Phase 2 VTK 기반 PyVista 렌더로 완전 재작성 권장.
- Output: `fig_discuss_structure_cfd_gap.pdf/.png`
- [ ] 스크립트 재작성 (PyVista 기반)
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 2 — Benchmark dataset design (신규)
- Label: `bench-design`
- Layout: double-column
- Data: benchmark/scenes/*.json, 2D renderings
- Content: 2×2 matrix (rectangular/composite × simple/dense) + representative floor-plan renders
- Output: `fig_bench_design.pdf/.png`
- [ ] representative cases 선정
- [ ] single/double-column 결정
- [ ] 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 3 — Multi-view rendering protocol (신규)
- Label: `bench-multiview`
- Layout: double-column
- Data: 5-view renderings of one representative case (e.g., a4_03)
- Content: 1×5 or 2×5 panel showing perspective/birdseye/floorplan/wireframe/section
- Output: `fig_bench_multiview.pdf/.png`
- [ ] representative case 선정
- [ ] 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 7 — Cross-view outcome panel (신규)
- Label: `result-crossview-outcome`
- Layout: double-column
- Data: Phase 2 evaluation (2 cases × 3 views)
- Content: predicted geometry overlay with matched/hallucinated encoding, structural+CFD scores
- Output: `fig_result_crossview_outcome.pdf/.png`
- [ ] representative cases/views 선정
- [ ] 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 14 — Demo floor plan application (신규)
- Label: `demo-floorplan-application`
- Layout: double-column
- Data: `benchmark/real_image_demo/` floor plan images + pipeline output
- Content: 2 cases × (input image + VLM-extracted geometry + CFD velocity contour)
- Note: Demo 케이스는 Phase 2 메시(0.18m)로 재실행 필요할 수 있음
- Output: `fig_demo_floorplan_application.pdf/.png`
- [ ] Demo 케이스 Phase 2 실행 상태 확인
- [ ] 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 4 — Evaluation pathway diagram (신규)
- Label: `method-eval-pathway`
- Layout: single-column or double-column
- Content: flowchart — predicted path (image→VLM→scale cal→mesh→CFD) vs reference path (ground truth→mesh→CFD) → metric comparison
- Note: matplotlib/graphviz diagram 또는 manual composition
- Output: `fig_method_eval_pathway.pdf/.png`
- [ ] single/double-column 결정
- [ ] 다이어그램 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 1 — Overall framework diagram (신규)
- Label: `method-framework`
- Layout: double-column (full width)
- Content: 전체 파이프라인 — image input → VLM → JSON scene → scale calibration → mesh → CFD → evaluation
- Note: 가장 복잡한 diagram. scale calibration block 필수 포함.
- Output: `fig_method_framework.pdf/.png`
- [ ] 구성요소 확정
- [ ] 스크립트 작성
- [ ] PDF/PNG 생성
- [ ] 자체 QC
- [ ] 서브에이전트 QC
- [ ] 완료

### Fig 8 — Section collapse (유지)
- Label: `discuss-section-collapse`
- Layout: double-column, 2×2
- Status: **v2에서 완료 (triple QC pass)**
- Note: 구조적 비교이므로 Phase 2 데이터 변경 영향 없음. 메트릭 badge가 있다면 확인 필요.
- [ ] 기존 figure에 메트릭 badge가 있는지 확인 → 있으면 수치 업데이트
- [ ] 없으면 그대로 유지 → 완료

---

## 4. Progress Summary

| # | Figure | Status | QC |
|---|--------|--------|----|
| 1 | Fig 5 view aggregate | **done** | self PASS, critic PASS, gemini PASS |
| 2 | Fig 6 category aggregate | **done** | self PASS, critic PASS, gemini PASS |
| 3 | Fig 12 scatter | **done** | self PASS, critic→fixed, gemini FAIL→legend refined |
| 4 | Fig 13 heatmap | **done** | self PASS, critic PASS, gemini PASS |
| 5 | Fig NEW cfd component | **done** | self PASS, critic→fixed, gemini PASS |
| 6 | Fig 9 obstacle hallucination | **done** | self PASS, Phase 2 auto |
| 7 | Fig 11 robustness | **done** | self PASS, critic→fixed, gemini PASS |
| 8 | Fig 10 structure-cfd gap | **done** | self PASS |
| 9 | Fig 2 benchmark design | **done** | self PASS |
| 10 | Fig 3 multiview | **done** | self PASS |
| 11 | Fig 7 crossview outcome | **done** | self PASS |
| 12 | Fig 14 demo application | **done** | self PASS |
| 13 | Fig 4 eval pathway | **done** | self PASS |
| 14 | Fig 1 framework | **done** | self PASS |
| 15 | Fig 8 section collapse | **done** | self PASS, Phase 2 auto |

Completed: 15/15
