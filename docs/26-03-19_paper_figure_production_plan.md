# Paper Figure Production Plan

> 작성일: 2026-03-19
> 프로젝트: OpenFOAM Image-to-CFD paper
> 목적: 논문용 고품질 벡터 PDF figure 후보 12종의 실제 생산 가능 여부를 점검하고, 즉시 생산 가능한 항목부터 자율적으로 제작/QC/체크오프하기 위한 작업 기준을 정리한다.

## 공통 원칙

- 최종 납품 형식은 기본적으로 `PDF + PNG` 동시 생성으로 한다.
- 가능하면 **벡터 PDF**를 우선한다. (matplotlib/pdf backend, SVG→PDF 등)
- 논문 본문에는 **벡터 PDF를 사용**하고, 이후 **LaTeX → PDF 빌드 경로**로 포함하는 것을 기본 원칙으로 한다.
- figure는 기본적으로 **2-column journal 호환** 스타일을 따른다.
- 각 figure마다 반드시 **1단 너비(single-column)** 로 갈지, **2단 너비(double-column)** 로 갈지를 먼저 결정하고 기록한다.
- subfigure가 있는 경우, 각 패널 라벨 `(a), (b), ...``, subfigure caption 처리 방식, 메인 caption과의 역할 분리를 먼저 설계한다.
- 모든 PNG 산출물은 **최소 600 dpi 이상**으로 저장한다.
- PDF와 PNG는 가능한 한 **동일 레이아웃 / 동일 여백 / 동일 annotation** 을 유지해야 한다.
- 축 라벨, 범례, tick, annotation, panel label의 폰트 크기는 논문 본문 삽입 후에도 읽히는 수준인지 반드시 검토한다.
- 기본적으로 single-column 삽입을 우선 고려하되, 패널 수가 많거나 비교 메시지가 강하게 깨지는 경우에만 double-column figure를 허용한다.
- figure 생성 직후에는 반드시 **visual QC**를 거친다.
- visual QC는 최소 3단계로 한다.
  1. 에이전트 자체 시각 검토(이미지/패널 레이아웃/라벨/범례/가독성/잘림/왜곡)
  2. 서브에이전트 기반 별도 시각 QC
  3. Gemini CLI 기반 외부 시각 QC
- **QC가 통과되기 전에는 figure 체크박스를 완료 처리하지 않는다.**
- 단순 생성 성공과 논문 품질 통과를 구분한다.

## Figure layout / typography 기준

- 기본 목표는 저널 2단 편집 기준에서 **재조정 없이 바로 LaTeX에 삽입 가능한 상태**다.
- 각 figure마다 아래 항목을 production log 또는 QC log에 명시한다.
  - intended width: `single-column` 또는 `double-column`
  - panel layout: 예) `1x1`, `1x2`, `2x2`, `2x3`
  - subfigure labels 사용 여부
  - main caption과 subcaption 분리 방식
  - export size (inch 또는 mm)
  - PNG dpi
  - PDF backend / vector 보존 여부
- 가능한 경우 글꼴은 figure 간 일관성을 유지한다.
- 축/범례/annotation이 많아 single-column에서 무너지면, 억지 축소 대신 double-column으로 승격한다.
- subfigure caption은 장황하게 쓰지 말고, 패널 내부 라벨과 메인 caption이 자연스럽게 이어지도록 설계한다.
- 색상, 선 두께, 마커 크기, 여백은 single-column 인쇄/축소 후에도 구분 가능해야 한다.

## 생산 상태 표기

- [ ] 미착수
- [-] 제작 중
- [x] 제작 + visual QC 완료
- [~] 제작 가능하지만 추가 조립/스크립트 보강 필요
- [!] 현재는 플레이스홀더 또는 추가 자산 필요

---

## Figure 1. Overall framework diagram (PaperBanana placeholder)

**의도**
- 논문 전체 프레임워크를 한 장으로 보여주는 메인 기술 구성도

**현재 생산 가능 여부**
- [!] 최종본은 아직 아님
- 이유: 사용자가 명시했듯 이 항목은 `paperbanana` 스킬로 별도 제작 예정

**현재 활용 가능한 자산**
- `results/indoor_pipeline_g31/indoor_pipeline_3d_comparison_v4.pdf`
- `results/indoor_pipeline_g33/indoor_pipeline_3d_comparison.pdf`
- `scripts/render_indoor_pipeline_3d.py`

**판정**
- 최종 figure는 아직 placeholder
- 단, paperbanana 제작 전 임시 설명용 레이아웃은 필요 시 생성 가능

**체크리스트**
- [ ] PaperBanana용 텍스트 spec 정리
- [ ] 기술 단계 박스/화살표/오류복구 단계 포함 여부 확정
- [ ] 최종 벡터 PDF 제작
- [ ] 자체 visual QC
- [ ] Gemini CLI visual QC

---

## Figure 2. Benchmark dataset design and difficulty matrix

**의도**
- A1/A2/A3/A4 2×2 난이도 체계를 대표 예시와 함께 보여줌

**현재 생산 가능 여부**
- [~] 즉시 생산 가능에 가까움

**근거 자산**
- `benchmark/manifests/scene_manifest.json`
- `benchmark/renderings/`
- `scripts/render_benchmark_views.py`
- category별 대표 case가 이미 확보됨

**필요 작업**
- category당 representative thumbnail 1개씩 선택
- 2×2 matrix 레이아웃 조립
- rectangular/composite, simple/dense 축 라벨링

**판정**
- 추가 실험 없이 제작 가능
- 새 조립 스크립트가 필요할 가능성 높음

**체크리스트**
- [ ] representative cases 선정 (A1/A2/A3/A4)
- [ ] 2×2 matrix 조립 스크립트 작성
- [ ] PDF/PNG 출력
- [ ] 자체 visual QC
- [ ] Gemini CLI visual QC

---

## Figure 3. Multi-view rendering protocol for each benchmark case

**의도**
- 하나의 장면이 5개 입력 뷰로 어떻게 변환되는지 설명

**현재 생산 가능 여부**
- [~] 즉시 생산 가능에 가까움

**근거 자산**
- `benchmark/renderings/` 내 5-view bundle
- `scripts/render_benchmark_views.py`
- `docs/26-03-18_benchmark_rendering_feasibility.md`

**필요 작업**
- 대표 rectangular case 1개, composite case 1개를 고를지 여부 결정
- 5-view 패널 정렬 및 공통 캡션/라벨 부여

**판정**
- 매우 생산 가능성이 높음

**체크리스트**
- [ ] representative case 선정
- [ ] 5-view 패널 조립
- [ ] 필요 시 rectangular/composite 2-row 버전 제작
- [ ] PDF/PNG 출력
- [ ] 자체 visual QC
- [ ] Gemini CLI visual QC

---

## Figure 4. Reference benchmark pipeline and artifact bundle

**의도**
- rule-based scene → reference CFD → render bundle → evaluation manifest의 reference path를 설명
- VLM 비개입 정답 경로라는 점을 강조

**현재 생산 가능 여부**
- [~] 생산 가능하나 개념도 조립이 필요

**근거 자산**
- `docs/26-03-17_image_to_cfd_paper_plan.md`
- `docs/26-03-18_benchmark_dataset_card.md`
- `benchmark/manifests/*`
- `scripts/build_benchmark_dataset_card.py`
- `scripts/verify_benchmark_bundle.py`

**필요 작업**
- artifact tree 요약
- reference-only generation path를 다이어그램으로 도식화
- 아마도 matplotlib + text box 기반 벡터 figure 제작 필요

**판정**
- 새 figure script 작성이 필요하지만 실험 자산은 충분함

**체크리스트**
- [ ] pipeline 박스 구성 확정
- [ ] artifact tree/summary 삽입 방식 결정
- [ ] 벡터 도식 스크립트 작성
- [ ] PDF/PNG 출력
- [ ] 자체 visual QC
- [ ] Gemini CLI visual QC

---

## Figure 5. Aggregate performance across input views

**의도**
- view별 structural/CFD 성능 차이를 headline figure로 제시

**현재 생산 가능 여부**
- [x] 제작 + visual QC 완료

**근거 자산**
- `benchmark/manifests/evaluation_aggregate_summary.json`
- `docs/26-03-19_cli_eval_aggregate_results.md`

**이미 확보된 핵심 메시지**
- floorplan strongest
- section weakest
- perspective less stable than expected

**필요 작업**
- plot script 작성 또는 aggregate script 확장
- single-column / double-column version 중 선택

**판정**
- 제작 완료. 현재 버전은 manuscript 삽입 가능한 production candidate로 승인함.

**산출물 / 로그**
- script: `scripts/paper_figures/make_figure5_view_aggregate.py`
- outputs:
  - `results/paper_figures/figure5_view_aggregate_performance.pdf`
  - `results/paper_figures/figure5_view_aggregate_performance.png`
- QC log: `docs/figure_qc/26-03-19_figure5_view_aggregate_qc.md`

**체크리스트**
- [x] structural / CFD / opening-wall match 중 표시 지표 확정
- [x] plotting script 작성
- [x] PDF/PNG 출력
- [x] 자체 visual QC
- [x] 서브에이전트 visual QC
- [x] Gemini CLI visual QC

---

## Figure 6. Aggregate performance across benchmark categories

**의도**
- A1/A2/A3/A4 category별 성능/해석 차이 제시

**현재 생산 가능 여부**
- [x] 제작 + visual QC 완료

**근거 자산**
- `benchmark/manifests/evaluation_aggregate_summary.json`
- `docs/26-03-19_cli_eval_aggregate_results.md`
- `docs/26-03-19_cli_eval_case_comparison_100of100.md`

**이미 확보된 핵심 메시지**
- A1 positive control
- A2 opening/blockage 민감
- A3 structure-vs-CFD decoupling
- A4 dense composite structure-physics gap

**판정**
- 제작 완료. 현재 버전은 manuscript 삽입 가능한 production candidate로 승인함.

**산출물 / 로그**
- script: `scripts/paper_figures/make_figure6_category_aggregate.py`
- outputs:
  - `results/paper_figures/figure6_category_aggregate_performance.pdf`
  - `results/paper_figures/figure6_category_aggregate_performance.png`
- QC log: `docs/figure_qc/26-03-19_figure6_category_aggregate_qc.md`

**체크리스트**
- [x] structural / CFD / opening-wall match 구성 확정
- [x] plotting script 작성
- [x] PDF/PNG 출력
- [x] 자체 visual QC
- [x] 서브에이전트 visual QC
- [x] Gemini CLI visual QC

---

## Figure 7. Representative success and failure cases across views

**의도**
- 같은 benchmark 구조라도 input-view에 따라 reconstruction/CFD 결과가 달라짐을 보여줌

**현재 생산 가능 여부**
- [~] 생산 가능

**근거 자산**
- `results/eval_bench_*/*` 아래 `comparison_1x2.png`, `indoor_pipeline_3d_comparison.pdf/png`, panel 이미지들
- per-task evaluation summary와 case notes 다수

**추천 대표 후보**
- strong case: `bench_a1_03`, `bench_a3_03`
- weak/stress case: `bench_a3_04`, `bench_a4_05`

**필요 작업**
- 어떤 케이스/뷰를 대표로 쓸지 선정
- 2×N panel 또는 4-panel 구성
- metric annotation을 과하지 않게 삽입

**판정**
- 즉시 조립 가능

**체크리스트**
- [ ] representative strong/weak case 확정
- [ ] source artifact 수집
- [ ] panel composition script 작성
- [ ] PDF/PNG 출력
- [ ] 자체 visual QC
- [ ] Gemini CLI visual QC

---

## Figure 8. Failure mode: composite-room collapse under section view

**의도**
- section view가 composite topology에 취약하다는 핵심 failure mode 시각화

**현재 생산 가능 여부**
- [x] 제작 + visual QC 완료

**근거 자산**
- `bench_a3_04/section`
- `bench_a4_05/section`
- 관련 aggregate note와 per-case notes
- 대응하는 reference/predicted scene JSON 및 evaluation summary 존재

**판정**
- 제작 완료. 현재 버전은 manuscript 삽입 가능한 production candidate로 승인함.

**산출물 / 로그**
- script: `scripts/paper_figures/make_figure8_section_collapse.py`
- outputs:
  - `results/paper_figures/figure8_section_view_composite_collapse.pdf`
  - `results/paper_figures/figure8_section_view_composite_collapse.png`
- QC log: `docs/figure_qc/26-03-19_figure8_section_collapse_qc.md`

**체크리스트**
- [x] 대표 case 2개 비교 선택 (`bench_a3_04/section`, `bench_a4_05/section`)
- [x] GT vs predicted geometry 패널 확정
- [x] CFD comparison 포함 여부 결정 (full CFD panel은 생략하고 structural/CFD penalty를 annotation으로 반영)
- [x] collapse annotation 추가
- [x] PDF/PNG 출력
- [x] 자체 visual QC
- [x] 서브에이전트 visual QC
- [x] Gemini CLI visual QC

---

## Figure 9. Failure mode: obstacle hallucination with limited CFD penalty

**의도**
- obstacle hallucination이 있어도 opening/topology가 맞으면 CFD가 유지될 수 있음을 보여줌

**현재 생산 가능 여부**
- [x] 제작 + visual QC 완료

**근거 자산**
- `bench_a3_01`, `bench_a3_02`, `bench_a3_05`
- aggregate tag: `opening_topology_preserved_with_obstacle_hallucination`
- 관련 evaluation outputs와 case notes 존재

**판정**
- 제작 완료. 현재 버전은 manuscript 삽입 가능한 production candidate로 승인함.

**산출물 / 로그**
- script: `scripts/paper_figures/make_figure9_obstacle_hallucination.py`
- outputs:
  - `results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.pdf`
  - `results/paper_figures/figure9_obstacle_hallucination_limited_cfd_penalty.png`
- QC log: `docs/figure_qc/26-03-19_figure9_obstacle_hallucination_qc.md`

**체크리스트**
- [x] 대표 A3 사례 확정 (`bench_a3_01/floorplan`, `bench_a3_05/floorplan`)
- [x] obstacle count mismatch와 CFD similarity를 함께 보여줄 패널 구성 설계
- [x] PDF/PNG 출력
- [x] 자체 visual QC
- [x] 서브에이전트 visual QC
- [x] Gemini CLI visual QC

---

## Figure 10. Dense composite regime: structural correctness vs CFD fidelity gap

**의도**
- dense composite에서는 구조 점수가 높아도 CFD fidelity가 낮을 수 있음을 보여줌

**현재 생산 가능 여부**
- [x] 생산 완료 (QC 통과)

**근거 자산**
- `bench_a4_01`, `bench_a4_02`, `bench_a4_04`, `bench_a4_05`
- aggregate tag: `dense_composite_structure_physics_gap`
- `bench_a4_04/floorplan` 등 대표 사례 확보

**판정**
- discussion용 핵심 figure로 생산 가능

**체크리스트**
- [x] 대표 A4 사례 확정
- [x] structural score high vs CFD moderate 사례 패널 구성
- [x] PDF/PNG 출력
- [x] 자체 visual QC
- [x] 서브에이전트 visual QC
- [x] Gemini CLI visual QC

---

## Figure 11. Robustness and stabilization pathway

**의도**
- repair salvage, mesh fallback, solver escalation, warning metadata를 요약

**현재 생산 가능 여부**
- [~] 생산 가능

**근거 자산**
- `benchmark/manifests/evaluation_aggregate_summary.json`
- `docs/26-03-18_frozen20_failure_signals.md`
- `docs/26-03-19_cli_eval_aggregate_results.md`

**핵심 메시지**
- 100/100 success라도 stabilization path는 중요한 engineering contribution

**판정**
- 정량+방법론 figure로 생산 가능

**체크리스트**
- [ ] 포함 항목 선정 (repair, fallback, ultra_robust, warnings)
- [ ] flowchart 또는 bar-summary 중 형식 결정
- [ ] plotting/diagram script 작성
- [ ] PDF/PNG 출력
- [ ] 자체 visual QC
- [ ] Gemini CLI visual QC

---

## Figure 12. Wild-image demonstration

**의도**
- synthetic benchmark 밖 실제 이미지에서의 pipeline feasibility를 정성적으로 시연

**현재 생산 가능 여부**
- [~] 제한적 생산 가능

**근거 자산**
- `results/full_pipeline_photo_run/comparison_1x2.png`
- `results/one_command_photo_run2/indoor_pipeline_3d_comparison.pdf`
- `results/indoor_pipeline_g31/indoor_pipeline_3d_comparison_v4.pdf`
- `results/indoor_pipeline_g33/indoor_pipeline_3d_comparison.pdf`
- photo-domain 관련 memory note 존재

**제약**
- benchmark만큼 정량적으로 닫힌 축은 아님
- 논문 본문에서는 qualitative demo로 다루는 편이 안전

**판정**
- appendix / late-results / qualitative demo로는 생산 가능
- 본문 핵심 정량 figure보다는 우선순위가 낮음

**체크리스트**
- [ ] 사용할 wild-image 사례 확정
- [ ] input / predicted scene / 결과 시각화 패널 선정
- [ ] 과도한 claim 없이 qualitative caption 정리
- [ ] PDF/PNG 출력
- [ ] 자체 visual QC
- [ ] Gemini CLI visual QC

---

# 종합 판정

## 즉시 자율 생산 우선순위 높음

1. Figure 5 — Aggregate performance across input views
2. Figure 6 — Aggregate performance across benchmark categories
3. Figure 8 — Composite-room collapse under section view
4. Figure 10 — Dense composite structure-vs-CFD gap
5. Figure 3 — Multi-view rendering protocol
6. Figure 7 — Representative success/failure cases across views
7. Figure 2 — Benchmark dataset design and difficulty matrix
8. Figure 11 — Robustness and stabilization pathway

## 보강 후 진행

10. Figure 4 — Reference benchmark pipeline and artifact bundle
11. Figure 12 — Wild-image demonstration

## placeholder 유지

12. Figure 1 — Overall framework diagram (PaperBanana)

---

# 자율 figure 생산 프로토콜 (cron용 요약)

각 figure에 대해 아래 순서를 따른다.

1. 후보 figure 1개 선택
2. source artifact / metrics / case를 확정
3. figure composition script 또는 plotting script 작성
4. `PDF + PNG` 출력 생성
5. 에이전트 자체 visual inspection 수행
6. 별도 서브에이전트 visual QC 수행
7. Gemini CLI 기반 visual QC 수행
8. QC 통과 시에만 해당 figure 체크박스를 `[x]`로 갱신
9. QC 실패 시에는 수정 후 재생성하고, 실패 원인을 문서에 짧게 남김

## QC 통과 기준

- 글자 크기와 축 라벨이 논문 삽입 시 읽힘
- 패널 정렬이 깨지지 않음
- 색상/선/범례가 중복되거나 과밀하지 않음
- 핵심 주장(예: section collapse, A3 hallucination-vs-CFD decoupling)이 그림만 보고도 읽힘
- PNG와 PDF 사이에 레이아웃 차이/잘림 없음
- CFD contour/geometry panel에 왜곡, 비정상 여백, 잘린 annotation이 없음

## 기록 위치 제안

- figure outputs: `results/paper_figures/`
- figure scripts: `scripts/paper_figures/`
- QC logs: `docs/figure_qc/`
- Google Drive delivery folder: `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/`
- 상태 관리 문서: 이 문서 자체를 갱신

## Google Drive 전달 규칙

- QC 3단계(자체 / 서브에이전트 / Gemini CLI)를 모두 통과한 figure만 Google Drive로 복사한다.
- Google Drive에는 **QC 통과본만** 들어가야 하며, 실패본/중간본/초안은 올리지 않는다.
- figure별로 가능하면 아래 구조를 유지한다.
  - `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig05_view_aggregate/`
  - `/mnt/c/Users/User/GoogleDrive/OpenFOAM_paper_figures_qc_passed/Fig06_category_aggregate/`
  - ...
- 각 figure 폴더에는 기본적으로 아래를 함께 둔다.
  - 최종 `PDF`
  - 최종 `PNG` (600 dpi 이상)
  - 짧은 `README.md` 또는 `qc_summary.md` (figure 의도, 배치 방식, QC 통과 메모)
- production plan 문서에서 체크 완료(`[x]`) 처리되기 전에는 Google Drive로 복사하지 않는다.
