# Benchmark Autopilot TODO

> 토픽: SurfClaw / topic 2242
> 모드: 사용자 의사결정이 꼭 필요할 때만 질문하고, 그 외에는 자율 진행
> Cron: `openfoam-benchmark-autopilot-topic2242`
> Job ID: `543fca81-0332-442a-beec-e8ef06e5e73d`
> 마지막 초기화: 2026-03-18

## 운영 원칙

- 실제로 사용자 판단이 필요한 경우가 아니면 자율적으로 계속 진행한다.
- 막히는 의사결정이 필요하면 먼저 cron job을 비활성화한 뒤 사용자에게 질문한다.
- 상태 기록은 짧고 행동 지향적으로 유지한다.
- 진행 보고와 질문은 한국어로 한다.

## 현재 목표

OpenFOAM Image-to-CFD 논문을 위한 benchmark dataset 파이프라인을 구축하고 검증한다.

1. rule-based benchmark scene 생성
2. reference CFD 생성
3. 결과 시각화
4. 3D→2D benchmark 입력 이미지 렌더링 정리
5. success/failure 분석 및 안정화
6. image-conditioned evaluation 준비 및 실행

## 상태 표기

- [ ] 대기
- [x] 완료
- [-] 진행중
- [!] 실패 / 수정 필요
- [~] 보류

## Active TODO

### A. Benchmark scene generation
- [x] `scripts/generate_benchmark_scenes.py` 구현
- [x] 초기 8개 benchmark scene 생성 (A1/A2/A3/A4 × 2)
- [x] 생성된 JSON scene 전체를 `validate_indoor_scene.py`로 검증
- [x] simple category에서 필요 시 obstacle 최소 1개를 보장하는 옵션 추가
- [x] scene별 category/seed/geometry manifest 추가
- [x] frozen scene set을 8개에서 12개로 확장 (A1/A2/A3/A4 × 3)
- [x] frozen scene set을 12개에서 16개로 확장 (A1/A2/A3/A4 × 4)
- [x] frozen scene set을 16개에서 20개로 확장 (A1/A2/A3/A4 × 5)

### B. Reference CFD pipeline
- [x] scene-JSON 기반 benchmark 실행에 stabilized pipeline entrypoint 재사용
- [x] batch wrapper 구현: `scripts/run_benchmark_reference_batch.py`
- [x] `run_indoor_stabilized.py`에서 `--end-time`이 case generation에 전달되도록 수정
- [x] rectangular benchmark scene이 CFD + 시각화까지 end-to-end 도달함을 확인
- [x] composite benchmark scene이 CFD + 시각화까지 end-to-end 도달함을 확인
- [x] 초기 8-scene pilot의 first-pass success rate 개선
- [x] 새로 추가된 `*_03` 4개 scene을 reference CFD batch flow로 실행
- [x] `*_04` tranche 실행
- [x] `*_05` tranche 실행

### C. Composite-room compatibility
- [x] validator와 Gmsh 경로에 `room.blocks` 지원 추가
- [x] image-based composite-room generation feasibility 확인
- [x] `room.blocks`에 대한 legacy repair-stage 비호환성 원인 파악
- [x] `run_indoor_stabilized.py`에서 repair-stage crash가 composite run을 막지 않도록 수정
- [x] `visualize_indoor_case.py`를 composite room에 맞게 수정
- [x] `render_indoor_pipeline_3d.py`를 composite room에 맞게 수정
- [x] repaired scene을 쓴 경우 solver가 실제 사용한 scene JSON 기준으로 렌더링되도록 보정
- [x] `repair_indoor_scene.py`에 composite-room 지원을 추가해 fallback/skip 의존도를 낮춤

### D2. Evaluation scaffolding
- [x] frozen-12 × 5-view bundle용 `benchmark/evaluations/` scaffold 생성
- [x] 60개 image-conditioned evaluation task용 manifest/symlink 생성
- [x] frozen-16 × 5-view bundle(80 tasks)로 scaffold 갱신
- [x] frozen-20 × 5-view bundle(100 tasks)로 scaffold 갱신
- [x] scaffolded evaluation task 하나를 end-to-end로 실행하고 `evaluation_summary.json`을 쓰는 runner 추가
- [x] geometry-only metric에서 normalized-grid VTK 비교 기반 CFD/result-side metric까지 확장
- [x] current cron shell에서 Gemini CLI cached auth를 활용하도록 image-conditioned path를 우회 복구 (`generate_indoor_scene_with_gemini.py`가 CLI multimodal `@image` 입력 지원, API env alias `GOOGLE_API_KEY`도 허용)
- [x] CLI backend로 image-conditioned smoke test 2개 재실행 성공: `bench_a1_01/perspective`, `bench_a1_01/section`
- [x] 대표 rectangular case `bench_a1_01`에 대해 5-view(`perspective`, `birdseye`, `floorplan`, `wireframe`, `section`) complete sweep 실행 및 per-view metric 확보
- [x] 대표 composite case `bench_a3_04`에 대해 5-view CLI complete sweep 실행 및 rectangular baseline과 1차 비교 신호 확보
- [x] obstacle-dense composite stress case `bench_a4_03`에 대해 5-view CLI complete sweep 실행 및 composite difficulty 신호 확장
- [x] laminar-fallback 계열 composite stress case `bench_a4_02`에 대해 5-view CLI complete sweep 실행 및 repair/solver escalation 신호 확보
- [x] 직사각형 dense-obstacle hard case `bench_a2_03`에 대해 5-view CLI complete sweep 실행 및 rectangular-vs-composite stress 비교축 확보
- [x] reference solver-stress 대표 case `bench_a1_04`에 대해 5-view CLI complete sweep 실행 및 prediction-side stress 완화 신호 확보
- [x] 직사각형 중간 난도 case `bench_a2_01`에 대해 5-view CLI complete sweep 실행 및 hard-case 편향 완화용 control 확보
- [x] composite 초기 pilot 계열 case `bench_a3_01`에 대해 5-view CLI complete sweep 실행 및 empty composite hallucination 신호 확보
- [x] rectangular multi-obstacle control case `bench_a2_02`에 대해 5-view CLI complete sweep 실행 및 opening-wall 민감도 신호 확보
- [x] rectangular easy-tail case `bench_a1_05`에 대해 5-view CLI complete sweep 실행 및 late rectangular coverage 확보
- [x] rectangular multi-obstacle tail case `bench_a2_04`에 대해 5-view CLI complete sweep 실행 및 repair-salvage 신호 확보
- [x] rectangular dense-tail control case `bench_a2_05`에 대해 5-view CLI complete sweep 실행 및 late rectangular coverage 균형 확보
- [x] `bench_a2_04/perspective` face-touching obstacle pair로 인한 repaired-scene meshing failure를 `repair_indoor_scene.py` clearance hardening으로 복구
- [x] empty composite companion case `bench_a3_02`에 대해 5-view CLI complete sweep 실행 및 empty-composite control 신호 확장
- [x] dense composite control case `bench_a4_01`에 대해 5-view CLI complete sweep 실행 및 A4 layout-fragile 신호 확보
- [x] late composite tail `bench_a3_05`에 대해 5-view CLI complete sweep 실행 및 tail positive-control / hallucination-vs-CFD 신호 확보
- [x] dense composite tail control `bench_a4_04`에 대해 5-view CLI complete sweep 실행 및 high-structure / mid-CFD 분리 신호 확보
- [x] dense composite tail stress-control `bench_a4_05`에 대해 5-view CLI complete sweep 실행 및 section composite-collapse / repair-sidecar 경고 신호 확보
- [x] frozen-20 image-conditioned CLI evaluation 100 tasks 전체를 100/100 success로 완료
- [x] 생성 직후 invalid composite scene도 salvage 가능하도록 `run_indoor_stabilized.py` / `repair_indoor_scene.py` hardening
- [x] 소규모 반복 실행용 `scripts/run_benchmark_evaluation_batch.py` 추가
- [~] API backend 자체는 현재 cron shell에서 여전히 env 미설정(`GEMINI_API_KEY`, `GOOGLE_API_KEY` 모두 없음)이지만, benchmark 진행은 CLI backend로 더 이상 막히지 않음

### D. Failure analysis / stabilization
- [x] 초기 8-scene pilot set batch 실행
- [x] 초기 success/failure 현황 요약
- [x] `bench_a1_02` 실패 원인 분석
- [x] `bench_a3_01` 실패 원인 분석
- [x] `bench_a3_02` 실패 원인 분석
- [x] generator / repair / meshing / solver 중 어디를 고쳐야 하는지 판단
- [x] 수정 후 실패 pilot scene 재실행
- [x] 재실행 후 pilot success-rate summary 갱신
- [x] 20-scene tranche에서도 새로운 stabilization regression이 없는지 점검

### E. Rendering / dataset feasibility
- [x] 성공한 benchmark run에서 3D comparison rendering이 생성됨을 확인
- [x] dataset case용 canonical rendering output structure 정의
- [x] benchmark scene용 자동 2D rendering export 경로 추가
- [x] benchmark view 타입 분리: `perspective`, `birdseye`, `floorplan`, `wireframe`
- [x] 논문 계획과 맞도록 canonical `section` view export 추가
- [x] 하나의 case에서 reference CFD 산출물과 benchmark input image asset이 함께 생성됨을 확인
- [x] 새로 추가된 `*_03` subset에서 integrated batch flow 확인
- [x] frozen-12 rendering manifest를 전체 `section` 포함 상태로 갱신
- [x] `*_05` tranche 이후 aggregate rendering manifest 갱신

### F. Documentation / records
- [x] generator 진행 문서 작성: `docs/26-03-18_benchmark_scene_generator.md`
- [x] 이 autopilot todo 문서 생성
- [x] pilot benchmark status note 작성 및 success/failure table 기록
- [x] major batch rerun 이후 manifest 갱신 유지
- [x] frozen-12 subset 후속 note 작성
- [x] frozen-16 subset 및 manifest-refresh fix 후속 note 작성
- [x] frozen-20 subset 및 scaffold refresh 후속 note 작성
- [x] normalized-grid CFD benchmark metric 관련 note 작성
- [x] frozen-20 reference bundle robustness / failure-signals note 작성
- [x] frozen-20 CLI evaluation 100/100 completion note 작성
- [x] frozen-20 CLI evaluation aggregate summary 생성 (`benchmark/manifests/evaluation_aggregate_summary.json`, `docs/26-03-19_cli_eval_aggregate_results.md`)
- [x] aggregate manifest에 해석용 태그 계층 추가 (opening/topology preserved with hallucination, blockage/layout failure, dense composite structure-vs-physics gap, non-blocking repair-sidecar warning 등)
- [x] 안정적인 체크포인트에서 프로젝트 저장소에 의미 있는 커밋 수행 (push는 사용자 요청 시만)

## 현재 benchmark 스냅샷

복구된 초기 8-scene pilot summary:
- 성공: `bench_a1_01`, `bench_a1_02`, `bench_a2_01`, `bench_a2_02`, `bench_a3_01`, `bench_a3_02`, `bench_a4_01`, `bench_a4_02`
- 성공률: 8/8 reference CFD

현재 frozen-set 상태:
- Frozen scene subset: **20 scenes** (`A1/A2/A3/A4 × 5`)
- Reference CFD 상태: **20/20 성공**
- Benchmark input-view export: `perspective`, `birdseye`, `floorplan`, `wireframe`, `section` 기준 **20/20 성공**
- Evaluation scaffold: **100 tasks** (`20 scenes × 5 views`)
- Image-conditioned evaluation 상태: **100 success**, **0 pending**
  - complete case sweep 완료: `bench_a1_01`, `bench_a1_02`, `bench_a1_03`, `bench_a1_04`, `bench_a1_05`, `bench_a2_01`, `bench_a2_02`, `bench_a2_03`, `bench_a2_04`, `bench_a2_05`, `bench_a3_01`, `bench_a3_02`, `bench_a3_03`, `bench_a3_04`, `bench_a3_05`, `bench_a4_01`, `bench_a4_02`, `bench_a4_03`, `bench_a4_04`, `bench_a4_05` 모두 5-view 전부 성공 (CLI backend)
  - API backend env는 현재 shell에서 여전히 비어 있으나, Gemini CLI cached auth 경로만으로 frozen-20 evaluation 전체를 완료함
- 실제 evaluation task는 geometry metric과 normalized-grid CFD metric (`cfd_metrics.json`)을 함께 기록함

## 주의해서 볼 케이스 / 신호

- `bench_a2_03`, `bench_a4_03`는 mesh 민감도가 상대적으로 큰 대표 케이스였다.
  - frozen-12 recovery 시점에는 `mesh_size=0.35`에서 실패 후 `0.25` fallback이 필요했다.
- `bench_a1_04`는 `0.35`에서 풀렸지만 `ultra_robust` preset이 필요했다.
- `bench_a2_01`은 새 rectangular mid-control로 추가되었고, 5-view 전부 성공했지만 `perspective`만 `mesh_size=0.35` import/checkMesh timeout 후 `0.25` fallback이 필요했다.
- `bench_a4_02`는 reference bundle에서는 `laminar_fallback` 대표 케이스였고, image-conditioned CLI sweep에서는 5-view 전부 성공했지만 `wireframe`이 repaired scene + `ultra_robust`까지 올라가는 강한 stress signal을 보였다.
- `bench_a4_01`은 dense composite control로 새롭게 확보되었고, 5-view 전부 `original + robust`에서 직접 성공했음에도 avg CFD가 `0.4124`로 낮았다.
  - `wireframe`은 `0.6718`로 가장 강했지만 `floorplan`/`perspective`는 각각 `0.2935`, `0.2904`에 머물러 dense composite에서는 floorplan 우세가 쉽게 깨질 수 있음을 보여준다.
  - `perspective`는 최종적으로 성공했지만 stabilization summary 안에는 비차단 repair sidecar failure가 남아 bookkeeping/regression 관점에서 추적 가치가 있다.
- `bench_a3_03`는 composite 중간 난도 positive-control로 새롭게 확보되었다.
  - 5-view 전부 `original + robust`에서 직접 성공했고, `floorplan`/`perspective`는 구조적으로 거의 정답에 가깝게 복원됐다.
  - 다만 `section`은 opening-wall mismatch로 CFD score가 크게 떨어졌고, `wireframe`은 원본 성공과 별개로 repair sidecar 오류 로그가 남아 bookkeeping 관점에서 추적 가치가 있다.
- `bench_a3_01`은 empty composite pilot 계열 case인데, 5-view 모두 composite room kind 자체는 복원했지만 reference obstacle이 0개인데도 모든 view가 obstacle 3개를 hallucinate했다.
  - 그럼에도 `birdseye`/`floorplan`/`section`은 CFD score가 각각 `0.7315` / `0.6965` / `0.5855`로 유지되어, empty composite에서는 obstacle-count mismatch보다 opening/topology fidelity가 더 중요할 수 있다는 신호를 준다.
  - `perspective`만 `ultra_robust`로 한 단계 상승했고, `section`은 outlet wall mismatch가 남았다.
- `bench_a3_02`는 empty composite companion control로, `perspective`는 obstacle 0개 / opening wall 일치 / `structural_score=1.0`까지 복원했지만 `floorplan`은 3-obstacle hallucination과 함께 `ultra_robust`까지 상승했다.
  - 반면 `wireframe`/`section`은 hallucinated obstacle 3개에도 CFD score가 각각 `0.6606`, `0.6248`로 유지되어 empty composite 축에서 hallucination burden과 실제 CFD penalty를 계속 분리해서 읽어야 함을 강화한다.
- `bench_a2_02`는 rectangular multi-obstacle control로, 5-view 모두 `original + robust`에서 직접 성공했지만 floorplan 외 4개 view가 opening wall을 하나 이상 틀렸다.
  - 그 결과 구조 점수는 `0.50~0.75` 범위이고 CFD 평균도 `0.3568`로 낮아, rectangular multi-obstacle case에서는 opening-wall fidelity가 실제 CFD 성능을 꽤 강하게 좌우한다.
- 새 `bench_a1_05`, `bench_a2_05`는 전 view가 `mesh_size=0.35`에서 direct success였고, `bench_a2_04`도 5-view 전부 성공으로 정리됐지만 `perspective`는 repaired-scene salvage가 필요했다.
- `bench_a2_04/perspective`는 original scene의 face-touching obstacle pair 때문에 Gmsh overlapping-facets 오류가 났고, repair clearance hardening 뒤 `repaired + robust`로 회수됐다.
- `bench_a3_05`는 late composite tail positive-control로, 전 view가 composite room kind와 opening wall을 유지했음에도 `birdseye`/`floorplan`/`wireframe`은 obstacle 1 → 3 hallucination을 보였다. 그럼에도 `floorplan` CFD가 `0.6930`까지 올라가 A3 계열에서 opening/topology fidelity가 obstacle exactness보다 더 중요하다는 신호를 강화한다.
- `bench_a4_04`는 dense composite tail control로, 5-view 모두 `original + robust + 0.35`에서 direct success였고 obstacle count `3`도 전부 맞췄다. 다만 `floorplan`의 구조 점수는 `1.0`인데 avg CFD는 `0.4664`에 머물러 dense composite에서는 고득점 structural recovery와 CFD fidelity가 분리될 수 있음을 보여준다.
- `bench_a4_05`는 dense composite tail stress/control로, 5-view 모두 task-level success였지만 `perspective`/`birdseye`에 non-blocking repair sidecar error가 남았고 `section`은 composite room 붕괴(`room_kind_match = false`), opening-wall mismatch, obstacle 4 → 3 축소와 함께 이번 tail 최저 점수(`structural_score = 0.4167`, `cfd_score = 0.3311`)를 기록했다.

## 해석 메모

- 초기 pilot 실패들은 benchmark category 자체의 문제라기보다 meshing/repair compatibility regression에 가까웠다.
- batch reference runner는 이제 reference-CFD 경로와 input-image 경로를 함께 만든다.
- 부분 rerun이 있으면 `reference_batch_summary.json`과 evaluation scaffold가 frozen artifact set을 계속 정확히 반영하는지 갱신이 필요하다.

## 다음 기본 액션

1. 100/100 completion 기준 aggregate summary와 comparison note를 유지하고, 필요하면 논문 Results/Discussion 본문 초안으로 바로 옮길 수 있게 케이스별 해석 문장을 더 압축한다.
2. [완료] composite case opening/topology-sensitive signal, rectangular multi-obstacle blockage/layout failure, empty composite hallucinated-obstacle burden, dense composite structure-vs-physics gap, non-blocking repair-sidecar warning을 aggregate manifest/summary 태그로 반영했다.
3. 다음 자동 액션은 이 태그/집계 결과를 논문용 Results/Discussion 초안 문장 또는 dataset card/README 인용 문장으로 더 압축하는 것이다.
4. API backend env(`GEMINI_API_KEY`/`GOOGLE_API_KEY`)는 별도 유지 과제로 남기되, benchmark core 결과는 CLI backend 기준으로 이미 완결되었으므로 이후 작업은 분석/문서화 중심으로 전환한다.
5. 현재 100/100 completion 상태의 로컬 commit checkpoint를 깔끔하게 유지한다.
6. 이후 meshing/solver 변경이 생기면 stress subset (`bench_a2_01`, `bench_a4_02`, `bench_a2_03`, `bench_a4_03`, `bench_a1_04`)에 더해 easy positive-controls (`bench_a1_02`, `bench_a1_03`, `bench_a1_05`), empty composite controls (`bench_a3_01`, `bench_a3_02`, `bench_a3_05`), dense-composite control (`bench_a4_01`, `bench_a4_04`), dense-composite tail stress/control (`bench_a4_05`), repaired-salvage control (`bench_a2_04`), rectangular opening-control (`bench_a2_02`)도 참고 신호로 본다.
