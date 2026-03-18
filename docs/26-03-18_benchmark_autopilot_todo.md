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
- Image-conditioned evaluation 상태: **2 success**, **98 pending**
  - 성공 smoke tests: `bench_a1_01/perspective`, `bench_a1_01/section` (CLI backend)
  - API backend env는 현재 shell에서 여전히 비어 있으나, Gemini CLI cached auth 경로로 평가 진행 가능
- 실제 evaluation task는 geometry metric과 normalized-grid CFD metric (`cfd_metrics.json`)을 함께 기록함

## 주의해서 볼 케이스 / 신호

- `bench_a2_03`, `bench_a4_03`는 mesh 민감도가 상대적으로 큰 대표 케이스였다.
  - frozen-12 recovery 시점에는 `mesh_size=0.35`에서 실패 후 `0.25` fallback이 필요했다.
- `bench_a1_04`는 `0.35`에서 풀렸지만 `ultra_robust` preset이 필요했다.
- `bench_a4_02`는 `laminar_fallback`으로 성공한 대표 케이스다.
- 반면 새 `*_05` scene들은 전부 `mesh_size=0.35 + robust`에서 직접 성공했다.

## 해석 메모

- 초기 pilot 실패들은 benchmark category 자체의 문제라기보다 meshing/repair compatibility regression에 가까웠다.
- batch reference runner는 이제 reference-CFD 경로와 input-image 경로를 함께 만든다.
- 부분 rerun이 있으면 `reference_batch_summary.json`과 evaluation scaffold가 frozen artifact set을 계속 정확히 반영하는지 갱신이 필요하다.

## 다음 기본 액션

1. CLI backend로 frozen-20 image-conditioned evaluation을 소규모 배치(예: case 1개 × 5 views 또는 category별 대표 subset) 실행해 view별 failure/success 신호를 모은다.
2. perspective/section 첫 성공 2건의 `cfd_metrics.json`을 기준으로, 논문용으로 slice-specific metric이나 scalar-field 추가 지표가 필요한지 판단한다.
3. API backend env(`GEMINI_API_KEY`/`GOOGLE_API_KEY`)는 별도 유지 과제로 남기되, benchmark 진행 자체는 CLI backend 기준으로 계속 전진한다.
4. frozen-20 benchmark + CLI-eval 복구 상태에서 깔끔한 로컬 commit checkpoint를 유지한다.
5. 이후 meshing/solver 변경이 생기면 stress subset (`bench_a2_01`, `bench_a4_02`, `bench_a2_03`, `bench_a4_03`, `bench_a1_04`)로 빠른 regression check를 수행한다.
