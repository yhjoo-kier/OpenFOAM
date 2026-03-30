# CLI-backed benchmark evaluation case sweep — `bench_a1_04`

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark 중 solver-stress 성격이 강했던 직사각형 case `bench_a1_04`에 대해 5-view image-conditioned evaluation을 모두 실행

## 목적

`bench_a1_04`는 reference bundle 기준으로 `mesh_size=0.35`에서 **`ultra_robust` preset까지 올라갔던 대표 직사각형 stress case**였다.
이번 sweep의 목적은 이 stress가 image-conditioned 경로에서도 그대로 드러나는지, 아니면 생성된 geometry 단순화 때문에 다른 failure mode로 바뀌는지 확인하는 것이다.

## 실행 경로

- runner: `scripts/run_benchmark_evaluation_batch.py`
- backend: Gemini CLI cached auth
- model: `gemini-3-flash-preview`
- task filter: `bench_a1_04`의 5개 view

실행 예시:

```bash
python3 scripts/run_benchmark_evaluation_batch.py \
  --cases bench_a1_04 \
  --backend cli \
  --model gemini-3-flash-preview \
  --statuses pending \
  --summary-out benchmark/manifests/evaluation_batch_summary_a1_04.json
```

## 결과 요약

`bench_a1_04`의 5개 view가 모두 **end-to-end 성공**했다.

| view | structural_score | cfd_score | solver path | 메모 |
|---|---:|---:|---|---|
| perspective | 0.875 | 0.410 | `robust` | 구조는 가장 좋았지만 CFD는 중간 수준 |
| birdseye | 0.750 | 0.444 | `conservative` | 5-view 중 유일하게 `conservative`까지 상승 |
| floorplan | 0.750 | 0.435 | `robust` | 안정적으로 수렴 |
| wireframe | 0.750 | 0.447 | `robust` | 이번 sweep 최고 CFD score |
| section | 0.750 | 0.407 | `robust` | 구조/CFD 모두 중간 |

모든 view가 **원본 predicted scene 기준으로 성공**했고, repair scene 의존은 없었다.
다만 mesh risk는 전 view에서 계속 높게 나타났고(`maxNonOrtho ≈ 86.6~87.7`, `aspectRatio ≈ 37.5~61.7`),
`birdseye`는 실제로 `laminar_fallback` 다음 `conservative`에서 수렴했다.

## 1차 해석

1. **reference solver-stress가 image-conditioned path에서는 약화되었다.**
   - reference bundle에서는 `bench_a1_04`가 `ultra_robust` 대표 케이스였다.
   - 하지만 이번 predicted scenes는 5개 전부 성공했고, 가장 어려운 `birdseye`도 `conservative`에서 멈췄다.
   - 즉 이 case의 reference-side 어려움은 VLM이 복원한 형상이 더 단순해지면서 일부 완화된 것으로 보인다.

2. **이 case의 주 failure mode는 “solver failure”보다 “geometry simplification with moderate CFD loss”에 가깝다.**
   - 구조 점수는 0.75~0.875로 무난하지만 CFD score는 0.41~0.45 수준에 머문다.
   - 완전히 깨지지는 않지만, 원래 stress를 만들던 세부 형상/개구부/장애물 조건이 예측 과정에서 부드러워졌을 가능성이 높다.

3. **view 간 차이는 크지 않다.**
   - `wireframe`이 최고 CFD score(≈0.447), `perspective`가 최고 structural score(0.875)를 보였지만,
     전체적으로는 한 view가 압도적으로 우세하다고 보기는 어렵다.
   - `bench_a2_03`처럼 특정 view(`floorplan`) 하나만 유난히 강한 패턴과는 다르다.

4. **stress-case 해석 시 reference hardness와 evaluation hardness를 구분해야 한다.**
   - `bench_a1_04`는 reference path에서는 solver-hard였지만,
     image-conditioned path에서는 오히려 “stress를 만드는 세부 geometry를 잃어버리는” 쪽으로 바뀌었다.
   - 반대로 `bench_a4_02`는 evaluation path에서도 repair/solver escalation 신호가 그대로 살아 있었다.

## 산출물

- batch summary:
  - `benchmark/manifests/evaluation_batch_summary_a1_04.json`
- per-view evaluation summaries:
  - `benchmark/evaluations/bench_a1_04/*/evaluation_summary.json`
- predicted scenes / cases / CFD metrics:
  - `benchmark/evaluations/bench_a1_04/*/predicted_scene.json`
  - `benchmark/evaluations/bench_a1_04/*/predicted_case`
  - `benchmark/evaluations/bench_a1_04/*/cfd_metrics.json`

## 현재까지의 비교 포인트

- `bench_a2_03`: obstacle-dense rectangular case → solver는 버티지만 obstacle layout fidelity가 크게 무너짐
- `bench_a1_04`: reference solver-stress case → predicted geometry가 단순화되며 solver stress는 일부 완화
- `bench_a3_04`, `bench_a4_02`, `bench_a4_03`: composite cases → topology / opening / repair / solver stress가 더 직접적으로 드러남

즉, 현재 frozen-20 benchmark의 hard cases는 적어도 다음 세 부류로 나뉜다.

1. **dense-obstacle inference hard** (`bench_a2_03`)
2. **reference-solver-stress but prediction-simplification hard** (`bench_a1_04`)
3. **composite topology + repair hard** (`bench_a3_04`, `bench_a4_02`, `bench_a4_03`)

## 다음 권장 액션

1. 지금까지 complete sweep이 끝난 6개 case를 묶어 **case-level aggregate summary**를 만든다.
2. 구조 metric 보강 후보를 두 갈래로 나눠 검토한다.
   - composite용: opening/topology-sensitive 신호
   - rectangular dense-obstacle용: occupancy / blockage-sensitive 신호
3. 다음 실행 후보는 아직 미평가인 쉬운/중간 case 일부(`bench_a1_02`, `bench_a2_01`, `bench_a3_03` 등)를 소규모로 추가해,
   hard-case 편향 없이 category별 baseline 평균을 만들 수 있게 한다.
