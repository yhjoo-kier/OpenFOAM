# Scale-Hint Variant Smoke Comparison (single-anchor vs dual-anchor)

> 작성일: 2026-03-19
> 범위: 4-case smoke subset (`bench_a1_01`, `bench_a2_03`, `bench_a3_03`, `bench_a4_03`) × 5 views = 20 tasks
> 비교 대상:
> - baseline: `benchmark/evaluations` (`no_scale_hint_baseline`)
> - single-anchor: `benchmark/evaluations_scale_hinted` (`scale_hinted_longest_horizontal_span_v1`)
> - dual-anchor: `benchmark/evaluations_scale_hinted_span_height` (`scale_hinted_longest_span_plus_height_v1`)

## 1. 목적

single-anchor smoke 결과에서 `Lx` scale mismatch는 완화되었지만 CFD score는 소폭 악화되었기 때문에,
`ceiling height + longest horizontal span` 보강안이 CFD degradation을 줄이는지 확인했다.

## 2. completion status

- single-anchor smoke: **20/20 success**
- dual-anchor smoke: **20/20 success**
- baseline smoke counterpart: **20/20 success**

즉, 두 scale-hinted variant 모두 task completion 자체는 baseline과 동등하게 유지되었다.

## 3. baseline 대비 mean delta

### 3.1 single-anchor (`longest horizontal span only`) - baseline
- structural score: **+0.0299**
- CFD score: **-0.0276**
- room bbox relative error
  - `Lx`: **-0.2583**
  - `Ly`: **+0.0428**
  - `Lz`: **+0.0537**
- room volume relative error: **-0.0030**
- room-kind match count delta: **0**
- opening-wall match count delta: **-1**

### 3.2 dual-anchor (`longest span + ceiling height`) - baseline
- structural score: **+0.0449**
- CFD score: **-0.0380**
- room bbox relative error
  - `Lx`: **-0.2733**
  - `Ly`: **+0.0375**
  - `Lz`: **-0.0704**
- room volume relative error: **-0.0884**
- room-kind match count delta: **0**
- opening-wall match count delta: **-2**

## 4. dual-anchor vs single-anchor delta

### dual-anchor - single-anchor
- structural score: **+0.0150**
- CFD score: **-0.0104**
- room bbox relative error
  - `Lx`: **-0.0149**
  - `Ly`: **-0.0052**
  - `Lz`: **-0.1241**
- room volume relative error: **-0.0854**
- room-kind match count delta: **0**
- opening-wall match count delta: **-1**

## 5. 해석

### 긍정 신호
1. dual-anchor는 single-anchor보다 **`Lz` 오차를 실제로 줄인다**.
2. dual-anchor는 single-anchor보다 **structural score도 소폭 높다**.
3. 두 variant 모두 **`Lx` absolute scale ambiguity 완화**에는 유효하다.

### 부정 신호
1. dual-anchor는 **CFD score를 회복하지 못했다**.
   - baseline 대비 하락폭이 single-anchor보다 더 크다.
2. dual-anchor는 **opening-wall fidelity도 더 나빠졌다**.
   - smoke 20-task 기준 single-anchor 대비 1건, baseline 대비 2건 더 줄었다.
3. dual-anchor는 **volume regularization이 더 강해졌지만**, 그 자체가 physics fidelity 개선으로 이어지지 않았다.

## 6. 운영 결론

현재 smoke evidence 기준으로는 다음이 가장 타당하다.

1. **dual-anchor를 frozen-20 full rerun의 main setting으로 채택하지 않는다.**
2. **single-anchor도 아직 main-setting full rerun으로 바로 승격하지 않는다.**
3. 따라서 **지금은 full frozen-20 scale-hinted rerun을 시작하지 않는다.**
4. 다음 단계는 full rerun이 아니라, 아래 둘 중 하나로 scale hint 설계를 한 번 더 다듬는 것이다.
   - hint wording 재설계 (opening/layout fidelity 보호 문구 강화)
   - span-only / span+height 외의 lighter regularization variant 검토

즉, 현 시점 full rerun scope 결정은 다음처럼 정리된다:

> **결정: frozen-20 full rerun 보류. 현재 20-task smoke evidence만으로는 scale-hinted setting을 main evaluation으로 승격하기 부족하다.**

## 7. 논문/figure 해석에 대한 즉시 영향

- 기존 baseline 100/100 결과는 계속 **no-scale-hint baseline artifact** 로 유지한다.
- scale-hinted setting은 방법론적으로 유효한 pivot이지만, **현재까지의 prompt/hint formulation은 geometry-scale 개선 대비 CFD 일관 개선이 부족**하다.
- 따라서 figure 재개 여부는 아직 이르며, scale-hinted main-setting이 실제로 정리된 뒤에만 다시 판단한다.
