# Scale-Hint Variant Smoke Comparison (layout-protected span-only redesign)

> 작성일: 2026-03-19  
> 범위: 4-case smoke subset (`bench_a1_01`, `bench_a2_03`, `bench_a3_03`, `bench_a4_03`) × 5 views = 20 tasks  
> 비교 대상:
> - baseline: `benchmark/evaluations` (`no_scale_hint_baseline`)
> - single-anchor: `benchmark/evaluations_scale_hinted` (`scale_hinted_longest_horizontal_span_v1`)
> - dual-anchor: `benchmark/evaluations_scale_hinted_span_height` (`scale_hinted_longest_span_plus_height_v1`)
> - layout-protected: `benchmark/evaluations_scale_hinted_layout_protected` (`scale_hinted_longest_span_layout_protected_v1`)

## 1. 목적

single-anchor / dual-anchor smoke에서 공통적으로 보였던 문제는 다음이었다.

- `Lx` absolute scale ambiguity는 줄어든다.
- 하지만 CFD score는 baseline 대비 하락한다.
- 특히 composite / opening-wall fidelity 보호가 충분하지 않다.

따라서 span-only는 유지하되, prompt wording을 **layout/topology 보호 우선**으로 다시 설계한
`layout-protected` variant를 smoke 수준에서 재검증했다.

## 2. completion status

- baseline smoke counterpart: **20/20 success**
- single-anchor smoke: **20/20 success**
- dual-anchor smoke: **20/20 success**
- layout-protected smoke: **20/20 success**

즉, 새 wording redesign도 smoke task completion 자체는 유지한다.

## 3. baseline 대비 mean delta

### 3.1 layout-protected - baseline
- structural score: **+0.0241**
- CFD score: **-0.0164**
- room bbox relative error
  - `Lx`: **-0.2733**
  - `Ly`: **+0.0454**
  - `Lz`: **+0.0591**
- room volume relative error: **-0.0183**
- room-kind match count delta: **0**
- opening-wall match count delta: **-1**

## 4. layout-protected vs earlier variants

### 4.1 layout-protected - single-anchor
- structural score: **-0.0058**
- CFD score: **+0.0113**
- room bbox relative error
  - `Lx`: **-0.0149**
  - `Ly`: **+0.0026**
  - `Lz`: **+0.0054**
- room volume relative error: **-0.0153**
- room-kind match count delta: **0**
- opening-wall match count delta: **0**

### 4.2 layout-protected - dual-anchor
- structural score: **-0.0208**
- CFD score: **+0.0216**
- room bbox relative error
  - `Lx`: **+0.0000**
  - `Ly`: **+0.0078**
  - `Lz`: **+0.1295**
- room volume relative error: **+0.0700**
- room-kind match count delta: **0**
- opening-wall match count delta: **+1**

## 5. 해석

### 긍정 신호
1. **layout-protected는 baseline 대비 `Lx` scale mismatch를 크게 줄인다.**
   - dual-anchor와 거의 같은 수준으로 longest-span anchor는 유지된다.
2. **CFD score degradation이 earlier variants보다 작다.**
   - baseline 대비 하락폭이 single-anchor(-0.0276), dual-anchor(-0.0380)보다 작다.
3. **opening-wall fidelity가 dual-anchor 대비 1 task 회복된다.**
   - dual-anchor에서 과하게 regularize되던 부분이 약간 완화된 것으로 볼 수 있다.

### 남아 있는 문제
1. **그래도 baseline 대비 CFD mean이 아직 음수다.**
   - 즉, prompt wording redesign만으로는 physics-side degradation을 완전히 제거하지 못했다.
2. **`Ly` / `Lz` regularity는 여전히 mixed다.**
   - longest-span anchor는 잘 걸리지만, non-anchored dimension까지 자동으로 좋아지지는 않는다.
3. **structural gain도 single-anchor 대비 오히려 약간 후퇴했다.**
   - layout 보존을 강화하면서 scale regularization 효과 일부를 trade-off한 셈이다.

## 6. 운영 결론

현재 smoke evidence 기준 운영 판단은 다음과 같다.

1. **layout-protected는 현재까지 본 3개 scale-hint variant 중 가장 균형 잡힌 후보**다.
   - `Lx` 개선 유지
   - CFD penalty 축소
   - opening-wall 손상 완화
2. 하지만 **baseline 대비 CFD degradation이 완전히 사라지지 않았으므로, frozen-20 full rerun main setting으로 즉시 승격하지는 않는다.**
3. 따라서 현재 full rerun scope는 여전히 다음처럼 유지한다.

> **결정: frozen-20 full rerun 보류 유지.**  
> layout-protected는 next-best redesign candidate로 채택하되, full evaluation에 올리기 전 한 번 더 small redesign 또는 targeted robustness check가 필요하다.

## 7. immediate implication

- baseline 100/100 artifact는 계속 `no-scale-hint baseline`으로 유지한다.
- scale-hinted 경로는 이제 다음 3개 variant evidence를 확보했다.
  - single-anchor
  - dual-anchor
  - layout-protected span-only
- 현 시점 figure 재개는 여전히 이르다.
- 다음 단계는 full rerun이 아니라, **layout-protected를 기준 후보로 삼아 추가 미세조정 또는 targeted follow-up smoke**를 설계하는 것이다.
