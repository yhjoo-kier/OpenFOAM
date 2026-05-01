# Scale-Hint Post-View-Guarded Gate Note

> 작성일: 2026-03-20  
> 범위: `layout-protected span-only` → `view-guarded` smoke 이후, frozen-20 full rerun gate를 다시 점검하고 다음 wording redesign 방향을 좁히기 위한 메모

## 1. 한 줄 결론

`view-guarded`는 **원래 겨냥했던 부작용 축(perspective / section / A4 dense composite)** 에서는 실제로 완화 효과를 보였지만, **overall smoke mean CFD를 baseline 이상으로 끌어올리지는 못했다.**

따라서 현재 운영 결론은 그대로다.

> **frozen-20 full rerun은 계속 보류한다.**

다만 이번 smoke로 다음 redesign의 타깃은 더 좁아졌다.

---

## 2. 이미 확인된 것

### 유지되는 강점
- `Lx` anchor 효과는 유지된다.
  - view-guarded overall mean `Lx` relative error: **0.0004**
- perspective / section / A4 dense composite에서 `layout-protected` 대비 개선 신호가 있다.
  - perspective mean CFD: **+0.0639** vs layout-protected
  - section mean CFD: **+0.0334** vs layout-protected
  - A4 mean CFD: **+0.0346** vs layout-protected

### 아직 못 푼 것
- overall mean CFD: **0.4926**
- layout-protected overall mean CFD: **0.5017**
- baseline 대비 view-guarded overall CFD delta: **-0.0254**
- opening-wall match rate: **0.70** (baseline 완전 회복 아님)

즉, 현재 문제는 더 이상 “scale hint가 전반적으로 틀렸나?”가 아니라,

> **어느 view/case에서는 view-guarded가 분명 맞는 방향인데, 다른 축에서 남는 penalty가 전체 평균을 끌어내린다**

로 보는 편이 더 정확하다.

---

## 3. per-task readout: 현재 overall mean을 끌어내리는 잔여 병목

아래는 4-case × 5-view smoke에서 `view-guarded - baseline` 기준 CFD delta가 음수로 남은 대표 축이다.

### A1 / wireframe
- baseline: **0.5997**
- layout-protected: **0.6255**
- view-guarded: **0.4225**
- delta vs baseline: **-0.1772**
- 해석:
  - 원래 잘 되던 쉬운 rectangular wireframe에 대해, 추가 guard가 오히려 불필요한 regularization side-effect를 만든 신호다.
  - 즉, `view-guarded` 문구가 모든 view에서 uniformly 더 좋은 것은 아니다.

### A4 / floorplan
- baseline: **0.5835**
- layout-protected: **0.6024**
- view-guarded: **0.4559**
- delta vs baseline: **-0.1276**
- 해석:
  - dense composite 보정 문구가 floorplan처럼 원래 layout-dominant인 view에도 과하게 개입했을 가능성이 있다.
  - A4에서 필요한 것은 dense-scene obstacle simplification이지, floorplan anchor 자체를 약하게 만드는 것이 아닐 수 있다.

### A3 / perspective
- baseline: **0.6317**
- layout-protected: **0.6495**
- view-guarded: **0.5594**
- delta vs baseline: **-0.0723**
- 해석:
  - perspective guard는 A4 perspective에는 유효했지만, light composite perspective(A3)에서는 오히려 geometry flexibility를 너무 줄였을 수 있다.
  - 즉, perspective guard도 **room kind / density** 에 따라 강도를 달리할 필요가 있다.

### A1 / section
- baseline: **0.4427**
- layout-protected: **0.3395**
- view-guarded: **0.3745**
- delta vs baseline: **-0.0682**
- opening-wall: baseline `True` → layout `False` → view `False`
- 해석:
  - section 전용 guard는 layout-protected 대비 개선되었지만, opening-wall fidelity를 baseline 수준으로 복구하지는 못했다.
  - section은 여전히 별도 최약축이다.

### A4 / birdseye
- baseline: **0.5280**
- layout-protected: **0.5432**
- view-guarded: **0.5012**
- delta vs baseline: **-0.0269**
- 해석:
  - A4 dense-scene wording이 birdseye에도 약간 과개입했을 가능성이 있다.

---

## 4. 반대로, view-guarded가 실제로 맞았던 곳

다음 축에서는 `view-guarded`가 “다음 방향이 맞다”는 근거를 준다.

### A4 / perspective
- baseline: **0.5287**
- layout-protected: **0.4541**
- view-guarded: **0.5626**
- delta vs baseline: **+0.0340**

### A4 / section
- baseline: **0.5014**
- layout-protected: **0.3801**
- view-guarded: **0.4942**
- delta vs baseline: **-0.0072**
- 사실상 baseline 근처까지 회복

### A4 / wireframe
- baseline: **0.4858**
- layout-protected: **0.3829**
- view-guarded: **0.5216**
- delta vs baseline: **+0.0359**

### A1 / perspective
- baseline: **0.6470**
- layout-protected: **0.3851**
- view-guarded: **0.6216**
- delta vs baseline: **-0.0254**
- 완전 회복은 아니지만, layout-protected의 perspective failure는 대부분 완화됨

즉, 이번 refinement의 가치는 분명하다.

> **문제는 guard 자체가 틀린 게 아니라, guard의 적용 범위와 강도가 아직 너무 균일하다는 점이다.**

---

## 5. 다음 redesign 가설

현재 가장 타당한 다음 가설은 아래와 같다.

### 가설: next variant는 "stronger hint"가 아니라 "view-/density-conditional guard weighting" 이어야 한다.

즉,
- `longest horizontal span only`
- `layout-protected`
- `opening/topology first`

의 큰 틀은 유지하되,

### (A) floorplan / birdseye
- 이미 layout signal이 강한 view에서는 추가 guard를 **최소화**한다.
- 특히 A4 dense composite에서도 floorplan은 anchor 자체를 약하게 만들기보다,
  - opening wall identity,
  - connected-room topology,
  - obstacle clearance
  만 남기고 과한 caution 문구는 줄이는 편이 더 나을 가능성이 있다.

### (B) perspective
- hidden depth guard는 유지하되,
- composite/light-composite에서는 과도하게 box-regularized 되지 않도록 완화한다.
- 즉 perspective guard는 `dense composite`와 `non-dense composite/rectangular`를 분기하는 편이 타당하다.

### (C) section
- section은 여전히 opening-wall fidelity 최우선이 맞다.
- 다만 현재 문구만으로는 baseline opening-wall을 회복하지 못했다.
- 따라서 section은 다음 variant에서도 별도 stress lane으로 유지하고,
  - opening wall identity before scale
  - do not relocate openings across walls
  를 더 직접적으로 못 박되,
  - unseen height suppression은 지금보다 더 간결하게 줄이는 편이 나을 수 있다.

### (D) wireframe
- wireframe은 원래 layout-protected에서 크게 나쁘지 않았다.
- 따라서 `view-guarded`의 추가 문구를 wireframe에까지 적극 주입하는 것은 오히려 과개입일 수 있다.
- 다음 variant에서는 wireframe을 layout-protected에 더 가깝게 되돌리는 것이 합리적이다.

---

## 6. 운영 제안

full rerun 전, 다음 실험은 여전히 **작은 smoke subset** 이 적절하다.

### 권장 working direction
- 임시 이름: `scale_hinted_longest_span_guard_weighted_v1`

### 핵심 원칙
1. stronger hint 추가 금지
2. `layout-protected span-only` 유지
3. `view-guarded`의 좋은 부분은 남기되,
   - floorplan / birdseye / wireframe에는 과한 guard를 줄이고
   - perspective / section에는 조건부로만 강화
   - dense composite와 non-dense case를 더 분리

### 최소 smoke 체크 대상
- `bench_a1_01 / wireframe`
- `bench_a1_01 / section`
- `bench_a3_03 / perspective`
- `bench_a4_03 / floorplan`
- `bench_a4_03 / perspective`
- `bench_a4_03 / section`

즉, 다음 smoke는 단순 4-case full 20-task 재반복보다,

> **현재 residual penalty를 만드는 6-task stress lane**

에 먼저 집중하는 편이 더 효율적이다.

---

## 7. 최종 메모

지금까지의 evidence를 실무적으로 요약하면 다음과 같다.

1. scale-hint pivot은 유지한다.
2. frozen-20 full rerun은 아직 아니다.
3. `view-guarded`는 실패가 아니라, **다음 redesign을 더 좁혀 준 중간 성공**이다.
4. 다음 한 수는 stronger hint가 아니라, **guard weighting / selective application 정리** 쪽이 더 타당하다.
