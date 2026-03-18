# CLI evaluation comparison snapshot — 30/100 milestone

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 중 complete 5-view sweep이 끝난 6개 case 비교

## 포함 case

- `bench_a1_01` — rectangular baseline
- `bench_a1_04` — reference solver-stress rectangular case
- `bench_a2_03` — dense-obstacle rectangular stress case
- `bench_a3_04` — composite representative case
- `bench_a4_02` — laminar-fallback composite stress case
- `bench_a4_03` — obstacle-dense composite stress case

현재 milestone:
- evaluation success **30 / 100**
- pending **70 / 100**

## case-level 평균 비교

| case | 성격 | avg structural | avg CFD | 핵심 신호 |
|---|---|---:|---:|---|
| `bench_a1_01` | rectangular baseline | 0.825 | 0.538 | 구조적으로 가장 안정적인 baseline |
| `bench_a1_04` | rectangular / reference solver-stress | 0.775 | 0.428 | 예측 경로에서는 solver stress가 일부 완화 |
| `bench_a2_03` | rectangular / dense-obstacle | 0.650 | 0.420 | obstacle layout fidelity가 주요 병목 |
| `bench_a3_04` | composite representative | 0.383 | 0.453 | structure-vs-CFD 괴리가 가장 큼 |
| `bench_a4_02` | composite / laminar-fallback | 0.738 | 0.334 | repair/solver escalation 신호가 직접 남음 |
| `bench_a4_03` | composite / obstacle-dense | 0.795 | 0.525 | composite hard case 중 가장 건강한 복원 |

## view-level 평균 비교

6개 case 전체에서 view별 평균은 다음과 같다.

| view | avg structural | avg CFD | 메모 |
|---|---:|---:|---|
| perspective | 0.660 | 0.444 | baseline에서는 강하지만 hard case에서 편차 큼 |
| birdseye | 0.660 | 0.448 | perspective와 비슷한 평균 |
| floorplan | 0.719 | 0.525 | 현재까지 가장 높은 평균 CFD score |
| wireframe | 0.690 | 0.433 | 중간 수준 |
| section | 0.743 | 0.398 | 구조는 높지만 CFD는 가장 낮음 |

## 현재까지 드러난 패턴

### 1. `floorplan`은 생각보다 강하다

- 6개 case 평균에서 `floorplan`이 **가장 높은 avg CFD score (0.525)** 를 보였다.
- 특히 `bench_a2_03`, `bench_a3_04`에서 `floorplan`이 다른 view보다 훨씬 좋은 CFD score를 냈다.
- 반면 구조적으로는 `section`이 더 높을 때도 있었기 때문에, **평면 정보가 opening/flow-path 추정에는 유리하지만 완전한 3D 복원과는 다른 축**일 가능성이 있다.

### 2. `section`은 structure-friendly, CFD-fragile 경향이 있다

- 전체 평균 structural score는 `section`이 가장 높다(0.743).
- 하지만 avg CFD score는 오히려 가장 낮다(0.398).
- `bench_a1_01`과 `bench_a4_03`에서도 비슷한 경향이 보였고, `bench_a3_04/section`은 composite room을 rectangular로 collapse하며 가장 취약했다.

### 3. hard case는 최소 세 부류로 나뉜다

1. **baseline-like**: `bench_a1_01`
   - room/opening/obstacle coarse structure가 안정적
2. **rectangular hard**
   - `bench_a2_03`: obstacle-dense inference hard
   - `bench_a1_04`: reference solver-hard였지만 prediction-side에서는 geometry simplification hard
3. **composite hard**
   - `bench_a3_04`: structure-vs-CFD 괴리가 큼
   - `bench_a4_02`: repair/solver-hard
   - `bench_a4_03`: composite지만 비교적 잘 복원됨

### 4. coarse structural score만으로는 설명이 부족하다

- `bench_a3_04`는 avg structural이 매우 낮음(0.383)에도 avg CFD가 `bench_a1_04`, `bench_a2_03`보다 약간 높다.
- `bench_a2_03/floorplan`은 obstacle IoU가 사실상 무너졌는데도 `cfd_score ≈ 0.655`가 나왔다.
- 따라서 현재 score 체계는 적어도 두 부분에서 추가 보강이 필요하다.
  - composite case용 **opening/topology-sensitive metric**
  - dense-obstacle case용 **occupancy/blockage-sensitive metric**

## case별 한 줄 판정

- `bench_a1_01`: baseline으로 계속 유지할 가치가 높다.
- `bench_a1_04`: reference stress case지만 predicted geometry가 stress를 완화하는지 보는 용도로 유용하다.
- `bench_a2_03`: obstacle-centric metric refinement의 핵심 regression case다.
- `bench_a3_04`: composite metric gap을 가장 잘 드러내는 대표 case다.
- `bench_a4_02`: repair + solver escalation lane의 핵심 regression case다.
- `bench_a4_03`: composite 성공 사례로서 positive control 역할을 한다.

## 다음 권장 액션

1. aggregate summary를 바탕으로 **metric-refinement 후보를 문서화**한다.
   - composite: opening/topology-sensitive
   - rectangular dense-obstacle: occupancy/blockage-sensitive
2. hard-case 편향을 줄이기 위해 다음 batch는 미평가 중간 난도 case(`bench_a1_02`, `bench_a2_01`, `bench_a3_03` 등) 중 하나를 선택한다.
3. 필요하면 이 6개 case를 우선 regression subset으로 묶어 이후 generator/repair/solver 변경 시 빠르게 재평가한다.
