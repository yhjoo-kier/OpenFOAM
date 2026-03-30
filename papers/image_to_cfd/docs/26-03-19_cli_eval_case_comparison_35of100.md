# CLI evaluation comparison snapshot — 35/100 milestone

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 중 complete 5-view sweep이 끝난 7개 case 비교

## 포함 case

- `bench_a1_01` — rectangular baseline
- `bench_a1_04` — reference solver-stress rectangular case
- `bench_a2_03` — dense-obstacle rectangular stress case
- `bench_a3_03` — mid-difficulty composite case
- `bench_a3_04` — composite representative case
- `bench_a4_02` — laminar-fallback composite stress case
- `bench_a4_03` — obstacle-dense composite stress case

현재 milestone:
- evaluation success **35 / 100**
- pending **65 / 100**

## case-level 평균 비교

| case | 성격 | avg structural | avg CFD | 핵심 신호 |
|---|---|---:|---:|---|
| `bench_a1_01` | rectangular baseline | 0.825 | 0.538 | 구조적으로 가장 안정적인 baseline |
| `bench_a1_04` | rectangular / reference solver-stress | 0.775 | 0.428 | 예측 경로에서는 solver stress가 일부 완화 |
| `bench_a2_03` | rectangular / dense-obstacle | 0.650 | 0.420 | obstacle layout fidelity가 주요 병목 |
| `bench_a3_03` | composite / mid-difficulty | 0.775 | 0.589 | composite positive-control로 활용 가능 |
| `bench_a3_04` | composite representative | 0.383 | 0.453 | structure-vs-CFD 괴리가 가장 큼 |
| `bench_a4_02` | composite / laminar-fallback | 0.738 | 0.334 | repair/solver escalation 신호가 직접 남음 |
| `bench_a4_03` | composite / obstacle-dense | 0.795 | 0.525 | composite hard case 중 가장 건강한 복원 |

## view-level 평균 비교

7개 case 전체에서 view별 평균은 다음과 같다.

| view | avg structural | avg CFD | 메모 |
|---|---:|---:|---|
| perspective | 0.708 | 0.471 | baseline/중간 난도 composite에서는 강함 |
| birdseye | 0.655 | 0.484 | obstacle 누락이 있어도 CFD는 의외로 유지되는 경우 존재 |
| floorplan | 0.759 | 0.554 | 여전히 가장 높은 평균 CFD score |
| wireframe | 0.681 | 0.452 | 중간 수준이나 repair-sidecar 잡음이 가끔 섞임 |
| section | 0.726 | 0.388 | 구조는 높지만 CFD는 가장 낮음 |

## 현재까지 드러난 패턴

### 1. composite도 hard case만 있는 것은 아니다

- `bench_a3_03`가 5/5 성공으로 추가되면서 composite 평가 축에 **중간 난도 positive-control**이 생겼다.
- `bench_a3_03`는 `bench_a3_04`보다 훨씬 건강했고, `bench_a4_02`보다도 평균 CFD가 크게 높았다.
- 따라서 composite category를 단순히 "어렵다"로 묶기보다, **복원 가능한 composite vs 구조-유동 괴리가 큰 composite**로 나눠 보는 편이 맞다.

### 2. `floorplan` 우세는 여전히 유지된다

- 35/100 milestone에서도 `floorplan`이 가장 높은 avg CFD score (`0.554`)를 유지했다.
- 이는 rectangular/composite를 가리지 않고 반복되고 있다.
- 현재 파이프라인에서는 시각적 realism보다 **배치/개구부/유동 경로에 대한 평면적 힌트**가 더 중요한 것으로 보인다.

### 3. `section`은 계속 structure-friendly, CFD-fragile 하다

- avg structural은 여전히 높지만 (`0.726`), avg CFD는 최하위 (`0.388`)다.
- `bench_a3_03`에서도 room kind는 맞췄지만 opening-wall mismatch 때문에 CFD가 크게 나빠졌다.
- 따라서 section은 “구조를 그럴듯하게 보이게” 만들 수는 있어도 실제 flow-path 복원에는 취약할 가능성이 높다.

### 4. obstacle metric과 CFD impact 사이의 간극이 계속 보인다

- `bench_a2_03`, `bench_a3_04`, `bench_a3_03/birdseye` 모두 obstacle fidelity가 약한데 CFD는 완전히 무너지지 않았다.
- 이는 현재 구조 점수가 obstacle IoU에 상대적으로 예민한 반면, 실제 유동은 더 coarse한 blockage/topology 정보에 좌우된다는 해석과 맞는다.
- dense-obstacle 케이스에는 **occupancy/blockage-sensitive metric**이 더 직접적일 수 있다.

### 5. repair bookkeeping은 failure 해석에서 분리해야 한다

- `bench_a4_02/wireframe`은 실제 repair + solver escalation stress case였다.
- 반면 `bench_a3_03/wireframe`은 원본이 바로 성공했는데도 repair sidecar 오류 로그가 남았다.
- 그래서 앞으로는 "repair attempted" 로그 자체와 실제 solver salvage 필요 여부를 구분해서 읽는 편이 안전하다.

## case별 한 줄 판정

- `bench_a1_01`: baseline으로 계속 유지할 가치가 높다.
- `bench_a1_04`: reference stress case지만 predicted geometry가 stress를 완화하는지 보는 용도로 유용하다.
- `bench_a2_03`: obstacle-centric metric refinement의 핵심 regression case다.
- `bench_a3_03`: composite positive-control로서 매우 중요하다.
- `bench_a3_04`: composite metric gap을 가장 잘 드러내는 대표 case다.
- `bench_a4_02`: repair + solver escalation lane의 핵심 regression case다.
- `bench_a4_03`: composite 성공 사례로서 positive hard-case control 역할을 한다.

## 다음 권장 액션

1. 다음 batch는 직사각형 중간 난도 case(`bench_a2_01` 또는 `bench_a1_02`)로 보내 hard-case 편향을 더 낮춘다.
2. metric-refinement 후보를 문서화한다.
   - composite: opening/topology-sensitive
   - dense-obstacle: occupancy/blockage-sensitive
3. 필요하면 이후 generator/repair/solver 변경 시 우선 재평가할 regression subset을 7개 case 기준으로 갱신한다.
