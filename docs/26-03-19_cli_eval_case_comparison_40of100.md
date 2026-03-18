# CLI evaluation comparison snapshot — 40/100 milestone

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 중 complete 5-view sweep이 끝난 8개 case 비교

## 포함 case

- `bench_a1_01` — rectangular baseline
- `bench_a1_04` — reference solver-stress rectangular case
- `bench_a2_01` — rectangular mid-difficulty case
- `bench_a2_03` — dense-obstacle rectangular stress case
- `bench_a3_03` — mid-difficulty composite case
- `bench_a3_04` — composite representative case
- `bench_a4_02` — laminar-fallback composite stress case
- `bench_a4_03` — obstacle-dense composite stress case

현재 milestone:
- evaluation success **40 / 100**
- pending **60 / 100**

## case-level 평균 비교

| case | 성격 | avg structural | avg CFD | 핵심 신호 |
|---|---|---:|---:|---|
| `bench_a1_01` | rectangular baseline | 0.825 | 0.538 | 구조적으로 가장 안정적인 baseline |
| `bench_a1_04` | rectangular / reference solver-stress | 0.775 | 0.428 | 예측 경로에서는 solver stress가 일부 완화 |
| `bench_a2_01` | rectangular / mid-difficulty | 0.805 | 0.474 | 5/5 direct success, 다만 obstacle 과분할과 perspective mesh 민감도 존재 |
| `bench_a2_03` | rectangular / dense-obstacle | 0.650 | 0.420 | obstacle layout fidelity가 주요 병목 |
| `bench_a3_03` | composite / mid-difficulty | 0.775 | 0.589 | composite positive-control로 활용 가능 |
| `bench_a3_04` | composite representative | 0.383 | 0.453 | structure-vs-CFD 괴리가 가장 큼 |
| `bench_a4_02` | composite / laminar-fallback | 0.738 | 0.334 | repair/solver escalation 신호가 직접 남음 |
| `bench_a4_03` | composite / obstacle-dense | 0.795 | 0.525 | composite hard case 중 가장 건강한 복원 |

## view-level 평균 비교

8개 case 전체에서 view별 평균은 다음과 같다.

| view | avg structural | avg CFD | 메모 |
|---|---:|---:|---|
| perspective | 0.713 | 0.469 | `bench_a2_01` 추가 후에도 depth-view 우위는 제한적 |
| birdseye | 0.691 | 0.487 | 구조는 좋아도 obstacle 과분할/누락이 계속 섞임 |
| floorplan | 0.770 | 0.559 | 여전히 가장 높은 평균 CFD score |
| wireframe | 0.689 | 0.444 | 구조 엣지는 보이지만 유동 복원은 중간 이하 |
| section | 0.726 | 0.392 | 구조는 높지만 CFD는 가장 낮은 패턴 유지 |

## 현재까지 드러난 패턴

### 1. hard-case 편향이 조금 더 줄어들었다

- 35/100 시점에는 rectangular 쪽 evidence가 `bench_a1_01` baseline과 `bench_a1_04`, `bench_a2_03` stress case에 치우쳐 있었다.
- `bench_a2_01`이 추가되면서 rectangular 쪽에도 **중간 난도 control**이 생겼다.
- 이제 rectangular/composite 모두 baseline–mid–stress 축을 조금 더 균형 있게 해석할 수 있다.

### 2. `floorplan` 우세는 여전히 매우 강하다

- 40/100 milestone에서도 `floorplan`이 가장 높은 avg CFD score (`0.559`)를 유지했다.
- `bench_a2_01`에서도 best CFD가 `floorplan`이었다.
- 현재 파이프라인은 photoreal depth cue보다 **평면 배치와 개구부 방향 힌트**를 더 잘 활용하는 것으로 보인다.

### 3. obstacle instance metric과 실제 CFD 영향의 간극이 반복된다

- `bench_a2_01`은 모든 view가 reference 2개 대신 obstacle 3개를 예측했다.
- 그런데 `floorplan`/`birdseye`는 CFD가 여전히 중간 이상으로 유지됐다.
- `bench_a2_03`, `bench_a3_03`, `bench_a3_04`에서 보이던 신호와 합치면, 현재 structural score는 obstacle IoU/instance fidelity에는 민감하지만 실제 blockage/topology 영향은 충분히 설명하지 못한다.

### 4. `section`은 계속 structure-friendly, CFD-fragile 하다

- 40/100 milestone에서도 `section` 평균 structural은 `0.726`으로 낮지 않다.
- 그러나 avg CFD는 계속 최하위 (`0.392`)다.
- `bench_a2_01/section`도 room height 과대추정 + opening wall mismatch로 이 패턴을 재확인했다.

### 5. 성공한 task 안에도 “경미한 stress” 정보가 숨어 있다

- `bench_a4_02/wireframe`은 repair + solver escalation이 필요했던 명시적 stress case였다.
- 새 `bench_a2_01/perspective`는 실패 없이 끝났지만 `mesh_size=0.35` import/checkMesh timeout 후 `0.25`로 내려가야 했다.
- 따라서 앞으로 regression subset은 단순 실패 기록뿐 아니라 **fallback이 필요했던 성공 사례**도 함께 보는 편이 낫다.

## case별 한 줄 판정

- `bench_a1_01`: baseline으로 계속 유지할 가치가 높다.
- `bench_a1_04`: reference stress case지만 predicted geometry가 stress를 완화하는지 보는 용도로 유용하다.
- `bench_a2_01`: rectangular mid-difficulty control로 중요하며, obstacle 과분할/mesh sensitivity 신호를 동시에 준다.
- `bench_a2_03`: obstacle-centric metric refinement의 핵심 regression case다.
- `bench_a3_03`: composite positive-control로서 매우 중요하다.
- `bench_a3_04`: composite metric gap을 가장 잘 드러내는 대표 case다.
- `bench_a4_02`: repair + solver escalation lane의 핵심 regression case다.
- `bench_a4_03`: composite 성공 사례로서 positive hard-case control 역할을 한다.

## 다음 권장 액션

1. 다음 batch는 더 쉬운 rectangular case(`bench_a1_02` 또는 `bench_a1_03`)로 보내 baseline–mid 비교축을 더 채운다.
2. metric-refinement 후보를 문서화한다.
   - composite: opening/topology-sensitive
   - rectangular multi-obstacle: occupancy/blockage-sensitive
3. regression subset 정의를 업데이트할 때 explicit failure뿐 아니라 fallback-success(`bench_a2_01/perspective`)도 포함할지 검토한다.
