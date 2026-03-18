# CLI evaluation comparison snapshot — 45/100 milestone

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 중 complete 5-view sweep이 끝난 9개 case 비교

## 포함 case

- `bench_a1_01` — rectangular baseline
- `bench_a1_02` — rectangular easy positive-control
- `bench_a1_04` — reference solver-stress rectangular case
- `bench_a2_01` — rectangular mid-difficulty case
- `bench_a2_03` — dense-obstacle rectangular stress case
- `bench_a3_03` — mid-difficulty composite case
- `bench_a3_04` — composite representative case
- `bench_a4_02` — laminar-fallback composite stress case
- `bench_a4_03` — obstacle-dense composite stress case

현재 milestone:
- evaluation success **45 / 100**
- pending **55 / 100**

## case-level 평균 비교

| case | 성격 | avg structural | avg CFD | 핵심 신호 |
|---|---|---:|---:|---|
| `bench_a1_01` | rectangular baseline | 0.825 | 0.538 | 구조적으로 가장 안정적인 baseline |
| `bench_a1_02` | rectangular / easy positive-control | 0.825 | 0.594 | 가장 깨끗한 easy-case CFD 신호, 다만 일부 view는 obstacle hallucination |
| `bench_a1_04` | rectangular / reference solver-stress | 0.775 | 0.428 | 예측 경로에서는 solver stress가 일부 완화 |
| `bench_a2_01` | rectangular / mid-difficulty | 0.805 | 0.474 | 5/5 direct success, 다만 obstacle 과분할과 perspective mesh 민감도 존재 |
| `bench_a2_03` | rectangular / dense-obstacle | 0.650 | 0.420 | obstacle layout fidelity가 주요 병목 |
| `bench_a3_03` | composite / mid-difficulty | 0.775 | 0.589 | composite positive-control로 활용 가능 |
| `bench_a3_04` | composite representative | 0.383 | 0.453 | structure-vs-CFD 괴리가 가장 큼 |
| `bench_a4_02` | composite / laminar-fallback | 0.738 | 0.334 | repair/solver escalation 신호가 직접 남음 |
| `bench_a4_03` | composite / obstacle-dense | 0.795 | 0.525 | composite hard case 중 가장 건강한 복원 |

## view-level 평균 비교

9개 case 전체에서 view별 평균은 다음과 같다.

| view | avg structural | avg CFD | 메모 |
|---|---:|---:|---|
| perspective | 0.718 | 0.486 | easy case에서도 obstacle hallucination이 반복됨 |
| birdseye | 0.726 | 0.509 | easy rectangular에서는 strongest positive-control signal |
| floorplan | 0.796 | 0.573 | 여전히 가장 높은 평균 CFD score |
| wireframe | 0.696 | 0.463 | edge cue는 주지만 object hallucination/mesh warning이 섞임 |
| section | 0.715 | 0.388 | 구조는 버티지만 CFD는 계속 최하위 |

## 현재까지 드러난 패턴

### 1. rectangular evidence가 이제 baseline–easy–mid–stress 축을 갖추기 시작했다

- `bench_a1_02`가 추가되면서 rectangular 쪽은
  - baseline (`bench_a1_01`)
  - easy positive-control (`bench_a1_02`)
  - mid-difficulty (`bench_a2_01`)
  - stress (`bench_a1_04`, `bench_a2_03`)
  로 해석 축이 더 분리됐다.
- 이제 rectangular category를 단순히 쉬움/어려움으로만 보지 않고, **solver stress vs perception ambiguity**로 나눠 읽을 수 있다.

### 2. `floorplan` 우세는 여전히 매우 강하고, easy case에서는 `birdseye`도 거의 동급이다

- `floorplan`이 전체 최고 avg CFD (`0.573`)를 유지했다.
- 새 `bench_a1_02`에서는 `birdseye`와 `floorplan`이 거의 동일한 최고 성능을 냈다.
- 즉, layout-preserving view가 가장 강한 흐름은 유지되지만, 쉬운 빈방 계열에서는 bird's-eye cue도 매우 경쟁력 있다.

### 3. easy rectangular에서도 obstacle hallucination은 사라지지 않는다

- `bench_a1_02` reference는 obstacle 0개인데, `perspective`/`wireframe`/`section`이 모두 3개 obstacle을 예측했다.
- 따라서 현재 image-conditioned path의 structural loss 중 일부는 benchmark 난도보다 **view-specific interpretation bias**에서 오는 것으로 봐야 한다.

### 4. `section`의 weakness는 이제 거의 일관된 법칙으로 보인다

- 45/100 시점에서도 `section` avg CFD는 최하위 (`0.388`)다.
- `bench_a1_02/section`처럼 room 자체는 쉬워도 opening wall mismatch와 partial-view ambiguity가 남는다.
- 그래서 section은 계속해서 “형상은 그럴듯, 유동은 취약” 패턴을 재생산하고 있다.

### 5. mesh-risk warning과 실제 실패를 분리해서 읽는 것이 더 중요해졌다

- `bench_a1_02` 5개 task는 모두 성공했지만 전부 `high` mesh-risk를 기록했다.
- 즉, 현재 risk flag는 failure predictor라기보다 conservative warning에 가깝다.
- 앞으로는 risk-level 자체보다 실제 fallback, repair, preset escalation 여부를 함께 봐야 regression 신호가 선명하다.

## case별 한 줄 판정

- `bench_a1_01`: baseline으로 계속 유지할 가치가 높다.
- `bench_a1_02`: easy rectangular positive-control로 중요하며, view-driven obstacle hallucination을 분리해서 보여준다.
- `bench_a1_04`: reference stress case지만 predicted geometry가 stress를 완화하는지 보는 용도로 유용하다.
- `bench_a2_01`: rectangular mid-difficulty control로 중요하며, obstacle 과분할/mesh sensitivity 신호를 동시에 준다.
- `bench_a2_03`: obstacle-centric metric refinement의 핵심 regression case다.
- `bench_a3_03`: composite positive-control로서 매우 중요하다.
- `bench_a3_04`: composite metric gap을 가장 잘 드러내는 대표 case다.
- `bench_a4_02`: repair + solver escalation lane의 핵심 regression case다.
- `bench_a4_03`: composite 성공 사례로서 positive hard-case control 역할을 한다.

## 다음 권장 액션

1. 다음 batch는 남아 있는 쉬운 rectangular case `bench_a1_03`로 보내 easy-case consistency를 더 확인한다.
2. metric-refinement 후보를 문서화한다.
   - composite: opening/topology-sensitive
   - rectangular multi-obstacle: occupancy/blockage-sensitive
   - easy empty-room: hallucinated-obstacle penalty를 별도 태깅할지 검토
3. 이후 meshing/solver 변경이 생기면 regression subset에 explicit failure뿐 아니라 fallback-success와 easy positive-control도 함께 포함할지 검토한다.
