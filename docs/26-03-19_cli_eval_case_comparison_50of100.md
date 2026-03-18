# CLI evaluation comparison snapshot — 50/100 milestone

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 중 complete 5-view sweep이 끝난 10개 case 비교

## 포함 case

- `bench_a1_01` — rectangular baseline
- `bench_a1_02` — rectangular easy positive-control
- `bench_a1_03` — rectangular easy positive-control / structural-vs-CFD gap case
- `bench_a1_04` — reference solver-stress rectangular case
- `bench_a2_01` — rectangular mid-difficulty case
- `bench_a2_03` — dense-obstacle rectangular stress case
- `bench_a3_03` — mid-difficulty composite case
- `bench_a3_04` — composite representative case
- `bench_a4_02` — laminar-fallback composite stress case
- `bench_a4_03` — obstacle-dense composite stress case

현재 milestone:
- evaluation success **50 / 100**
- pending **50 / 100**

## case-level 평균 비교

| case | 성격 | avg structural | avg CFD | 핵심 신호 |
|---|---|---:|---:|---|
| `bench_a1_01` | rectangular baseline | 0.825 | 0.538 | 구조적으로 가장 안정적인 baseline |
| `bench_a1_02` | rectangular / easy positive-control | 0.825 | 0.594 | clean easy-case signal, 일부 view만 obstacle hallucination |
| `bench_a1_03` | rectangular / easy positive-control | 0.750 | 0.696 | 구조 점수는 낮은데 CFD는 가장 강한 structural-vs-CFD gap case |
| `bench_a1_04` | rectangular / reference solver-stress | 0.775 | 0.428 | 예측 경로에서는 solver stress가 일부 완화 |
| `bench_a2_01` | rectangular / mid-difficulty | 0.805 | 0.474 | 5/5 direct success, 다만 obstacle 과분할과 perspective mesh 민감도 존재 |
| `bench_a2_03` | rectangular / dense-obstacle | 0.650 | 0.420 | obstacle layout fidelity가 주요 병목 |
| `bench_a3_03` | composite / mid-difficulty | 0.775 | 0.589 | composite positive-control로 활용 가능 |
| `bench_a3_04` | composite representative | 0.383 | 0.453 | structure-vs-CFD 괴리가 가장 큼 |
| `bench_a4_02` | composite / laminar-fallback | 0.738 | 0.334 | repair/solver escalation 신호가 직접 남음 |
| `bench_a4_03` | composite / obstacle-dense | 0.795 | 0.525 | composite hard case 중 가장 건강한 복원 |

## view-level 평균 비교

10개 case 전체에서 view별 평균은 다음과 같다.

| view | avg structural | avg CFD | 메모 |
|---|---:|---:|---|
| perspective | 0.721 | 0.504 | 이제 평균 CFD가 0.5를 넘었지만 obstacle hallucination은 여전 |
| birdseye | 0.728 | 0.532 | easy rectangular positive-control에서 매우 강함 |
| floorplan | 0.791 | 0.582 | 여전히 가장 높은 평균 CFD score |
| wireframe | 0.702 | 0.487 | object hallucination/mesh warning이 섞여도 CFD는 종종 유지 |
| section | 0.718 | 0.420 | 여전히 최하위지만 `bench_a1_03`로 인해 '항상 약함' 해석은 완화 |

## 현재까지 드러난 패턴

### 1. rectangular 축은 이제 baseline–easy–mid–stress coverage가 꽤 탄탄해졌다

- completed rectangular cases가 `bench_a1_01`, `bench_a1_02`, `bench_a1_03`, `bench_a1_04`, `bench_a2_01`, `bench_a2_03`까지 확장됐다.
- 덕분에 rectangular 쪽에서는
  - baseline
  - easy positive-controls
  - mid-difficulty control
  - solver/obstacle stress
  를 꽤 분리해서 읽을 수 있다.

### 2. `bench_a1_03`는 structural score만으로는 설명되지 않는 강한 positive signal이다

- 모든 view가 structural score `0.75`에 머물렀다.
- 그런데 avg CFD는 전체 최고(`0.696`)였다.
- 이는 current structural aggregate가 easy open-room cases에서 flow-relevant similarity를 충분히 대변하지 못할 수 있음을 보여준다.

### 3. `floorplan` 우세는 유지되지만 easy rectangular에서는 `birdseye`와 `section`도 반등할 수 있다

- 전체 평균에서는 여전히 `floorplan`이 최고 CFD (`0.582`)다.
- 하지만 `bench_a1_03`에서는 `birdseye`가 최고였고, `section`도 높은 CFD를 기록했다.
- 따라서 뷰 타입 우열은 case topology와 opening fidelity에 따라 달라지는 면이 있다.

### 4. obstacle hallucination은 easy case에서도 반복되며, 이 현상 자체가 benchmark 해석 포인트다

- `bench_a1_02`와 `bench_a1_03` 모두 일부 view가 존재하지 않는 obstacle을 여러 개 만들었다.
- 그런데 CFD는 유지되거나 오히려 강했다.
- 그래서 향후 metric refinement에서는 'hallucinated obstacle count'를 기록하되, 실제 blockage/topology 영향과 분리해 읽는 편이 좋다.

### 5. mesh-risk warning은 점점 더 '주의 표지'에 가깝게 보인다

- `bench_a1_02`, `bench_a1_03` 모두 5/5 direct success였지만 다수 run이 `high` risk를 유지했다.
- 즉, 현 시점의 risk flags는 failure predictor보다는 conservative mesh-quality warning 역할에 가깝다.
- regression subset에서는 단순 risk-level보다 fallback/repair/escalation 여부를 우선 신호로 보는 편이 낫다.

## case별 한 줄 판정

- `bench_a1_01`: baseline으로 계속 유지할 가치가 높다.
- `bench_a1_02`: easy rectangular positive-control로 중요하며, 일부 view의 hallucination을 분리해서 보여준다.
- `bench_a1_03`: 구조 점수 대비 CFD가 과도하게 강한 핵심 metric-gap case다.
- `bench_a1_04`: reference stress case지만 predicted geometry가 stress를 완화하는지 보는 용도로 유용하다.
- `bench_a2_01`: rectangular mid-difficulty control로 중요하며, obstacle 과분할/mesh sensitivity 신호를 동시에 준다.
- `bench_a2_03`: obstacle-centric metric refinement의 핵심 regression case다.
- `bench_a3_03`: composite positive-control로서 매우 중요하다.
- `bench_a3_04`: composite metric gap을 가장 잘 드러내는 대표 case다.
- `bench_a4_02`: repair + solver escalation lane의 핵심 regression case다.
- `bench_a4_03`: composite 성공 사례로서 positive hard-case control 역할을 한다.

## 다음 권장 액션

1. 다음 batch는 아직 미평가인 case들 중 coverage 확장 가치가 높은 `bench_a2_02` 또는 `bench_a3_01`로 보낸다.
2. metric-refinement 후보를 문서화한다.
   - composite: opening/topology-sensitive
   - rectangular multi-obstacle: occupancy/blockage-sensitive
   - easy empty-room: hallucinated-obstacle tagging + topology-preserving score 보조축
3. 50/100 milestone 기준으로 로컬 commit checkpoint를 남겨 이후 regression 비교 기준점으로 쓴다.
