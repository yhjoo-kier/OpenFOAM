# CLI evaluation comparison snapshot — 60/100 milestone

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 중 complete 5-view sweep이 끝난 12개 case 비교

## 포함 case

- `bench_a1_01` — rectangular baseline
- `bench_a1_02` — rectangular easy positive-control
- `bench_a1_03` — rectangular easy positive-control / structural-vs-CFD gap case
- `bench_a1_04` — reference solver-stress rectangular case
- `bench_a2_01` — rectangular mid-difficulty case
- `bench_a2_02` — rectangular multi-obstacle control
- `bench_a2_03` — dense-obstacle rectangular stress case
- `bench_a3_01` — empty composite hallucination case
- `bench_a3_03` — mid-difficulty composite case
- `bench_a3_04` — composite representative case
- `bench_a4_02` — laminar-fallback composite stress case
- `bench_a4_03` — obstacle-dense composite stress case

현재 milestone:
- evaluation success **60 / 100**
- pending **40 / 100**

## case-level 평균 비교

| case | 성격 | avg structural | avg CFD | 핵심 신호 |
|---|---|---:|---:|---|
| `bench_a1_01` | rectangular baseline | 0.825 | 0.538 | 구조적으로 가장 안정적인 baseline |
| `bench_a1_02` | rectangular / easy positive-control | 0.825 | 0.594 | clean easy-case signal, 일부 view만 obstacle hallucination |
| `bench_a1_03` | rectangular / easy positive-control | 0.750 | 0.696 | 구조 점수는 낮은데 CFD는 가장 강한 structural-vs-CFD gap case |
| `bench_a1_04` | rectangular / reference solver-stress | 0.775 | 0.428 | 예측 경로에서는 solver stress가 일부 완화 |
| `bench_a2_01` | rectangular / mid-difficulty | 0.805 | 0.474 | 5/5 direct success, 다만 obstacle 과분할과 perspective mesh 민감도 존재 |
| `bench_a2_02` | rectangular / multi-obstacle control | 0.650 | 0.357 | opening-wall mismatch가 CFD를 강하게 깎는 control case |
| `bench_a2_03` | rectangular / dense-obstacle | 0.650 | 0.420 | obstacle layout fidelity가 주요 병목 |
| `bench_a3_01` | composite / empty-room hallucination | 0.675 | 0.588 | obstacle 3개 hallucination에도 CFD는 surprisingly 강함 |
| `bench_a3_03` | composite / mid-difficulty | 0.775 | 0.589 | composite positive-control로 활용 가능 |
| `bench_a3_04` | composite representative | 0.383 | 0.453 | structure-vs-CFD 괴리가 가장 큼 |
| `bench_a4_02` | composite / laminar-fallback | 0.738 | 0.334 | repair/solver escalation 신호가 직접 남음 |
| `bench_a4_03` | composite / obstacle-dense | 0.795 | 0.525 | composite hard case 중 가장 건강한 복원 |

## view-level 평균 비교

12개 case 전체에서 view별 평균은 다음과 같다.

| view | avg structural | avg CFD | 메모 |
|---|---:|---:|---|
| perspective | 0.715 | 0.477 | 이제 rectangular/control case까지 넓어지며 평균 CFD가 소폭 하락 |
| birdseye | 0.711 | 0.532 | `bench_a3_01` 덕분에 composite empty-room에서 매우 강한 신호 확인 |
| floorplan | 0.784 | 0.586 | 여전히 가장 높은 평균 CFD score |
| wireframe | 0.699 | 0.478 | 구조 점수는 낮아도 종종 CFD는 버틴다 |
| section | 0.692 | 0.425 | opening-wall mismatch 누적 시 가장 취약한 축으로 남음 |

## 이번 60/100 milestone에서 새로 드러난 패턴

### 1. `bench_a3_01`은 empty composite에서 obstacle hallucination과 CFD 성능이 분리될 수 있음을 보여준다

- reference obstacle 수는 `0`인데, 5개 view 모두 `3`개 obstacle을 생성했다.
- 그런데 avg CFD는 `0.588`로 높고, `birdseye`는 `0.732`, `floorplan`도 `0.697`였다.
- 즉, empty composite case에서는 hallucinated obstacle count 자체보다 room topology와 opening fidelity가 더 중요할 수 있다.

### 2. `bench_a2_02`는 rectangular multi-obstacle control에서 opening-wall fidelity가 핵심이라는 점을 강화한다

- 5-view 모두 direct success였지만 avg CFD는 `0.357`에 그쳤다.
- `floorplan`만 opening wall을 정확히 맞췄고, 이 view만 CFD가 `0.515`까지 올라갔다.
- 반대로 나머지 4개 view는 opening-wall mismatch가 남으며 `0.30~0.34` 범위에 머물렀다.

### 3. 이제 “hallucination”과 “opening mismatch”를 같은 종류의 구조 오류로 읽으면 안 된다

- `bench_a3_01`에서는 obstacle hallucination이 커도 CFD가 유지됐다.
- `bench_a2_02`에서는 obstacle 수가 크게 틀리지 않아도 opening-wall mismatch 때문에 CFD가 낮았다.
- 그래서 향후 metric refinement에서는 적어도
  - hallucinated-obstacle burden
  - opening-wall correctness
  - room/topology correctness
  를 분리해서 봐야 한다.

### 4. `floorplan` 우세는 계속 유지되지만, 왜 강한지는 case마다 다르다

- `bench_a3_01`에서는 composite block recovery + opening fidelity가 강점이었다.
- `bench_a2_02`에서는 opening-wall correctness가 거의 유일한 차별점이었다.
- 즉 `floorplan`이 늘 강한 건 맞지만, underlying reason은 empty composite와 rectangular obstacle case에서 다르게 읽어야 한다.

## 해석 메모

- benchmark는 이제 “파이프라인이 깨지는가?”보다 “무엇을 잘못 이해했을 때 CFD가 실제로 얼마나 무너지는가?”를 묻는 단계로 넘어가고 있다.
- `bench_a3_01`과 `bench_a2_02`는 둘 다 전부 direct success였지만, 물리적 해석 포인트는 전혀 다르다.
- 이 두 case가 추가되면서 current evidence는 baseline / easy / mid / control / stress 축이 조금 더 균형을 갖추게 됐다.

## 다음 권장 액션

1. 다음 batch는 아직 미평가인 composite case `bench_a3_02` 또는 `bench_a4_01`로 보내 coverage를 더 넓힌다.
2. aggregate note에 보조 태그를 추가한다.
   - empty-room hallucination
   - opening-wall mismatch severity
   - solver escalation / repair use
3. 60/100 milestone 기준으로 로컬 commit checkpoint를 남겨 이후 metric/prompt 변경의 비교 기준점으로 쓴다.
