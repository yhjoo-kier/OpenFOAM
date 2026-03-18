# CLI evaluation tail completion — final composite sweeps to 100/100

> 작성일: 2026-03-19
> 범위: `bench_a3_05`, `bench_a4_04`, `bench_a4_05`

## Summary

Late composite tail 3개 case에 대한 5-view Gemini CLI sweep을 모두 마쳤다.
이로써 frozen-20 benchmark의 image-conditioned evaluation은 **100/100 success**로 닫혔다.

이번 tail completion에서 관찰된 핵심은 다음 세 가지다.

1. **A3 tail(`bench_a3_05`)은 opening-wall fidelity가 유지되면 obstacle hallucination이 있어도 CFD가 꽤 잘 버틴다.**
2. **A4 tail control(`bench_a4_04`)은 구조 복원은 매우 강하지만 CFD는 중간 수준에 머문다.**
3. **A4 tail stress/control(`bench_a4_05`)은 전 view task-level success였지만, `section`에서 composite 붕괴와 opening-wall mismatch가 다시 나타났다.**

---

## 1) `bench_a3_05` — late composite tail (A3)

평균 점수:
- 평균 structural score: **0.6750**
- 평균 CFD score: **0.5425**

| View | Structural | CFD | Solver path | 해석 메모 |
|---|---:|---:|---|---|
| perspective | 0.6250 | 0.5042 | original + `robust` @ `0.25` | 유일한 mesh-size fallback case |
| birdseye | 0.6250 | 0.5576 | original + `robust` @ `0.35` | obstacle 1 → 3 hallucination |
| floorplan | 0.7500 | 0.6930 | original + `robust` @ `0.35` | 이번 sweep 최고 CFD |
| wireframe | 0.6250 | 0.4399 | original + `robust` @ `0.35` | obstacle 1 → 3 hallucination |
| section | 0.7500 | 0.5175 | original + `robust` @ `0.35` | obstacle 수와 opening wall 모두 유지 |

### Key signals

- 모든 view가 **composite room kind**와 **opening-wall orientation**을 유지했다.
- 그럼에도 `birdseye` / `floorplan` / `wireframe`은 reference obstacle 1개를 3개로 과분할했다.
- 그런데 CFD는 `floorplan` 기준 **0.6930**까지 올라갔고, `birdseye`도 **0.5576**으로 준수했다.
- 즉 A3 composite tail에서도 **opening/topology fidelity가 obstacle instance exactness보다 더 중요**하다는 기존 신호가 다시 확인됐다.
- `perspective`만 `mesh_size=0.25` fallback이 필요해, 이 케이스는 “구조는 맞지만 meshing 민감도는 남아 있는” composite tail positive-control로 볼 수 있다.

---

## 2) `bench_a4_04` — dense composite tail control (A4)

평균 점수:
- 평균 structural score: **0.8167**
- 평균 CFD score: **0.4664**

| View | Structural | CFD | Solver path | 해석 메모 |
|---|---:|---:|---|---|
| perspective | 0.6250 | 0.3886 | original + `robust` @ `0.35` | 구조는 보통, CFD는 약함 |
| birdseye | 0.9167 | 0.4735 | original + `robust` @ `0.35` | 구조 강함 |
| floorplan | 1.0000 | 0.4596 | original + `robust` @ `0.35` | 구조는 완벽하지만 CFD는 중간 |
| wireframe | 0.9167 | 0.5017 | original + `robust` @ `0.35` | 이번 sweep 최고 CFD |
| section | 0.6250 | 0.5086 | original + `robust` @ `0.35` | opening wall mismatch 있음 |

### Key signals

- 다섯 view 모두 **original scene + robust + 0.35**에서 직접 끝났다. repair나 solver escalation이 필요하지 않았다.
- 모든 view가 reference obstacle count `3`을 그대로 맞췄고, `section`을 제외하면 opening wall도 맞췄다.
- 그런데 `floorplan`은 structural score가 **1.0**인데도 CFD는 **0.4596**에 그쳤다.
- 반대로 `section`은 opening wall mismatch가 있는데도 CFD가 **0.5086**으로 가장 높지는 않지만 강한 편이다.
- 따라서 dense composite tail에서는 **구조 점수 고득점 자체가 곧바로 높은 CFD fidelity를 보장하지 않는다**는 점이 더 분명해졌다.

---

## 3) `bench_a4_05` — dense composite tail stress/control (A4)

평균 점수:
- 평균 structural score: **0.6208**
- 평균 CFD score: **0.3947**

| View | Structural | CFD | Solver path | 해석 메모 |
|---|---:|---:|---|---|
| perspective | 0.7500 | 0.3797 | original + `robust` @ `0.35` | opening-wall mismatch + non-blocking repair sidecar error |
| birdseye | 0.5625 | 0.4280 | original + `robust` @ `0.35` | non-blocking repair sidecar error |
| floorplan | 0.7500 | 0.4564 | original + `robust` @ `0.35` | 이번 sweep 최고 CFD |
| wireframe | 0.6250 | 0.3782 | original + `robust` @ `0.35` | count는 유지했지만 fidelity는 보통 |
| section | 0.4167 | 0.3311 | original + `robust` @ `0.35` | composite → non-composite collapse, obstacle 4 → 3 |

### Key signals

- task-level로는 5/5 success지만, 세부 신호는 꽤 거칠다.
- `perspective`와 `birdseye`는 최종 성공과 별개로 **repair sidecar error**가 남았다.
- `perspective`는 obstacle count `4`를 유지했음에도 opening wall mismatch 때문에 구조/CFD가 모두 눌렸다.
- `section`은 이번 tail 중 가장 취약한 뷰로,
  - `room_kind_match = false`
  - opening wall mismatch
  - obstacle count `4 → 3`
  - structural score **0.4167**
  - CFD score **0.3311**
  로 마감됐다.
- 즉 A4 dense tail에서는 여전히 **section view가 composite topology 붕괴를 일으키는 대표 취약 뷰**다.

---

## Tail completion interpretation

이번 최종 3개 case는 frozen-20 CLI benchmark를 다음처럼 닫아준다.

- **A3 tail (`bench_a3_05`)**: composite positive-control with mild meshing sensitivity
- **A4 tail control (`bench_a4_04`)**: structurally strong but CFD-mid dense composite case
- **A4 tail stress/control (`bench_a4_05`)**: task-level success yet topology/opening fragility remains visible

즉 benchmark는 이제 단순 성공률 문제가 아니라,

1. obstacle hallucination이 실제 CFD를 얼마나 망치는지,
2. opening-wall mismatch가 언제 결정적으로 작동하는지,
3. dense composite에서 section view가 왜 반복적으로 약한지,
4. repair sidecar failure를 bookkeeping 상 별도 태그로 둘 필요가 있는지

를 읽어내는 해석 단계로 넘어갔다.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_bench_a3_05_a4_04_a4_05.json`
- Tail rerun summary: `benchmark/manifests/evaluation_batch_bench_a4_05_tail.json`
- Aggregate status: `benchmark/evaluations/summary.json`
