# CLI evaluation comparison snapshot — 70/100 milestone

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 중 complete 5-view sweep이 끝난 14개 case 비교

## 포함 case

- `bench_a1_01` — rectangular baseline
- `bench_a1_02` — rectangular easy positive-control
- `bench_a1_03` — rectangular easy positive-control / structural-vs-CFD gap case
- `bench_a1_04` — reference solver-stress rectangular case
- `bench_a2_01` — rectangular mid-difficulty case
- `bench_a2_02` — rectangular multi-obstacle control
- `bench_a2_03` — dense-obstacle rectangular stress case
- `bench_a3_01` — empty composite hallucination case
- `bench_a3_02` — empty composite control with one near-perfect perspective recovery
- `bench_a3_03` — mid-difficulty composite case
- `bench_a3_04` — composite representative case
- `bench_a4_01` — dense composite control / layout-fragile case
- `bench_a4_02` — laminar-fallback composite stress case
- `bench_a4_03` — obstacle-dense composite stress case

현재 milestone:
- evaluation success **70 / 100**
- pending **30 / 100**

## 이번 milestone에서 새로 확보된 핵심 신호

### 1. empty composite 축이 더 이상 단일 anecdote가 아니다

- `bench_a3_01`은 전 view obstacle hallucination에도 CFD가 강하게 유지되는 케이스였다.
- 새 `bench_a3_02`는 `perspective`가 거의 정답 수준(`structural_score=1.0`)으로 empty composite를 복원했고,
  반대로 `floorplan`은 `ultra_robust`까지 올라가며 더 약한 CFD를 보였다.
- 즉 empty composite에서는
  - hallucination burden
  - opening-wall correctness
  - view별 solver stress
  가 서로 독립적으로 움직일 수 있다.

### 2. A4 composite도 `stress`와 `control`이 분화되기 시작했다

- `bench_a4_03`은 dense composite hard case 중에서도 비교적 건강한 positive hard-case였다.
- 새 `bench_a4_01`은 모든 view가 original + `robust`에서 직접 끝났지만 avg CFD가 `0.412`로 낮았다.
- 따라서 A4는 이제
  - `bench_a4_01`: layout-fragile control
  - `bench_a4_02`: repair / solver-escalation stress
  - `bench_a4_03`: healthier positive hard-case
  의 3축으로 읽을 수 있다.

### 3. `floorplan` 우세는 전체 평균 사실이지만, case별로는 깨진다

- 전체적으로는 여전히 `floorplan`이 강한 view였지만,
- `bench_a3_02`에서는 `wireframe`이,
- `bench_a4_01`에서는 `wireframe`/`birdseye`가 `floorplan`보다 훨씬 강했다.
- 그래서 이후 논문 해석에서는 “floorplan best”를 단정적으로 쓰기보다,
  **easy rectangular + 일부 composite representative에서 강하지만 dense-composite ambiguity에는 취약할 수 있다** 쪽이 더 정확하다.

### 4. obstacle-count correctness만으로는 dense composite 품질을 설명할 수 없다

- `bench_a4_01`은 전 view가 reference obstacle count `3`을 그대로 맞췄다.
- 그런데 CFD는 `0.29`부터 `0.67`까지 넓게 벌어졌다.
- 반대로 empty composite에서는 obstacle hallucination이 커도 CFD가 버티는 경우가 있었다.
- 따라서 현재 benchmark 해석은 obstacle-instance metric만으로는 충분하지 않다.

### 5. bookkeeping 관점에서는 “비차단 sidecar failure”도 따로 태깅할 가치가 있다

- `bench_a4_01/perspective`는 최종적으로 original scene에서 성공했지만,
  stabilization summary 안에는 repair sidecar 실패 흔적이 남았다.
- 이건 benchmark failure는 아니지만, 파이프라인 위생 / regression 추적에는 분리 태그가 있는 편이 낫다.

## 현재 해석 메모

- benchmark는 이제 성공률 자체보다 **어떤 구조 오류가 실제 CFD를 얼마나 망가뜨리는가**를 해석하는 단계가 더 중요해졌다.
- empty composite 계열(`bench_a3_01`, `bench_a3_02`)은 hallucination-aware interpretation을 요구한다.
- dense composite 계열(`bench_a4_01`, `bench_a4_02`, `bench_a4_03`)은 opening/topology fidelity와 solver-stress bookkeeping을 함께 봐야 한다.
- rectangular 쪽은 이미 easy/mid/control/stress 축이 어느 정도 갖춰졌고, 남은 `bench_a1_05`, `bench_a2_04`, `bench_a2_05`를 채우면 카테고리 coverage도 거의 정리된다.

## 다음 권장 액션

1. 남은 rectangular coverage를 우선 채운다.
   - `bench_a1_05`
   - `bench_a2_04`
   - `bench_a2_05`
2. 이후 late composite tail (`bench_a3_05`, `bench_a4_04`, `bench_a4_05`)을 마저 평가해 frozen-20 CLI sweep completion을 향해 간다.
3. aggregate note / metric 보조 태그에 아래 두 항목을 추가 검토한다.
   - non-blocking repair-sidecar failure
   - hallucinated-obstacle burden vs opening-wall mismatch 분리
