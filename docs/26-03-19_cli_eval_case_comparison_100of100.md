# CLI evaluation comparison snapshot — 100/100 completion

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark의 image-conditioned CLI evaluation 전체 20개 case / 100개 task

## Completion summary

최종 상태:
- evaluation success **100 / 100**
- pending **0 / 100**
- blocked **0 / 100**
- running **0 / 100**

이제 frozen-20 benchmark는 아래 4축이 모두 닫혔다.

1. **scene JSON** — 20/20
2. **reference CFD** — 20/20
3. **5-view render bundle** — 20/20
4. **image-conditioned CLI evaluation** — 100/100

## Case-level average snapshot

| Category | Cases | Mean structural score | Mean CFD score | 해석 메모 |
|---|---:|---:|---:|---|
| A1 | 5 | 0.7850 | 0.5567 | 쉬운 rectangular 계열은 여전히 가장 안정적 |
| A2 | 5 | 0.6727 | 0.4438 | multi-obstacle rectangular에서 opening/blockage 민감도 큼 |
| A3 | 5 | 0.6567 | 0.5466 | composite인데도 CFD는 의외로 강한 편; hallucination burden 해석 필요 |
| A4 | 5 | 0.7124 | 0.4266 | dense composite는 구조 점수 대비 CFD가 눌리는 경향 |

## 이번 100/100 milestone에서 확정된 핵심 해석

### 1. benchmark는 더 이상 “성공률” 문제가 아니다

- reference bundle도 20/20,
- image-conditioned CLI eval도 100/100
으로 닫혔다.

즉 지금부터 중요한 것은 **실패율 개선**보다
**어떤 구조 오류가 실제 CFD fidelity를 얼마나 깎는지**를 읽는 일이다.

### 2. A3 composite는 obstacle hallucination이 있어도 CFD가 버틸 수 있다

- `bench_a3_01`, `bench_a3_02`, `bench_a3_05`에서 반복적으로 보였듯,
  obstacle count mismatch가 커도 opening / room topology가 맞으면 CFD가 꽤 유지된다.
- 특히 `bench_a3_05/floorplan`은 obstacle 1 → 3 hallucination인데도 `cfd_score ≈ 0.693`까지 올라갔다.

따라서 composite A3 계열에서는
**hallucinated-obstacle burden**과 **실제 CFD penalty**를 분리 태깅하는 쪽이 타당하다.

### 3. A4 dense composite는 “구조 복원”과 “CFD 복원”이 분리된다

- `bench_a4_04/floorplan`은 structural score `1.0`인데 CFD는 `0.4596`이다.
- `bench_a4_01`, `bench_a4_05`도 obstacle count와 room kind를 대체로 맞추지만 CFD는 중간 이하로 남는다.
- 반대로 `bench_a4_03`은 dense composite hard-case인데도 평균 CFD가 `0.5255`로 가장 건강하다.

즉 dense composite에서는
**obstacle count / room kind correctness만으로는 downstream CFD fidelity를 설명하기 어렵다.**

### 4. section view는 여전히 composite topology에 취약하다

- `bench_a3_03`, `bench_a3_04`, `bench_a4_05`에서 section view는 구조적으로 그럴듯해 보여도 opening-wall mismatch나 topology collapse로 CFD가 약해졌다.
- 특히 `bench_a4_05/section`은 composite 붕괴와 obstacle count 감소까지 동반했다.

논문 해석에서는 section view를 단순 “높이 정보가 있는 중간 난도 view”로 보기보다,
**composite topology에는 취약한 특수 뷰**로 다루는 편이 더 정확하다.

### 5. non-blocking repair sidecar failure는 별도 bookkeeping 가치가 있다

- `bench_a4_01/perspective`, `bench_a4_05/perspective`, `bench_a4_05/birdseye` 등은 최종 성공했지만 repair sidecar failure가 로그에 남았다.
- 이들은 benchmark failure는 아니지만, 파이프라인 위생과 regression 추적에는 중요하다.

따라서 후속 metric/manifest에는 최소한 아래 두 태그가 있으면 좋다.
- `repair_sidecar_warning`
- `hallucinated_obstacle_burden`

## 추천 후속 작업

1. **aggregate analysis note 고도화**
   - case별 평균 structural / CFD / opening-wall mismatch / obstacle hallucination burden을 한 표로 정리
2. **metric 보조 태그 추가 검토**
   - hallucinated-obstacle burden
   - non-blocking repair-sidecar warning
   - opening-wall mismatch severity
3. **stress subset 재확인 lane 정리**
   - `bench_a2_01`, `bench_a4_02`, `bench_a2_03`, `bench_a4_03`, `bench_a1_04`
   - + positive/control anchors: `bench_a1_02`, `bench_a1_03`, `bench_a1_05`, `bench_a3_01`, `bench_a3_02`, `bench_a4_01`, `bench_a2_02`, `bench_a2_04`
4. **논문용 결과 섹션 초안 준비**
   - category × view별 대표 signal과 failure mode를 Results/Discussion 문장으로 옮길 수 있는 상태다.

## Artifact paths

- Aggregate evaluation summary: `benchmark/evaluations/summary.json`
- Tail completion note: `docs/26-03-19_cli_eval_tail_completion_100of100.md`
- Final tail batch manifests:
  - `benchmark/manifests/evaluation_batch_bench_a3_05_a4_04_a4_05.json`
  - `benchmark/manifests/evaluation_batch_bench_a4_05_tail.json`
