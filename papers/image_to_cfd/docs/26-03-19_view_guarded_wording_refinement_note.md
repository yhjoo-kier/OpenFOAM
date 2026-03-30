# View-Guarded Wording Refinement Note

> 작성일: 2026-03-19  
> 범위: scale-hinted follow-up (`layout-protected span-only` → `view-guarded` wording refinement)

## 이번 턴에서 한 일

`layout-protected span-only`를 유지한 채, 다음 세 축에만 더 직접적으로 개입하는 wording refinement를 코드에 반영했다.

1. **perspective**
   - hidden depth / recessed volume 과잉 추론 억제
   - visible front/back ordering과 side-wall depth를 보수적으로 유지하라는 문구 추가

2. **section**
   - unseen ceiling height / off-cut topology / opening-wall relocation 과잉 보정 억제
   - cut-plane geometry를 먼저 믿고, regularized box로 밀어붙이지 말라는 문구 추가

3. **dense composite / A4**
   - obstacle detail보다 opening placement + connected-room topology 우선
   - uncertain dense scene에서는 speculative small boxes 대신 fewer larger solver-friendly obstacles를 선호하도록 강화

## 반영 위치

- `scripts/generate_indoor_scene_with_gemini.py`
  - 공통 prompt의 scale-hint 해석을 `soft anchor` 중심으로 완화
  - obstacle rules에 topology/opening 우선, dense-scene obstacle simplification 문구 추가

- `scripts/run_benchmark_evaluation_task.py`
  - `scale_hinted_longest_span_view_guarded_v1`의 `scale_hint.prompt_text`를 더 직접적인 guard 문구로 갱신
  - scenario builder가 `view`, `room_kind`, `obstacle_count`를 함께 사용하도록 확장
  - `perspective`, `section`, composite, dense-scene에 대해 조건부 wording을 다르게 주입하도록 수정
  - 기존 task scaffold에 남아 있던 이전 `prompt_text`도 rerun 시 자동 refresh 되도록 보정
  - 중복 `__main__` 꼬리도 함께 정리

## 현재 실행 상태

새 wording 기준 smoke subset 재실행은 완료됐다.

- evaluation root: `benchmark/evaluations_scale_hinted_view_guarded`
- smoke subset:
  - `bench_a1_01`
  - `bench_a2_03`
  - `bench_a3_03`
  - `bench_a4_03`
- 최종 상태: **20/20 success**
- 관련 산출물:
  - aggregate: `benchmark/manifests/evaluation_aggregate_summary_scale_hinted_view_guarded_smoke.json`
  - aggregate note: `docs/26-03-19_scale_hinted_view_guarded_smoke_aggregate.md`
  - comparison note: `docs/26-03-19_scale_hint_view_guarded_smoke_comparison.md`

추가 메모:
- 첫 `bench_a4_03/wireframe` 시도는 dense composite obstacle overlap으로 meshing 실패했다.
- stronger hint를 추가하지 않고, dense-scene obstacle clearance guard 한 줄만 더 넣은 뒤 단일 재시도로 recover했다.

## 현 시점 운영 판단

핵심 결론은 다음과 같다.

> `view-guarded`는 targeted side-effect mitigation에는 유효했지만, overall smoke 기준 baseline 대비 CFD degradation을 제거하지는 못했다.

좀 더 구체적으로는:
- overall mean `Lx` error 개선은 유지된다.
- perspective / section / A4 dense composite 축에서는 `layout-protected` 대비 CFD가 의미 있게 완화된다.
- 하지만 overall mean CFD는 여전히 baseline보다 낮다.

따라서 현재 상태는 **wording redesign 완료 + smoke 완료 + full rerun 보류 유지** 다.
