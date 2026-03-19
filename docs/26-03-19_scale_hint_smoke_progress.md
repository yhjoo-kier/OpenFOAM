# Scale-Hint Smoke Progress (complete, 20/20 tasks)

> 작성일: 2026-03-19
> 상태: 완료
> 대상 smoke subset: `bench_a1_01` (easy rectangular), `bench_a2_03` (rectangular multi-obstacle), `bench_a3_03` (light composite), `bench_a4_03` (dense composite)
> 범위: 4 cases × 5 views = 20 tasks
> 결과: **20/20 success**

## 1. 이번 턴에서 정리한 구현/운영 상태

### 코드/운영 결함 수정
- `scripts/run_benchmark_evaluation_batch.py`
  - `--evaluation-root` 인자를 추가하여 baseline root와 `evaluations_scale_hinted` root를 명시적으로 분리 실행 가능하게 수정
  - 기존의 미정의 변수(`evaluation_root`)와 중복 `if __name__ == "__main__"` / stray `Exit(main())` 문제를 정리
- `scripts/run_benchmark_evaluation_task.py`
  - evaluation index refresh 시 summary의 `evaluation_root`가 항상 baseline root를 가리키던 bookkeeping 오류 수정
- `scripts/build_benchmark_evaluation_aggregate.py`
  - 임의 evaluation root를 집계할 수 있게 정리
  - `build_payload()` 호출/정의 불일치 버그 수정
  - scale-sensitive 집계를 위해 room bbox relative error, room volume relative error, laminar fallback count도 같이 집계 가능하도록 확장

### 실행 명령
```bash
python3 scripts/run_benchmark_evaluation_batch.py \
  --evaluation-root benchmark/evaluations_scale_hinted \
  --cases bench_a1_01 bench_a2_03 bench_a3_03 bench_a4_03 \
  --statuses pending \
  --backend cli \
  --model gemini-3-flash-preview \
  --mesh-size 0.35 \
  --end-time 250 \
  --solver-timeout 600 \
  --summary-out benchmark/manifests/evaluation_batch_summary_scale_hinted_smoke.json
```

## 2. smoke subset completion status

- `bench_a1_01`: 5/5 views success
- `bench_a2_03`: 5/5 views success
- `bench_a3_03`: 5/5 views success
- `bench_a4_03`: 5/5 views success

즉, 선택한 smoke subset에서는 **baseline과 동일하게 task success 자체는 20/20** 이다.

## 3. baseline 대비 smoke summary (20 task complete)

비교 기준:
- baseline root: `benchmark/evaluations`
- scale-hinted root: `benchmark/evaluations_scale_hinted`

### mean delta (scale-hinted - baseline)
- structural score: **+0.0299**
- CFD score: **-0.0276**
- room bbox relative error:
  - `Lx`: **-0.2583**
  - `Ly`: **+0.0428**
  - `Lz`: **+0.0537**
- room volume relative error: **-0.0030**
- room-kind match count delta: **0**
- opening-wall match count delta: **-1**

### case-level mean delta
- `bench_a1_01`
  - structural **-0.0750**
  - CFD **-0.0024**
  - `Lx` rel error **-0.0797**
  - `Ly` rel error **+0.1697**
  - `Lz` rel error **+0.0085**
- `bench_a2_03`
  - structural **+0.0667**
  - CFD **+0.0183**
  - `Lx` rel error **-0.3761**
  - `Ly` rel error **+0.0640**
  - `Lz` rel error **+0.1442**
- `bench_a3_03`
  - structural **+0.0000**
  - CFD **-0.0815**
  - `Lx` rel error **-0.1938**
  - `Ly` rel error **-0.0772**
  - `Lz` rel error **-0.0020**
- `bench_a4_03`
  - structural **+0.1281**
  - CFD **-0.0449**
  - `Lx` rel error **-0.3837**
  - `Ly` rel error **+0.0145**
  - `Lz` rel error **+0.0641**

## 4. 현재 해석

### 긍정 신호
1. **absolute horizontal scale ambiguity는 실제로 줄어든다.**
   - smoke 전체 평균에서 `Lx` relative error가 **-0.2583** 개선됐다.
   - 특히 `A2`, `A4` representative case에서 longest horizontal span anchor 효과가 크다.

2. **task success는 손상되지 않았다.**
   - smoke subset에서는 baseline과 동일하게 **20/20 success** 를 유지했다.

3. **bookkeeping 분리는 이제 안정적이다.**
   - baseline `benchmark/evaluations/`와 scale-hinted `benchmark/evaluations_scale_hinted/`가 섞이지 않는다.
   - task/evaluation summary에 `setting`과 `scale_hint`가 같이 남는다.

### 부정/주의 신호
1. **longest-horizontal-span only hint는 아직 physics 개선을 보장하지 않는다.**
   - smoke 평균 CFD는 오히려 **-0.0276** 하락했다.
   - `A3`, `A4` representative case에서는 Lx scale은 좋아졌지만 CFD는 악화되는 방향이 나타났다.

2. **`Ly`, `Lz`는 같이 좋아지지 않는다.**
   - `Lx`는 분명히 개선되지만 `Ly`, `Lz` relative error는 소폭 악화되었다.
   - 즉, current hint는 “가로 span anchor”에는 유효하지만 전체 3D metric regularization으로는 불충분하다.

3. **opening/layout fidelity를 자동으로 보장하지 않는다.**
   - opening-wall match count는 smoke 전체에서 1건 감소했다.
   - structural score는 평균적으로 소폭 상승했지만, 그 자체가 CFD 개선으로 이어지지 않았다.

## 5. 잠정 결론

### 이미 확정 가능한 것
- benchmark dataset 자체는 그대로 유지해도 된다.
- scale-hinted setting을 별도 root로 운영하는 전략은 타당하다.
- `longest horizontal span` 힌트는 **room long-span scale mismatch 완화용 1차 anchor** 로는 유효하다.

### 아직 확정하면 안 되는 것
- `longest horizontal span only`를 최종 main setting으로 바로 고정하는 결정
- frozen-20 full rerun을 지금 즉시 시작하는 결정

현재 smoke 결과는 다음처럼 요약된다:

> **"geometry long-span scale은 꽤 좋아지지만, CFD fidelity는 아직 mixed-to-negative다."**

따라서 다음 단계는 full rerun보다 먼저,

1. `ceiling height + longest horizontal span` 2-anchor variant를 소규모로 추가 비교하고,
2. 그 결과가 현재 single-anchor보다 낫거나 최소한 CFD degradation을 줄이는지 본 뒤,
3. 그 다음 frozen-20 full rerun 범위를 확정하는 것이 합리적이다.

## 6. immediate next action
- `ceiling height + longest horizontal span` 보강안을 prompt/hint spec에 넣은 variant를 소규모 subset으로 추가 실행
- same 4-case smoke subset 또는 더 작은 2-case stress subset으로 single-anchor vs dual-anchor 비교
- dual-anchor가 CFD degradation을 줄이지 못하면, full rerun 없이 scale-hinted pivot을 재설계
