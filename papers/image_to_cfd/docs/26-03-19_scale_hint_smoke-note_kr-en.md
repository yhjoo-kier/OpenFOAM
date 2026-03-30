# Scale-Hint Smoke Note / 1st Pivot Check

작성일: 2026-03-19  
프로젝트: OpenFOAM Image-to-CFD benchmark  
토픽: SurfClaw topic 2242

## 1. 결정사항 (current working spec)

- **Baseline 보존**: 기존 `benchmark/evaluations`는 그대로 `no-scale-hint baseline` artifact로 유지한다.
- **Scale-hinted 분리 저장**: 새 실험은 `benchmark/evaluations_scale_hinted` 아래에 별도 scaffold/evaluation root로 저장한다.
- **Setting identifier**
  - baseline: `no_scale_hint_baseline`
  - hinted: `scale_hinted_longest_horizontal_span_v1`
- **Scale hint rule v1**: reference scene의 **longest horizontal span = max(Lx, Ly)** 를 absolute scale anchor로 사용한다.
- **Prompt wording**: “the longest horizontal span of the room is approximately X m” 형태로 주입한다.
- **Post-generation handling**: 현재는 Gemini prompt 주입 + task/evaluation bookkeeping은 동작하며, generator 쪽에는 hinted span 반영/기록 경로를 넣었다.

## 2. 코드 반영 범위

수정 파일(working tree 기준):
- `scripts/scaffold_benchmark_evaluations.py`
- `scripts/run_benchmark_evaluation_task.py`
- `scripts/run_indoor_stabilized.py`
- `scripts/generate_indoor_scene_with_gemini.py`
- `scripts/build_benchmark_evaluation_aggregate.py`

핵심 반영 내용:
- scale-hinted 전용 evaluation root scaffold 가능
- task/evaluation summary에 `setting`, `scale_hint`, `evaluation_root` 기록
- `run_indoor_stabilized.py`가 `--scale-hint`를 generation 단계로 전달
- `generate_indoor_scene_with_gemini.py`에 `--scale-hint` 인자/summary 경로 추가

## 3. Smoke subset (actual run status)

대표 subset 의도:
- easy rectangular: `bench_a1_01`
- rectangular multi-obstacle: `bench_a2_03`
- light composite: `bench_a3_03`
- dense composite: `bench_a4_03`

실제 완료된 smoke rerun:
- `bench_a1_01 / perspective`
- `bench_a1_01 / floorplan`
- `bench_a2_03 / perspective`

배치 loop는 중간에 종료되어 composite 2건은 아직 미실행.

## 4. Early comparison vs baseline

### bench_a1_01 / perspective
- **Lx error**: 개선 (`1.334 -> 0.666`)
- **Ly error**: 악화 (`0.140 -> 1.140`)
- **CFD score**: 소폭 악화 (`0.646986 -> 0.632667`)
- 해석: longest-span anchor는 잡히지만 aspect ratio drift가 남아 있음.

### bench_a1_01 / floorplan
- **structural score**: 개선 (`0.875 -> 1.000`)
- **CFD score**: 거의 동일/소폭 개선 (`0.392970 -> 0.395409`)
- **bbox error**는 mixed result
- 해석: floorplan 계열은 scale hint가 layout structure 쪽엔 도움 가능성.

### bench_a2_03 / perspective
- 최신 rerun 기준:
  - **Lx error**: 크게 개선 (`1.447 -> 0.003`)
  - **Ly error**: 악화 (`1.064 -> 1.864`)
  - **structural score**: 개선 (`0.625 -> 0.7083`)
  - **CFD score**: 소폭 악화 (`0.373465 -> 0.348134`)
- 해석: 힌트는 longest span 자체는 매우 잘 고정하지만, 나머지 span/개구 wall/orientation mismatch는 아직 잔존.

## 5. Current reading

현재 v1 scale hint는 다음 성격을 보임:
- **장점**: longest horizontal span mismatch는 실제로 줄일 수 있음.
- **한계**: one-number hint만으로는 **aspect ratio / opening wall / obstacle placement** 까지 자동 개선되지는 않음.
- 따라서 full rerun 전에는 최소한 아래를 더 확인해야 함:
  1. composite 2 case smoke 보완
  2. floorplan / perspective별 일관성 확인
  3. 필요 시 v2 hint 문구(예: longest span + not-both-spans warning, or span+height) 검토

## 6. Provisional full-rerun decision

- **아직 full frozen-20 rerun 확정 아님**
- 이유: geometry scale anchor는 유효하지만, CFD score improvement가 아직 일관적이지 않음.
- 현 단계 권고: **composite smoke 추가 후 rerun 범위 재결정**
