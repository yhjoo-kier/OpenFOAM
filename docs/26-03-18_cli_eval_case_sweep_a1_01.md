# CLI-backed benchmark evaluation case sweep — `bench_a1_01`

> 작성일: 2026-03-18
> 범위: frozen-20 benchmark 중 대표 A1 case 하나(`bench_a1_01`)에 대해 5-view image-conditioned evaluation을 모두 실행

## 목적

이전 cron turn에서 `perspective`, `section` 두 개의 smoke test 성공을 확인했다.
이번에는 같은 case에서 나머지 `birdseye`, `floorplan`, `wireframe`까지 모두 채워서,
**동일 형상 기준 view-type별 구조/CFD 신호 차이**를 빠르게 확인하는 첫 complete case sweep을 만든다.

## 실행 경로

- runner 추가: `scripts/run_benchmark_evaluation_batch.py`
- backend: Gemini CLI cached auth
- model: `gemini-3-flash-preview`
- task filter: `bench_a1_01`의 pending view 3개

실행 예시:

```bash
python3 scripts/run_benchmark_evaluation_batch.py \
  --cases bench_a1_01 \
  --backend cli \
  --model gemini-3-flash-preview \
  --statuses pending \
  --summary-out benchmark/manifests/evaluation_batch_summary_bench_a1_01.json
```

## 결과 요약

`bench_a1_01`의 5개 view가 모두 **end-to-end 성공**했다.

| view | structural_score | cfd_score | 메모 |
|---|---:|---:|---|
| perspective | 0.750 | 0.647 | 현재 A1 대표 case에서 가장 높은 CFD score |
| birdseye | 0.750 | 0.606 | perspective와 유사한 구조 점수, CFD는 약간 낮음 |
| wireframe | 0.750 | 0.600 | birdseye와 유사 |
| floorplan | 0.875 | 0.393 | 구조는 꽤 잘 맞지만 CFD score는 더 낮음 |
| section | 1.000 | 0.443 | 구조는 완전 일치, 하지만 CFD는 중간 수준 |

## 1차 해석

1. **view별 구조 추론과 CFD 유사도는 단순히 같이 움직이지 않는다.**
   - `section`은 구조적으로는 완전 복원했지만 `cfd_score`는 최고가 아니다.
   - 즉, 동일한 topology를 맞춰도 opening 세부 배치/장애물 위치/solver 민감도 차이 때문에 유동장 오차는 별개로 남을 수 있다.

2. **A1 대표 case에서는 perspective가 가장 실용적인 입력 신호로 보인다.**
   - structural score는 최고가 아니지만 현재 5-view 중 `cfd_score`가 가장 높다.
   - depth cue가 opening/obstacle 배치 추정에 유리했을 가능성이 있다.

3. **floorplan은 평면 구조 추정에는 강하지만 CFD fidelity는 별도 검증이 필요하다.**
   - 높이/깊이 ambiguity가 남아도 구조 metric 일부는 좋게 나올 수 있다.
   - 논문에서는 structure-only metric과 CFD metric을 함께 보여줘야 함을 다시 확인했다.

## 산출물

- batch summary:
  - `benchmark/manifests/evaluation_batch_summary_bench_a1_01.json`
- per-view evaluation summaries:
  - `benchmark/evaluations/bench_a1_01/*/evaluation_summary.json`
- predicted cases/results:
  - `benchmark/evaluations/bench_a1_01/*/predicted_case`
  - `benchmark/evaluations/bench_a1_01/*/predicted_results`

## 다음 권장 액션

1. **composite 대표 case 하나(A3 또는 A4)도 동일한 5-view sweep**을 실행해서,
   직사각형 case에서 보인 경향이 composite room에도 유지되는지 확인한다.
2. stress signal이 큰 case(`bench_a2_03`, `bench_a4_03`, `bench_a1_04`, `bench_a4_02`)는
   모든 view를 돌리기 전에 대표 view 1~2개로 먼저 실패/성공 신호를 수집한다.
3. case가 2~3개 누적되면 view-type별 평균 `structural_score`, `cfd_score`, success rate를
   집계하는 aggregate summary script를 추가한다.
