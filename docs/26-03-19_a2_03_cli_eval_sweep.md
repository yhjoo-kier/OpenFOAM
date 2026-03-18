# CLI-backed benchmark evaluation case sweep — `bench_a2_03`

> 작성일: 2026-03-19
> 범위: frozen-20 benchmark 중 직사각형 hard case `bench_a2_03`에 대해 5-view image-conditioned evaluation을 모두 실행

## 목적

이전까지는 `bench_a1_01`(직사각형 baseline)과 composite hard/representative cases(`bench_a3_04`, `bench_a4_02`, `bench_a4_03`) 위주로 신호를 모았다.
이번 sweep은 **직사각형이지만 obstacle-dense한 stress case**인 `bench_a2_03`를 같은 5-view 프레임에서 실행해,
composite difficulty와 rectangular difficulty를 분리해서 비교하기 위한 것이다.

## 실행 경로

- runner: `scripts/run_benchmark_evaluation_batch.py`
- backend: Gemini CLI cached auth
- model: `gemini-3-flash-preview`
- task filter: `bench_a2_03`의 5개 view

실행 예시:

```bash
python3 scripts/run_benchmark_evaluation_batch.py \
  --cases bench_a2_03 \
  --backend cli \
  --model gemini-3-flash-preview \
  --statuses pending \
  --summary-out benchmark/manifests/evaluation_batch_summary_a2_03.json
```

## 결과 요약

`bench_a2_03`의 5개 view가 모두 **end-to-end 성공**했다.

| view | structural_score | cfd_score | 메모 |
|---|---:|---:|---|
| perspective | 0.625 | 0.373 | 방은 직사각형으로 맞췄지만 obstacle 배치가 크게 붕괴 |
| birdseye | 0.625 | 0.356 | perspective와 거의 같은 실패 패턴 |
| floorplan | 0.750 | 0.655 | opening-wall은 맞추고 CFD도 가장 높았지만 obstacle IoU는 여전히 0 |
| wireframe | 0.625 | 0.355 | obstacle 인식 실패 + mesh risk 높음 |
| section | 0.625 | 0.359 | 단면도도 obstacle 복원에는 거의 도움을 못 줌 |

공통적으로 모든 view가 **solver는 `mesh_size=0.35 + robust`에서 원본 scene 그대로 성공**했다.
즉, 이 case의 핵심 난점은 현재로서는 solver 안정성보다 **image-conditioned obstacle reconstruction fidelity** 쪽에 더 가깝다.

## 세부 관찰

1. **모든 view에서 obstacle 개수는 3개로 맞췄지만, 실제 위치/크기/겹침 관계는 거의 복원하지 못했다.**
   - `floorplan`조차 obstacle match 기준 IoU는 전부 0이었다.
   - 구조 점수는 room kind / opening wall / obstacle count 같은 coarse signal 덕분에 0.75까지 올라가지만,
     실제 obstacle geometry fidelity는 훨씬 낮다.

2. **직사각형 hard case에서는 composite보다 failure mode가 더 단순하다.**
   - `bench_a3_04`, `bench_a4_02`, `bench_a4_03`는 composite room collapse, repair, solver escalation 같은 복합 failure가 섞였다.
   - 반면 `bench_a2_03`는 room kind 자체는 안정적으로 유지되며, 핵심 손실은 장애물 layout inference에 집중된다.

3. **`floorplan`은 rectangular obstacle-dense case에서 가장 유용한 입력 신호로 보인다.**
   - 다른 4개 view는 CFD score가 0.35~0.37 수준에 머물렀다.
   - `floorplan`만 `cfd_score ≈ 0.655`로 baseline `bench_a1_01/perspective`와 비슷한 수준까지 올라갔다.
   - 다만 obstacle match IoU가 0이므로, 현재 CFD metric은 obstacle topology mismatch를 충분히 강하게 벌주지 못할 가능성도 같이 시사한다.

4. **이 case는 evaluation metric refinement 후보이기도 하다.**
   - coarse structural score(0.625~0.750)와 실제 obstacle-layout fidelity 사이 괴리가 크다.
   - rectangular case에서도 opening/topology뿐 아니라 obstacle occupancy / passage blockage 민감도를 별도 보강할 필요가 있다.

## 산출물

- batch summary:
  - `benchmark/manifests/evaluation_batch_summary_a2_03.json`
- per-view evaluation summaries:
  - `benchmark/evaluations/bench_a2_03/*/evaluation_summary.json`
- predicted scenes / cases / CFD metrics:
  - `benchmark/evaluations/bench_a2_03/*/predicted_scene.json`
  - `benchmark/evaluations/bench_a2_03/*/predicted_case`
  - `benchmark/evaluations/bench_a2_03/*/cfd_metrics.json`

## 1차 해석

- **Composite hard cases**는 room topology collapse + repair/solver stress가 주요 이슈였다.
- **Rectangular hard case `bench_a2_03`**는 topology는 맞추되 dense obstacle layout을 놓치는 쪽이 주요 이슈다.
- 따라서 다음 분석 단계에서는 composite 전용 opening/topology-sensitive 지표와 함께,
  rectangular dense-obstacle case를 위한 **obstacle occupancy / blockage-sensitive metric**도 병행 검토하는 편이 낫다.

## 다음 권장 액션

1. `bench_a1_04` 같은 solver-stress 직사각형 case도 5-view sweep해서,
   `bench_a2_03`의 geometry-hard signal과 `bench_a1_04`의 solver-hard signal을 분리 비교한다.
2. 현재까지 완료된 5개 case(`bench_a1_01`, `bench_a2_03`, `bench_a3_04`, `bench_a4_02`, `bench_a4_03`)를 묶어,
   view별 평균 구조/CFD 점수와 failure mode 태그를 집계하는 aggregate note/script를 추가한다.
3. 구조 metric 보강 후보로 다음 둘을 우선 검토한다.
   - obstacle occupancy / blockage overlap
   - opening-connected flow-path / free-corridor sensitivity
