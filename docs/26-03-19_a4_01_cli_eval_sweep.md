# A4 Composite Control Case CLI Evaluation Sweep (`bench_a4_01`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the composite dense-obstacle case `bench_a4_01`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **70/100 success**.

Completed full-case sweeps so far:

- `bench_a1_01`: 5/5 success
- `bench_a1_02`: 5/5 success
- `bench_a1_03`: 5/5 success
- `bench_a1_04`: 5/5 success
- `bench_a2_01`: 5/5 success
- `bench_a2_02`: 5/5 success
- `bench_a2_03`: 5/5 success
- `bench_a3_01`: 5/5 success
- `bench_a3_02`: 5/5 success
- `bench_a3_03`: 5/5 success
- `bench_a3_04`: 5/5 success
- `bench_a4_01`: 5/5 success
- `bench_a4_02`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Solver path | Notes |
|---|---:|---:|---|---|
| perspective | 0.3750 | 0.2904 | original + `robust` at `0.35` | weakest view; opening wall mismatch remained, and a non-blocking repair sidecar error was logged |
| birdseye | 0.7500 | 0.4888 | original + `robust` at `0.35` | solid composite recovery with correct opening walls |
| floorplan | 0.3750 | 0.2935 | original + `robust` at `0.35` | unexpectedly weak for floorplan; room kind stayed composite but geometry fidelity collapsed |
| wireframe | 0.7500 | 0.6718 | original + `robust` at `0.35` | best CFD in this sweep and the cleanest dense-composite signal |
| section | 0.7083 | 0.3176 | original + `robust` at `0.35` | structure looked plausible, but opening mismatch still hurt CFD badly |

Case-level averages:

- avg structural score: **0.5917**
- avg CFD score: **0.4124**

## Key signals

1. **`bench_a4_01` is not another positive hard-case like `bench_a4_03`.**
   - All 5 views solved directly on the original scene.
   - But the case-level average stayed modest, especially because `perspective` and `floorplan` both underperformed.

2. **Dense composite rooms still split views sharply even when obstacle count is correct.**
   - Every view predicted the correct obstacle count (`3`) and preserved the composite room kind.
   - Yet CFD ranged from `0.2904` to `0.6718`, so count agreement alone says very little here.

3. **`wireframe` clearly outperformed `floorplan` on this case.**
   - `wireframe` reached the best CFD score of the sweep (`0.6718`).
   - `floorplan` dropped to `0.2935`, which is one of the clearest examples so far that top-down abstraction can fail when dense-composite geometry is laid out ambiguously.

4. **Opening-wall mismatch again behaves like a first-order failure mode.**
   - `perspective` and `section` both carried opening-wall mismatches and both had poor CFD agreement.
   - `birdseye`/`wireframe`, which kept opening walls correct, were much healthier.

5. **A non-blocking repair-sidecar error appeared on `perspective`.**
   - The original scene solved successfully, so the task still counts as a clean success.
   - But the stabilization summary recorded a failed repair attempt in the background path.
   - This mirrors earlier bookkeeping noise on successful tasks and is worth keeping as a regression signal even though it did not affect the final result.

## Interpretation

`bench_a4_01` is useful because it fills a gap between the healthier dense-composite success of `bench_a4_03`
and the more explicitly stressed `bench_a4_02` lane.
This case shows that even without repair dependence or preset escalation,
A4 scenes can still be semantically fragile at the image-to-geometry stage.

The important detail is that the weakness is not dominated by obstacle-count failure.
Instead, it looks more like a geometry-layout / opening-fidelity problem inside an otherwise valid composite-room prediction.
That makes `bench_a4_01` a good argument for keeping topology/opening-sensitive interpretation separate from raw instance matching.

## Suggested next steps

1. Keep `bench_a4_01` alongside `bench_a4_02` and `bench_a4_03` as a three-point A4 comparison set.
2. Track the non-blocking repair-sidecar failure on `perspective` as a bookkeeping/regression note rather than a benchmark failure.
3. Continue remaining coverage with pending rectangular (`bench_a1_05`, `bench_a2_04`, `bench_a2_05`) and late composite (`bench_a3_05`, `bench_a4_04`, `bench_a4_05`) cases.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_bench_a3_02_a4_01.json`
- Case root: `benchmark/evaluations/bench_a4_01/`
- Aggregate summary: `benchmark/evaluations/summary.json`
