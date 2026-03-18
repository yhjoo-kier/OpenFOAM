# A2 Rectangular Multi-Obstacle Control CLI Evaluation Sweep (`bench_a2_02`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the rectangular multi-obstacle control case `bench_a2_02`.
All five views succeeded end-to-end, lifting the frozen-20 image-conditioned evaluation scaffold to **60/100 success**.

Completed full-case sweeps so far:

- `bench_a1_01`: 5/5 success
- `bench_a1_02`: 5/5 success
- `bench_a1_03`: 5/5 success
- `bench_a1_04`: 5/5 success
- `bench_a2_01`: 5/5 success
- `bench_a2_02`: 5/5 success
- `bench_a2_03`: 5/5 success
- `bench_a3_01`: 5/5 success
- `bench_a3_03`: 5/5 success
- `bench_a3_04`: 5/5 success
- `bench_a4_02`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Solver path | Notes |
|---|---:|---:|---|---|
| perspective | 0.7500 | 0.3002 | original + `robust` | room/obstacle count reasonable, but one opening wall mismatched |
| birdseye | 0.6250 | 0.3358 | original + `robust` | extra obstacle hallucinated and opening wall mismatch remained |
| floorplan | 0.7500 | 0.5153 | original + `robust` | only view with full opening-wall match |
| wireframe | 0.6250 | 0.3170 | original + `robust` | extra obstacle + opening wall mismatch |
| section | 0.5000 | 0.3158 | original + `robust` | both opening walls wrong, weakest structural view |

Case-level averages:

- avg structural score: **0.6500**
- avg CFD score: **0.3568**

## Key signals

1. **This is a clean success-only case, but downstream quality is modest.**
   - All 5 views solved directly on `original + robust` with `mesh_size=0.35`.
   - No repair, fallback mesh, or solver escalation was needed.
   - Yet the average CFD score stayed relatively low (`0.3568`).

2. **Opening-wall fidelity looks like the main bottleneck here.**
   - `floorplan` was the only view with perfect opening-wall match.
   - It also had the best CFD score of the sweep (`0.5153`).
   - The other 4 views each missed at least one opening wall, and all stayed near `0.30~0.34` CFD.

3. **This case is a good rectangular counterweight to empty-room hallucination cases.**
   - The reference scene contains 2 obstacles.
   - `perspective` and `floorplan` recovered 2 obstacles, while `birdseye`, `wireframe`, and `section` over-predicted 3.
   - Unlike `bench_a3_01`, the extra-obstacle effect here is mixed with opening-wall mismatch, making it a better control for real blockage-sensitive layouts.

4. **`section` underperformed mostly through topology/BC misunderstanding, not solver instability.**
   - `section` had the lowest structural score (`0.5000`).
   - Its `wall_match_ratio` dropped to `0.0`.
   - But it still solved directly, reinforcing that some benchmark weaknesses are now interpretive/model-quality issues rather than pipeline-break issues.

5. **Rectangular multi-obstacle cases may need a more explicit opening-aware diagnostic.**
   - Current structural scores distinguish the views somewhat, but the CFD spread suggests that opening-wall correctness matters more than raw obstacle count in this case.
   - A blockage-only metric would miss that; an opening/topology tag would explain it much better.

## Interpretation

`bench_a2_02` is useful as a rectangular control that is neither trivial nor mesh-fragile.
It shows that the pipeline can remain operational while the actual physical fidelity remains mediocre.
That makes it a good “quality-not-stability” case for later prompt tuning or metric refinement.

In particular, this case argues that rectangular multi-obstacle evaluation should not be treated as only an occupancy problem.
The openings still dominate the final flow behavior.
So a future refined score probably needs separate visibility for:

- opening-wall correctness
- obstacle-count / blockage correctness
- room-size / topology correctness

instead of letting one structural average blur the diagnosis.

## Suggested next steps

1. Compare this case directly against `bench_a2_03` to separate moderate vs dense rectangular obstacle stress.
2. Add `bench_a2_02` to the control subset for prompt changes that target opening placement.
3. Continue coverage expansion with an unswept composite case (`bench_a3_02` or `bench_a4_01`) so the 60/100 milestone does not become too rectangular-heavy.

## Artifact paths

- Batch summary: `benchmark/manifests/eval_batch_bench_a2_02.json`
- Case root: `benchmark/evaluations/bench_a2_02/`
- Aggregate summary: `benchmark/evaluations/summary.json`
