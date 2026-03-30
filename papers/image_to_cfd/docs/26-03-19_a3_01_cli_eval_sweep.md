# A3 Empty-Composite Pilot Case CLI Evaluation Sweep (`bench_a3_01`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the empty-composite pilot case `bench_a3_01`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **55/100 success** at the time this case finished.

Completed full-case sweeps at this checkpoint:

- `bench_a1_01`: 5/5 success
- `bench_a1_02`: 5/5 success
- `bench_a1_03`: 5/5 success
- `bench_a1_04`: 5/5 success
- `bench_a2_01`: 5/5 success
- `bench_a2_03`: 5/5 success
- `bench_a3_01`: 5/5 success
- `bench_a3_03`: 5/5 success
- `bench_a3_04`: 5/5 success
- `bench_a4_02`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Solver path | Notes |
|---|---:|---:|---|---|
| perspective | 0.6250 | 0.3795 | original + `ultra_robust` | only view that needed escalation above `robust` |
| birdseye | 0.6250 | 0.7315 | original + `robust` | best CFD despite only partial block-level match |
| floorplan | 0.7500 | 0.6965 | original + `robust` | strongest structural recovery of the sweep |
| wireframe | 0.7500 | 0.5482 | original + `robust` | good topology recovery with lower downstream CFD |
| section | 0.6250 | 0.5855 | original + `robust` | outlet wall mismatch, but CFD still stayed healthy |

Case-level averages:

- avg structural score: **0.6750**
- avg CFD score: **0.5882**

## Key signals

1. **This case is empty in the reference geometry, but Gemini hallucinated obstacles in every view.**
   - The reference scene has `0` obstacles.
   - All 5 predicted scenes produced `3` obstacles.
   - Even so, the case still solved successfully in every view and several views retained strong CFD similarity.

2. **Obstacle-count mismatch is not the whole story for empty composite rooms.**
   - `birdseye` reached the best CFD score of the sweep (`0.7315`) while still hallucinating 3 obstacles.
   - `floorplan` also stayed strong (`0.6965`).
   - This suggests that for some empty composite layouts, opening placement and gross topology matter more than obstacle-count correctness.

3. **Composite room kind was preserved in all 5 views, but block fidelity varied sharply by view.**
   - All views correctly selected a composite room (`room_kind_match = True`).
   - `floorplan`, `wireframe`, and `section` achieved full room-block F1 (`1.0`).
   - `perspective` and `birdseye` only matched one block confidently (`room_block_match.f1 = 0.5`).

4. **`section` remained opening-sensitive, but less catastrophically than in some other composite cases.**
   - `section` had an outlet wall mismatch (`wall_match_ratio = 0.5`).
   - Even so, CFD only dropped to `0.5855`, which is still materially healthier than the weak-section behavior seen in harder composite cases.
   - That makes `bench_a3_01` a useful counterexample to the simple rule that section views are always poor for composite recovery.

5. **`perspective` looks like the mesh/stability outlier, not the geometry winner.**
   - It was the only view that escalated to `ultra_robust`.
   - It also produced the lowest CFD score (`0.3795`).
   - The case therefore adds a lightweight composite stress signal without becoming a full failure case.

## Interpretation

`bench_a3_01` is useful because it broadens the composite evidence beyond only mid-difficulty and stress-heavy examples.
It shows a slightly awkward regime: the model recognizes that the room is composite, places the openings mostly sensibly, but still invents flow obstacles that do not exist.

That mismatch hurts structural fidelity, but not always CFD fidelity.
So for empty or nearly empty composite rooms, the current aggregate structural score likely over-penalizes obstacle hallucination relative to actual flow-path preservation.

This case therefore strengthens the argument for tracking at least one auxiliary metric that separates:

- topology / room-kind correctness
- opening-wall correctness
- hallucinated-obstacle burden

rather than collapsing them too early into one structural aggregate.

## Suggested next steps

1. Pair this case with a rectangular multi-obstacle control (`bench_a2_02`) so empty-composite hallucination and real-obstacle recovery can be compared directly.
2. Add a lightweight hallucinated-obstacle tag to aggregate milestone notes.
3. Keep `bench_a3_01` in the control subset for future prompt/metric changes, because it is informative without being solver-fragile.

## Artifact paths

- Batch summary: `benchmark/manifests/eval_batch_bench_a3_01.json`
- Case root: `benchmark/evaluations/bench_a3_01/`
- Aggregate summary: `benchmark/evaluations/summary.json`
