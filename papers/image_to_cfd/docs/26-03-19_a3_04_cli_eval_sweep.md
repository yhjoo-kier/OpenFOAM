# A3 Composite Case CLI Evaluation Sweep (`bench_a3_04`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the composite benchmark case `bench_a3_04`.
All five image-conditioned tasks succeeded end-to-end, so the frozen-20 evaluation scaffold is now at **10/100 success** overall:

- `bench_a1_01` (rectangular baseline): 5/5 success
- `bench_a3_04` (composite representative): 5/5 success

This gives the benchmark its first side-by-side rectangular vs composite complete-case comparison under the same CLI backend path.

## `bench_a3_04` per-view results

| View | Structural score | CFD score | Solver stabilization note |
|---|---:|---:|---|
| perspective | 0.3750 | 0.4156 | `robust`, mesh 0.35 |
| birdseye | 0.3750 | 0.4169 | `robust`, mesh 0.35 |
| floorplan | 0.3750 | 0.6840 | `robust`, mesh 0.35 |
| wireframe | 0.3750 | 0.4123 | `robust`, mesh 0.35 |
| section | 0.4167 | 0.3352 | needed `ultra_robust` after `robust` FPE |

## Comparison against the rectangular baseline (`bench_a1_01`)

| View | `bench_a1_01` structural | `bench_a1_01` CFD | `bench_a3_04` structural | `bench_a3_04` CFD |
|---|---:|---:|---:|---:|
| perspective | 0.7500 | 0.6470 | 0.3750 | 0.4156 |
| birdseye | 0.7500 | 0.6061 | 0.3750 | 0.4169 |
| floorplan | 0.8750 | 0.3930 | 0.3750 | 0.6840 |
| wireframe | 0.7500 | 0.5997 | 0.3750 | 0.4123 |
| section | 1.0000 | 0.4427 | 0.4167 | 0.3352 |

## Interpretation

1. **Composite room inference is materially harder than rectangular inference.**
   - The rectangular baseline preserved room/opening structure much better across all views.
   - In `bench_a3_04`, the model often recognized that the room was composite, but the second block was poorly sized or misplaced.

2. **Floorplan is the most interesting composite-case outlier.**
   - `bench_a3_04/floorplan` kept a low structural score (`0.375`) but achieved the best CFD score in the composite sweep (`0.6840`).
   - The predicted room remained composite with the correct block count and correct opening type count, even though one room block missed the IoU threshold and opening walls were wrong.
   - This suggests the current structural score is still too harsh or too geometry-localized to explain downstream CFD similarity by itself.

3. **Section remains fragile for composite scenes.**
   - `bench_a3_04/section` collapsed the room from composite to rectangular, over-predicted obstacles (`1 -> 3`), and swapped opening walls from `north/south` to `east/west`.
   - It also needed a solver preset escalation from `robust` to `ultra_robust` after a floating-point failure.
   - Section views appear especially risky when the composite layout is not visually obvious from a single cut.

4. **Wireframe / perspective / birdseye behaved similarly here.**
   - All three landed near `structural_score ≈ 0.375` and `cfd_score ≈ 0.41`.
   - For this case, they did not separate as clearly as they did in the rectangular baseline.

## Practical next steps

1. Add or refine a metric that better captures **opening-orientation correctness** and **composite-room topology fidelity**.
   - The current aggregate structural score under-explains why `bench_a3_04/floorplan` produced decent CFD agreement.
2. Run the next full CLI sweep on a harder composite stress case, preferably `bench_a4_03` or `bench_a4_02`.
   - This will show whether the composite-vs-rectangular gap grows with obstacle density.
3. Consider an analysis table for the paper with both:
   - structure fidelity metrics
   - CFD similarity metrics
   - solver stabilization cost (`preset`, `attempt_count`, fallback use)

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_summary_a3_04.json`
- Case root: `benchmark/evaluations/bench_a3_04/`
- Aggregate summary: `benchmark/evaluations/summary.json`
