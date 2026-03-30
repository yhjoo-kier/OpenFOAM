# A3 Mid-Difficulty Composite Case CLI Evaluation Sweep (`bench_a3_03`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the mid-difficulty composite case `bench_a3_03`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **35/100 success**.

Completed full-case sweeps so far:

- `bench_a1_01`: 5/5 success
- `bench_a1_04`: 5/5 success
- `bench_a2_03`: 5/5 success
- `bench_a3_03`: 5/5 success
- `bench_a3_04`: 5/5 success
- `bench_a4_02`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Solver path | Notes |
|---|---:|---:|---|---|
| perspective | 1.0000 | 0.6317 | original + `robust` | strong geometry recovery with stable CFD |
| birdseye | 0.6250 | 0.6955 | original + `robust` | best non-floorplan CFD despite missing obstacle match |
| floorplan | 1.0000 | 0.7276 | original + `robust` | best overall view in this sweep |
| wireframe | 0.6250 | 0.5680 | original + `robust` | run succeeded directly, but repair sidecar failed noisily |
| section | 0.6250 | 0.3237 | original + `robust` | weakest CFD despite composite room recovery |

Case-level averages:

- avg structural score: **0.7750**
- avg CFD score: **0.5893**

## Key signals

1. **This is a healthier composite case than the previously swept hard composite cases.**
   - `bench_a3_03` sits clearly above `bench_a3_04` (`avg structural 0.383`, `avg CFD 0.453`) and `bench_a4_02` (`avg structural 0.738`, `avg CFD 0.334`).
   - It is therefore useful as a composite **mid-difficulty positive control**, not just another failure-heavy stress case.

2. **`floorplan` remains the strongest downstream view even for composite rooms.**
   - `floorplan` achieved perfect structural score and the best CFD score of the sweep (`0.7276`).
   - `perspective` also reconstructed the composite geometry well, but still underperformed `floorplan` on CFD.
   - This continues the running pattern that layout-preserving 2D views are surprisingly competitive for the final flow field.

3. **`section` still looks structurally plausible but CFD-fragile.**
   - `section` preserved the composite room kind, but only matched one opening wall and missed the obstacle match entirely.
   - Its CFD score dropped to `0.3237`, much lower than the other four views.
   - This reinforces the existing concern that section views can overstate structure recovery while under-delivering on flow-path fidelity.

4. **Obstacle-aware structure scoring is still under-expressive for some views.**
   - `birdseye` and `wireframe` both scored only `0.625` structurally because obstacle matching collapsed.
   - Even so, their CFD scores stayed relatively healthy (`0.6955`, `0.5680`).
   - That again suggests the current obstacle metric is stricter than the downstream CFD impact in some single-obstacle composite cases.

5. **Wireframe produced a repair-sidecar hygiene issue, but not a task failure.**
   - The original generated scene solved directly with `robust` + RAS.
   - However, the optional repair sidecar still threw an error before being discarded.
   - This is not blocking benchmark progress, but it is a small bookkeeping/regression signal worth watching if repair logging is later used as a failure proxy.

## Interpretation

`bench_a3_03` helps reduce the current hard-case bias in the evaluation evidence.
So far the benchmark had several strong signals from composite stress cases (`bench_a3_04`, `bench_a4_02`, `bench_a4_03`), but less evidence for whether moderate composite layouts are routinely recoverable.

This sweep suggests the answer is **yes, often**, at least through the Gemini CLI path:

- all 5 views solved on the original scene
- no solver escalation above `robust` was needed
- composite room type was recovered for all views
- `floorplan` and `perspective` both reached near-identity structural recovery

The remaining weakness is not raw pipeline stability but **which structural errors actually matter for CFD**.
A missed obstacle or one wrong opening wall can hurt the geometry score a lot, yet the CFD score does not always collapse proportionally.
That keeps supporting the need for more topology/opening/blockage-sensitive interpretation rather than relying on a single coarse structural aggregate.

## Suggested next steps

1. Add one rectangular mid-difficulty case next (`bench_a2_01` or `bench_a1_02`) to balance the new composite mid-case evidence.
2. Fold `bench_a3_03` into the aggregate comparison note so the benchmark is no longer dominated by only baseline + hard-case narratives.
3. If repair bookkeeping becomes analytically important, clean up the wireframe repair-sidecar false-alarm path so "repair attempted" does not read like a real failure when the original run already succeeded.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_summary_a3_03.json`
- Case root: `benchmark/evaluations/bench_a3_03/`
- Aggregate summary: `benchmark/evaluations/summary.json`
