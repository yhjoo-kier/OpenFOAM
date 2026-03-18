# A4 Laminar-Fallback Composite Case CLI Evaluation Sweep (`bench_a4_02`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the laminar-fallback composite stress case `bench_a4_02`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **20/100 success**.

Completed full-case sweeps so far:

- `bench_a1_01`: 5/5 success
- `bench_a3_04`: 5/5 success
- `bench_a4_02`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Solver path | Notes |
|---|---:|---:|---|---|
| perspective | 0.6250 | 0.2875 | original + `robust` | weakest overall view in this sweep |
| birdseye | 0.7500 | 0.3391 | original + `robust` | stable but mesh-risk-high |
| floorplan | 0.8125 | 0.4023 | original + `robust` | best structure and best CFD in this sweep |
| wireframe | 0.7500 | 0.2963 | repaired + `ultra_robust` | original scene invalid, repair salvaged run |
| section | 0.7500 | 0.3447 | original + `robust` | moderate geometry / CFD agreement |

## Key signals

1. **The laminar-fallback reference family is still evaluable through the CLI path.**
   - The reference benchmark case `bench_a4_02` had been notable as a laminar-fallback success during dataset construction.
   - In image-conditioned evaluation, however, 4 of 5 views solved directly with `robust` + RAS at mesh size `0.35`.
   - This means the old reference-side fragility does not automatically transfer to the predicted-scene evaluation path.

2. **`floorplan` again gives the strongest downstream CFD agreement for a composite case.**
   - This now matches the same qualitative pattern seen in `bench_a3_04` and `bench_a4_03`.
   - The current benchmark appears to reward preserving large-scale room/opening layout more than fine perspective realism.

3. **`wireframe` exposed a different salvage path from `bench_a4_03`.**
   - Here the initial generation already violated composite-room containment (`obstacles[2] must lie fully inside one room block of the composite room`).
   - Meshing failed at all original-scene mesh sizes (`0.35`, `0.25`, `0.18`).
   - The repaired scene then meshed successfully, but `robust` still hit an FPE and required `ultra_robust` to finish.
   - So this case combines both geometry repair pressure and solver-stability pressure in a single evaluation task.

4. **High mesh-risk can still end in clean task success.**
   - Every successful view remained in the `high` mesh-risk bucket.
   - Even so, the current stabilization ladder recovered all five views without manual intervention.

## Interpretation

Compared with `bench_a4_03`, this case looks less faithful physically:

- structure scores are lower overall (`0.625` to `0.8125` vs. `0.7083` to `0.9167` for `bench_a4_03`)
- CFD scores are also lower overall (`0.2875` to `0.4023` vs. `0.4858` to `0.5835` for `bench_a4_03`)

This supports the idea that not all dense composite scenes are equally recoverable, and that `bench_a4_02` remains a harder stress case despite full task-level success.

The wireframe failure mode is especially useful for the paper because it shows a realistic chain:

**image ambiguity → invalid composite geometry → repair → turbulence instability → ultra-robust recovery**

That is a stronger benchmark stress signal than a simple one-shot fail/succeed label.

## Suggested next steps

1. Add a lightweight composite-topology / opening-layout analysis pass across `bench_a3_04`, `bench_a4_02`, and `bench_a4_03`.
2. Extend the sweep set to one rectangular hard case (`bench_a2_03` or `bench_a1_04`) to compare composite stress against known mesh-sensitive rectangular scenes.
3. Cut a local commit once the docs/TODO snapshot is refreshed around the 20/100 evaluation milestone.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_summary_a4_02.json`
- Case root: `benchmark/evaluations/bench_a4_02/`
- Aggregate summary: `benchmark/evaluations/summary.json`
