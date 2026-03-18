# A4 Composite Stress Case CLI Evaluation Sweep (`bench_a4_03`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the obstacle-dense composite stress case `bench_a4_03`.
After the sweep and one immediate rerun of the initially failed perspective task, all five views now succeed end-to-end.

This raises the frozen-20 image-conditioned evaluation scaffold to **15/100 success**:

- `bench_a1_01`: 5/5 success
- `bench_a3_04`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Notes |
|---|---:|---:|---|
| perspective | 0.7083 | 0.5287 | initial attempt exposed invalid composite-origin scene sample; rerun succeeded |
| birdseye | 0.7083 | 0.5280 | stable on `robust` |
| floorplan | 0.7500 | 0.5835 | best CFD score in this sweep |
| wireframe | 0.8929 | 0.4858 | strongest structure fidelity among non-section views |
| section | 0.9167 | 0.5014 | best structural score overall |

All successful runs finished on mesh size `0.35` with the `robust` preset; no solver preset escalation or mesh fallback was needed once the generation-stage issue was bypassed.

## Key signal

Unlike the earlier `bench_a3_04` composite representative case, `bench_a4_03` shows a much healthier relation between structure fidelity and downstream CFD agreement:

- room kind remained composite in all 5 views
- opening-wall match stayed correct in all 5 views
- structure scores stayed relatively high (`0.7083` to `0.9167`)
- CFD scores were moderate but consistent (`0.4858` to `0.5835`)

This suggests the CLI path can handle at least some dense composite scenes more reliably than the first A3 representative case implied.

## Generation/repair pipeline signal

The first `perspective` batch attempt failed **before meshing/solver**, during generation validation, because the sampled scene JSON used a composite second block with negative `y` origin and invalid downstream geometry placement.

To harden this failure class:

- `scripts/repair_indoor_scene.py` was updated to normalize composite room origins back to the `room_min_corner` convention before applying the usual repair steps.
- `scripts/run_indoor_stabilized.py` was updated so a nonzero generation exit does not immediately abort when a scene JSON was still written; this preserves a salvage path for future invalid-but-repairable generations.

A direct repair smoke test on the previously invalid `eval_bench_a4_03_perspective.json` then became valid after an origin shift of `y += 2.0`.

## Interpretation

1. **Composite difficulty is not uniform.**
   - `bench_a3_04` remained structurally fragile.
   - `bench_a4_03` is harder in obstacle density, yet more recoverable across views.

2. **Floorplan is still the strongest CFD view in both recent composite sweeps.**
   - This pattern now appears in both `bench_a3_04` and `bench_a4_03`.
   - It is becoming plausible that floorplan preserves layout/opening information that matters more to the current CFD metric than some richer-looking perspective views.

3. **Section is no longer universally bad for composite scenes.**
   - It was the weakest path in `bench_a3_04`.
   - It became the strongest structure view in `bench_a4_03`.
   - Section-view reliability likely depends heavily on whether the cut actually exposes the decisive composite corner/leg geometry.

## Suggested next steps

1. Compare `bench_a3_04` vs `bench_a4_03` more explicitly to isolate what made A3 structurally brittle.
2. Run another stress case, likely `bench_a4_02`, to see whether the current hardening also helps the laminar-fallback family.
3. Add a lightweight analysis metric for **room-topology preservation** so the paper can discuss composite-room inference more precisely than with the current aggregate structural score alone.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_summary_a4_03.json`
- Case root: `benchmark/evaluations/bench_a4_03/`
- Aggregate summary: `benchmark/evaluations/summary.json`
