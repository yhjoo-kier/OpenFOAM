# A1 Easy Rectangular Case CLI Evaluation Sweep (`bench_a1_02`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the easy rectangular case `bench_a1_02`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **45/100 success**.

Completed full-case sweeps so far:

- `bench_a1_01`: 5/5 success
- `bench_a1_02`: 5/5 success
- `bench_a1_04`: 5/5 success
- `bench_a2_01`: 5/5 success
- `bench_a2_03`: 5/5 success
- `bench_a3_03`: 5/5 success
- `bench_a3_04`: 5/5 success
- `bench_a4_02`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Solver path | Notes |
|---|---:|---:|---|---|
| perspective | 0.7500 | 0.6266 | original + `robust` at `0.35` | room recovered, but 3 spurious obstacles appeared |
| birdseye | 1.0000 | 0.6873 | original + `robust` at `0.35` | best CFD in this sweep; openings and empty-room structure recovered cleanly |
| floorplan | 1.0000 | 0.6860 | original + `robust` at `0.35` | nearly tied with birdseye; strongest layout fidelity |
| wireframe | 0.7500 | 0.6089 | original + `robust` at `0.35` | also added 3 spurious obstacles; mesh skewness warning remained |
| section | 0.6250 | 0.3593 | original + `robust` at `0.35` | opening wall mismatch plus hallucinated obstacles hurt CFD most |

Case-level averages:

- avg structural score: **0.8250**
- avg CFD score: **0.5936**

## Key signals

1. **`bench_a1_02` is now the cleanest easy rectangular positive-control after `bench_a1_01`.**
   - All 5 views finished directly on the original scene with `mesh_size=0.35 + robust`.
   - No repair path, fallback mesh-size reduction, or preset escalation was needed.

2. **For an empty-room reference, `birdseye` and `floorplan` recover the geometry almost perfectly.**
   - Both views reached structural score `1.0` and CFD score about `0.687`.
   - This reinforces the current interpretation that layout-centric cues remain the most reliable input for downstream flow agreement.

3. **Photoreal / edge-emphasized views still hallucinate obstacles even in an easy case.**
   - `perspective`, `wireframe`, and `section` all predicted **3 obstacles** for a reference case with **0 obstacles**.
   - The pipeline still solved successfully, but the structural penalty here is a clean example of view-induced object hallucination rather than true benchmark difficulty.

4. **`section` remains the most CFD-fragile view even when the room itself is simple.**
   - The generated room stayed rectangular, but one opening wall flipped and extra obstacles were introduced.
   - That dropped the CFD score to `0.3593`, far below the other four views.

5. **High mesh-risk warnings do not necessarily imply task-level failure on simple cases.**
   - All five successful runs still reported `high` mesh-risk scores, mainly from non-orthogonality / aspect-ratio warnings.
   - So this case is useful as a reminder that current risk flags are conservative and should be interpreted alongside actual solver outcomes.

## Interpretation

`bench_a1_02` helps separate two different phenomena that had been partially entangled in earlier sweeps:

- true benchmark difficulty
- view-dependent hallucination / ambiguity in the image-conditioned path

This case is geometrically easy and operationally stable, yet several views still invent obstacles.
That means some structural-score loss is not coming from solver fragility or composite-room complexity,
but from how the multimodal generation step interprets visually rich or partial cues.

At the same time, the CFD scores stayed relatively strong for `perspective` and `wireframe`, which again suggests
that obstacle-instance mismatch can overstate the downstream physical error when the main room topology and opening directions are preserved.

## Suggested next steps

1. Run the next easy rectangular case (`bench_a1_03`) to see whether the empty/simple-room positive-control signal stays consistent across another seed.
2. Refresh the aggregate comparison note to the **45/100** milestone so the rectangular evidence includes two easy positives plus mid/stress cases.
3. Keep `bench_a1_02` in the rectangular control subset because it isolates view-driven hallucination with minimal solver-side confounders.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_bench_a1_02.json`
- Case root: `benchmark/evaluations/bench_a1_02/`
- Aggregate summary: `benchmark/evaluations/summary.json`
