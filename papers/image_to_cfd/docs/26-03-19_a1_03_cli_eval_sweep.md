# A1 Easy Rectangular Case CLI Evaluation Sweep (`bench_a1_03`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the easy rectangular case `bench_a1_03`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **50/100 success**.

Completed full-case sweeps so far:

- `bench_a1_01`: 5/5 success
- `bench_a1_02`: 5/5 success
- `bench_a1_03`: 5/5 success
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
| perspective | 0.7500 | 0.6686 | original + `robust` at `0.35` | opening walls preserved; 3 spurious obstacles remained |
| birdseye | 0.7500 | 0.7356 | original + `robust` at `0.35` | best CFD in this sweep despite extra obstacles |
| floorplan | 0.7500 | 0.6607 | original + `robust` at `0.35` | 1 obstacle hallucinated; still strong layout fidelity |
| wireframe | 0.7500 | 0.7064 | original + `robust` at `0.35` | high CFD despite 3 hallucinated obstacles and skewness warning |
| section | 0.7500 | 0.7092 | original + `robust` at `0.35` | unlike prior section cases, opening walls stayed correct and CFD remained strong |

Case-level averages:

- avg structural score: **0.7500**
- avg CFD score: **0.6961**

## Key signals

1. **`bench_a1_03` became the strongest easy-case CFD positive-control so far.**
   - All 5 views finished directly on the original scene with `mesh_size=0.35 + robust`.
   - The average CFD score (`0.6961`) is currently higher than every previously completed case.

2. **This case cleanly exposes a structural-metric vs CFD gap.**
   - Every view only scored `0.75` structurally, mainly because hallucinated obstacles were introduced.
   - Yet every view still produced strong CFD agreement (`0.66`–`0.74`).
   - That is a direct signal that current object-level structural scoring can be harsher than the real flow penalty in easy open-room geometries.

3. **`section` is not universally weak; its failures seem topology/opening dependent.**
   - Here `section` reached CFD `0.7092`, one of the best scores in the whole sweep.
   - That contrasts sharply with earlier section cases where opening-wall mismatch dominated.
   - So the section-view weakness is better read as a conditional failure mode, not a universal rule.

4. **Bird's-eye cue is extremely competitive on easy rectangular rooms.**
   - `birdseye` produced the best CFD score of the sweep (`0.7356`).
   - Together with `bench_a1_02`, this suggests that easy rectangular scenes may not require strict floorplan abstraction to obtain strong downstream flow agreement.

5. **High mesh-risk warnings still coexist with fully healthy outcomes.**
   - Several views kept `high` risk flags from non-orthogonality / skewness thresholds.
   - Nevertheless, all five runs solved directly without repair or fallback.
   - This further supports separating conservative mesh warnings from true operational instability.

## Interpretation

`bench_a1_03` sharpens the main metric question of the current benchmark pipeline.
The multimodal generator still hallucinates obstacles, so structure scores stay capped,
but the resulting predicted rooms preserve enough global topology and opening placement that the CFD fields remain close to reference.

That means easy rectangular rooms are now giving two complementary signals:

- **`bench_a1_02`**: some views are nearly perfect structurally, while `section` still breaks CFD when openings drift
- **`bench_a1_03`**: structure stays imperfect across all views, but CFD remains strong anyway

Together, they argue for adding a more flow-relevant interpretation layer rather than relying only on object-instance matching.

## Suggested next steps

1. Refresh the aggregate comparison note to the **50/100** milestone, explicitly calling out the structural-vs-CFD gap now visible in `bench_a1_03`.
2. Keep `bench_a1_03` as an easy positive-control when testing future metric refinements, especially any occupancy/blockage-aware scoring.
3. Move the next batch toward still-pending rectangular or composite cases depending on whether the priority is coverage expansion or metric analysis.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_bench_a1_03.json`
- Case root: `benchmark/evaluations/bench_a1_03/`
- Aggregate summary: `benchmark/evaluations/summary.json`
