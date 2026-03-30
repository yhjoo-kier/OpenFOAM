# A3 Empty Composite Case CLI Evaluation Sweep (`bench_a3_02`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the empty composite case `bench_a3_02`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **65/100 success**.

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
- `bench_a4_02`: 5/5 success
- `bench_a4_03`: 5/5 success

## Per-view results

| View | Structural score | CFD score | Solver path | Notes |
|---|---:|---:|---|---|
| perspective | 1.0000 | 0.4913 | original + `robust` at `0.35` | cleanest recovery; composite room + zero obstacles + opening walls all matched |
| birdseye | 0.7500 | 0.5538 | original + `robust` at `0.35` | obstacle hallucination appeared, but topology/openings stayed correct |
| floorplan | 0.7500 | 0.4703 | original + `ultra_robust` at `0.35` | only view needing preset escalation; still solved on original scene |
| wireframe | 0.7500 | 0.6606 | original + `robust` at `0.35` | best CFD in this sweep despite 3 hallucinated obstacles |
| section | 0.6250 | 0.6248 | original + `robust` at `0.35` | opening wall mismatch appeared, but CFD stayed surprisingly strong |

Case-level averages:

- avg structural score: **0.7750**
- avg CFD score: **0.5601**

## Key signals

1. **`bench_a3_02` is a second strong empty-composite control next to `bench_a3_01`, but it behaves differently.**
   - All 5 views solved directly on the original scene with `mesh_size=0.35`.
   - Unlike `bench_a3_01`, `perspective` here recovered the empty composite almost perfectly.

2. **Obstacle hallucination again does not map cleanly to CFD degradation.**
   - `birdseye`, `floorplan`, `wireframe`, and `section` all hallucinated 3 obstacles for a reference case with 0 obstacles.
   - Even so, `wireframe` and `section` still produced CFD scores around `0.66` and `0.62`.

3. **`floorplan` is not always the strongest composite view.**
   - This time `floorplan` needed `ultra_robust` and delivered the lowest CFD score of the sweep.
   - `wireframe` was the best CFD view instead, which is useful counter-evidence against a too-simple “floorplan always wins” story.

4. **Opening-wall mismatch is again more consequential than obstacle-count fidelity.**
   - `section` dropped structurally because an opening wall flipped.
   - But even there the CFD penalty stayed milder than the structural penalty, reinforcing the need to separate error families rather than collapsing them into one score.

5. **The empty-composite lane now has two complementary controls.**
   - `bench_a3_01`: strong CFD despite broad obstacle hallucination across all views.
   - `bench_a3_02`: one nearly perfect recovery (`perspective`) plus one solver-escalated floorplan.
   - Together they make the empty-composite story much less anecdotal.

## Interpretation

`bench_a3_02` strengthens the idea that empty composite rooms are a special regime in this benchmark.
The model often preserves the global L-shaped topology and the opening pair correctly,
while still inventing obstacles that do not always ruin the downstream flow field.

That means composite evaluation should not treat all structure errors as interchangeable.
At minimum, the benchmark should continue distinguishing:

- room-topology correctness
- opening-wall correctness
- hallucinated-obstacle burden
- solver stress / preset escalation

## Suggested next steps

1. Keep `bench_a3_02` in the empty-composite control subset together with `bench_a3_01`.
2. Use the pair to test any future hallucination-aware metric refinement.
3. Continue the remaining coverage expansion with pending rectangular/A4 cases.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_bench_a3_02_a4_01.json`
- Case root: `benchmark/evaluations/bench_a3_02/`
- Aggregate summary: `benchmark/evaluations/summary.json`
