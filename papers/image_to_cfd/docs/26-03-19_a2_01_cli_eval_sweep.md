# A2 Mid-Difficulty Rectangular Case CLI Evaluation Sweep (`bench_a2_01`)

> Date: 2026-03-19

## Summary

A full 5-view Gemini CLI evaluation sweep was completed for the rectangular mid-difficulty case `bench_a2_01`.
All five views succeeded end-to-end, raising the frozen-20 image-conditioned evaluation scaffold to **40/100 success**.

Completed full-case sweeps so far:

- `bench_a1_01`: 5/5 success
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
| perspective | 0.7500 | 0.4543 | original + `robust` at `0.25` after `0.35` import/checkMesh timeout | longest path in this sweep; obstacle/layout fidelity weakest |
| birdseye | 0.9500 | 0.5083 | original + `robust` | strongest structure recovery; opening walls correct |
| floorplan | 0.8500 | 0.5937 | original + `robust` | best overall CFD in this sweep |
| wireframe | 0.7500 | 0.3908 | original + `robust` | full success, but worst CFD among the five views |
| section | 0.7250 | 0.4216 | original + `robust` | opening wall mismatch plus inflated room height |

Case-level averages:

- avg structural score: **0.8050**
- avg CFD score: **0.4737**

## Key signals

1. **`bench_a2_01` is a useful rectangular mid-difficulty control, not a pure stress case.**
   - All 5 views succeeded without repair or solver-preset escalation beyond `robust`.
   - Unlike `bench_a2_03`, the main issue here was not outright pipeline fragility but view-dependent geometry fidelity.

2. **`floorplan` again gives the best downstream CFD signal.**
   - `floorplan` reached the best CFD score of the sweep (`0.5937`) despite imperfect obstacle recovery.
   - This reinforces the running pattern that layout-preserving views outperform more visually rich views for final flow agreement.

3. **Obstacle-count overprediction is now a repeated rectangular pattern.**
   - Every view predicted **3 obstacles** for a reference case with **2 obstacles**.
   - `birdseye` still matched the two true obstacles reasonably well, but `perspective` and `wireframe` collapsed to zero obstacle precision/recall.
   - This strengthens the idea that dense/medium obstacle cases need a blockage- or occupancy-sensitive interpretation layer beyond raw object matching.

4. **`perspective` remained the only real stabilization warning in the sweep.**
   - The first `mesh_size=0.35` attempt timed out during import/checkMesh.
   - The retry at `0.25` solved successfully on the original scene, so this was not a repair-path failure.
   - That makes `bench_a2_01/perspective` a good lightweight regression case for mesh-import sensitivity even though the task finished successfully.

5. **`section` is still CFD-fragile for reasons that structural score only partly captures.**
   - The generated room height was inflated (`Lz=4.0` vs reference `2.733`) and one opening wall flipped (`north -> east`).
   - The resulting CFD score fell to `0.4216`, with even the mean velocity direction cosine turning negative.
   - This keeps the earlier interpretation intact: section images often look structurally serviceable while still misrouting flow.

## Interpretation

`bench_a2_01` reduces the earlier hard-case bias on the rectangular side.
Compared with `bench_a2_03`, this case is easier and cleaner operationally, yet it still exposes a meaningful interpretation gap:
CFD can remain moderate even when obstacle-level matching is weak or oversegmented.

That matters because the benchmark evidence is no longer just “easy baseline vs hard failures.”
We now have a rectangular mid-difficulty case showing that:

- the CLI path can stably finish 5/5 without repair
- mesh sensitivity can still appear locally (`perspective` at `0.35`)
- the strongest view remains `floorplan`
- obstacle-instance metrics still look harsher than the downstream CFD penalty in several views

## Suggested next steps

1. Add one easier rectangular case next (`bench_a1_02` or `bench_a1_03`) to separate “mid-difficulty obstacle ambiguity” from “baseline geometry recovery.”
2. Refresh the aggregate comparison note to the **40/100** milestone so rectangular evidence is no longer centered only on baseline + stress.
3. Consider whether the current structural aggregate should expose a separate occupancy/blockage component for rectangular multi-obstacle rooms.

## Artifact paths

- Batch summary: `benchmark/manifests/evaluation_batch_summary_a2_01.json`
- Case root: `benchmark/evaluations/bench_a2_01/`
- Aggregate summary: `benchmark/evaluations/summary.json`
