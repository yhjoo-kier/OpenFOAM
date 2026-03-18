# Frozen 12-Scene Benchmark Reference Status

> Date: 2026-03-18

## Summary

- Frozen benchmark subset size: 12 scenes (`A1/A2/A3/A4 × 3`)
- Reference CFD success: **12/12**
- Benchmark input-view export: **12/12** (`perspective`, `birdseye`, `floorplan`, `wireframe`, `section`)
- 3D comparison rendering bundle: **12/12** (`comparison_1x2.png`, `indoor_pipeline_3d_comparison.png`)

This pass extended the previously recovered 8-scene pilot set to the full 12-scene frozen subset and recovered the new `*_03` cases.

## What changed in this pass

1. `run_indoor_stabilized.py` was hardened for reruns:
   - remove any pre-existing case directory before recreating a case
   - remove stale `constant/polyMesh` before `gmshToFoam`
   - add bounded timeouts for mesh import / `checkMesh`
2. `run_benchmark_reference_batch.py` gained non-destructive bookkeeping modes:
   - `--collect-only`
   - `--skip-existing-success`
3. `render_benchmark_views.py` was extended with a canonical `section` renderer aligned to the opening pair when possible.
4. Batch manifests were refreshed from the live artifact state.

## New `*_03` subset outcome

| Case | Cat. | Room | Obst. | Success | Preset | Mesh | Notes |
|---|---|---|---:|---|---|---:|---|
| `bench_a1_03` | A1 | rectangular | 1 | yes | robust | 0.35 | completed directly |
| `bench_a2_03` | A2 | rectangular | 3 | yes | robust | 0.25 | initial `0.35` path timed out around `checkMesh`, recovered at `0.25` |
| `bench_a3_03` | A3 | composite | 1 | yes | robust | 0.35 | completed directly |
| `bench_a4_03` | A4 | composite | 3 | yes | robust | 0.25 | initial `0.35` path timed out around `checkMesh`, recovered at `0.25` |

## Interpretation

- The frozen 12-scene subset is now operational enough to serve as the first benchmark slice for the paper workflow.
- Composite-room scenes are no longer blocked at the visualization / 3D rendering layer.
- The benchmark input bundle now covers the paper-plan 5-view matrix (`perspective`, `birdseye`, `floorplan`, `wireframe`, `section`) for all 12 scenes.
- The remaining instability is not category-level failure; it is mesh-quality / runtime-robustness sensitivity for some harder cases.
- `bench_a2_03` and `bench_a4_03` still carry elevated mesh-risk signals (`defaultFaces`, high non-orthogonality, high aspect ratio / skewness), but they are now solvable through the existing stabilized fallback path.

## Artifact / manifest state

Updated manifests:
- `benchmark/manifests/reference_batch_summary.json`
- `benchmark/manifests/frozen12_reference_status.json`
- `benchmark/manifests/frozen12_reference_status.csv`

Canonical artifact roots:
- `benchmark/reference_cfd/<case>`
- `benchmark/visualizations/<case>`
- `benchmark/renderings/<case>/<view>/...`

## Recommended next actions

1. Take a local commit checkpoint for the 12/12 frozen subset recovery.
2. Expand the frozen generator output toward the next tranche (e.g. 16 or 20 scenes).
3. Start `benchmark/evaluations/` scaffolding so image-conditioned pipeline outputs can be compared directly against the new reference set.
4. Revisit `repair_indoor_scene.py` to replace the remaining composite fallback/skip behavior with native composite repair logic.
