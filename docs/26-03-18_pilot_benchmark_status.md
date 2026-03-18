# Pilot Benchmark Reference Status

> Date: 2026-03-18

## Summary

- Pilot set size: 8 scenes
- Reference CFD success: 8/8
- Root cause fixed in this pass: `scene_to_gmsh.py` called Boolean subtraction with an empty obstacle list, which broke the simple 0-obstacle pilot scenes.
- Compatibility improvement: `repair_indoor_scene.py` now handles composite `room.blocks` opening normalization using actual exposed wall segments.

## Current pilot table

| Case | Cat. | Room | Obst. | Success | Preset | Mesh | 3D bundle | Notes |
|---|---|---|---:|---|---|---:|---|---|
| bench_a1_01 | A1 | rectangular | 1 | yes | robust | 0.35 | partial | pilot baseline retained |
| bench_a1_02 | A1 | rectangular | 0 | yes | robust | 0.35 | partial | 0-obstacle path recovered |
| bench_a2_01 | A2 | rectangular | 2 | yes | conservative | 0.35 | partial | pilot baseline retained |
| bench_a2_02 | A2 | rectangular | 2 | yes | robust | 0.35 | partial | pilot baseline retained |
| bench_a3_01 | A3 | composite | 0 | yes | robust | 0.35 | partial | 0-obstacle path recovered |
| bench_a3_02 | A3 | composite | 0 | yes | robust | 0.35 | partial | 0-obstacle path recovered |
| bench_a4_01 | A4 | composite | 3 | yes | robust | 0.35 | partial | pilot baseline retained |
| bench_a4_02 | A4 | composite | 4 | yes | laminar_fallback | 0.35 | partial | pilot baseline retained |

## Artifact structure now in use

- `benchmark/scenes/*.json` — frozen benchmark scene definitions
- `benchmark/manifests/scene_manifest.json` — category/seed/geometry manifest
- `benchmark/manifests/pilot_reference_status.{json,csv}` — aggregated pilot reference status
- `benchmark/reference_cfd/<case>` — symlink to canonical OpenFOAM case directory under `cases/<case>`
- `benchmark/visualizations/<case>` — symlink to canonical visualization/result directory under `results/<case>`

## Rendering feasibility status

- All 8 pilot scenes now have reusable per-case 3D rendering artifacts in `results/<case>/` and via `benchmark/visualizations/<case>/`.
- Current reusable bundle: `comparison_1x2.png`, `indoor_pipeline_3d_comparison.png/.pdf`, `panel_geometry_3d.png`, `panel_flow_3d.png`.
- Explicit benchmark input views are now also exported under `benchmark/renderings/<case>/<view>/<case>_<view>.png`.
- Canonical input-view set currently used: `perspective`, `birdseye`, `floorplan`, `wireframe`.
- `run_benchmark_reference_batch.py` now performs reference CFD execution and benchmark input-view export in one batch flow.

## Next recommended actions

1. Finish the newly launched `*_03` reference runs and inspect whether the 12-scene frozen subset remains fully stable.
2. Write a follow-up status note for the 12-scene subset once the new four runs finish.
3. If the 12-scene set is stable, take a local commit checkpoint and then expand toward the 12-20 scene target range.
