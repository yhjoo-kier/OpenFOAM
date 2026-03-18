# Benchmark Dataset Card

> Date: 2026-03-18

## Summary

- Frozen benchmark bundle: **20 scenes**
- Reference CFD success: **20/20**
- Render bundles: **20**
- Evaluation scaffold: **100 tasks** across **20 cases**
- Stress subset: **5 cases**

## Bundle composition

### Categories

- A1: 5
- A2: 5
- A3: 5
- A4: 5

### Room kinds

- composite: 10
- rectangular: 10

### Obstacle histogram

- 0 obstacles: 3
- 1 obstacles: 7
- 2 obstacles: 3
- 3 obstacles: 5
- 4 obstacles: 2

### View coverage

- Views: perspective, birdseye, floorplan, wireframe, section
- Section axis coverage: {'x': 11, 'y': 9}

Per-category view counts:

- A1: {'birdseye': 5, 'floorplan': 5, 'perspective': 5, 'section': 5, 'wireframe': 5}
- A2: {'birdseye': 5, 'floorplan': 5, 'perspective': 5, 'section': 5, 'wireframe': 5}
- A3: {'birdseye': 5, 'floorplan': 5, 'perspective': 5, 'section': 5, 'wireframe': 5}
- A4: {'birdseye': 5, 'floorplan': 5, 'perspective': 5, 'section': 5, 'wireframe': 5}

## Stabilization profile

- Successful preset counts: {'conservative': 1, 'laminar_fallback': 1, 'robust': 17, 'ultra_robust': 1}
- Successful mode counts: {'RAS': 19, 'laminar': 1}
- Successful mesh size counts: {'0.25': 2, '0.35': 18}
- Stress cases: bench_a1_04, bench_a2_01, bench_a2_03, bench_a4_02, bench_a4_03

## Evaluation profile

- Status counts: {'blocked': 1, 'pending': 99}
- Note: Blocked tasks in the current cron shell are backend-environment blocks, not benchmark reference failures.

## Integrity / release notes

- Integrity summary: {'scene_count': 20, 'reference_case_count': 20, 'reference_success_count': 20, 'render_bundle_count': 20, 'evaluation_case_count': 20, 'evaluation_task_count': 100, 'expected_evaluation_task_count': 100, 'hard_issue_count': 0, 'soft_stress_signal_count': 11}
- Soft stress signal counts: {'stress_case_multi_attempt': 5, 'stress_case_preset_fallback': 3, 'stress_case_mesh_fallback': 2, 'stress_case_mode_fallback': 1}
- Missing reference cases: []
- Missing render cases: []
- Missing evaluation cases: []

## Artifact paths

- `scene_manifest`: `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/scene_manifest.json`
- `reference_summary`: `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/reference_batch_summary.json`
- `renderings_manifest`: `/home/yhjoo/projects/OpenFOAM/benchmark/renderings/renderings_manifest.json`
- `evaluation_summary`: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations/summary.json`
- `integrity_summary`: `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/dataset_integrity_summary.json`

## Notes

- This dataset card is generated from the benchmark manifests rather than maintained manually.
- It is intended to be a stable paper/release snapshot for the current frozen benchmark bundle.
- Re-run this script after future tranche expansions or partial-rerun manifest refreshes.
