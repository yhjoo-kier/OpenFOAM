# Frozen 20-Scene Benchmark Reference Status

> Date: 2026-03-18

## Summary

- Frozen benchmark subset size: 20 scenes (`A1/A2/A3/A4 × 5`)
- Reference CFD success: **20/20**
- Benchmark input-view export: **20/20** (`perspective`, `birdseye`, `floorplan`, `wireframe`, `section`)
- Evaluation scaffold size: **100 tasks** (`20 scenes × 5 views`)

This pass extended the previously recovered frozen-16 subset by adding the `*_05` tranche and then refreshing the aggregate manifests/scaffold back to full-set state.

## New `*_05` tranche outcome

| Case | Cat. | Room | Obst. | Success | Preset | Mesh | Notes |
|---|---|---|---:|---|---|---:|---|
| `bench_a1_05` | A1 | rectangular | 1 | yes | `robust` | 0.35 | solved directly on first attempt |
| `bench_a2_05` | A2 | rectangular | 3 | yes | `robust` | 0.35 | solved directly on first attempt |
| `bench_a3_05` | A3 | composite | 1 | yes | `robust` | 0.35 | solved directly on first attempt |
| `bench_a4_05` | A4 | composite | 4 | yes | `robust` | 0.35 | solved directly on first attempt |

## What changed in this pass

### 1. Frozen-set expansion

`scripts/generate_benchmark_scenes.py` was used to add one new scene per category with:

- `--start-index 5`
- `--count 1`
- `--min-simple-obstacles 1`

This keeps the simple categories (`A1`, `A3`) from regressing toward obstacle-free variants in later frozen tranches.

### 2. Reference/path stability check at the next tranche

All four newly added scenes converged at the initial `mesh_size=0.35` using the `robust` preset.

Interpretation:
- the current stabilization ladder did **not** need mesh fallback for this tranche,
- the next expansion tranche is still plausible without immediate generator redesign,
- recent failures remain concentrated in older harder variants (`bench_a2_03`, `bench_a4_03`) rather than worsening systematically with later seeds.

### 3. Aggregate bookkeeping refresh

Because `reference_batch_summary.json` is batch-oriented, a targeted 4-case rerun temporarily leaves the aggregate manifest incomplete for downstream tools.

After the `*_05` runs, the full-set aggregate state was refreshed with:

- `python3 scripts/run_benchmark_reference_batch.py --collect-only`
- `python3 scripts/scaffold_benchmark_evaluations.py`

This restored the canonical dataset state to:
- `benchmark/manifests/scene_manifest.json`: 20 scenes
- `benchmark/manifests/reference_batch_summary.json`: 20 reference cases
- `benchmark/renderings/renderings_manifest.json`: 20 rendering bundles
- `benchmark/evaluations/summary.json`: 100 evaluation tasks

## Interpretation

- The reference benchmark path is now stable across a **20-scene frozen bundle**.
- Composite-room cases remain compatible through scene generation, meshing, CFD, rendering, and evaluation scaffolding.
- The main blocker is still not the reference pipeline; it is the lack of Gemini backend credentials in the current cron shell for image-conditioned evaluation runs.
- The bookkeeping refresh step is now an important part of incremental tranche expansion and should continue after partial reruns.

## Recommended next actions

1. Extend evaluation summaries from geometry-only metrics toward CFD/result-side metrics.
2. Re-run scaffolded image-conditioned tasks once Gemini backend credentials are present in the runtime environment.
3. Take a clean local commit checkpoint for the frozen-20 benchmark state.
4. If expansion continues beyond 20 scenes, monitor whether first-attempt `0.35 + robust` success starts to degrade again.
