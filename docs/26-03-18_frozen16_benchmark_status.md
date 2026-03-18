# Frozen 16-Scene Benchmark Reference Status

> Date: 2026-03-18

## Summary

- Frozen benchmark subset size: 16 scenes (`A1/A2/A3/A4 × 4`)
- Reference CFD success: **16/16**
- Benchmark input-view export: **16/16** (`perspective`, `birdseye`, `floorplan`, `wireframe`, `section`)
- Evaluation scaffold size: **80 tasks** (`16 scenes × 5 views`)

This pass extended the previously recovered frozen-12 subset by adding the `*_04` tranche and also fixed an aggregate-manifest bookkeeping regression discovered during that expansion.

## New `*_04` tranche outcome

| Case | Cat. | Room | Obst. | Success | Preset | Mesh | Notes |
|---|---|---|---:|---|---|---:|---|
| `bench_a1_04` | A1 | rectangular | 1 | yes | `ultra_robust` | 0.35 | solved at initial mesh; required the strongest preset among the new four |
| `bench_a2_04` | A2 | rectangular | 2 | yes | `robust` | 0.35 | solved directly |
| `bench_a3_04` | A3 | composite | 1 | yes | `robust` | 0.35 | solved directly |
| `bench_a4_04` | A4 | composite | 3 | yes | `robust` | 0.35 | solved directly |

## What was fixed in this pass

### 1. Aggregate rendering manifest refresh

`render_benchmark_views.py` previously rewrote `benchmark/renderings/renderings_manifest.json` using only the most recently rendered batch.

That was harmless for one-shot runs, but it broke downstream tooling after incremental benchmark expansion because the global rendering manifest no longer represented the full frozen set.

The renderer now rebuilds the aggregate rendering manifest by scanning all per-case render manifests under:

```text
benchmark/renderings/<case>/manifest.json
```

### 2. Evaluation scaffolding decoupled from frozen-12-only manifest naming

`scripts/scaffold_benchmark_evaluations.py` originally hard-coded:

```text
benchmark/manifests/frozen12_reference_status.json
```

It now accepts manifest paths via CLI arguments and defaults to the generic aggregate reference manifest:

```text
benchmark/manifests/reference_batch_summary.json
```

It also preserves existing task state when the scaffold is refreshed, rather than resetting all tasks to `pending` blindly.

## Dataset state after refresh

- `benchmark/manifests/scene_manifest.json` now contains **16** scene entries.
- `benchmark/manifests/reference_batch_summary.json` now reflects the full 16-scene artifact set after `--collect-only` refresh.
- `benchmark/renderings/renderings_manifest.json` now reflects the full 16-scene rendering bundle.
- `benchmark/evaluations/summary.json` now reports **80** tasks with the prior failed smoke-test record preserved.

## Interpretation

- The benchmark bundle is no longer just 12-scene-stable; it is now 16-scene-stable for the reference path.
- The main remaining blocker has shifted from reference CFD generation to image-conditioned evaluation backend availability (`GEMINI_API_KEY` absent in the current cron shell for API-mode runs).
- The evaluation scaffold now distinguishes a true benchmark failure from a runtime-environment block: the smoke-tested `bench_a1_01/perspective` task is recorded as `blocked`, not `failed`.
- Composite-room repair fidelity also improved in this pass: `repair_indoor_scene.py` now infers/persists a repair opening axis instead of blindly forcing west↔east normalization.
- The manifest-refresh fix matters because future frozen-set expansion would otherwise silently corrupt aggregate bookkeeping.

## Recommended next actions

1. Take a local commit checkpoint for the 16-scene + manifest-refresh state.
2. Add richer evaluation metrics beyond bookkeeping-only summaries.
3. Resume image-conditioned evaluation once the Gemini backend environment is available in the runtime shell.
4. Consider the next expansion tranche to 20 scenes while watching preset/mesh fallback rates.
