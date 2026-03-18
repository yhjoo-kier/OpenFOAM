# Benchmark Evaluation Scaffold (Frozen-12 × 5 views)

> Date: 2026-03-18

## Summary

The benchmark bundle now has a canonical `benchmark/evaluations/` scaffold so image-conditioned pipeline runs can be organized against the frozen reference set without ad-hoc folder creation.

Current scaffold size:
- Cases: 12
- Views per case: 5 (`perspective`, `birdseye`, `floorplan`, `wireframe`, `section`)
- Evaluation tasks: 60 total

## What was added

### 1. Evaluation scaffold script

New helper:
- `scripts/scaffold_benchmark_evaluations.py`

It reads:
- `benchmark/manifests/scene_manifest.json`
- `benchmark/manifests/frozen12_reference_status.json`
- `benchmark/renderings/renderings_manifest.json`

and materializes:
- `benchmark/evaluations/summary.json`
- `benchmark/evaluations/manifest.json`
- `benchmark/evaluations/<case>/manifest.json`
- `benchmark/evaluations/<case>/<view>/task.json`

### 2. Canonical per-task structure

Each evaluation task directory now contains symlinked reference inputs/artifacts:

```text
benchmark/evaluations/<case>/<view>/
├── input.png
├── reference_scene.json
├── reference_case
├── reference_results
└── task.json
```

Expected future outputs are pre-declared in `task.json`:
- `predicted_scene.json`
- `predicted_case/`
- `predicted_results/`
- `evaluation_summary.json`

## Why this helps

This removes ambiguity about where future benchmark-evaluation runs should write artifacts.

It also makes it easier to:
1. sweep all 60 benchmark tasks systematically,
2. attach run status later (`pending` / `running` / `success` / `failed`),
3. compare predicted geometry/CFD outputs against the frozen reference bundle.

## Small pipeline fix bundled with this checkpoint

`run_indoor_stabilized.py` now uses the **actual scene JSON that succeeded** (`original` or `repaired`) when exporting the 2D comparison and 3D rendering artifacts.

This avoids visualization drift in repaired-scene recoveries, where a successful fallback run could previously still be rendered using the original unrepaired scene JSON.

## Runner added

New helper:
- `scripts/run_benchmark_evaluation_task.py`

What it does for one scaffold task:
- reads `benchmark/evaluations/<case>/<view>/task.json`
- marks task status `running` → `success` / `failed`
- calls the existing `scripts/run_indoor_stabilized.py` image-conditioned pipeline
- writes/captures scaffold outputs when available:
  - `predicted_scene.json`
  - `predicted_case/` (symlink to the canonical generated case)
  - `predicted_results/` (symlink to the canonical generated results)
  - `evaluation_summary.json`
- refreshes aggregate files:
  - `benchmark/evaluations/summary.json`
  - `benchmark/evaluations/manifest.json`
  - `benchmark/evaluations/<case>/manifest.json`

The first summary is intentionally lightweight: task/run bookkeeping, pipeline return code, stabilization summary passthrough, and a small scene-structure comparison (`room_blocks`, obstacle count, opening count/types).

## Smoke-test note

Smoke test attempted on:
- `bench_a1_01 / perspective`

Command:
- `python3 scripts/run_benchmark_evaluation_task.py --case-name bench_a1_01 --view perspective --end-time 200`

Outcome:
- runner wiring worked and task/aggregate manifests updated correctly,
- but the actual image-conditioned generation failed immediately because `GEMINI_API_KEY` was not set in the current shell, so the task was recorded as `failed` with the captured backend error in `evaluation_summary.json`.
