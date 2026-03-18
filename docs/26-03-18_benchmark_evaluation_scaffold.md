# Benchmark Evaluation Scaffold (Frozen-20 × 5 views)

> Date: 2026-03-18

## Summary

The benchmark bundle now has a canonical `benchmark/evaluations/` scaffold so image-conditioned pipeline runs can be organized against the frozen reference set without ad-hoc folder creation.

Current scaffold size:
- Cases: 20
- Views per case: 5 (`perspective`, `birdseye`, `floorplan`, `wireframe`, `section`)
- Evaluation tasks: 100 total
- Current task status snapshot: 1 `blocked`, 99 `pending`

## What was added

### 1. Evaluation scaffold script

New helper:
- `scripts/scaffold_benchmark_evaluations.py`

It reads:
- `benchmark/manifests/scene_manifest.json`
- `benchmark/manifests/reference_batch_summary.json` (or another compatible aggregate reference-status manifest)
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
1. sweep all 100 benchmark tasks systematically,
2. attach run status later (`pending` / `running` / `success` / `failed` / `blocked`),
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
- the actual image-conditioned generation is still blocked in this cron shell because `GEMINI_API_KEY` is not set,
- `scripts/run_benchmark_evaluation_task.py` now preflights backend availability and records this as `status=blocked` (not a false benchmark failure) with the captured backend reason in `evaluation_summary.json`.

## Follow-up improvement in this pass

Two useful follow-up upgrades landed around the scaffold:

1. A native composite-room repair upgrade was completed in `scripts/repair_indoor_scene.py`.
   - instead of always forcing repaired scenes onto a west↔east opening pair, the repair stage now infers the better repair axis (`x` or `y`) from the original opening layout plus exposed outer-wall availability,
   - preserves the inlet wall when it already lies on the chosen axis,
   - records the selected axis in `meta.repair_opening_axis` and the repair info payload.

2. `scripts/run_benchmark_evaluation_task.py` now emits richer geometry-comparison metrics even before CFD-result comparison is added.
   - this means future image-conditioned runs can be judged immediately on scene-structure quality, not just pass/fail bookkeeping,
   - and blocked runs remain cleanly distinguishable from true benchmark failures.

3. A first CFD/result-side comparison helper has now been added as well:
   - `scripts/compute_benchmark_cfd_metrics.py`
   - it samples reference/predicted VTK volumes on a shared normalized room-coordinate lattice,
   - and `run_benchmark_evaluation_task.py` writes the result to `benchmark/evaluations/<case>/<view>/cfd_metrics.json` after successful runs.

The composite-room repair upgrade reduces geometry drift for composite benchmark cases whose stable repair is better aligned with a south↔north corridor.
