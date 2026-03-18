# Benchmark CFD Comparison Metrics (normalized-grid prototype)

> Date: 2026-03-18

## Summary

The benchmark evaluation scaffold now has a first result-side comparison path for CFD outputs, not just geometry.

New helper:
- `scripts/compute_benchmark_cfd_metrics.py`

This compares a reference case and a predicted case by sampling both VTK volumes on the same **normalized room-coordinate lattice**.

## Why this approach

Direct cell-to-cell comparison is brittle because:
- predicted geometry may differ from the reference,
- meshes are not aligned,
- composite rooms can have different local block topology even when the high-level layout is similar.

Instead, the new comparator:
1. reads the reference and predicted room bounds from scene JSON,
2. builds a shared normalized lattice in `(x/Lx, y/Ly, z/Lz)` coordinates,
3. samples both VTK fields at corresponding normalized positions,
4. compares only points that are valid in both fluid domains.

This gives a coarse but robust room-relative similarity measure for benchmark scoring.

## Current metrics

For the overlap region of valid sampled points, the script reports:

- overlap coverage
  - `ratio_vs_total`
  - `ratio_vs_union`
- velocity vector metrics
  - L2 MAE / RMSE
  - relative RMSE vs reference RMS
  - mean direction cosine
- velocity magnitude metrics
  - MAE / RMSE
  - relative RMSE
- pressure metrics
  - MAE / RMSE
  - relative RMSE
- aggregate `cfd_score`
  - currently a simple average of overlap similarity, velocity-magnitude similarity, direction similarity, and pressure similarity when available

## Smoke tests

### Identity check: rectangular benchmark case

Command:

```bash
docker run --rm \
  -v /home/yhjoo/projects/OpenFOAM:/home/yhjoo/projects/OpenFOAM \
  -w /home/yhjoo/projects/OpenFOAM \
  openfoam-pipeline-local:latest \
  python3 scripts/compute_benchmark_cfd_metrics.py \
  --reference-scene benchmark/scenes/a1_01.json \
  --reference-case cases/bench_a1_01 \
  --predicted-scene benchmark/scenes/a1_01.json \
  --predicted-case cases/bench_a1_01 \
  --output benchmark/generated_a1_01_selfcheck_cfd_metrics.json
```

Outcome:
- overlap ratio vs union: `1.0`
- velocity / pressure errors: `0.0`
- aggregate `cfd_score`: `1.0`

### Identity check: composite benchmark case

Command:

```bash
docker run --rm \
  -v /home/yhjoo/projects/OpenFOAM:/home/yhjoo/projects/OpenFOAM \
  -w /home/yhjoo/projects/OpenFOAM \
  openfoam-pipeline-local:latest \
  python3 scripts/compute_benchmark_cfd_metrics.py \
  --reference-scene benchmark/scenes/a3_04.json \
  --reference-case cases/bench_a3_04 \
  --predicted-scene benchmark/scenes/a3_04.json \
  --predicted-case cases/bench_a3_04 \
  --output benchmark/generated_a3_04_selfcheck_cfd_metrics.json
```

Outcome:
- overlap ratio vs union: `1.0`
- velocity / pressure errors: `0.0`
- aggregate `cfd_score`: `1.0`

### Non-identity sanity check

`bench_a1_01` vs `bench_a1_03` produced a clearly lower score:
- overlap ratio vs union: `0.92963`
- velocity-magnitude relative RMSE: `0.905871`
- pressure relative RMSE: `0.727587`
- aggregate `cfd_score`: `0.513593`

So the metric is not just returning trivially high scores.

## Integration status

`run_benchmark_evaluation_task.py` now attempts to compute these CFD metrics automatically after a successful image-conditioned run and writes:

```text
benchmark/evaluations/<case>/<view>/cfd_metrics.json
```

If the comparison helper fails, the evaluation summary keeps the failure payload instead of crashing the benchmark bookkeeping.

## Limits / next step

This is intentionally a first coarse comparison layer.

Still worth improving later:
1. add plane-specific metrics aligned with paper figures (`floor`, `section`, longitudinal slice),
2. consider temperature/scalar comparisons if thermal runs are introduced,
3. calibrate the aggregate `cfd_score` weighting before manuscript reporting.
