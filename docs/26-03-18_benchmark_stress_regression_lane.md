# Benchmark Stress-Subset Regression Lane

> Date: 2026-03-18

## Why this was added

The frozen benchmark bundle is currently integrity-clean (`20/20` reference success, `20/20` render bundles, `100/100` evaluation tasks), so the next useful step was not another blind full rerun.

Instead, the pipeline now has a dedicated **stress-subset regression lane** for the five known hard cases already identified by `benchmark/manifests/dataset_integrity_summary.json`:

- `bench_a1_04`
- `bench_a2_01`
- `bench_a2_03`
- `bench_a4_02`
- `bench_a4_03`

These are the cases that previously signaled:
- multi-attempt recovery,
- preset fallback,
- mesh fallback, and/or
- laminar fallback.

## New helper

Added:

```text
scripts/run_benchmark_stress_subset.py
```

This script:
1. reads `benchmark/manifests/dataset_integrity_summary.json`,
2. extracts the current `stress_cases`,
3. maps them back to canonical benchmark scene JSON files,
4. reruns only that subset through `scripts/run_benchmark_reference_batch.py`, and
5. can optionally refresh the full frozen-set manifests and re-run bundle verification.

## Typical usage

### 1. Quick subset bookkeeping / selection smoke test

```bash
python3 scripts/run_benchmark_stress_subset.py --collect-only --verify-after
```

### 2. Real regression rerun after meshing / solver changes

```bash
python3 scripts/run_benchmark_stress_subset.py --verify-after
```

### 3. Subset rerun + restore full frozen-set aggregate manifests

```bash
python3 scripts/run_benchmark_stress_subset.py \
  --refresh-full-manifests \
  --verify-after
```

## Why this matters for the paper work

This is a practical bridge between two needs:

1. **Engineering need**
   - future meshing/solver/pipeline edits should be checked against the hardest known layouts first,
   - but rerunning the entire frozen set every time is slower and obscures where regressions start.

2. **Paper need**
   - the benchmark now has an explicit notion of **robustness-signaling scenes**, not just binary success.
   - that supports later reporting of stabilization cost / difficulty concentration rather than only aggregate success rate.

## Current status

- The canonical stress subset is now machine-readable via the integrity summary.
- A dedicated regression entrypoint now exists to reuse the existing indoor CFD batch pipeline on just that subset.
- A collect-only smoke test succeeded with:

```bash
python3 scripts/run_benchmark_stress_subset.py --collect-only --refresh-full-manifests --verify-after
```

  and re-confirmed the expected 5 stress cases while restoring the full frozen-set aggregates and a clean verifier result.
- The helper also writes a run summary to:

```text
benchmark/manifests/stress_subset_run_summary.json
```

## Recommended operating pattern

After any future change to:
- `scene_to_gmsh.py`
- `repair_indoor_scene.py`
- `run_indoor_stabilized.py`
- solver presets / stabilization ladders
- render/reference manifest plumbing

run the stress subset first, then refresh/verify the full bundle only if the subset remains clean.
