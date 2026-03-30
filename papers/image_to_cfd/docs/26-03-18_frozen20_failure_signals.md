# Frozen 20-Scene Benchmark Failure Signals / Robustness Notes

> Date: 2026-03-18

## Summary

The frozen benchmark bundle is currently **20/20 successful** on the reference path, but not all scenes are equally easy. The aggregate manifest now provides a clear picture of which cases needed extra stabilization beyond the default `0.35 + robust + RAS` lane.

Quick breakdown:
- Total scenes: **20**
- Solved on first attempt: **15**
- Needed more than one attempt: **5**
- Needed mesh fallback from `0.35` to `0.25`: **2**
- Needed a non-`robust` preset: **3**
- Needed laminar fallback: **1**

## Aggregate stabilization profile

### Preset distribution

- `robust`: 17
- `conservative`: 1
- `ultra_robust`: 1
- `laminar_fallback`: 1

### Solver mode distribution

- `RAS`: 19
- `laminar`: 1

### Successful mesh size distribution

- `0.35`: 18
- `0.25`: 2

### Attempt-count distribution

- `1 attempt`: 15
- `2 attempts`: 3
- `3 attempts`: 1
- `4 attempts`: 1

## High-signal cases

| Case | Category | Attempts | Final preset | Final mode | Final mesh | Signal |
|---|---|---:|---|---|---:|---|
| `bench_a2_01` | A2 | 4 | `conservative` | `RAS` | 0.35 | hardest preset-ladder case in the current frozen set |
| `bench_a4_02` | A4 | 3 | `laminar_fallback` | `laminar` | 0.35 | only case needing laminar fallback; re-run after composite repair-path fix confirmed the old repair error was stale metadata, not an active compatibility blocker |
| `bench_a1_04` | A1 | 2 | `ultra_robust` | `RAS` | 0.35 | simple geometry but still sensitive to stabilization settings |
| `bench_a2_03` | A2 | 2 | `robust` | `RAS` | 0.25 | mesh-sensitivity case; required finer mesh fallback |
| `bench_a4_03` | A4 | 2 | `robust` | `RAS` | 0.25 | composite + obstacle-dense mesh-sensitivity case |

## Interpretation

1. **Most of the benchmark is now routine.**
   - 15/20 scenes finish on the first attempt with the default stabilized lane.
   - This is good enough to treat the reference dataset as operational rather than experimental.

2. **The remaining difficulty is concentrated, not random.**
   - `bench_a2_03` and `bench_a4_03` remain the clearest mesh-sensitivity pair.
   - `bench_a2_01` and `bench_a4_02` are preset/solver-stability outliers rather than generator-schema failures.

3. **Composite rooms are no longer the main failure story by themselves.**
   - Both rectangular and composite families appear in the high-signal list.
   - The current bottleneck has shifted from room-schema compatibility toward solver robustness / mesh quality for a handful of layouts.
   - A targeted March 18 rerun of `bench_a4_02` verified that `repair_indoor_scene.py` now repairs/validates that composite scene cleanly; the previously recorded repair failure was just stale pre-fix metadata in an older stabilization summary.

4. **For paper reporting, success-rate alone will hide useful nuance.**
   - A small “stabilization cost” table (attempt count, final preset, mesh fallback use) would better communicate benchmark difficulty than a binary success flag.

## Practical follow-up ideas

### Near-term

- Preserve these five scenes as the canonical **stress subset** for future pipeline changes.
- When modifying meshing or solver presets, re-run at least:
  - `bench_a2_01`
  - `bench_a4_02`
  - `bench_a2_03`
  - `bench_a4_03`

### For the paper

Potential lightweight robustness descriptors to report per scene:
- first-attempt success (yes/no)
- successful preset tier
- successful mesh size
- attempt count
- whether laminar fallback was required

### For future expansion beyond 20 scenes

Watch for these regression patterns:
- more cases dropping from `0.35` to `0.25`
- more cases escaping RAS into laminar fallback
- simple-category scenes (`A1`, `A3`) unexpectedly needing `ultra_robust` or worse

If those rates grow materially in later tranches, the next fix should likely target meshing quality heuristics rather than scene-generation diversity.
