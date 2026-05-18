# Paint-booth dynamic surrogate v1 command-history demo

## Objective

Use the calibrated metric-level dynamic surrogate v1 to test a practical blower-command history before launching additional transient CFD.

The purpose is not to replace validation CFD, but to check whether the current dynamic surrogate gives useful control-oriented insight for multi-step operating scenarios.

## Inputs

- Branch: `feat/paint-booth-transient-poc`
- Surrogate script: `scripts/simulate_paint_booth_dynamic_metrics.py`
- Steady library: `cases/paint_booth_l2_u_sweep_v0/post_steady_metric_library.csv`
- Output directory:
  - `cases/paint_booth_l2_dynamic_metric_surrogate_v1_command_demo/`
- Command history CSV:
  - `cases/paint_booth_l2_dynamic_metric_surrogate_v1_command_demo/command_history_multistep.csv`
- Prediction CSV:
  - `cases/paint_booth_l2_dynamic_metric_surrogate_v1_command_demo/predicted_metrics_multistep.csv`
- Prediction summary JSON:
  - `cases/paint_booth_l2_dynamic_metric_surrogate_v1_command_demo/predicted_metrics_multistep_summary.json`
- Analysis JSON:
  - `cases/paint_booth_l2_dynamic_metric_surrogate_v1_command_demo/analysis_multistep.json`
- Figure:
  - `docs/figures/26-05-03_paint_booth_dynamic_surrogate_command_history_demo.svg`

Generated `cases/paint_booth_*` outputs are intentionally ignored by Git. The documentation and figure preserve the result summary.

## Command scenario

The scenario uses `dt = 0.25 s`, `0 <= t <= 130 s`.

To preserve the initial field convention used in transient validation, the command changes at the first sample after each listed switch time. For example, `t = 10.00 s` is still the pre-change state, and `t = 10.25 s` is the first `5.45 m/s` sample.

| Interval | Command |
|---|---:|
| 0.00--10.00 s | 4.36 m/s |
| 10.25--40.00 s | 5.45 m/s |
| 40.25--70.00 s | 4.36 m/s |
| 70.25--100.00 s | 6.00 m/s |
| 100.25--130.00 s | 4.36 m/s |

## Execution command

```bash
OUT=cases/paint_booth_l2_dynamic_metric_surrogate_v1_command_demo
python3 scripts/simulate_paint_booth_dynamic_metrics.py \
  --steady-library cases/paint_booth_l2_u_sweep_v0/post_steady_metric_library.csv \
  --command-csv "$OUT/command_history_multistep.csv" \
  --output-csv "$OUT/predicted_metrics_multistep.csv" \
  --output-json "$OUT/predicted_metrics_multistep_summary.json" \
  > "$OUT/log.simulate_multistep" 2>&1
```

## Key selected samples

| time [s] | U [m/s] | floor outflow [m3/s] | work-zone Uz [m/s] | near-car Uz [m/s] | near-car reverse fraction | filter Uz [m/s] |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 4.36 | 7.8431 | -0.20001 | -0.14350 | 0.23268 | -0.32863 |
| 10.00 | 4.36 | 7.8431 | -0.20001 | -0.14350 | 0.23268 | -0.32863 |
| 10.25 | 5.45 | 9.3927 | -0.24788 | -0.18880 | 0.16334 | -0.39357 |
| 10.50 | 5.45 | 9.7175 | -0.25623 | -0.19287 | 0.18602 | -0.40718 |
| 12.00 | 5.45 | 9.8037 | -0.25616 | -0.18815 | 0.23732 | -0.41079 |
| 20.00 | 5.45 | 9.8037 | -0.25472 | -0.18371 | 0.23176 | -0.41079 |
| 40.00 | 5.45 | 9.8037 | -0.25471 | -0.18366 | 0.22948 | -0.41079 |
| 40.25 | 4.36 | 8.2541 | -0.20708 | -0.14492 | 0.27138 | -0.34585 |
| 41.00 | 4.36 | 7.8469 | -0.19928 | -0.13988 | 0.24061 | -0.32879 |
| 50.00 | 4.36 | 7.8431 | -0.20001 | -0.14346 | 0.23012 | -0.32863 |
| 70.00 | 4.36 | 7.8431 | -0.20001 | -0.14350 | 0.23260 | -0.32863 |
| 70.25 | 6.00 | 10.1747 | -0.27228 | -0.21195 | 0.16323 | -0.42633 |
| 71.00 | 6.00 | 10.7874 | -0.28668 | -0.21595 | 0.21468 | -0.45200 |
| 80.00 | 6.00 | 10.7931 | -0.28261 | -0.20425 | 0.23105 | -0.45224 |
| 100.00 | 6.00 | 10.7931 | -0.28260 | -0.20417 | 0.22853 | -0.45224 |
| 100.25 | 4.36 | 8.4614 | -0.21068 | -0.14564 | 0.27046 | -0.35454 |
| 101.00 | 4.36 | 7.8488 | -0.19891 | -0.13803 | 0.23977 | -0.32887 |
| 110.00 | 4.36 | 7.8431 | -0.20001 | -0.14344 | 0.22985 | -0.32863 |
| 130.00 | 4.36 | 7.8431 | -0.20001 | -0.14350 | 0.23258 | -0.32863 |

## Window extrema

### Baseline 4.36 m/s, 0--10 s

- Work-zone `Uz_mean`: constant `-0.20001 m/s`
- Near-car `Uz_mean`: constant `-0.14350 m/s`
- Near-car reverse fraction: constant `0.23268`
- Filter-layer `Uz_mean`: constant `-0.32863 m/s`

### Upward 4.36 -> 5.45 m/s, 10--40 s

- Work-zone `Uz_mean`
  - start: `-0.20001 m/s`
  - end: `-0.25471 m/s`
  - most negative: `-0.25753 m/s @ 10.75 s`
- Near-car `Uz_mean`
  - start: `-0.14350 m/s`
  - end: `-0.18366 m/s`
  - most negative: `-0.19287 m/s @ 10.50 s`
- Near-car reverse fraction
  - start: `0.23268`
  - minimum: `0.16334 @ 10.25 s`
  - rebound maximum: `0.24138 @ 13.00 s`
  - end: `0.22948`
- Filter-layer `Uz_mean`
  - start: `-0.32863 m/s`
  - end: `-0.41079 m/s`
  - most negative: `-0.41079 m/s @ 15.75 s`

### Downward 5.45 -> 4.36 m/s, 40--70 s

- Work-zone `Uz_mean`
  - start: `-0.25471 m/s`
  - end: `-0.20001 m/s`
  - positive overshoot/weakest downdraft: `-0.19906 m/s @ 40.75 s`
- Near-car `Uz_mean`
  - start: `-0.18366 m/s`
  - end: `-0.14350 m/s`
  - weakest downdraft: `-0.13955 m/s @ 40.75 s`
- Near-car reverse fraction
  - start: `0.22948`
  - fast maximum: `0.27138 @ 40.25 s`
  - slow minimum: `0.22383 @ 43.25 s`
  - end: `0.23260`
- Filter-layer `Uz_mean`
  - start: `-0.41079 m/s`
  - end: `-0.32863 m/s`

### Upward 4.36 -> 6.00 m/s, 70--100 s

- Work-zone `Uz_mean`
  - start: `-0.20001 m/s`
  - end: `-0.28260 m/s`
  - most negative: `-0.28685 m/s @ 70.75 s`
- Near-car `Uz_mean`
  - start: `-0.14350 m/s`
  - end: `-0.20417 m/s`
  - most negative: `-0.21808 m/s @ 70.50 s`
- Near-car reverse fraction
  - start: `0.23260`
  - minimum: `0.16323 @ 70.25 s`
  - rebound maximum: `0.24102 @ 73.00 s`
  - end: `0.22853`
- Filter-layer `Uz_mean`
  - start: `-0.32863 m/s`
  - end: `-0.45224 m/s`
  - most negative: `-0.45224 m/s @ 76.00 s`

### Downward 6.00 -> 4.36 m/s, 100--130 s

- Work-zone `Uz_mean`
  - start: `-0.28260 m/s`
  - end: `-0.20001 m/s`
  - weakest downdraft: `-0.19858 m/s @ 100.75 s`
- Near-car `Uz_mean`
  - start: `-0.20417 m/s`
  - end: `-0.14350 m/s`
  - weakest downdraft: `-0.13753 m/s @ 100.75 s`
- Near-car reverse fraction
  - start: `0.22853`
  - fast maximum: `0.27046 @ 100.25 s`
  - slow minimum: `0.22320 @ 103.25 s`
  - end: `0.23258`
- Filter-layer `Uz_mean`
  - start: `-0.45224 m/s`
  - end: `-0.32863 m/s`

## Threshold-style control observation

For work-zone downdraft magnitude, using `-work_zone_Uz_mean`:

| Threshold | Predicted intervals satisfying threshold |
|---:|---|
| >= 0.20 m/s | `0.00--40.25 s`, `43.25--100.25 s`, `103.50--130.00 s` |
| >= 0.22 m/s | `10.25--40.00 s`, `70.25--100.00 s` |
| >= 0.24 m/s | `10.25--40.00 s`, `70.25--100.00 s` |
| >= 0.25 m/s | `10.50--40.00 s`, `70.25--100.00 s` |

This means a command above roughly `5.45 m/s` keeps the work-zone average downdraft above `0.24--0.25 m/s` almost immediately in the surrogate, but dropping back to `4.36 m/s` creates a short weak-downdraft overshoot around `0.75 s` after the switch.

## Reverse-fraction settling observation

Near-car reverse fraction remains the slow/wake-sensitive metric. Approximate 5% settling time relative to each 30 s segment endpoint:

| Switch | Approx. reverse-fraction 5% settling after switch |
|---|---:|
| 4.36 -> 5.45 | 22.00 s |
| 5.45 -> 4.36 | 23.25 s |
| 4.36 -> 6.00 | 22.00 s |
| 6.00 -> 4.36 | 22.75 s |

The initial mean-flow response is sub-second, but wake/reverse metrics carry a memory on the order of `20--25 s` under this v1 model.

## Interpretation

1. **Mean-flow metrics respond fast.**
   - Floor outflow, work-zone mean `Uz`, near-car mean `Uz`, and filter-layer `Uz` reach near-new operating levels within roughly `0.25--1 s` after a step.

2. **Mean-flow overshoot exists at command changes.**
   - Upward commands briefly overshoot toward stronger downdraft.
   - Downward commands briefly overshoot toward weaker downdraft.
   - The predicted weakest near-car downdraft after `6.00 -> 4.36 m/s` is `-0.13753 m/s` at `100.75 s`, weaker than the final `4.36 m/s` value `-0.14350 m/s`.

3. **Reverse fraction is the limiting dynamic state.**
   - Upward steps cause an immediate reverse-fraction dip, then rebound.
   - Downward steps cause an immediate reverse-fraction spike, then undershoot/recovery.
   - The slow recovery is why short command cycling could create history-dependent wake conditions even if average downdraft has already settled.

4. **Command frequency matters.**
   - This demo used long 30 s holds after each step. The reverse-fraction memory is still visible for about `22--23 s`, so shorter cycling periods would likely produce accumulated wake-state memory.

## Recommended next step

Before running more CFD, extend this into a small command-design sweep:

- same high/low commands but shorter hold times, e.g. `5 s`, `10 s`, `15 s`, `30 s`;
- ramp-limited command instead of pure step;
- compute objective metrics such as:
  - time below work-zone downdraft threshold;
  - peak near-car reverse fraction;
  - integral of reverse fraction above baseline;
  - command energy proxy or total flow cost.

This will identify which scenarios are worth validating with expensive transient CFD.
