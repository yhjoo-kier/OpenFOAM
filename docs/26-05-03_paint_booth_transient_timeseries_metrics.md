# Paint-booth transient step-response time-series metrics

## Objective

Use the completed 60 s transient step-response case to extract metric time series and estimate first response-time numbers for the dynamic quasi-steady surrogate baseline.

This is the first practical use of the literature/methodology note:

```text
steady CFD library + transient step-response identification
```

## Case

- Case: `cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s`
- Input: supply velocity step from initial steady `4.36 m/s` field to boundary `5.45 m/s` at `t = 0`
- Solver: `pimpleFoam`
- Written times: `0, 2, 4, ..., 60 s`
- VTK conversion: `foamToVTK` for all written times

## New script

Created:

```text
scripts/postprocess_paint_booth_transient_timeseries.py
```

Outputs:

```text
cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s/post_transient_metrics_timeseries.csv
cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s/post_transient_metrics_timeseries.json
cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s/post_transient_metrics_timeseries_summary.json
```

Important implementation note:

- OpenFOAM v1912 `foamToVTK` directory suffixes are write indices, e.g. `_170`, `_342`, not physical times.
- The script therefore parses `VTK/<case>.vtm.series` to recover the authoritative physical time mapping.

## Commands

Convert all written times:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD/cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s" \
  openfoam-pipeline-local:latest \
  bash -lc 'foamToVTK > log.foamToVTK_all_times 2>&1'
```

Extract time-series metrics:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 scripts/postprocess_paint_booth_transient_timeseries.py cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s > cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s/log.postprocess_timeseries 2>&1'
```

## Extracted metrics

For every written time, the script extracts:

```text
time_s
supply_inflow_m3s
floor_outflow_m3s
net_out_minus_in_m3s
relative_imbalance_out_minus_in
work_zone_Uz_mean
near_car_Uz_mean
near_car_reverse_fraction_Uz_positive
filter_layer_panel_Uz_mean
filter_pressure_drop_proxy_m2_per_s2
```

## Time-series summary

Number of extracted times:

```text
31
```

Time range:

```text
0 s to 60 s
```

Selected early values:

```text
t = 0 s
  work_zone_Uz_mean                  = -0.2000069320
  near_car_Uz_mean                   = -0.1434996426
  near_car_reverse_fraction          =  0.2326842387
  filter_layer_panel_Uz_mean         = -0.3286292255

t = 2 s
  work_zone_Uz_mean                  = -0.2599384189
  near_car_Uz_mean                   = -0.1855722070
  near_car_reverse_fraction          =  0.2335184332
  filter_layer_panel_Uz_mean         = -0.4038213193

t = 10 s
  work_zone_Uz_mean                  = -0.2561358809
  near_car_Uz_mean                   = -0.1792097241
  near_car_reverse_fraction          =  0.2304680205
  filter_layer_panel_Uz_mean         = -0.4038751125
```

Final values at `t = 60 s`:

```text
work_zone_Uz_mean                  = -0.2561563551
near_car_Uz_mean                   = -0.1791233122
near_car_reverse_fraction          =  0.2306298791
filter_layer_panel_Uz_mean         = -0.4038849771
```

## Response-time estimates

The script estimates 63%, 90%, and 95% crossing times using the normalized response:

```text
r(t) = [y(t) - y(0)] / [y(60) - y(0)]
```

Because the write interval is `2 s`, sub-2 s crossing values are linear interpolations between `t = 0` and `t = 2 s`; they should be interpreted as coarse estimates only.

### Work-zone `Uz_mean`

```text
initial: -0.2000069320
final:   -0.2561563551
delta:   -0.0561494231

t63: 1.18 s
t90: 1.69 s
t95: 1.78 s
settling within 10% band: 2 s
settling within 5% band:  6 s
overshoot beyond final:   6.7% of observed delta
```

Interpretation:

- Work-zone mean downdraft reacts very quickly to this imposed boundary step.
- There is a small early overshoot: the `t = 2 s` value is slightly more negative than the `t = 60 s` value.

### Near-car `Uz_mean`

```text
initial: -0.1434996426
final:   -0.1791233122
delta:   -0.0356236696

t63: 1.07 s
t90: 1.52 s
t95: 1.61 s
settling within 10% band: 6 s
settling within 5% band:  6 s
overshoot beyond final:   18.1% of observed delta
```

Interpretation:

- The near-car mean also reacts quickly but has a stronger early overshoot than the work-zone mean.
- This suggests the near-car averaged flow metric is not purely first order.

### Filter-layer panel `Uz_mean`

```text
initial: -0.3286292255
final:   -0.4038849771
delta:   -0.0752557516

t63: 1.27 s
t90: 1.80 s
t95: 1.90 s
settling within 10% band: 2 s
settling within 5% band:  2 s
overshoot beyond final:   0%
```

Interpretation:

- The filter-layer mean response is nearly immediate at this write resolution.
- A shorter write interval would be needed to resolve the true time constant.

### Near-car reverse-flow fraction

```text
initial: 0.2326842387
final:   0.2306298791
delta:  -0.0020543596

t63: 9.17 s
t90: 9.67 s
t95: 9.76 s
settling within 10% band: 22 s
settling within 5% band:  26 s
```

Interpretation:

- Reverse-flow fraction behaves differently from mean velocity metrics.
- It has non-monotonic behavior and a small final net change, so normalized overshoot/undershoot values are sensitive.
- This metric is likely a better indicator of wake/recirculation adjustment than region-mean `Uz` alone.

## Main conclusion

For the imposed `4.36 -> 5.45 m/s` supply velocity boundary step:

1. Mean flow metrics (`work_zone_Uz_mean`, `near_car_Uz_mean`, `filter_layer_panel_Uz_mean`) respond within the first few seconds at the current `2 s` sampling resolution.
2. The reverse-flow fraction has a slower and more non-monotonic response, settling on the order of `20--30 s` by the current metric definition.
3. A single first-order time constant can probably approximate mean downdraft metrics, but not the reverse-flow/wake metric.
4. For a control surrogate, the first baseline should therefore use:

```text
fast lag for mean downdraft metrics
slower/nonlinear or region-specific lag for reverse-flow / wake-sensitive metrics
```

## Implication for dynamic quasi-steady surrogate

A practical first scalar surrogate can be:

```text
y_dyn(k+1) = y_dyn(k) + alpha_y [ y_ss(q_k) - y_dyn(k) ]
alpha_y = 1 - exp(-dt / tau_y)
```

Candidate initial `tau_y` values from this case:

```text
work_zone_Uz_mean:          tau ≈ 1--2 s, with small overshoot
near_car_Uz_mean:           tau ≈ 1--2 s, with stronger overshoot
filter_layer_panel_Uz_mean: tau <~ 2 s, unresolved by 2 s write interval
near_car_reverse_fraction:  not well represented by single first-order tau; use settling 20--30 s or fit a second-order/nonlinear model later
```

## Recommended next cases

To make the dynamic model useful beyond one step, run at least:

```text
5.45 -> 4.36 m/s  downward step
4.36 -> 6.54 m/s  larger upward step
2.18 -> 4.36 m/s  low-to-nominal step
```

The downward step is especially important because recirculation/reverse-flow metrics may be asymmetric.
