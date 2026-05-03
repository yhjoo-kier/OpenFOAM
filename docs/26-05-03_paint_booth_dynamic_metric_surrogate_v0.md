# Paint-booth dynamic metric surrogate baseline v0

## Objective

Build the first control-oriented dynamic quasi-steady surrogate using the existing steady CFD metric library and transient step-response data.

This baseline is intentionally metric-level rather than full-field:

```text
U_supply(t) -> steady metric interpolation -> dynamic lag/overshoot model -> metric(t)
```

It is meant to answer whether the current CFD data are enough to create a useful low-cost dynamic surrogate before launching more transient CFD cases.

## Inputs

### Steady metric library

Generated from the 9-point L2 steady sweep:

```text
cases/paint_booth_l2_u_sweep_v0/post_steady_metric_library.csv
```

Source metrics:

```text
cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u*/post_plenum_metrics.json
```

Supply points:

```text
2.18, 2.72, 3.27, 3.81, 4.36, 4.91, 5.45, 6.00, 6.54 m/s
```

Columns:

```text
U_supply_mps
case_dir
supply_inflow_m3s
floor_outflow_m3s
relative_mass_imbalance
work_zone_Uz_mean
near_car_Uz_mean
near_car_reverse_fraction_Uz_positive
filter_layer_panel_Uz_mean
filter_pressure_drop_proxy_m2_per_s2
```

### Transient validation cases

```text
4.36 -> 5.45 m/s, 10 s, writeInterval = 0.25 s
4.36 -> 5.45 m/s, 60 s, writeInterval = 2 s
5.45 -> 4.36 m/s, 60 s, writeInterval = 2 s
```

## Scripts added

```text
scripts/build_paint_booth_steady_metric_library.py
scripts/simulate_paint_booth_dynamic_metrics.py
```

`build_paint_booth_steady_metric_library.py` extracts the steady JSON metrics into a compact CSV table.

`simulate_paint_booth_dynamic_metrics.py`:

1. reads the steady metric library;
2. linearly interpolates metric targets by `U_supply_mps`;
3. applies metric-specific first-order lag;
4. adds an optional decaying overshoot state for upward mean-downdraft steps;
5. compares predictions against transient CFD time-series CSV files.

## Model equations

Steady target:

```text
y_ss(t) = interp_U[y_steady(U_supply(t))]
```

First-order lag:

```text
alpha = 1 - exp(-dt / tau_y)
y_base[k+1] = y_base[k] + alpha * (y_ss[k] - y_base[k])
```

Optional overshoot correction for step changes:

```text
z[k_step+] += beta_y(direction) * (y_ss,new - y_ss,old)
z[k+1] = z[k] exp(-dt / tau_decay,y)
y_pred[k] = y_base[k] + z[k]
```

The overshoot term is currently only enabled for upward mean `Uz` channels, calibrated from the high-resolution `4.36 -> 5.45 m/s` run.

## Default v0 parameters

```text
supply_inflow_m3s: tau = 0.02 s
floor_outflow_m3s: tau = 0.16 s
filter_layer_panel_Uz_mean: tau = 0.16 s
filter_pressure_drop_proxy_m2_per_s2: tau = 0.16 s
work_zone_Uz_mean: tau = 0.16 s
near_car_Uz_mean: tau = 0.14 s
near_car_reverse_fraction_Uz_positive: tau = 8.0 s
```

Overshoot:

```text
work_zone_Uz_mean:
  upward overshoot fraction = 0.10
  decay tau = 1.5 s

near_car_Uz_mean:
  upward overshoot fraction = 0.34
  decay tau = 1.8 s
```

## Validation outputs

Generated folder:

```text
cases/paint_booth_l2_dynamic_metric_surrogate_v0/
```

Files:

```text
u4p36_to_u5p45_10s_pred.csv
u4p36_to_u5p45_10s_pred_summary.json
u4p36_to_u5p45_60s_pred.csv
u4p36_to_u5p45_60s_pred_summary.json
u5p45_to_u4p36_60s_pred.csv
u5p45_to_u4p36_60s_pred_summary.json
```

## Validation error summary

### 4.36 -> 5.45 m/s, 10 s high-resolution

```text
work_zone_Uz_mean:
  RMSE = 0.00363 m/s
  MAE  = 0.00297 m/s
  max  = 0.01408 m/s

near_car_Uz_mean:
  RMSE = 0.00346 m/s
  MAE  = 0.00329 m/s
  max  = 0.00449 m/s

near_car_reverse_fraction_Uz_positive:
  RMSE = 0.01801
  MAE  = 0.01026
  max  = 0.07942

filter_layer_panel_Uz_mean:
  RMSE = 0.00685 m/s
  MAE  = 0.00674 m/s
  max  = 0.00960 m/s

floor_outflow_m3s:
  RMSE = 0.06723 m3/s
  MAE  = 0.01986 m3/s
  max  = 0.41718 m3/s
```

### 4.36 -> 5.45 m/s, 60 s

```text
work_zone_Uz_mean:
  RMSE = 0.00173 m/s
  MAE  = 0.00159 m/s
  max  = 0.00407 m/s

near_car_Uz_mean:
  RMSE = 0.00431 m/s
  MAE  = 0.00421 m/s
  max  = 0.00453 m/s

near_car_reverse_fraction_Uz_positive:
  RMSE = 0.00283
  MAE  = 0.00159
  max  = 0.01288

filter_layer_panel_Uz_mean:
  RMSE = 0.00679 m/s
  MAE  = 0.00668 m/s
  max  = 0.00697 m/s

floor_outflow_m3s:
  RMSE = 0.00752 m3/s
  MAE  = 0.00740 m3/s
  max  = 0.00785 m3/s
```

### 5.45 -> 4.36 m/s, 60 s

```text
work_zone_Uz_mean:
  RMSE = 0.00329 m/s
  MAE  = 0.00318 m/s
  max  = 0.00553 m/s

near_car_Uz_mean:
  RMSE = 0.00151 m/s
  MAE  = 0.00145 m/s
  max  = 0.00164 m/s

near_car_reverse_fraction_Uz_positive:
  RMSE = 0.00743
  MAE  = 0.00713
  max  = 0.01260

filter_layer_panel_Uz_mean:
  RMSE = 0.00531 m/s
  MAE  = 0.00523 m/s
  max  = 0.00541 m/s

floor_outflow_m3s:
  RMSE = 0.00612 m3/s
  MAE  = 0.00602 m3/s
  max  = 0.00657 m3/s
```

## Selected comparison: high-resolution upward step

```text
t = 0.00 s:
  work-zone Uz: pred -0.20001, CFD -0.20001
  near-car Uz:  pred -0.14350, CFD -0.14350
  reverse frac: pred  0.23268, CFD  0.23268
  filter Uz:    pred -0.32863, CFD -0.32863

t = 0.25 s:
  work-zone Uz: pred -0.24788, CFD -0.26196
  near-car Uz:  pred -0.18880, CFD -0.19137
  reverse frac: pred  0.23258, CFD  0.15316
  filter Uz:    pred -0.39357, CFD -0.40317

t = 2.00 s:
  work-zone Uz: pred -0.25616, CFD -0.25995
  near-car Uz:  pred -0.18815, CFD -0.18561
  reverse frac: pred  0.23196, CFD  0.23344
  filter Uz:    pred -0.41079, CFD -0.40381

t = 10.0 s:
  work-zone Uz: pred -0.25472, CFD -0.25614
  near-car Uz:  pred -0.18371, CFD -0.17921
  reverse frac: pred  0.23034, CFD  0.23049
  filter Uz:    pred -0.41079, CFD -0.40388
```

## Interpretation

The v0 surrogate is already useful for mean-flow metrics:

- Work-zone and near-car mean `Uz` errors are small, typically a few `0.001 m/s` to `0.004 m/s`.
- The high-resolution upward run confirms that an overshoot correction is needed for mean downdraft channels.
- The 60 s upward/downward validation is reasonably good for mean metrics.

The main limitation is still the reverse-flow fraction:

- The high-resolution run shows a large non-monotonic dip at `t = 0.25 s` that the simple slow-lag model does not reproduce.
- The scalar initial-final change in reverse fraction is small, while transient excursions are large.
- This metric should not be modeled as a simple first-order response.

The filter-layer metric also shows a systematic offset because the transient endpoint does not exactly match the existing steady library endpoint. This may reflect solver/settings differences between steady and transient cases or the fact that the transient has not fully matched the exact steady reference for every local metric.

## Recommended next step

Run a paired high-resolution downward step:

```text
5.45 -> 4.36 m/s
endTime = 10 s
writeInterval = 0.25 s
parallel: 4 MPI ranks
```

Reason:

- The upward high-resolution run revealed an early overshoot that was invisible in 2 s output.
- The downward 60 s case has only 2 s output, so its early response may also hide under/overshoot.
- Direction-specific overshoot/undershoot parameters are needed before improving v1 of the surrogate.

After that, update the surrogate from v0 to v1:

```text
mean-Uz channels:
  direction-specific overshoot/decay

reverse-flow channel:
  empirical response kernel or two-state wake model

filter/through-flow channel:
  fast lag plus endpoint calibration check
```
