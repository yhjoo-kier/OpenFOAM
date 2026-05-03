# Paint-booth dynamic metric surrogate v1

## Objective

Update the metric-level dynamic quasi-steady surrogate after obtaining paired high-resolution transient data:

```text
upward:   4.36 -> 5.45 m/s, 10 s, writeInterval = 0.25 s
downward: 5.45 -> 4.36 m/s, 10 s, writeInterval = 0.25 s
```

The v0 surrogate used steady metric interpolation, metric-specific first-order lags, and upward-only mean-`Uz` overshoot correction. The high-resolution downward case showed that downward steps also have early overshoot and that near-car reverse fraction is strongly non-monotonic in both directions.

## Script updated

```text
scripts/simulate_paint_booth_dynamic_metrics.py
```

The default model was updated from:

```text
metric_level_dynamic_quasi_steady_v0
```

to:

```text
metric_level_dynamic_quasi_steady_v1
```

## Model structure

Steady target interpolation:

```text
y_ss(t) = interp_U[y_steady(U_supply(t))]
```

First-order base channel:

```text
alpha = 1 - exp(-dt / tau_y)
y_base[k+1] = y_base[k] + alpha [y_ss[k] - y_base[k]]
```

Direction-specific mean-flow overshoot:

```text
z_y[k_step+] += beta_y(direction) [y_ss,new - y_ss,old]
z_y[k+1] = z_y[k] exp(-dt / tau_decay_y(direction))
y_pred = y_base + z_y
```

Reverse-fraction empirical impulse kernel:

```text
r_pred = r_base + sum_i a_i(direction) exp(-t_since_step / tau_i(direction))
```

This is not a full physics model of wake topology. It is a compact empirical response kernel to reproduce the fast dip/rise and slower recovery seen in high-resolution CFD.

## Default v1 parameters

Base lags:

```text
supply_inflow_m3s: tau = 0.02 s
floor_outflow_m3s: tau = 0.16 s
filter_layer_panel_Uz_mean: tau = 0.16 s
filter_pressure_drop_proxy_m2_per_s2: tau = 0.16 s
work_zone_Uz_mean: tau = 0.16 s
near_car_Uz_mean: tau = 0.14 s
near_car_reverse_fraction_Uz_positive: tau = 8.0 s
```

Mean `Uz` overshoot:

```text
work_zone_Uz_mean:
  upward beta = 0.10, tau_decay = 1.5 s
  downward beta = 0.14, tau_decay = 0.45 s

near_car_Uz_mean:
  upward beta = 0.34, tau_decay = 1.8 s
  downward beta = 0.15, tau_decay = 2.0 s
```

Reverse fraction kernels:

```text
upward 4.36 -> 5.45:
  a1 = -0.140, tau1 = 0.90 s
  a2 = +0.040, tau2 = 3.00 s

downward 5.45 -> 4.36:
  a1 = +0.080, tau1 = 0.90 s
  a2 = -0.020, tau2 = 4.00 s
```

## Validation outputs

Generated folder:

```text
cases/paint_booth_l2_dynamic_metric_surrogate_v1/
```

Files:

```text
u4p36_to_u5p45_10s_pred.csv
u4p36_to_u5p45_10s_pred_summary.json
u5p45_to_u4p36_10s_pred.csv
u5p45_to_u4p36_10s_pred_summary.json
u4p36_to_u5p45_60s_pred.csv
u4p36_to_u5p45_60s_pred_summary.json
u5p45_to_u4p36_60s_pred.csv
u5p45_to_u4p36_60s_pred_summary.json
```

## v1 validation errors

### Upward high-resolution: `4.36 -> 5.45`, 10 s

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
  RMSE = 0.00357
  MAE  = 0.00264
  max  = 0.01018

filter_layer_panel_Uz_mean:
  RMSE = 0.00685 m/s
  MAE  = 0.00674 m/s

floor_outflow_m3s:
  RMSE = 0.06723 m3/s
  MAE  = 0.01986 m3/s
```

### Downward high-resolution: `5.45 -> 4.36`, 10 s

```text
work_zone_Uz_mean:
  RMSE = 0.00466 m/s
  MAE  = 0.00438 m/s
  max  = 0.01107 m/s

near_car_Uz_mean:
  RMSE = 0.00226 m/s
  MAE  = 0.00164 m/s
  max  = 0.00859 m/s

near_car_reverse_fraction_Uz_positive:
  RMSE = 0.00589
  MAE  = 0.00546
  max  = 0.00803

filter_layer_panel_Uz_mean:
  RMSE = 0.00639 m/s
  MAE  = 0.00576 m/s

floor_outflow_m3s:
  RMSE = 0.06456 m3/s
  MAE  = 0.01762 m3/s
```

### Upward 60 s

```text
work_zone_Uz_mean: RMSE = 0.00173 m/s
near_car_Uz_mean:  RMSE = 0.00431 m/s
reverse fraction:  RMSE = 0.00145
filter-layer Uz:   RMSE = 0.00679 m/s
floor outflow:     RMSE = 0.00752 m3/s
```

### Downward 60 s

```text
work_zone_Uz_mean: RMSE = 0.00330 m/s
near_car_Uz_mean:  RMSE = 0.00164 m/s
reverse fraction:  RMSE = 0.00672
filter-layer Uz:   RMSE = 0.00531 m/s
floor outflow:     RMSE = 0.00612 m3/s
```

## v0 -> v1 comparison

The main improvement is reverse fraction.

```text
up10 reverse RMSE:
  v0 = 0.01801
  v1 = 0.00357
  improvement = 80.1%

down10 reverse RMSE:
  v0 = 0.01352
  v1 = 0.00589
  improvement = 56.4%

up60 reverse RMSE:
  v0 = 0.00283
  v1 = 0.00145
  improvement = 48.8%

down60 reverse RMSE:
  v0 = 0.00743
  v1 = 0.00672
  improvement = 9.6%
```

Mean `Uz` changes were smaller:

```text
down10 work-zone Uz RMSE:
  v0 = 0.00494 m/s
  v1 = 0.00466 m/s
  improvement = 5.6%

down10 near-car Uz RMSE:
  v0 = 0.00254 m/s
  v1 = 0.00226 m/s
  improvement = 11.0%
```

The 60 s downward mean-`Uz` errors changed slightly worse because v1 is calibrated to resolve the early 0.25 s transient, while 60 s validation is coarse at 2 s output. This is acceptable for the present purpose.

## Reverse fraction shape check

### Upward high-resolution reverse fraction

```text
time    v1 pred    CFD
0.00    0.23268    0.23268
0.25    0.16334    0.15316
0.50    0.18602    0.17905
1.00    0.21487    0.20749
2.00    0.23732    0.23344
4.00    0.24029    0.24423
10.0    0.23176    0.23049
```

### Downward high-resolution reverse fraction

```text
time    v1 pred    CFD
0.00    0.22940    0.22940
0.25    0.27131    0.27375
0.50    0.25785    0.26223
1.00    0.24054    0.24502
2.00    0.22666    0.22769
4.00    0.22427    0.21812
10.0    0.23010    0.22397
```

The v1 reverse kernel now captures the direction-specific qualitative shape:

```text
upward:   early dip -> recovery/overshoot -> final relaxation
downward: early rise -> decay/undershoot -> recovery
```

## Interpretation

The current surrogate is now a practical control-oriented metric baseline:

- Through-flow/filter metrics are treated as very fast channels.
- Mean downdraft metrics are represented by fast lag plus direction-specific overshoot.
- Reverse fraction is represented by a compact empirical wake-response kernel.

This still is not a full-field surrogate and should not replace CFD for topology changes. But it is now good enough to test command histories, controller logic, and data requirements before running many more transient CFD cases.

## Limitations

1. Parameters are calibrated only on one paired step around `4.36 <-> 5.45 m/s`.
2. The reverse kernel is empirical and may not extrapolate to low-flow or high-flow regimes.
3. The steady library has a known low-flow kink around `2.72..3.27 m/s`.
4. Filter-layer endpoint offsets remain between transient and steady references, likely from solver/time-history differences.
5. Full-field phase lag is not modeled; this is metric-level only.

## Recommended next step

Use v1 to simulate a few realistic supply command histories, e.g. blower command sequences with holds and reversals. Then select the most informative additional CFD cases by disagreement/uncertainty:

```text
1. larger step: 4.36 -> 6.54 m/s
2. low-flow step: 2.18 -> 4.36 m/s
3. kink region: 2.72 <-> 3.27 m/s
```

Before more CFD, it is also reasonable to package the current transient/surrogate scripts and docs in a commit, while leaving generated case artifacts untracked.
