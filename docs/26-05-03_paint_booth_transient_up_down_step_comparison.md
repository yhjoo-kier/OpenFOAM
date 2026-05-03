# Paint-booth transient up/down step-response comparison

## Objective

Compare the first upward and downward supply-velocity step-response cases to check whether the dynamic quasi-steady surrogate can use symmetric response times or needs direction-dependent dynamics.

Cases:

- Upward step: `4.36 -> 5.45 m/s`, case `u4p36_to_u5p45_60s`
- Downward step: `5.45 -> 4.36 m/s`, case `u5p45_to_u4p36_60s`

Both cases use:

- steady latest-time field at `t = 200` copied to transient `0/`
- `pimpleFoam`
- `endTime = 60 s`
- `writeInterval = 2 s`
- `maxCo = 1.0`

## Downward step solver verification

Case:

```text
cases/paint_booth_l2_transient_step_poc/u5p45_to_u4p36_60s
```

`checkMesh -constant`:

```text
points: 265059
faces: 712407
internal faces: 679819
cells: 224769
Mesh OK.
```

Written time directories:

```text
0, 2, 4, ..., 58, 60
```

`pimpleFoam` reached final time:

```text
Time = 60
End
```

Final solver state:

```text
Courant Number mean: 0.0818932
Courant Number max:  0.999621--0.999623
deltaT:              0.0145985
ExecutionTime:       3361.68 s
ClockTime:           3362 s
```

Final continuity:

```text
time step continuity errors:
  sum local  = 7.82551e-08
  global     = 2.32549e-10
  cumulative = 7.71072e-06
```

`foamToVTK -latestTime` succeeded for `t = 60`.

## Downward latest-time metrics at 60 s

Output:

```text
cases/paint_booth_l2_transient_step_poc/u5p45_to_u4p36_60s/post_plenum_metrics.json
```

Key values:

```text
Supply inflow:              7.8480005522 m3/s
Floor outflow:              7.8492978178 m3/s
Relative mass imbalance:    1.6530e-4  (~0.0165%)
Work-zone Uz mean:         -0.2030571401 m/s
Near-car Uz mean:          -0.1418604702 m/s
Near-car reverse fraction:  0.2253632481
Filter-layer Uz mean:      -0.3232269287 m/s
```

These values are close to the `4.36 m/s` steady reference for global flow and mean downdraft metrics, but the reverse-flow fraction differs more noticeably from the initial/target steady values and should be treated as a wake-sensitive metric.

## Time-series response extraction

Time-series outputs for the downward case:

```text
cases/paint_booth_l2_transient_step_poc/u5p45_to_u4p36_60s/post_transient_metrics_timeseries.csv
cases/paint_booth_l2_transient_step_poc/u5p45_to_u4p36_60s/post_transient_metrics_timeseries.json
cases/paint_booth_l2_transient_step_poc/u5p45_to_u4p36_60s/post_transient_metrics_timeseries_summary.json
```

The same script was used for the upward case:

```text
scripts/postprocess_paint_booth_transient_timeseries.py
```

The response fraction is computed as:

```text
r(t) = [y(t) - y(0)] / [y(60) - y(0)]
```

and the script estimates `t63`, `t90`, `t95`, 10%/5% settling time, and overshoot relative to the observed final value.

Because the write interval is `2 s`, sub-2 s response times are linear interpolations between `t = 0` and `t = 2 s`; they should be treated as coarse estimates.

## Upward vs downward response summary

### Work-zone `Uz_mean`

Upward `4.36 -> 5.45`:

```text
t63 = 1.18 s
t90 = 1.69 s
t95 = 1.78 s
5% settling = 6 s
overshoot beyond final = 6.7% of observed delta
```

Downward `5.45 -> 4.36`:

```text
t63 = 1.32 s
t90 = 1.87 s
t95 = 1.98 s
5% settling = 2 s
overshoot beyond final ≈ 0%
```

Interpretation:

- Mean work-zone downdraft responds quickly in both directions.
- Upward response has small overshoot and longer 5% settling.
- Downward response is more monotonic for this metric.

### Near-car `Uz_mean`

Upward `4.36 -> 5.45`:

```text
t63 = 1.07 s
t90 = 1.52 s
t95 = 1.61 s
5% settling = 6 s
overshoot beyond final = 18.1% of observed delta
```

Downward `5.45 -> 4.36`:

```text
t63 = 1.36 s
t90 = 1.93 s
t95 = 5.49 s
5% settling = 6 s
overshoot beyond final ≈ 0%
```

Interpretation:

- Near-car mean velocity is direction-dependent.
- Upward step reaches the final mean quickly but overshoots strongly.
- Downward step has a slower approach to 95%, but little overshoot.

### Filter-layer panel `Uz_mean`

Upward `4.36 -> 5.45`:

```text
t63 = 1.27 s
t90 = 1.80 s
t95 = 1.90 s
5% settling = 2 s
overshoot beyond final = 0%
```

Downward `5.45 -> 4.36`:

```text
t63 = 1.27 s
t90 = 1.80 s
t95 = 1.90 s
5% settling = 2 s
overshoot beyond final ≈ 0%
```

Interpretation:

- Filter-layer mean response is nearly symmetric at the current sampling resolution.
- The true time constant may be below the current 2 s write interval.

### Near-car reverse-flow fraction

Upward `4.36 -> 5.45`:

```text
t63 = 9.17 s
t90 = 9.67 s
t95 = 9.76 s
5% settling = 26 s
overshoot beyond final = 43.3% of observed delta
```

Downward `5.45 -> 4.36`:

```text
t63 = 2.18 s
t90 = 2.41 s
t95 = 2.45 s
5% settling = 30 s
overshoot beyond final = 180.2% of observed delta
```

Interpretation:

- Reverse-flow fraction is clearly not a simple first-order metric.
- The net initial-to-final change is small, so normalized overshoot is sensitive.
- Both directions require much longer settling (`~26--30 s`) than mean velocity metrics.
- This metric likely reflects wake/recirculation topology adjustment and should receive a separate slower/nonlinear model in the dynamic surrogate.

## Main conclusion

For the current paint-booth L2 case and the `4.36 <-> 5.45 m/s` step pair:

1. Mean downdraft metrics are fast:

```text
work-zone / filter-layer / near-car mean Uz: O(1--6 s)
```

2. Wake-sensitive reverse-flow fraction is slow and non-monotonic:

```text
near-car reverse fraction: O(20--30 s) settling
```

3. Upward and downward steps are not perfectly symmetric:

```text
upward: stronger overshoot in near-car mean Uz
downward: near-car mean reaches 95% more slowly but with little overshoot
reverse-flow: non-monotonic in both directions
```

4. A single global first-order lag is not sufficient for all outputs.

Recommended first surrogate structure:

```text
mean downdraft metrics:
  direction-dependent first-order lag, tau ≈ 1--2 s, plus optional overshoot correction

filter-layer mean:
  very fast lag, tau <~ 2 s; high-resolution short run needed to resolve

reverse-flow / wake-sensitive metric:
  separate slow settling model, nonlinear/second-order or empirical step-response kernel
```

## Implication for next simulations

The immediate simulation priority is no longer simply "does transient run?". It is to characterize response kernels over more of the flow range.

Recommended next cases:

```text
2.18 -> 4.36 m/s  low-to-nominal upward step
4.36 -> 6.54 m/s  nominal-to-high upward step
6.54 -> 4.36 m/s  high-to-nominal downward step
```

Also useful:

```text
4.36 -> 5.45 m/s short high-resolution run
endTime = 10 s, writeInterval = 0.25--0.5 s
```

Purpose of the high-resolution run:

- resolve true time constants for mean metrics that currently cross 63--95% before the first written `t = 2 s` point.

## Files generated

Downward case:

```text
cases/paint_booth_l2_transient_step_poc/u5p45_to_u4p36_60s/
```

Main logs:

```text
log.checkMesh
log.pimpleFoam
log.foamToVTK
log.foamToVTK_all_times
log.postprocess
log.postprocess_timeseries
```

Metric outputs:

```text
post_plenum_metrics.json
post_transient_metrics_timeseries.csv
post_transient_metrics_timeseries.json
post_transient_metrics_timeseries_summary.json
```
