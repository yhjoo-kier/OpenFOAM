# Paint-booth downward high-resolution transient and surrogate check

## Objective

Run the paired high-resolution downward step response:

```text
5.45 -> 4.36 m/s
endTime = 10 s
writeInterval = 0.25 s
```

This complements the upward high-resolution case and checks whether early downward overshoot/undershoot was hidden by the previous 2 s output interval.

## Case

```text
case: cases/paint_booth_l2_transient_step_poc/u5p45_to_u4p36_10s_hi_res_parallel
source steady case: cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u5p45
source latest time: 200
solver: pimpleFoam
initial condition: latest steady fields copied to 0/
target supply velocity: 4.36 m/s
endTime: 10 s
deltaT: 0.02 s
adjustTimeStep: yes
maxCo: 1.0
maxDeltaT: 0.1 s
writeInterval: 0.25 s
```

## Parallel execution

The first 4-rank attempt used:

```foam
method scotch;
```

but failed at `decomposePar` because this Docker image does not have the scotch decomposition library loaded:

```text
Attempted to use <scotch> without the scotchDecomp library loaded.
Please install <scotch> and ensure libscotch.so is in LD_LIBRARY_PATH.
```

The decomposition was changed to `hierarchical`:

```foam
numberOfSubdomains 4;
method hierarchical;

hierarchicalCoeffs
{
    n     (2 2 1);
    delta 0.001;
    order xyz;
}
```

Dry-run decomposition quality:

```text
Finished decomposition into 4 domains
Cells min: 56192
Cells max: 56193
median: 56192
balance: 100.002%
```

Final run sequence:

```text
checkMesh -constant
decomposePar -force
mpirun --allow-run-as-root -np 4 pimpleFoam -parallel
reconstructPar
foamToVTK
```

Wrapper log:

```text
logs/26-05-03_paint_booth_transient_u5p45_to_u4p36_10s_hi_res_parallel_retry_hierarchical.log
```

Runtime:

```text
start: Sun May  3 22:30:28 KST 2026
end:   Sun May  3 22:38:05 KST 2026
wall:  about 457 s including checkMesh/decompose/solve/reconstruct/VTK
```

The parallel `pimpleFoam` log itself reported:

```text
ExecutionTime = 329.45 s
ClockTime = 330 s
```

For comparison, the upward 10 s high-resolution serial run took about `735 s` solver clock time. The 4-rank hierarchical run therefore reduced solver wall time by about `2.2x`, despite reconstruction/VTK overhead.

## Execution verification

Mesh:

```text
points: 265059
faces: 712407
internal faces: 679819
cells: 224769
Mesh OK.
```

Solver final status:

```text
Time = 10
End
Courant mean = 0.0823466
Courant max = 0.999404
deltaT = 0.0147059
cumulative continuity = 3.90494e-06
```

Reconstruction:

```text
reconstructPar completed through Time = 10
End
```

VTK/time-series:

```text
n_times = 41
time range = 0..10 s
output times = 0, 0.25, 0.5, ..., 9.75, 10
```

Processor directories retained as generated artifacts:

```text
processor0, processor1, processor2, processor3
```

## Latest metrics at 10 s

```text
Supply inflow:             7.8480005522 m3/s
Floor outflow:             7.8492905568 m3/s
Relative mass imbalance:   1.6437e-4
Work-zone Uz mean:        -0.2033843398 m/s
Near-car Uz mean:         -0.1423479170 m/s
Near-car reverse fraction: 0.2239687737
Filter-layer Uz mean:     -0.3232219219 m/s
Filter pressure-drop proxy: 55.6286432743 m2/s2
```

## Response summary from 0.25 s output

```text
floor_outflow_m3s:
  initial = 9.8036927569
  final_observed = 7.8492905568
  t63 = 0.158 s
  t90 = 0.225 s
  t95 = 0.238 s
  settling_time_5pct = 0.25 s

work_zone_Uz_mean:
  initial = -0.2547136545
  final_observed = -0.2033843398
  t63 = 0.138 s
  t90 = 0.197 s
  t95 = 0.208 s
  settling_time_10pct = 0.5 s
  settling_time_5pct = 0.75 s
  overshoot_beyond_final_fraction = 14.4%

near_car_Uz_mean:
  initial = -0.1836552769
  final_observed = -0.1423479170
  t63 = 0.138 s
  t90 = 0.196 s
  t95 = 0.207 s
  settling_time_10pct = 0.5 s
  settling_time_5pct = 4.75 s
  overshoot_beyond_final_fraction = 14.6%

near_car_reverse_fraction_Uz_positive:
  initial = 0.2293972633
  final_observed = 0.2239687737
  t63 = 2.167 s
  t90 = 2.330 s
  t95 = 2.366 s
  settling_time_10pct = 9.5 s
  settling_time_5pct = 9.75 s
  overshoot_beyond_final_fraction = 113.8%
  undershoot_before_initial_fraction = 817.1%

filter_layer_panel_Uz_mean:
  initial = -0.4107867181
  final_observed = -0.3232219219
  t63 = 0.159 s
  t90 = 0.226 s
  t95 = 0.239 s
  settling_time_5pct = 0.25 s
```

For reverse fraction, normalized overshoot/undershoot percentages are large because the net initial-final change is small while the transient excursion is large. The qualitative conclusion is the important part: reverse/wake behavior is strongly non-monotonic and not first-order.

## Selected time samples

```text
t=0.00:
  work Uz = -0.25471
  near-car Uz = -0.18366
  reverse = 0.22940
  filter Uz = -0.41079
  floor outflow = 9.80369

t=0.25:
  work Uz = -0.19601
  near-car Uz = -0.13632
  reverse = 0.27375
  filter Uz = -0.32365
  floor outflow = 7.85019

t=0.50:
  work Uz = -0.19919
  near-car Uz = -0.13928
  reverse = 0.26223
  filter Uz = -0.32357
  floor outflow = 7.85013

t=1.00:
  work Uz = -0.20245
  near-car Uz = -0.14232
  reverse = 0.24502
  filter Uz = -0.32341
  floor outflow = 7.85001

t=2.00:
  work Uz = -0.20506
  near-car Uz = -0.14467
  reverse = 0.22769
  filter Uz = -0.32329
  floor outflow = 7.84967

t=4.00:
  work Uz = -0.20548
  near-car Uz = -0.14482
  reverse = 0.21812
  filter Uz = -0.32325
  floor outflow = 7.84930

t=10.0:
  work Uz = -0.20338
  near-car Uz = -0.14235
  reverse = 0.22397
  filter Uz = -0.32322
  floor outflow = 7.84929
```

## Surrogate v0 comparison

Prediction files:

```text
cases/paint_booth_l2_dynamic_metric_surrogate_v0/u5p45_to_u4p36_10s_hi_res_pred.csv
cases/paint_booth_l2_dynamic_metric_surrogate_v0/u5p45_to_u4p36_10s_hi_res_pred_summary.json
```

Error against the high-resolution downward CFD:

```text
floor_outflow_m3s:
  RMSE = 0.06456 m3/s
  MAE  = 0.01762 m3/s
  max  = 0.40386 m3/s

work_zone_Uz_mean:
  RMSE = 0.00494 m/s
  MAE  = 0.00446 m/s
  max  = 0.01546 m/s

near_car_Uz_mean:
  RMSE = 0.00254 m/s
  MAE  = 0.00131 m/s
  max  = 0.01391 m/s

near_car_reverse_fraction_Uz_positive:
  RMSE = 0.01352
  MAE  = 0.01134
  max  = 0.04425

filter_layer_panel_Uz_mean:
  RMSE = 0.00639 m/s
  MAE  = 0.00576 m/s
  max  = 0.02220 m/s
```

## Interpretation

The paired high-resolution runs now show direction-specific early transients:

### Upward `4.36 -> 5.45`

```text
work-zone mean Uz:
  first t95 ~0.215 s
  5% settling ~4.25 s
  overshoot ~10.4%

near-car mean Uz:
  first t95 ~0.177 s
  5% settling ~5.5 s
  overshoot ~34.0%
```

### Downward `5.45 -> 4.36`

```text
work-zone mean Uz:
  first t95 ~0.208 s
  5% settling ~0.75 s
  overshoot ~14.4%

near-car mean Uz:
  first t95 ~0.207 s
  5% settling ~4.75 s
  overshoot ~14.6%
```

Both directions have sub-second initial response, but the overshoot magnitude and relaxation differ by region and direction.

The reverse-flow fraction is clearly not a simple lag:

- upward high-res had a strong early dip then recovery;
- downward high-res has a strong early rise, then decay below final, then recovery;
- 10 s is barely enough for this metric to settle within 5% of the 10 s observed final value;
- normalized overshoot/undershoot is sensitive because the initial-final difference is small.

## Recommended surrogate v1 update

1. Add direction-specific overshoot parameters for mean `Uz` channels:

```text
work_zone_Uz_mean:
  upward beta ~0.10, decay ~1.5 s
  downward beta ~0.14, faster decay ~0.5..1.0 s

near_car_Uz_mean:
  upward beta ~0.34, decay ~1.8 s
  downward beta ~0.15, decay ~2..4 s
```

2. Replace reverse fraction slow first-order lag with a two-state or empirical kernel model:

```text
reverse_fraction = slow baseline + fast signed excursion + slower recovery
```

3. Keep filter/floor through-flow as a fast channel:

```text
tau ~0.16 s
settling <=0.25 s at this output resolution
```

4. For future larger transient CFD, use 4-rank hierarchical decomposition unless the Docker image is rebuilt with scotch support.
