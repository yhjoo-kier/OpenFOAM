# Paint-booth transient high-resolution initial response

## Objective

Refine the early-time response estimate for the `4.36 -> 5.45 m/s` supply-velocity step. The previous 60 s case used `writeInterval = 2 s`, which under-resolved the fast mean-flow response. This case uses a 10 s window with 0.25 s output spacing.

## Case

```text
case: cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_10s_hi_res
source steady case: cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u4p36
source latest time: 200
solver: pimpleFoam
initial condition: latest steady fields copied to 0/
target supply velocity: 5.45 m/s
endTime: 10 s
deltaT: 0.02 s
adjustTimeStep: yes
maxCo: 1.0
maxDeltaT: 0.1 s
writeInterval: 0.25 s
```

## Execution verification

Run sequence:

```text
checkMesh -constant
pimpleFoam
foamToVTK
```

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
ExecutionTime = 734.72 s
ClockTime = 735 s
Courant mean = 0.0838617
Courant max = 1.00466
deltaT = 0.0119048
cumulative continuity = -5.02976e-06
```

VTK/time-series:

```text
n_times = 41
time range = 0..10 s
output times = 0, 0.25, 0.5, ..., 9.75, 10
```

## Latest metrics at 10 s

```text
Supply inflow:             9.8100000465 m3/s
Floor outflow:             9.8114827460 m3/s
Relative mass imbalance:   1.5114e-4
Work-zone Uz mean:        -0.2561351359 m/s
Near-car Uz mean:         -0.1792134792 m/s
Near-car reverse fraction: 0.2304866965
Filter-layer Uz mean:     -0.4038751423 m/s
Filter pressure-drop proxy: 86.8039348125 m2/s2
```

## Response summary from 0.25 s output

```text
work_zone_Uz_mean:
  initial = -0.2000069320
  final_observed = -0.2561351359
  t63 = 0.143 s
  t90 = 0.204 s
  t95 = 0.215 s
  settling_time_10pct = 0.5 s
  settling_time_5pct = 4.25 s
  overshoot = 10.4%

near_car_Uz_mean:
  initial = -0.1434996426
  final_observed = -0.1792134792
  t63 = 0.118 s
  t90 = 0.168 s
  t95 = 0.177 s
  settling_time_10pct = 4.25 s
  settling_time_5pct = 5.5 s
  overshoot = 34.0%

filter_layer_panel_Uz_mean:
  initial = -0.3286292255
  final_observed = -0.4038751423
  t63 = 0.160 s
  t90 = 0.227 s
  t95 = 0.240 s
  settling_time_10pct = 0.25 s
  settling_time_5pct = 0.25 s
  overshoot = 0%

floor_outflow_m3s:
  initial = 7.8430873630
  final_observed = 9.8114827460
  t63 = 0.158 s
  t90 = 0.225 s
  t95 = 0.238 s
  settling_time_10pct = 0.25 s
  settling_time_5pct = 0.25 s

near_car_reverse_fraction_Uz_positive:
  initial = 0.2326842387
  final_observed = 0.2304866965
  t63/t90/t95 are not physically meaningful here because the net initial-final change is very small and the signal is non-monotonic.
  settling_time_10pct = 9.75 s
  settling_time_5pct = 10.0 s within this 10 s window
```

Selected time samples:

```text
t=0.00: work Uz=-0.2000, near-car Uz=-0.1435, reverse=0.2327, filter Uz=-0.3286
t=0.25: work Uz=-0.2620, near-car Uz=-0.1914, reverse=0.1532, filter Uz=-0.4032
t=0.50: work Uz=-0.2607, near-car Uz=-0.1895, reverse=0.1791, filter Uz=-0.4034
t=1.00: work Uz=-0.2600, near-car Uz=-0.1874, reverse=0.2075, filter Uz=-0.4037
t=2.00: work Uz=-0.2599, near-car Uz=-0.1856, reverse=0.2334, filter Uz=-0.4038
t=4.00: work Uz=-0.2591, near-car Uz=-0.1830, reverse=0.2442, filter Uz=-0.4039
t=5.50: work Uz=-0.2576, near-car Uz=-0.1810, reverse=0.2394, filter Uz=-0.4039
t=10.0: work Uz=-0.2561, near-car Uz=-0.1792, reverse=0.2305, filter Uz=-0.4039
```

## Interpretation

The high-resolution run shows that the first approach to the new flow-rate condition is much faster than the 2 s output case suggested:

- Mass/flow-through and filter-layer mean response are essentially complete by `0.25 s` at this output resolution.
- Work-zone and near-car mean `Uz` reach their first 95% crossing around `0.18..0.22 s`, but then relax back from an overshoot, so practical 5% settling remains several seconds.
- Near-car mean `Uz` has a stronger overshoot than the 2 s run indicated, because the peak occurs at `t = 0.25 s` and was missed by the coarser output.
- Reverse-flow fraction remains a wake-sensitive/non-monotonic signal. A first-order lag fit is not reliable for this metric.

For a dynamic quasi-steady surrogate, split the dynamics into at least two channels:

```text
fast through-flow/filter channel:
  characteristic response <~0.25 s at current resolution

mean work-zone / near-car downdraft channel:
  first crossing ~0.1..0.25 s, but overshoot/relaxation settling ~4..6 s

wake/reverse-flow channel:
  non-monotonic and slower; retain separate empirical response/settling model
```

## CPU/parallel note

The run above was serial. WSL and Docker expose 14 logical CPUs, so the one-core behavior is not a WSL resource-allocation limit. It is because OpenFOAM solvers run serial unless the case is decomposed and launched with MPI, e.g.:

```bash
decomposePar -force
mpirun --allow-run-as-root -np 4 pimpleFoam -parallel
reconstructPar
foamToVTK
```

Future larger transient cases should use 4 MPI ranks unless serial reproducibility is specifically desired.
