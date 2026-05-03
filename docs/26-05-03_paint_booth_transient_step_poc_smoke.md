# Paint-booth transient step PoC smoke result

Date: 2026-05-03

## Objective

Verify that the current L2 steady paint-booth baseline can be converted to a transient `pimpleFoam` case and advanced stably for a short step-response smoke test.

## Case

- Source steady case: `cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u4p36`
- Transient smoke case: `cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke`
- Case-prep script: `scripts/create_paint_booth_transient_poc.py`
- Initial field: latest steady result at `U_supply = 4.36 m/s`, copied into `0/`
- Transient boundary from `t = 0`: `supplyInlet U = (0 0 -5.45) m/s`
- Solver: `pimpleFoam`
- Turbulence: `kOmegaSST`

## Generated smoke settings

```text
endTime = 5 s
deltaT = 0.02 s
adjustTimeStep = yes
maxCo = 1.0
maxDeltaT = 0.1 s
writeControl = adjustableRunTime
writeInterval = 1 s
```

## Commands

Generate case:

```bash
python3 scripts/create_paint_booth_transient_poc.py \
  --source-case cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u4p36 \
  --case-dir cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke \
  --target-supply-velocity 5.45 \
  --end-time 5 \
  --delta-t 0.02 \
  --max-delta-t 0.1 \
  --write-interval 1 \
  --force
```

Run smoke test:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD/cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke" \
  openfoam-pipeline-local:latest \
  bash -lc 'checkMesh -constant > log.checkMesh 2>&1 && pimpleFoam > log.pimpleFoam 2>&1 && foamToVTK -latestTime > log.foamToVTK 2>&1'
```

Latest-time metrics:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 scripts/postprocess_paint_booth_plenum_metrics.py cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke'
```

## Verification results

### Runtime availability

`pimpleFoam` exists in the Docker image:

```text
/usr/bin/pimpleFoam
```

### Mesh

`checkMesh -constant` passed:

```text
points:           265059
faces:            712407
internal faces:   679819
cells:            224769
Mesh OK.
```

### Solver

`pimpleFoam` reached `End` at `t = 5 s`.

Final log tail indicators:

```text
Courant Number mean: 0.0825073 max: 0.992863
deltaT = 0.0117647
Time = 5
...
time step continuity errors : sum local = 1.50814e-09, global = -8.92306e-11, cumulative = -4.55735e-06
End
```

Written time directories:

```text
0, 1, 2, 3, 4, 5
```

Latest VTK export exists for `t = 5`.

### Latest-time scalar metrics at 5 s

From `post_plenum_metrics.json`:

```text
supply inflow:        9.8100000465 m3/s
relative imbalance:   0.0001185902  (~0.0119%)
work-zone Uz mean:   -0.2580942512 m/s
near-car Uz mean:    -0.1815650165 m/s
near-car reverse:     0.241592689   (~24.16%)
```

For reference, the 5.45 m/s steady case had approximately:

```text
work-zone Uz mean: -0.2547 m/s
near-car Uz mean:  -0.1837 m/s
```

So by 5 s, the transient latest-time region metrics are already close to the corresponding steady 5.45 m/s values. This suggests the response time may be short for these coarse region-averaged quantities, but this should be confirmed using time-series metrics over all written times.

## Performance observation

The 5 s smoke run took about:

```text
ExecutionTime ≈ 471 s
```

With the observed adaptive `deltaT ≈ 0.01176 s`, a 60 s case may take roughly 1.5 hours on this WSL/Docker setup. Long transient runs should therefore be launched as background jobs with completion notification.

## Decision

The transient PoC smoke test passes.

The next appropriate step is a longer `60 s` step-response run with the same setup and `writeInterval = 2 s`, followed by time-series postprocessing over written times to estimate response time metrics.
