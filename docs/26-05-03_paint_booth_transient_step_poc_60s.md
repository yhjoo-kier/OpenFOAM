# Paint-booth L2 transient step-response PoC: 60 s run

## Objective

Validate a first longer transient paint-booth step-response case after the 5 s `pimpleFoam` smoke test.

The case reuses the converged steady L2 field at `U_supply = 4.36 m/s` as the transient initial condition and applies a step change of the supply inlet boundary to `U_supply = 5.45 m/s` from `t = 0`.

## Case

- Source steady case: `cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u4p36`
- Source latest time copied to transient `0/`: `200`
- Transient case: `cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s`
- Solver: `pimpleFoam`
- Turbulence model: `kOmegaSST`
- Mesh cells: `224,769`
- Step input:
  - initial field: steady `U_supply = 4.36 m/s`
  - transient boundary target: `supplyInlet U = (0 0 -5.45) m/s`

## Runtime settings

- `endTime = 60 s`
- `deltaT = 0.02 s`
- `adjustTimeStep = yes`
- `maxCo = 1.0`
- `maxDeltaT = 0.1 s`
- `writeControl = adjustableRunTime`
- `writeInterval = 2 s`

## Commands

Case generation:

```bash
python3 scripts/create_paint_booth_transient_poc.py \
  --source-case cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u4p36 \
  --case-dir cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s \
  --target-supply-velocity 5.45 \
  --end-time 60 \
  --delta-t 0.02 \
  --max-delta-t 0.1 \
  --write-interval 2 \
  --force
```

Run:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD/cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s" \
  openfoam-pipeline-local:latest \
  bash -lc 'checkMesh -constant > log.checkMesh 2>&1 && pimpleFoam > log.pimpleFoam 2>&1 && foamToVTK -latestTime > log.foamToVTK 2>&1'
```

Latest-time metric postprocess:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 scripts/postprocess_paint_booth_plenum_metrics.py cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s > cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s/log.postprocess 2>&1'
```

## Solver verification

`checkMesh -constant`:

```text
points: 265059
faces: 712407
internal faces: 679819
cells: 224769
Mesh OK.
```

Time directories written:

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
Courant Number mean: 0.0828168 max: 0.997358
deltaT = 0.0117647
ExecutionTime = 4861.65 s  ClockTime = 4863 s
```

Final continuity:

```text
time step continuity errors:
  sum local  = 8.07759e-08
  global     = 5.31194e-10
  cumulative = -4.46874e-06
```

`foamToVTK -latestTime` exported the final time:

```text
Create mesh for time = 60
Internal: VTK/u4p36_to_u5p45_60s_5093/internal.vtu
Boundary: supplyInlet, plenumTopWall, floorExhaust, frontWall, rearWall, sideWalls, carBody
End
```

## Latest-time metrics at 60 s

Output file:

```text
cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_60s/post_plenum_metrics.json
```

Key scalar metrics:

- Supply inflow: `9.8100000465 m3/s`
- Floor outflow: `9.8113594594 m3/s`
- Relative mass imbalance: `1.3857e-4` (`0.0139%`)
- Work-zone `Uz_mean`: `-0.2561563551 m/s`
- Near-car `Uz_mean`: `-0.1791233122 m/s`
- Near-car reverse-flow fraction (`Uz > 0`): `0.2306298791`
- Filter-layer panel `Uz_mean`: `-0.4038849771 m/s`
- Filter pressure-drop proxy: `86.8264074326 m2/s2`

## Comparison against reference cases

Reference values from the existing steady L2 U-sweep:

### Initial steady reference: `U_supply = 4.36 m/s`

- Supply inflow: `7.8480005522 m3/s`
- Work-zone `Uz_mean`: `-0.2000069320 m/s`
- Near-car `Uz_mean`: `-0.1434996426 m/s`
- Near-car reverse-flow fraction: `0.2326842387`
- Filter-layer panel `Uz_mean`: `-0.3286292255 m/s`

### Target steady reference: `U_supply = 5.45 m/s`

- Supply inflow: `9.8100000465 m3/s`
- Work-zone `Uz_mean`: `-0.2547136545 m/s`
- Near-car `Uz_mean`: `-0.1836552769 m/s`
- Near-car reverse-flow fraction: `0.2293972633`
- Filter-layer panel `Uz_mean`: `-0.4107867181 m/s`

### 60 s transient latest-time interpretation

At `t = 60 s`, the global and regional scalar metrics are close to the target `U_supply = 5.45 m/s` steady reference:

- Supply inflow exactly matches the imposed target inlet flow.
- Work-zone `Uz_mean` is slightly more negative than the steady target by about `0.00144 m/s`.
- Near-car `Uz_mean` is slightly less negative than the steady target by about `0.00453 m/s`.
- Near-car reverse-flow fraction is close to the steady target (`0.2306` vs `0.2294`).
- Filter-layer `Uz_mean` is close but slightly less negative than the steady target (`-0.4039` vs `-0.4108 m/s`).

This supports the conclusion that the 60 s step-response run is numerically stable and has moved close to the new steady operating point in the currently tracked scalar regions.

## Caveats

- This document reports latest-time metrics only. The case wrote fields every 2 s, so the next required step is time-series postprocessing over all written times.
- `foamToVTK -latestTime` exported only the final time for this run. A time-series exporter should either:
  - run `foamToVTK` for selected written times, or
  - parse OpenFOAM fields directly and compute metrics without converting every time to VTK.
- The first transient result establishes numerical feasibility; it is not yet a full transient dataset specification.

## Recommended next steps

1. Implement a transient time-series metric exporter:
   - supply/floor flow rates vs time
   - work-zone `Uz_mean(t)`
   - near-car `Uz_mean(t)`
   - reverse-flow fraction vs time
   - response-time estimates: 63%, 90%, 95% of final change
2. Decide transient dataset settings from the response curve:
   - write interval
   - history-window length
   - target prediction horizon
3. Prepare a closed-network HPC handoff package:
   - OpenFOAM case tarball
   - Apptainer/Singularity-first run scripts
   - Slurm/PBS templates
   - postprocess scripts
   - verification checklist
4. Only after time-series validation, expand to additional step/ramp/sinusoidal operating profiles.
