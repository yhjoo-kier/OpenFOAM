# Paint Booth Plenum-Filter Case 생성 및 1차 검증

Date: 2026-05-02 KST

## 1. 목적

사용자 피드백을 반영하여 단순 ceiling uniform inlet 대신, 상부 duct/supply jet가 plenum으로 들어오고 filter/저항층을 통과해 작업공간으로 내려오는 구조를 1차 구현했다.

핵심 변경:

- 기존 높이 3.0 m domain을 3.8 m로 확장
- z = 3.05 ~ 3.8 m 상부 plenum 추가
- z = 2.95 ~ 3.05 m thin porous filter layer 추가
- top central supply opening에서 고속 jet 유입
- floor 전체 pressure outlet
- smooth single-shell car body 유지

## 2. 생성 파일

- Generator: `scripts/create_paint_booth_plenum_filter_case.py`
- Case: `cases/paint_booth_plenum_filter_full_floor_v035/`
- Geometry render: `docs/figures/26-05-02_paint_booth_plenum_filter_geometry.png`
- Velocity mid-plane render: `docs/figures/26-05-02_paint_booth_plenum_filter_velocity_midplane.png`
- Metric JSON: `cases/paint_booth_plenum_filter_full_floor_v035/post_plenum_metrics.json`

## 3. Geometry / domain

### Domain

- x: -1.5 ~ 6.5 m
- y: -2.0 ~ 2.0 m
- z: 0.0 ~ 3.8 m
- Overall size: 8.0 × 4.0 × 3.8 m

### Work zone

- z: 0.0 ~ 2.95 m

### Filter layer

- z: 2.95 ~ 3.05 m
- Thickness: 0.10 m
- Current 1st implementation: full cross-section porous layer
- cellZone name: `filterZone`
- Selected cells after `topoSet`: 1,088 cells

### Plenum

- z: 3.05 ~ 3.8 m
- Height: 0.75 m

### Supply opening

- top central opening at z = 3.8 m
- x: 1.75 ~ 3.25 m
- y: -0.60 ~ 0.60 m
- Area: 1.8 m²
- U = (0 0 -4.36) m/s
- Estimated flow rate: 7.848 m³/s

This corresponds to approximately 0.35 m/s over a 22.4 m² reference filter panel, but the current implementation uses a full cross-section porous layer for numerical robustness.

## 4. Boundary conditions

### `supplyInlet`

- OpenFOAM patch type: `patch`
- `U`: `fixedValue uniform (0 0 -4.36)`
- `p`: `zeroGradient`
- `k`, `omega`: fixedValue based on 5% turbulence intensity and 0.15 m length scale

### `floorExhaust`

- OpenFOAM patch type: `patch`
- `U`: `pressureInletOutletVelocity`
- `p`: `fixedValue uniform 0`

### Walls

- `plenumTopWall`
- `frontWall`
- `rearWall`
- `sideWalls`
- `carBody`

Boundary condition:

- `U`: `noSlip`
- `p`: `zeroGradient`
- `k`: `kqRWallFunction`
- `omega`: `omegaWallFunction`
- `nut`: `nutkWallFunction`

## 5. Porous filter model

Filter layer is represented with OpenFOAM `fvOptions`:

- type: `explicitPorositySource`
- model: `DarcyForchheimer`
- selection: `cellZone filterZone`
- Darcy coefficient: d = (0 0 0)
- Forchheimer coefficient: f = (6800 6800 6800)

This coefficient corresponds to the earlier pressure-drop scale estimate for about 50 Pa at 0.35 m/s across a 0.10 m layer, interpreted through a Forchheimer-type loss model.

OpenFOAM log confirmed fvOptions activation:

```text
Creating finite volume options from "system/fvOptions"
Selecting finite volume options type explicitPorositySource
Source: filterPorosity
selecting cells using cellZone filterZone
Porosity region filterPorosity:
selecting model: DarcyForchheimer
creating porous zone: filterZone
```

## 6. Mesh and solver validation

Validation command sequence:

```bash
surfaceCheck constant/triSurface/simplified_car_shell.stl
blockMesh
snappyHexMesh -overwrite
checkMesh
topoSet
simpleFoam
foamToVTK -latestTime
```

### surfaceCheck

- Surface closed: yes
- All edges connected to two faces: yes
- Number of unconnected parts: 1
- Number of zones with consistent normal: 1

### checkMesh

- Points: 32,373
- Faces: 77,551
- Internal faces: 68,113
- Cells: 23,046
- Domain bounding box: `(-1.5 -2 0)` to `(6.5 2 3.8)`
- Max aspect ratio: 5.00002
- Max non-orthogonality: 45.1403
- Average non-orthogonality: 11.8256
- Max skewness: 2.63588
- Result: `Mesh OK.`

### topoSet

- Created cellZoneSet: `filterZone`
- `filterZone` cell count: 1,088

### simpleFoam

- Completed to time step 200
- Exit code: 0
- Final cumulative continuity error: -0.0293959
- Result: `End`

The first plenum-filter case is therefore numerically runnable in the current Docker/OpenFOAM environment.

## 7. First flow-field metrics

Metrics were computed from `VTK/paint_booth_plenum_filter_full_floor_v035_200/internal.vtu`.

### Near-car region

Region:

- x: -0.25 ~ 4.75 m
- y: -1.25 ~ 1.25 m
- z: 0.10 ~ 2.00 m

Results:

- Cells: 13,866
- Mean |U|: 0.272 m/s
- 95th percentile |U|: 0.427 m/s
- Mean Uz: -0.220 m/s
- Uz std: 0.132 m/s
- Uz min/max: -0.450 / 0.034 m/s
- Reverse-flow fraction, Uz > 0: 5.94%
- Low-speed fraction, |U| < 0.05 m/s: 10.76%

### Roof region

- Mean |U|: 0.306 m/s
- Mean Uz: -0.216 m/s
- Reverse-flow fraction: 0.0%
- Low-speed fraction: 0.22%

### Filter-below region

Region:

- z: 2.75 ~ 2.95 m

Results:

- Mean |U|: 0.457 m/s
- Mean Uz: -0.457 m/s
- Uz std: 0.042 m/s
- Reverse-flow fraction: 0.0%

### Plenum region

Region:

- z: 3.05 ~ 3.75 m

Results:

- Mean |U|: 1.198 m/s
- 95th percentile |U|: 3.205 m/s
- Mean Uz: -0.295 m/s
- Uz std: 0.835 m/s
- Uz min/max: -4.107 / 1.140 m/s
- Reverse-flow fraction: 38.1%

Interpretation:

- Plenum contains strong jet/recirculating motion, as intended.
- Work zone below the filter is much lower-speed and primarily downward.
- Near-car downdraft is lower than the direct filter-below value due to car blockage and full-floor exhaust redistribution.

## 8. Important limitations of the current 1st implementation

This is a successful **numerical prototype**, not yet a final physical filter model.

Current simplifications:

1. Filter is modeled as a full cross-section thin porous layer, not yet a central 7.0 × 3.2 m panel with solid perimeter frame.
2. The porous coefficient is an engineering first guess; it should be calibrated with target pressure drop/flow rate.
3. The current mesh is still coarse around the 0.10 m filter layer.
4. `simpleFoam` convergence completed, but cumulative continuity error is larger than the earlier skeleton case.
5. No mass-flow integral check has yet been added through supply, filter, and floor patches.
6. The floor exhaust is still full-floor, not grating/slot based.

## 9. Recommended next work

### Immediate next step

Add functionObjects or post-processing scripts for:

- supply mass flow rate
- floor exhaust mass flow rate
- filter pressure drop
- filter-normal velocity uniformity
- work-zone Uz uniformity

### Model refinement

Then improve the filter geometry:

1. Split z = 3.0 plane into central filter panel and solid ceiling frame.
2. Use either:
   - central porous volume only plus blocked perimeter, or
   - baffle/porous jump patch if practical.
3. Run resistance sweep:
   - f or equivalent Δp: 30, 50, 80, 120 Pa target
4. Compare against the idealized direct downdraft case.

## 10. Conclusion

A first plenum-filter paint booth case was implemented and validated:

- `cases/paint_booth_plenum_filter_full_floor_v035/`
- Total domain height expanded to 3.8 m
- Top central supply jet introduced
- Thin porous filter layer added via `fvOptions`
- Full floor exhaust used
- `blockMesh`, `snappyHexMesh`, `checkMesh`, `topoSet`, `simpleFoam`, and `foamToVTK` completed successfully

This confirms that the project can now move from a simple ceiling-inlet idealization toward a more realistic **supply jet → plenum → filter resistance → downdraft work zone** modeling workflow.
