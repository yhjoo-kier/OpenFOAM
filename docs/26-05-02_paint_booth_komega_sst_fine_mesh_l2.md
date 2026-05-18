# Paint Booth CFD: k-omega SST Fine/L2 경계층 메쉬 결과

## 목적

L1에서 `kOmegaSST` + carBody prism layer가 안정적으로 동작하고 carBody y+가 대체로 1 근처로 맞춰짐을 확인했다. 이번 단계에서는 사용자의 요청에 따라 유속/유량 조건 재보정보다는 **fine mesh/L2 해상도 상승**을 진행했다.

현재 case는 아직 실제 도장부스 유량 정보가 확정되지 않은 상태의 CFD prototype이므로, 목적은 calibration이 아니라 다음과 같다.

1. 차체 주변 공간 해상도 상승
2. carBody surface refinement 상승
3. k-omega SST에 맞는 near-wall y+ 유지
4. mesh/solver 안정성 확인
5. L1 대비 grid trend 확인

## Case

- Case: `cases/paint_booth_panel_frame_komega_l2_fine_yplus/`
- Geometry/model: central filter panel + high-resistance sealed frame
- Supply velocity: `4.36 m/s`
- Filter Forchheimer coefficient: `10000`
- Solver: `simpleFoam`
- Turbulence model: `RASModel kOmegaSST`
- Runtime: Docker image `openfoam-pipeline-local:latest`

## 생성 조건

```bash
python3 scripts/create_paint_booth_plenum_filter_case.py \
  --central-filter-panel-frame \
  --case-dir cases/paint_booth_panel_frame_komega_l2_fine_yplus \
  --supply-velocity 4.36 \
  --filter-forchheimer 10000 \
  --cell-size 0.125 \
  --filter-z-cells 6 \
  --car-refinement-min 3 \
  --car-refinement-max 4 \
  --add-layers \
  --n-surface-layers 5 \
  --absolute-layer-sizes \
  --final-layer-thickness 0.004 \
  --min-layer-thickness 0.0003 \
  --expansion-ratio 1.2 \
  --force
```

L1 대비 주요 변경:

- base cell size: `0.18 m -> 0.125 m`
- filter z cells: `4 -> 6`
- carBody refinement: `2/3 -> 3/4`
- carBody prism layers: `4 -> 5`
- final layer thickness: `0.005 m -> 0.004 m`

## 검증 workflow

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD/cases/paint_booth_panel_frame_komega_l2_fine_yplus" \
  openfoam-pipeline-local:latest \
  bash -lc '
    surfaceCheck constant/triSurface/simplified_car_shell.stl > log.surfaceCheck 2>&1
    blockMesh > log.blockMesh 2>&1
    snappyHexMesh -overwrite > log.snappyHexMesh 2>&1
    checkMesh > log.checkMesh_snappy 2>&1
    topoSet > log.topoSet 2>&1
    timeout 600 simpleFoam > log.simpleFoam 2>&1
    simpleFoam -postProcess -func yPlus -latestTime > log.simpleFoam_post_yPlus 2>&1
    foamToVTK -latestTime > log.foamToVTK 2>&1
  '
```

추가 후처리:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc '
    python3 scripts/postprocess_paint_booth_plenum_metrics.py cases/paint_booth_panel_frame_komega_l2_fine_yplus
    python3 scripts/postprocess_mesh_yplus_summary.py cases/paint_booth_panel_frame_komega_l2_fine_yplus
  '
```

## Mesh 결과

`checkMesh`:

- Points: `265,059`
- Faces: `712,407`
- Internal faces: `679,819`
- Cells: `224,769`
- Max aspect ratio: `16.2478` OK
- Average non-orthogonality: `7.7488`
- Max skewness: `1.1476` OK
- Result: `Mesh OK`

Cell count 비교:

- L0 coarse: 약 `22,953`
- L1 layered: `87,522`
- L2 fine layered: `224,769`

L2는 L0 대비 약 `9.8배`, L1 대비 약 `2.6배`의 cell count를 가진다.

## Boundary layer 결과

`snappyHexMesh` carBody layer summary:

- carBody faces at layering: `20,692`
- Average layers added: `4.98 / 5`
- Overall layer thickness: `0.0143 m`
- Layer coverage: `99.8 %`
- Added layer cells: `103,050`

해석:

- 차체 표면 대부분에 5개 prism layer가 안정적으로 삽입됐다.
- L1의 4-layer coverage `99.6 %`보다 약간 개선됐다.
- 현재 설정은 carBody 중심 layer 전략이며, booth walls에는 아직 prism layer를 넣지 않았다.

## y+ 결과

carBody y+ 통계:

- min: `0.0445`
- p05: `0.135`
- median: `0.529`
- mean: `0.708`
- p95: `1.737`
- p99: `2.327`
- max: `23.496`

분포:

- `y+ < 1`: `74.6 %`
- `1 <= y+ <= 5`: `25.3 %`
- `5 < y+ <= 30`: `0.039 %`
- `y+ > 30`: `0.0 %`

해석:

- carBody는 거의 전 영역이 `y+ < 5`이며, median이 `0.53`이다.
- k-omega SST의 low-y+/wall-resolved 쪽으로는 충분히 fine한 편이다.
- L1보다 y+가 더 낮아졌으므로, wall-resolved 관점에서는 안정적이다.
- 계산비용을 줄이고 싶다면 L2에서 first/near-wall layer를 약간 키워 median y+를 1 근처로 올릴 수 있지만, 현재 목적이 fine mesh라면 이 설정은 보수적이고 적절하다.

OpenFOAM y+ log:

```text
patch plenumTopWall y+ : min = 200.096, max = 1217.74, average = 713.203
patch frontWall y+ : min = 101.2, max = 408.01, average = 215.039
patch rearWall y+ : min = 101.269, max = 407.988, average = 210.701
patch sideWalls y+ : min = 108.194, max = 2088.07, average = 492.835
patch carBody y+ : min = 0.0444715, max = 23.4964, average = 0.707875
```

주의:

- booth wall y+는 높다.
- 현재 물리적 관심이 차체 주변 downdraft/recirculation이면 carBody layer 우선 전략이 맞다.
- booth wall shear나 벽면 박리까지 보려면 wall patch에도 layer를 추가해야 한다.

## Cell-size 통계

VTK cell volume의 등가 cube size `volume^(1/3)` 기준이다.

### Near-car region

- Cells: `160,634`
- Mean: `0.0253 m`
- Median: `0.0155 m`
- p95: `0.0630 m`
- Min/Max: `0.00544 / 0.126 m`

### carBody surface face scale

- Faces: `20,692`
- Surface area: `20.59 m²`
- Mean `sqrt(face area)`: `0.0312 m`
- Median `sqrt(face area)`: `0.0306 m`
- p95 `sqrt(face area)`: `0.0381 m`

L1 대비 carBody surface face median은 `0.0443 m -> 0.0306 m`로 감소했다.

### Filter panel

- Cells: `8,736` from `topoSet filterZone`
- Sealed frame cells: `3,552`
- Filter panel equivalent cell size median: 약 `0.062 m` 수준

## Solver 결과

`simpleFoam`:

- Final time: `200`
- Exit: `End`
- Runtime: 약 `153 s`
- Continuity final log 기준 안정적으로 종료

마지막 로그:

```text
GAMG:  Solving for p, Initial residual = 6.0956e-06, Final residual = 6.08465e-07, No Iterations 2
time step continuity errors : sum local = 6.78126e-06, global = 3.28156e-07, cumulative = -0.152575
smoothSolver:  Solving for omega, Initial residual = 2.56615e-07, Final residual = 5.11186e-09, No Iterations 2
smoothSolver:  Solving for k, Initial residual = 1.29765e-05, Final residual = 3.68202e-07, No Iterations 2
ExecutionTime = 152.9 s  ClockTime = 151 s
End
```

## Flow metric

후처리: `scripts/postprocess_paint_booth_plenum_metrics.py`

- Supply inflow: `7.8480 m³/s`
- Floor exhaust outflow: `7.8431 m³/s`
- Relative mass imbalance: `-0.0626 %`
- Filter pressure drop proxy, rho=1.2 kg/m³: `66.84 Pa`
- Filter-below mean Uz: `-0.847 m/s`
- Filter-below Uz CV: `0.0367`
- Work-zone mean Uz: `-0.200 m/s`
- Work-zone reverse-flow fraction: `9.83 %`
- Near-car mean Uz: `-0.1435 m/s`
- Near-car reverse-flow fraction: `23.27 %`

L1과 비교하면 L2 metric은 대체로 비슷한 방향을 유지한다.

- L1 pressure drop proxy: `61.94 Pa`
- L2 pressure drop proxy: `66.84 Pa`
- L1 work-zone mean Uz: `-0.204 m/s`
- L2 work-zone mean Uz: `-0.200 m/s`
- L1 near-car mean Uz: `-0.148 m/s`
- L2 near-car mean Uz: `-0.1435 m/s`
- L1 near-car reverse fraction: `24.5 %`
- L2 near-car reverse fraction: `23.3 %`

즉, L1에서 이미 coarse L0와 크게 달라졌고, L2는 L1과 비교해 더 안정적인 grid trend를 보인다.

## 현재 결론

이번 L2 fine mesh는 다음 조건을 만족한다.

1. `kOmegaSST` solver 정상 완료
2. `checkMesh`: Mesh OK
3. carBody prism layer coverage `99.8 %`
4. carBody y+ median `0.53`, p95 `1.74`
5. carBody y+ > 30 없음
6. L1 대비 cell count 약 `2.6배`, L0 대비 약 `9.8배`
7. L1과 flow metric trend가 비교적 일관됨

따라서 `cases/paint_booth_panel_frame_komega_l2_fine_yplus/`는 현재 프로젝트의 **fine mesh kOmegaSST baseline**으로 사용할 수 있다.

## 다음 권장 작업

유량/유속 조건이 아직 확정되지 않았으므로 재보정 sweep보다는 다음 품질 개선이 더 유의미하다.

1. **L1/L2 grid comparison figure 생성**
   - velocity mid-plane
   - near-car Uz field
   - y+ map on carBody
2. **booth wall layer 적용 여부 결정**
   - 차체 중심이면 현 상태 유지
   - booth wall recirculation/벽면 shear까지 보려면 wall layer 추가
3. **geometry fidelity 개선**
   - sealed frame을 high-resistance porous zone에서 실제 baffle/wall로 전환
   - floor exhaust를 full-floor outlet에서 grating/slot exhaust로 전환
4. **장기적으로 L3/grid-independence 검토**
   - 현재 L2는 fine baseline이지만 production-grade grid independence를 말하려면 L3 또는 localized refinement study가 필요하다.
