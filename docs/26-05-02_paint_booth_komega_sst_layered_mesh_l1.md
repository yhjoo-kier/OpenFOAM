# Paint Booth CFD: k-omega SST 기반 L1 경계층 메쉬 도입 결과

## 목적

추천 조건(`U_supply = 4.36 m/s`, `filter Forchheimer = 10000`)을 기준으로 기존 coarse mesh(L0)보다 높은 수준의 메쉬를 만들고, 차체 표면에 prism boundary layer를 도입해 `kOmegaSST` 난류모델에 맞는 near-wall 해상도를 확인했다.

현재 목표는 production CFD 최종 메쉬가 아니라, 다음 단계 grid refinement/y+ study를 위한 **검증된 L1 경계층 baseline**을 확보하는 것이다.

## 생성/검증 case

- Case: `cases/paint_booth_panel_frame_komega_l1_yplus/`
- Geometry/model: central filter panel + high-resistance sealed frame
- Supply velocity: `4.36 m/s`
- Filter Forchheimer coefficient: `10000`
- Solver: `simpleFoam`
- Turbulence model: `RASModel kOmegaSST`
- Runtime: Docker image `openfoam-pipeline-local:latest`

## 주요 generator 변경

`scripts/create_paint_booth_plenum_filter_case.py`에 다음 옵션을 추가/검증했다.

- Base mesh refinement:
  - `--cell-size`
  - `--filter-z-cells`
  - `--car-refinement-min`
  - `--car-refinement-max`
- Boundary layer mesh:
  - `--add-layers`
  - `--n-surface-layers`
  - `--expansion-ratio`
  - `--final-layer-thickness`
  - `--min-layer-thickness`
  - `--absolute-layer-sizes`
- 난류모델:
  - `constant/turbulenceProperties`를 `kOmegaSST`로 작성
  - `0/omega` 및 `omegaWallFunction` 포함

이번 L1 case 생성 명령:

```bash
python3 scripts/create_paint_booth_plenum_filter_case.py \
  --central-filter-panel-frame \
  --case-dir cases/paint_booth_panel_frame_komega_l1_yplus \
  --supply-velocity 4.36 \
  --filter-forchheimer 10000 \
  --cell-size 0.18 \
  --filter-z-cells 4 \
  --car-refinement-min 2 \
  --car-refinement-max 3 \
  --add-layers \
  --n-surface-layers 4 \
  --absolute-layer-sizes \
  --final-layer-thickness 0.005 \
  --min-layer-thickness 0.0005 \
  --expansion-ratio 1.2 \
  --force
```

검증 명령 요약:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD/cases/paint_booth_panel_frame_komega_l1_yplus" \
  openfoam-pipeline-local:latest \
  bash -lc '
    surfaceCheck constant/triSurface/simplified_car_shell.stl > log.surfaceCheck 2>&1
    blockMesh > log.blockMesh 2>&1
    snappyHexMesh -overwrite > log.snappyHexMesh 2>&1
    checkMesh > log.checkMesh_snappy 2>&1
    topoSet > log.topoSet 2>&1
    timeout 300 simpleFoam > log.simpleFoam 2>&1
    simpleFoam -postProcess -func yPlus -latestTime > log.simpleFoam_post_yPlus 2>&1
    foamToVTK -latestTime > log.foamToVTK 2>&1
  '
```

## Mesh 결과

`checkMesh` 결과:

- Points: `107,144`
- Faces: `280,975`
- Internal faces: `265,158`
- Cells: `87,522`
- Max aspect ratio: `16.259` OK
- Average non-orthogonality: `8.8557`
- Max skewness: `1.52721` OK
- Result: `Mesh OK`

기존 L0 추천 case와 비교하면 cell 수는 약 `22,953 -> 87,522`로 약 `3.8배` 증가했다.

## Boundary layer 결과

`snappyHexMesh` layer summary, carBody 기준:

- carBody faces at layering: `9,935`
- Average layers added: `3.97 / 4`
- Overall layer thickness: `0.0155 m`
- Layer coverage: `99.6 %`

즉, 이번 설정에서는 차체 표면 대부분에 4개 prism layer가 안정적으로 들어갔다.

## y+ 결과

`simpleFoam -postProcess -func yPlus -latestTime` 및 VTK 기반 `carBody` y+ 통계:

- carBody y+ min: `0.139`
- carBody y+ p05: `0.378`
- carBody y+ median: `1.300`
- carBody y+ mean: `1.611`
- carBody y+ p95: `3.526`
- carBody y+ p99: `3.922`
- carBody y+ max: `38.846`
- Fraction y+ < 1: `32.7 %`
- Fraction 1 <= y+ <= 5: `67.2 %`
- Fraction 5 < y+ <= 30: `0.08 %`
- Fraction y+ > 30: `0.02 %`

해석:

- 차체 표면은 대부분 `y+ < 5`이며 median이 약 `1.3`이라 `kOmegaSST`의 low-y+/integrated near-wall 해석 기준으로 상당히 적절하다.
- 소수의 outlier만 `y+ > 30`이다. 현재 simplified car shell의 sharp/curvature/mesh transition 부근일 가능성이 높으며, 다음 L2에서 local refinement 또는 layer 품질 개선 대상으로 보면 된다.
- 이번 layer는 **carBody에만 적용**했다. booth wall들은 아직 prism layer가 없어서 wall-function 영역의 높은 y+가 나타난다.

OpenFOAM wall y+ log:

```text
patch plenumTopWall y+ : min = 385.91, max = 1589.58, average = 1027.56
patch frontWall y+ : min = 179.22, max = 619.288, average = 309.229
patch rearWall y+ : min = 177.598, max = 618.97, average = 305.074
patch sideWalls y+ : min = 213.971, max = 2429.6, average = 743.522
patch carBody y+ : min = 0.13921, max = 38.8462, average = 1.61059
```

## Cell-size 통계

VTK cell volume의 등가 cube size `volume^(1/3)` 기준:

### Near-car region

- Cells: `66,426`
- Mean: `0.0366 m`
- Median: `0.0223 m`
- p95: `0.0907 m`
- Min/Max: `0.00851 / 0.1828 m`

### Car close box

- Cells: `52,710`
- Mean: `0.0325 m`
- Median: `0.0217 m`
- p95: `0.0899 m`

### Filter panel

- Cells: `2,888`
- Mean: `0.0918 m`
- Median: `0.0918 m`
- p95: `0.0928 m`

### Car body surface face scale

- Faces: `9,935`
- Surface area: `20.59 m²`
- Mean `sqrt(face area)`: `0.0448 m`
- Median `sqrt(face area)`: `0.0443 m`
- p95 `sqrt(face area)`: `0.0555 m`

## Flow metric 변화

후처리: `scripts/postprocess_paint_booth_plenum_metrics.py`

- Supply inflow: `7.8480 m³/s`
- Floor exhaust outflow: `7.8328 m³/s`
- Relative mass imbalance: `-0.194 %`
- Filter pressure drop proxy, rho=1.2 kg/m³: `61.94 Pa`
- Filter-below mean Uz: `-0.867 m/s`
- Filter-below Uz CV: `0.0324`
- Work-zone mean Uz: `-0.204 m/s`
- Work-zone reverse-flow fraction: `11.3 %`
- Near-car mean Uz: `-0.148 m/s`
- Near-car reverse-flow fraction: `24.5 %`

L0 추천 case와 비교하면:

- Pressure drop proxy는 약 `51.0 -> 61.9 Pa`로 증가
- Work-zone mean Uz는 `-0.285 -> -0.204 m/s`로 변화
- Near-car mean Uz는 `-0.240 -> -0.148 m/s`로 변화
- Near-car reverse fraction은 `7.8 % -> 24.5 %`로 증가

이 변화는 mesh/layer 도입으로 차체 주변 near-wall 및 wake/recirculation이 더 해상되면서 coarse mesh의 과도하게 매끈한 downdraft trend가 달라진 것으로 해석된다. 따라서 기존 calibration sweep의 최적점은 L1/L2 mesh에서 재보정해야 한다.

## 추가 후처리 스크립트

`scripts/postprocess_mesh_yplus_summary.py`를 추가했다.

역할:

- `checkMesh` 주요 항목 파싱
- `snappyHexMesh` carBody layer coverage 파싱
- VTK 내부 field에서 region별 cell equivalent size 계산
- VTK boundary `carBody`에서 face size 및 y+ 통계 계산
- output: `mesh_yplus_summary.json`

사용 예:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 scripts/postprocess_mesh_yplus_summary.py cases/paint_booth_panel_frame_komega_l1_yplus'
```

## 결론

이번 L1 case는 다음 기준을 만족한다.

1. `kOmegaSST` 기반 solver 정상 완료
2. `checkMesh`: Mesh OK
3. 차체 표면 prism layer coverage 약 `99.6 %`
4. carBody y+ median 약 `1.3`, p95 약 `3.5`
5. 차체 주변 cell scale이 L0보다 뚜렷하게 개선됨

따라서 이 case는 **kOmegaSST + 차체 경계층 도입의 첫 번째 유효 L1 baseline**으로 사용할 수 있다.

## 다음 권장 작업

1. **L1 조건에서 재보정 sweep**
   - L0 calibration이 L1에서는 유동 metric이 변했으므로 `U_supply`와 `filter Forchheimer`를 다시 조정해야 한다.
   - 우선 후보:
     - `U_supply = 4.36, 5.0, 5.5 m/s`
     - `filter Forchheimer = 6800, 10000, 15000`
2. **L2 fine mesh**
   - `cell-size = 0.12~0.15 m`
   - filter z cells `5~6`
   - car refinement `3/4`
   - carBody layer는 현재 y+ 기준을 유지하도록 first/overall layer thickness 조정
3. **booth wall layer 정책 결정**
   - 차체 표면 도장 품질 중심이면 현재처럼 carBody 우선
   - booth wall shear/벽면 recirculation도 중요하면 벽면 layer도 추가하고 wall y+를 별도로 관리
4. **near-car reverse-flow 저감 설계 탐색**
   - floor exhaust 분포, side exhaust/grating, sealed frame 실제 baffle화, plenum diffuser 구조를 조정해야 한다.
