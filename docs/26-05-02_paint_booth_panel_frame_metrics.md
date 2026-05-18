# Paint Booth Plenum/Filter 후처리 Metric 및 Panel-Frame v036 검증

## 목적

이 문서는 `paint_booth_plenum_filter_full_floor_v035` 이후 제안한 다음 작업을 수행한 기록이다.

1. plenum/filter case의 후처리 metric 자동화
2. 기존 full cross-section porous filter를 중앙 filter panel + 외곽 sealed frame 구조로 개선
3. 개선 case의 OpenFOAM mesh/solver/VTK/postprocess smoke test 검증

## 생성/수정 파일

```text
scripts/postprocess_paint_booth_plenum_metrics.py
scripts/create_paint_booth_plenum_filter_case.py
cases/paint_booth_plenum_filter_panel_frame_v036/
docs/figures/26-05-02_paint_booth_panel_frame_geometry.png
docs/figures/26-05-02_paint_booth_panel_frame_velocity_midplane.png
docs/26-05-02_paint_booth_panel_frame_metrics.md
```

## 후처리 자동화 스크립트

신규 스크립트:

```text
scripts/postprocess_paint_booth_plenum_metrics.py
```

역할:

- `foamToVTK -latestTime` 결과 자동 탐색
- `internal.vtu` cell data 기반 region metric 계산
- `supplyInlet`, `floorExhaust` patch VTP 기반 patch flow rate 계산
- filter 상/하부 slab 평균 압력차 proxy 계산
- filter 아래 / work-zone / near-car / plenum 영역 속도 통계 계산
- 결과를 case 내부 JSON으로 저장

출력 예:

```text
cases/<case>/post_plenum_metrics.json
```

주의:

- OpenFOAM incompressible `p`는 kinematic pressure, 단위는 `m2/s2`이다.
- Pa 단위 압력으로 해석하려면 `rho`를 곱해야 한다.
- 내부 virtual plane flux는 VTK slice/cell-data 기반 proxy이다. 정확한 보존형 internal face flux는 향후 faceZone/sampled-surface 기반 OpenFOAM functionObject로 개선하는 것이 좋다.

## 기존 v035 case 재후처리 결과

대상:

```text
cases/paint_booth_plenum_filter_full_floor_v035/
```

Patch flow:

```text
supplyInlet inflow:   7.8480 m3/s
floorExhaust outflow: 7.8405 m3/s
relative imbalance:  -0.096 %
```

Pressure drop proxy:

```text
filter slab delta p: 14.544 m2/s2
rho = 1.2 kg/m3 가정 시: 약 17.45 Pa
```

주요 region metric:

```text
filter_below:
  mean |U| = 0.463 m/s
  mean Uz  = -0.463 m/s
  Uz CV    = 0.106
  reverse flow fraction = 0.0 %

work_zone:
  mean |U| = 0.299 m/s
  mean Uz  = -0.251 m/s
  Uz CV    = 0.427
  reverse flow fraction = 1.03 %

near_car:
  mean |U| = 0.272 m/s
  mean Uz  = -0.220 m/s
  Uz CV    = 0.601
  reverse flow fraction = 5.94 %
```

## Panel-frame v036 case

신규 case:

```text
cases/paint_booth_plenum_filter_panel_frame_v036/
```

생성 명령:

```bash
python3 scripts/create_paint_booth_plenum_filter_case.py \
  --central-filter-panel-frame \
  --force
```

### 구조 변경

기존 v035:

```text
filter layer 전체 cross-section이 porous layer
x = -1.5 ~ 6.5 m
y = -2.0 ~ 2.0 m
z =  2.95 ~ 3.05 m
```

v036:

```text
central porous filter panel:
  x = -1.0 ~ 6.0 m
  y = -1.6 ~ 1.6 m
  z =  2.95 ~ 3.05 m
  area = 22.4 m2

sealed frame zone:
  central panel 외곽 filter-layer 영역
  high-resistance Darcy-Forchheimer porous zone으로 수치적 차폐
```

현재 sealed frame은 실제 내부 baffle/wall face가 아니라, 매우 큰 저항을 가진 porous cellZone으로 근사했다.

```text
filterZone Forchheimer coefficient:      6.8e3
sealedFrameZone Forchheimer coefficient: 5.0e7
```

이는 blockMesh 내부면을 baffle/wall로 재구성하기 전의 안정적인 1차 구현이다.

## v036 검증 결과

실행 순서:

```bash
surfaceCheck
blockMesh
snappyHexMesh -overwrite
checkMesh
topoSet
simpleFoam
foamToVTK -latestTime
postprocess_paint_booth_plenum_metrics.py
```

### checkMesh

```text
points:         32,290
faces:          77,266
internal faces: 67,756
cells:          22,953
Max aspect ratio: 5.00003 OK
Max skewness:     2.63588 OK
Mesh OK
```

### topoSet

```text
filterZone size:      728 cells
sealedFrameZone size: 360 cells
```

### fvOptions 적용 확인

`simpleFoam` log에서 두 porous region이 모두 생성됨을 확인했다.

```text
Porosity region filterPorosity:
  creating porous zone: filterZone

Porosity region sealedFramePorosity:
  creating porous zone: sealedFrameZone
```

### solver

```text
Time = 200
End
```

### VTK export

```text
foamToVTK -latestTime: passed
VTK/paint_booth_plenum_filter_panel_frame_v036_200/internal.vtu 생성 확인
```

## v036 후처리 metric

Patch flow:

```text
supplyInlet inflow:   7.8480 m3/s
floorExhaust outflow: 7.8257 m3/s
relative imbalance:  -0.284 %
```

Pressure drop proxy:

```text
filter slab delta p: 29.076 m2/s2
rho = 1.2 kg/m3 가정 시: 약 34.89 Pa
```

주요 region metric:

```text
filter_below:
  mean |U| = 0.654 m/s
  mean Uz  = -0.652 m/s
  Uz CV    = 0.054
  reverse flow fraction = 0.0 %

work_zone:
  mean |U| = 0.349 m/s
  mean Uz  = -0.291 m/s
  Uz CV    = 0.447
  reverse flow fraction = 2.00 %

near_car:
  mean |U| = 0.309 m/s
  mean Uz  = -0.247 m/s
  Uz CV    = 0.635
  reverse flow fraction = 7.76 %

plenum:
  mean |U| = 1.321 m/s
  mean Uz  = -0.517 m/s
  reverse flow fraction = 32.23 %
```

## v035 vs v036 해석

v036은 중앙 panel로 filter 통과 영역을 제한하고 외곽 frame을 고저항 영역으로 막았기 때문에, 같은 supply 유량에서 filter panel 부근의 downdraft가 강해졌다.

비교:

```text
v035 filter_below mean Uz: -0.463 m/s
v036 filter_below mean Uz: -0.652 m/s

v035 pressure drop proxy: 약 17.45 Pa, rho=1.2 기준
v036 pressure drop proxy: 약 34.89 Pa, rho=1.2 기준
```

이는 구조적으로 자연스러운 변화다. 외곽 frame을 막으면 유동이 중앙 filter panel에 더 집중된다.

다만 목표 평균 downdraft `0.35 m/s`에 비해 현재 filter_below proxy는 과대이다. 다음 단계에서는 아래 두 요소를 함께 calibration해야 한다.

1. supply flow rate / inlet velocity
2. filter 및 sealed-frame resistance coefficient

## 이미지 산출물

Geometry overview:

```text
docs/figures/26-05-02_paint_booth_panel_frame_geometry.png
```

Velocity mid-plane:

```text
docs/figures/26-05-02_paint_booth_panel_frame_velocity_midplane.png
```

## 현재 한계

1. sealed frame은 실제 wall/baffle이 아니라 high-resistance porous cellZone 근사이다.
2. 내부 filter plane flux는 VTK slice 기반 proxy이며, 정확한 보존형 face flux는 아니다.
3. filter pressure drop은 slab 평균 압력차 proxy이다.
4. 목표 `0.35 m/s`에 맞는 pressure-drop calibration은 아직 수행하지 않았다.
5. mesh refinement/grid study 및 y+ 평가는 아직 수행하지 않았다.
6. floor grating/slot 구조는 아직 full-floor pressure outlet이다.

## 다음 추천 작업

1. v036 기준으로 supply velocity sweep 수행
   - 예: `U_supply = 2.3, 3.0, 3.5, 4.36 m/s`
   - 목적: work-zone / near-car downdraft 목표 범위 탐색

2. filter resistance sweep 수행
   - 예: `f = 6.8e3, 1.0e4, 2.0e4, 4.0e4`
   - 목적: pressure drop과 uniformity trade-off 확인

3. internal faceZone 또는 sampled-surface 기반 정확한 filter flux/pressure-drop postprocess 구현

4. sealed frame을 실제 baffle/wall face로 구현하는 v037 검토

5. floor exhaust를 full floor에서 grating/slot variant로 확장
