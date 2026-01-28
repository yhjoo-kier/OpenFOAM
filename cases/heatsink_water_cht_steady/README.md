# Water-Cooled Heatsink CHT Analysis (Steady-State)

수냉식 히트싱크의 Conjugate Heat Transfer(CHT) 정상상태 해석 케이스입니다.

## 형상

```
        [Outlet Tube]
              |
    +---------+----------+
    |  o o o o o o o o   |  <- 8x8 Pin Array
    |  o o o o o o o o   |
    |  o o o o o o o o   |
    |  o o o o o o o o   |     Fluid Channel
    |  o o o o o o o o   |     (Water)
    |  o o o o o o o o   |
    |  o o o o o o o o   |
    |  o o o o o o o o   |
    +---------+----------+
              |
        [Inlet Tube]

    =====================
       Base Plate (Al)
    =====================
         Heat Source
```

| 항목 | 치수 |
|------|------|
| 베이스 플레이트 | 80 x 83 x 2 mm |
| 핀 배열 | 8 x 8 (64개) |
| 핀 직경 | 6 mm |
| 핀 높이 | 15 mm |
| 핀 피치 | 10 mm |
| 냉각수 튜브 직경 | 5 mm |
| 튜브 길이 | 25 mm |

## 물성치

### Fluid (Water)

| 물성 | 값 |
|------|-----|
| 밀도 (rho) | 998 kg/m³ |
| 비열 (Cp) | 4182 J/kg·K |
| 점성계수 (mu) | 0.001 Pa·s |
| 프란틀 수 (Pr) | 7 |
| 열전도도 (k) | ~0.6 W/m·K |

### Solid (Aluminum)

| 물성 | 값 |
|------|-----|
| 밀도 (rho) | 2700 kg/m³ |
| 비열 (Cp) | 900 J/kg·K |
| 열전도도 (kappa) | 200 W/m·K |

## 경계조건

### Fluid Region

| Patch | U | p_rgh | T |
|-------|---|-------|---|
| inlet | pressureInletOutletVelocity | totalPressure (10 Pa) | fixedValue (293.15 K) |
| outlet | pressureInletOutletVelocity | fixedValue (0 Pa) | inletOutlet |
| fluid_to_solid | noSlip | fixedFluxPressure | turbulentTemperatureCoupledBaffleMixed |
| defaultFaces | noSlip | fixedFluxPressure | zeroGradient |

### Solid Region

| Patch | T |
|-------|---|
| heat_source | externalWallHeatFluxTemperature (60 kW/m²) |
| solid_to_fluid | turbulentTemperatureCoupledBaffleMixed |
| defaultFaces | zeroGradient |

## 솔버 설정

- **Application**: chtMultiRegionSimpleFoam
- **Algorithm**: SIMPLE
- **Iterations**: 2000
- **Turbulence**: Laminar

## 실행 방법

### 처음부터 실행

```bash
cd /workspaces/OpenFOAM
./run_water_case_steady.sh
```

### 기존 케이스 재실행

```bash
cd /workspaces/OpenFOAM/cases/heatsink_water_cht_steady

# 이전 결과 삭제
rm -rf [1-9]* VTK log.*

# 시뮬레이션 실행
chtMultiRegionSimpleFoam

# VTK 변환
foamToVTK -region fluid
foamToVTK -region solid
```

## 결과 확인

### ParaView

```bash
# Fluid 온도장/유속장
paraview VTK/fluid/heatsink_water_cht_steady_2000.vtm

# Solid 온도장
paraview VTK/solid/heatsink_water_cht_steady_2000.vtm
```

### 대표 결과 (60 kW/m²)

| 항목 | 값 |
|------|-----|
| Fluid 입구 온도 | 293 K (20°C) |
| Fluid 최대 온도 | 335 K (62°C) |
| Fluid 온도 상승 | ~42 K |
| Solid 온도 범위 | 327 ~ 335 K |

## 파일 구조

```
heatsink_water_cht_steady/
├── 0/                      # 초기 조건
│   ├── fluid/              # U, p, p_rgh, T
│   └── solid/              # T, p
├── constant/               # 물성치
│   ├── fluid/              # thermophysicalProperties, turbulenceProperties
│   ├── solid/              # thermophysicalProperties
│   ├── regionProperties    # 영역 정의
│   └── g                   # 중력 벡터
├── system/                 # 솔버 설정
│   ├── fluid/              # fvSchemes, fvSolution
│   ├── solid/              # fvSchemes, fvSolution
│   ├── controlDict         # 실행 제어
│   └── fvSchemes/fvSolution
└── README.md
```

## Grid Study Configuration

이 케이스는 Grid Independence Study 프레임워크와 호환됩니다.

### 설정 값

| 항목 | 값 | 설명 |
|------|-----|------|
| metric_patch | `heat_source` | 베이스 플레이트 바닥면 |
| metric_field | `T` | 평균 온도 |
| metric_region | `solid` | 고체 영역 |
| convergence_threshold | 1.0% | 권장 수렴 기준 |

### 실행 예시

```python
import sys
sys.path.insert(0, 'src')
from grid_study import run_grid_study

analysis = run_grid_study(
    base_case="cases/heatsink_water_cht_steady",
    adaptive=True,
    threshold=1.0,
    max_cells=1_000_000,
    metric_patch="heat_source",
    metric_field="T",
    metric_region="solid",
)
```

### 검증 결과 (2026-01-26)

| Level | Cells | T_base_avg [K] | Δ [%] | Status |
|-------|-------|----------------|-------|--------|
| L1_coarse | 28,145 | 379.44 | - | - |
| L2_medium | 72,921 | 355.88 | 6.62 | FAIL |
| L3_adaptive | 158,576 | 326.72 | 8.93 | FAIL |
| L4_adaptive | 352,157 | 311.77 | 4.80 | FAIL |
| L5_adaptive | 742,826 | 308.79 | 0.96 | **PASS** |

- Richardson extrapolated: **307.95 K**
- Recommended level: **L4_adaptive** (352k cells)

## 주의사항

1. **메시 생성**: 이 케이스는 `scripts/generate_mesh.py`로 파라메트릭 메시를 생성합니다.
2. **열 결합**: `turbulentTemperatureCoupledBaffleMixed` BC로 유체-고체 인터페이스 열전달을 처리합니다.
3. **Heat Flux 조절**: `0/solid/T`의 `heat_source` 패치에서 `q uniform` 값을 변경하여 열유속을 조절할 수 있습니다.
4. **Grid Study**: `mesh_factor` 파라미터로 격자 조밀도 조절 (1.0=레퍼런스, 0.7=조밀, 2.0=성김)
