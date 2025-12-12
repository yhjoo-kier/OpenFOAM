# topopt_rect_step_flow

## Overview

FEniCSx 기반 위상최적화(`fenics_toptim/topopt/results_rect/solid.step`)로 생성된 **3D 고체 형상**을
직사각 채널에 배치된 **장애물(obstacle)**로 간주하고, OpenFOAM에서 유동해석을 수행하는 케이스입니다.

- **Geometry unit**: STEP 원본은 **mm** → 메시 생성 시 **1e-3 스케일(mm→m)** 적용
- **Solver**: `simpleFoam` (steady, incompressible), `laminar`
- **2D-equivalent**: `front/back = symmetryPlane`

## Geometry

- **Channel**: \(x \in [0, L], y \in [0, H], z \in [0, T]\)
  - \(L = 2.2\,\mathrm{mm}\)
  - \(H = 0.41\,\mathrm{mm}\)
  - \(T = 10.0\,\mathrm{mm}\)
- **Obstacle**: `geometry/solid.step`

## Boundary Conditions

| Boundary | Type | Value |
|----------|------|-------|
| inlet | `fixedValue` | 균일 유입: \(U_{mean} = 0.3\,\mathrm{m/s}\) |
| outlet | p=`fixedValue` | p = 0 |
| topWall/bottomWall/obstacle | wall | `noSlip` |
| front/back | symmetryPlane | 2D 등가 |

> **Note**: 원래 inlet은 포물선 프로파일(`codedFixedValue`: \(u_x = 6U_{mean}(y/H)(1-y/H)\))로 설계되었으나,
> root 사용자 환경에서 보안 제한으로 `codedFixedValue`가 동작하지 않아 균일 유입(`fixedValue`)으로 대체되었습니다.

## How to Run

### 1) Mesh generation (Gmsh)

```bash
cd cases/topopt_rect_step_flow

# STEP(mm) -> (scale 1e-3) -> mesh.msh
python3 scripts/generate_mesh.py --output mesh.msh

# Convert to OpenFOAM polyMesh
gmshToFoam mesh.msh
```

### 2) Fix patch types (after gmshToFoam)

`gmshToFoam`는 기본적으로 모든 경계를 `patch`로 가져오는 경우가 많습니다.
아래 스크립트로 `wall/symmetryPlane/patch` 타입을 자동 보정합니다.

```bash
python3 scripts/setup_openfoam.py
```

### 3) Run solver

```bash
checkMesh
simpleFoam | tee log.simpleFoam
```

### 4) Post-processing

```bash
python3 scripts/post_process.py
python3 scripts/visualize_results.py
```

## Outputs

- `postProcessing/`: inlet/outlet 평균압, 유량(phi) 기록 (system/controlDict의 functionObjects)
- `results_summary.txt`: Δp, 유량 요약 (`scripts/post_process.py`)
- `images/geometry.png`, `images/velocity_field.png`, `images/pressure_field.png`

## Notes

- 유체 물성(점성)은 `constant/transportProperties`의 `nu`를 수정하세요.
- 메시는 약 100만 개의 tetrahedra 셀로 구성됩니다.

## Troubleshooting

### 1) `codedFixedValue` 보안 오류 (root 사용자)

**증상**: simpleFoam 실행 시 다음 오류 발생
```
FOAM FATAL IO ERROR:
This code should not be executed by someone with administrator rights for security reasons.
```

**원인**: OpenFOAM의 `codedFixedValue`는 동적 라이브러리를 생성하므로, root 사용자에서 보안상 차단됩니다.

**해결**: `0/U`의 inlet 경계를 `fixedValue`로 변경:
```cpp
inlet
{
    type            fixedValue;
    value           uniform (0.3 0 0);
}
```

### 2) `div((nuEff*dev2(T(grad(U)))))` 스킴 누락

**증상**: simpleFoam 실행 시 `Entry 'div(...)' not found` 오류

**해결**: `system/fvSchemes`의 `divSchemes`에 추가:
```cpp
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;  // 추가
}
```

### 3) 경계면(inlet/outlet 등) 미검출

**증상**: gmshToFoam 후 `obstacle`만 인식되고 inlet/outlet 등이 `defaultFaces`로 분류됨

**원인**: `generate_mesh.py`의 평면 감지 tolerance가 OCC geometry 오차(~1e-7)보다 작았음

**해결**: `generate_mesh.py`에서 tolerance를 증가 (이미 수정됨):
```python
# 기존: tol = max(1e-12, 1e-6 * min(Lm, Hm, Tm))
tol = max(1e-6, 1e-4 * min(Lm, Hm, Tm))  # 수정됨
```

### 4) checkMesh에서 2개 region 감지

**증상**: `Number of regions: 2` 경고 발생

**원인**: Boolean cut 후 obstacle 내부에 작은 solid region이 남음 (정상 동작)

**해결**: 무시해도 됨. "Mesh OK" 메시지가 출력되면 시뮬레이션 가능.

## Example Results

- **Inlet 평균 압력**: ~508 Pa
- **Outlet 평균 압력**: 0 Pa
- **압력 강하 (Δp)**: ~508 Pa
- **유량**: 1.23×10⁻⁶ m³/s
- **수렴**: 약 400 iterations


