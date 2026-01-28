# OpenFOAM Case Creation Guide

## Standard Case Structure
새로운 CFD 케이스 생성 시 아래 구조를 따를 것:

```
cases/{case_name}/
├── 0/                          # 초기/경계 조건
│   ├── U                       # 속도장
│   ├── p                       # 압력장
│   └── T                       # 온도장 (열해석 시)
│
├── constant/                   # 물성치 및 메시
│   ├── polyMesh/               # OpenFOAM 메시
│   ├── transportProperties     # 유체 물성
│   └── turbulenceProperties    # 난류 모델
│
├── system/                     # 솔버 설정
│   ├── controlDict             # 실행 제어
│   ├── fvSchemes               # 이산화 방법
│   └── fvSolution              # 솔버 설정
│
├── scripts/                    # Python 스크립트 (필수)
│   ├── generate_mesh.py        # 메시 생성 (Gmsh 등)
│   ├── setup_openfoam.py       # OpenFOAM 딕셔너리 설정
│   ├── post_process.py         # 결과 후처리
│   └── visualize_results.py    # 시각화
│
├── images/                     # 대표 이미지 (필수, 4-5개)
│   ├── geometry.png            # 형상/메시
│   ├── velocity_field.png      # 속도장
│   ├── temperature_field.png   # 온도장
│   └── results_summary.png     # 결과 요약
│
├── README.md                   # 케이스 설명 (필수)
└── results_summary.txt         # 수치 결과 요약 (선택)
```

## Naming Conventions

### Case Name
- 소문자, 언더스코어 사용: `staggered_finned_tube_cht`
- 의미 있는 이름: 형상_해석유형

### Image Files
- 내용을 명확히 설명하는 이름
- 번호 prefix 불필요
- 예: `velocity_field.png`, `temperature_contour.png`

## README.md Template
```markdown
# {Case Name}

## Overview
케이스에 대한 간단한 설명

## Geometry
- 주요 치수
- 특징

## Boundary Conditions
| Boundary | Type | Value |
|----------|------|-------|
| inlet    | ...  | ...   |

## Results
주요 결과 요약

## How to Run
실행 방법
```

## CHT (Conjugate Heat Transfer) Cases
CHT 케이스의 경우 다중 영역 구조:
```
0/
├── fluid/
└── solid/
constant/
├── fluid/
└── solid/
system/
├── fluid/
└── solid/
```

## Scripts Guidelines

### generate_mesh.py
- Gmsh API 사용 권장
- 메시 품질 체크 포함
- 출력: `mesh.msh`
- **Grid Study 호환 필수** (아래 섹션 참고)

### setup_openfoam.py
- OpenFOAM 딕셔너리 자동 생성
- 영역 분할, 경계 조건 설정

### post_process.py
- 정량적 결과 추출
- 출력: `results_summary.txt`

### visualize_results.py
- PyVista 또는 matplotlib 사용
- 출력: `images/` 폴더에 PNG 저장

---

## Grid Study 호환 케이스 작성 가이드

새 케이스는 **처음부터 Grid Study 파이프라인을 고려**하여 작성해야 합니다.

### 필수 요구사항

#### 1. generate_mesh.py 시그니처

`mesh_factor` 파라미터를 **반드시** 지원해야 합니다:

```python
#!/usr/bin/env python3
"""Parametric mesh generation for {case_name}."""

import gmsh
from pathlib import Path


def generate_mesh(
    output_path: Path,
    mesh_factor: float = 1.0,           # 필수: 격자 크기 승수
    bl_first_height: float = 0.0005,    # 선택: 경계층 첫 번째 셀 높이
    bl_growth_ratio: float = 1.2,       # 선택: 경계층 성장률
    bl_num_layers: int = 0,             # 선택: 경계층 레이어 수
) -> dict:
    """
    Generate parametric mesh.

    Args:
        output_path: Output .msh file path
        mesh_factor: Mesh size multiplier (2.0=coarse, 1.0=reference, 0.5=fine)
        bl_first_height: First cell height for boundary layer
        bl_growth_ratio: Growth ratio for boundary layer
        bl_num_layers: Number of boundary layer cells (0=disabled)

    Returns:
        dict with mesh info: {"num_cells": int, "num_nodes": int, ...}
    """
    # 기준 격자 크기 (레퍼런스에서 적절한 값)
    BASE_LC = 0.002  # [m]

    # mesh_factor 적용
    lc = BASE_LC * mesh_factor

    gmsh.initialize()
    gmsh.model.add("mesh")

    # ... 형상 정의 ...

    # 격자 크기 설정
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

    # 경계층 설정 (bl_num_layers > 0일 때)
    if bl_num_layers > 0:
        # Distance + Threshold field 사용
        dist_field = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist_field, "SurfacesList", interface_surfaces)

        thresh_field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
        gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", bl_first_height)
        gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", lc)
        gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", 0)
        gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", bl_first_height * bl_num_layers)
        gmsh.model.mesh.field.setAsBackgroundMesh(thresh_field)

    # 메시 생성
    gmsh.model.mesh.generate(3)

    # OpenFOAM 호환 포맷으로 저장
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # 필수!
    gmsh.write(str(output_path))

    # 통계 반환
    num_nodes = len(gmsh.model.mesh.getNodes()[0])
    num_cells = len(gmsh.model.mesh.getElementsByType(4)[0])  # tetrahedra

    gmsh.finalize()

    return {
        "num_cells": num_cells,
        "num_nodes": num_nodes,
        "mesh_factor": mesh_factor,
    }


if __name__ == "__main__":
    # 레퍼런스 메시 생성 (테스트용)
    generate_mesh(Path("mesh.msh"), mesh_factor=1.0)
```

#### 2. mesh_factor 의미 규약

| mesh_factor | 의미 | 예상 셀 수 비율 |
|-------------|------|-----------------|
| 2.0 | 매우 성긴 격자 | ~12.5% |
| 1.5 | 성긴 격자 | ~30% |
| **1.0** | **레퍼런스 (기준)** | **100%** |
| 0.7 | 조밀한 격자 | ~300% |
| 0.5 | 매우 조밀한 격자 | ~800% |

**규약**: `mesh_factor = 1.0`이 레퍼런스 해석에 사용된 격자

#### 3. 모니터링 지표 정의

Grid Study에서 사용할 **주요 모니터링 지표**를 README.md에 명시:

```markdown
## Grid Study Configuration

| 항목 | 값 | 설명 |
|------|-----|------|
| metric_patch | `heat_source` | 모니터링할 경계 패치 |
| metric_field | `T` | 모니터링할 물리량 |
| metric_region | `solid` | 다중 영역 시 대상 영역 |
| convergence_threshold | 1.0% | 권장 수렴 기준 |
```

#### 4. 케이스 검증 체크리스트

새 케이스 작성 완료 전 확인:

- [ ] `generate_mesh.py`가 `mesh_factor` 파라미터 지원
- [ ] `mesh_factor=1.0`으로 레퍼런스 해석 성공
- [ ] `mesh_factor=2.0` (성긴 격자)에서도 수렴
- [ ] `mesh_factor=0.7` (조밀한 격자)에서도 수렴
- [ ] Gmsh 출력 포맷이 2.2 버전
- [ ] README.md에 Grid Study 설정 명시

### Grid Study 실행 예시

케이스가 준비되면:

```python
import sys
sys.path.insert(0, 'src')
from grid_study import run_grid_study

analysis = run_grid_study(
    base_case="cases/new_case",
    adaptive=True,
    threshold=1.0,
    max_cells=1_000_000,
    metric_patch="heat_source",  # README에 명시된 값
    metric_field="T",
    metric_region="solid",
)
```

### 지원되지 않는 케이스

현재 Grid Study 프레임워크가 지원하지 않는 케이스:

| 유형 | 이유 | 대안 |
|------|------|------|
| blockMesh 전용 | 파라메트릭 메시 생성 불가 | Gmsh로 변환 |
| 비정상 해석 | 시간 적분 비용 | 정상상태 근사 |
| 외부 메시 파일 | mesh_factor 적용 불가 | Gmsh 스크립트 작성 |

## Git Rules for Cases
### Include
- `0/`, `constant/`, `system/` (초기 설정)
- `scripts/`, `images/`, `README.md`

### Exclude (자동)
- 시간 디렉토리 (`100/`, `1000/` 등)
- `VTK/`, `postProcessing/`
- `*.msh`, `log.*`
