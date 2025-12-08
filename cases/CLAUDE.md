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

### setup_openfoam.py
- OpenFOAM 딕셔너리 자동 생성
- 영역 분할, 경계 조건 설정

### post_process.py
- 정량적 결과 추출
- 출력: `results_summary.txt`

### visualize_results.py
- PyVista 또는 matplotlib 사용
- 출력: `images/` 폴더에 PNG 저장

## Git Rules for Cases
### Include
- `0/`, `constant/`, `system/` (초기 설정)
- `scripts/`, `images/`, `README.md`

### Exclude (자동)
- 시간 디렉토리 (`100/`, `1000/` 등)
- `VTK/`, `postProcessing/`
- `*.msh`, `log.*`
