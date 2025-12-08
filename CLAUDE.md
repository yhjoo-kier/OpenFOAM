# OpenFOAM CFD Analysis Project

## Project Overview
OpenFOAM 기반 CFD 시뮬레이션 프로젝트. 열유체 해석 케이스들을 포함하며, 단순한 예제부터 복잡한 CHT(Conjugate Heat Transfer) 해석까지 다양한 난이도의 케이스를 제공.

## Project Structure
```
/workspaces/OpenFOAM/
├── cases/                      # 모든 시뮬레이션 케이스
│   ├── channelFlow/            # 기본 채널 유동
│   ├── heatsink_flow/          # 히트싱크 냉각 유동
│   └── staggered_finned_tube_cht/  # 핀튜브 CHT 해석
├── scripts/                    # 공용 스크립트 (있는 경우)
├── src/                        # 소스 코드
└── CLAUDE.md                   # 이 파일
```

## Case Structure Convention
각 케이스는 self-contained 구조를 따름:
```
cases/{case_name}/
├── 0/                          # OpenFOAM 초기 조건
├── constant/                   # 물성치, 메시 설정
├── system/                     # 솔버 설정
├── scripts/                    # Python 스크립트 모음
│   ├── generate_mesh.py        # 메시 생성
│   ├── setup_openfoam.py       # OpenFOAM 설정
│   ├── post_process.py         # 후처리
│   └── visualize_results.py    # 시각화
├── images/                     # 대표 결과 이미지 (4-5개)
├── README.md                   # 케이스 설명
└── results_summary.txt         # 결과 요약 (선택)
```

## Common Commands

### Mesh Generation
```bash
# blockMesh 사용
blockMesh

# Gmsh 메시 변환
gmshToFoam mesh.msh
```

### Solvers
```bash
# 단순 유동
simpleFoam

# CHT (Conjugate Heat Transfer)
chtMultiRegionSimpleFoam

# 다중 영역 설정
splitMeshRegions -cellZones -overwrite
```

### Post-processing
```bash
# VTK 변환 (ParaView용)
foamToVTK

# 특정 영역만
foamToVTK -region fluid
foamToVTK -region solid
```

## Git Rules

### Tracked (Git에 포함)
- `0/` - 초기 조건
- `constant/` - 물성치, 메시 설정
- `system/` - 솔버 설정
- `scripts/` - Python 스크립트
- `images/` - 대표 이미지 (4-5개)
- `README.md` - 케이스 설명
- `results_summary.txt` - 결과 요약

### Ignored (Git에서 제외)
- 시간 디렉토리 (`[1-9]*/`, `[1-9][0-9]*/` 등)
- `VTK/` - ParaView 데이터
- `postProcessing/` - 후처리 시계열 데이터
- `processor*/` - 병렬 분해 데이터
- `*.msh` - 생성된 메시 파일
- `log.*` - 솔버 로그
- `__pycache__/` - Python 캐시

## Code Style
- Python: PEP 8 준수
- OpenFOAM 딕셔너리: 표준 OpenFOAM 포맷

## Adding New Cases
새 케이스 생성 시 `cases/CLAUDE.md`의 템플릿을 참고할 것.
