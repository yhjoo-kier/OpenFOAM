# OpenFOAM CFD Analysis Project

## Project Overview
OpenFOAM 기반 CFD 시뮬레이션 프로젝트. 열유체 해석 케이스들을 포함하며, 단순한 예제부터 복잡한 CHT(Conjugate Heat Transfer) 해석까지 다양한 난이도의 케이스를 제공.

## Project Structure
```
/workspaces/OpenFOAM/
├── cases/                      # 모든 시뮬레이션 케이스
│   ├── channelFlow/            # 기본 채널 유동
│   ├── heatsink_flow/          # 히트싱크 냉각 유동
│   ├── heatsink_water_cht_steady/  # 수냉 히트싱크 CHT (정상상태)
│   └── staggered_finned_tube_cht/  # 핀튜브 CHT 해석
├── scripts/                    # 공용 스크립트
├── src/                        # 소스 코드
│   ├── geometry/               # 형상 생성 스크립트
│   ├── visualization/          # 시각화 스크립트
│   └── grid_study/             # 격자 독립성 검증 프레임워크
├── .claude/commands/           # Claude Code skills
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

## Grid Study Framework

격자 독립성 검증을 위한 자동화 프레임워크 (`src/grid_study/`).

### Adaptive Mode (권장)

수렴할 때까지 자동으로 격자를 조밀화:

```python
import sys
sys.path.insert(0, 'src')

from grid_study import run_grid_study

# 한 줄로 실행 (Adaptive 모드)
analysis = run_grid_study(
    base_case="cases/heatsink_water_cht_steady",
    adaptive=True,           # 수렴까지 자동 조밀화
    threshold=1.0,           # 1% 수렴 기준
    max_cells=1_000_000,     # 최대 셀 수 제한
)

# 결과 확인
print(f"수렴: {analysis.is_converged}")
print(f"추천 레벨: {analysis.recommended_level}")
print(f"종료 이유: {analysis.stop_reason}")
```

### 종료 조건

| 조건 | 설명 |
|------|------|
| `converged` | Δ < threshold (성공) |
| `max_cells_exceeded` | 셀 수 초과 |
| `max_levels_reached` | 레벨 수 도달 |
| `max_runtime_exceeded` | 시간 초과 |

### Standard Mode

미리 정의된 레벨로 실행:

```python
from grid_study import GridStudyConfig, GridStudy, MeshLevel

config = GridStudyConfig(
    base_case_path=Path("cases/heatsink_water_cht_steady"),
    mesh_levels=[
        MeshLevel("L1", mesh_factor=2.0, bl_first_height=0.001, bl_growth_ratio=1.2, bl_num_layers=0),
        MeshLevel("L2", mesh_factor=1.0, bl_first_height=0.0005, bl_growth_ratio=1.2, bl_num_layers=0),
    ],
)
study = GridStudy(config)
analysis = study.run()
```

### CLI
```bash
python -m grid_study run cases/heatsink_water_cht_steady --adaptive -t 1.0
```

### Skill
- `/grid-study` - 격자 독립성 검증 실행

## Papers

논문별 작업물은 `papers/` 아래에 self-contained 구조로 관리. 각 논문 폴더의 `CLAUDE.md` 참조.

| 논문 | 폴더 | 상태 |
|------|------|------|
| Image-to-CFD (VLM 기반 실내 CFD 자동화) | `papers/image_to_cfd/` | B&E desk reject → 재투고 예정 |
| CFD Visual QA | `papers/cfd_visual_qa/` | 작업 중 |

### 공통 자산 (루트)
- `benchmark/` — 벤치마크 데이터셋
- `cases/` — OpenFOAM 시뮬레이션 케이스
- `scripts/run_indoor_stabilized.py` — Image-to-CFD 프레임워크 엔트리포인트
- `src/` — 소스 코드 (geometry, visualization, grid_study)

### 공통 인용 관리
- BibTeX: `C:\Vaults\Research\Zotero_YHJoo.bib` (Zotero BetterBibTeX 자동 내보내기)
- 검색: PaperSearch (`~/Projects/PaperSearch`)
- 초안: `[저자 연도 키워드]` → 확정: `[@citekey]`

## Adding New Cases
새 케이스 생성 시 `cases/CLAUDE.md`의 템플릿을 참고할 것.
