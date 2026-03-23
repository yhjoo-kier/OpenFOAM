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

## Image-to-CFD Paper (진행 중)

2D 이미지로부터 VLM(Gemini)을 이용해 3D 실내 형상을 추상화하고 CFD까지 자동화하는 프레임워크의 논문화 작업.

- **계획 문서**: `docs/26-03-17_image_to_cfd_paper_plan.md`
- **프레임워크 엔트리포인트**: `scripts/run_indoor_stabilized.py`
- **핵심 접근**: rule-based 벤치마크 데이터셋 (형상 난이도 × 뷰 타입) + reference CFD 정답 확보 → 프레임워크 정량 평가
- **데이터셋 위치**: `benchmark/` — rule-based 형상 생성, reference CFD, 2D 렌더링, 평가 결과 일체

### Paper Data Convention (Critical — 반드시 준수)

논문 결과는 **preset-matched** 데이터만 사용한다. Solver robustness preset마다 inlet velocity가 다르므로, reference와 predicted 케이스가 동일한 preset(=동일 inlet velocity)으로 수렴한 pair만 비교해야 한다.

**⚠️ 사용 금지:**
- `cfd_metrics.json` — preset 미매칭 상태의 CFD 비교 결과 (inlet velocity 불일치)
- `evaluation_aggregate_summary_posthoc_scaled_longest_span.json`의 원본 CFD score (구 버전)

**✅ 사용해야 하는 파일:**
- `cfd_metrics_matched.json` — 각 eval dir 내, preset-matched reference 기준 CFD 비교
- `benchmark/manifests/evaluation_aggregate_summary_phase2.json` — 매칭된 집계 (최신)
- `benchmark/manifests/evaluation_statistics_phase2.json` — 매칭된 통계 (최신)
- `benchmark/manifests/evaluation_statistics_matched.json` — 매칭 전용 통계

**Reference 케이스 경로:**
- Preset-matched: `cases/phase2_ref_{scene}_preset_{preset}/` (논문용)
- 원본 Phase 2: `cases/phase2_ref_{scene}/` (참고용만, 논문에 사용 금지)

**Predicted 케이스 경로:**
- `cases/phase2_pred_bench_{case}_{view}/`

**최종 수치 (2026-03-23):**
- Structural: 0.781 ± 0.151, CFD agreement: 0.477 ± 0.158 (n=97)
- Best view: floorplan (struct 0.884, CFD 0.572)

**배경:** Solver preset별 inlet velocity가 다름 (robust=0.02, ultra_robust=0.005, laminar=0.01 m/s). Phase 2 초기에는 ref/pred가 다른 preset으로 수렴하여 velocity magnitude similarity ≈ 0이었음. `--force-preset` 플래그 도입 + 26개 reference 재해석으로 해결 (2026-03-23).

## Paper Writing Convention (Markdown → LaTeX → PDF)

### Figure
- 폰트: Arial (불가 시 Liberation Sans → DejaVu Sans fallback, QC log에 기록)
- 파일 형식: PDF (벡터, 기본) + PNG (600 dpi 이상)
- 네이밍: `fig_섹션키워드_설명.pdf` (숫자 번호 사용 금지, LaTeX에서 자동 번호 할당)
- 키워드 예시: bench (benchmark), method (methodology), result (results), discuss (discussion), demo (application demo)
- 본문 참조: `[Fig:label]` → LaTeX 변환 시 `Fig.~\ref{fig:label}`

### Table
- 캡션+레이블을 표 상단에 배치: `[Table:키워드-설명]`
- 숫자 번호 사용 금지
- LaTeX 변환 시 `Table~\ref{tab:키워드-설명}`

### Equation
- `\tag{Eq:키워드-변수}` 형식, 숫자 번호 금지
- LaTeX 변환 시 `Eq.~\ref{eq:키워드-변수}`

### Algorithm
- `[Alg:키워드-설명]` 형식 레이블+캡션 상단 배치
- algorithm2e 또는 algorithmicx 패키지 사용

### Citation
- 초안 작성 시: `[저자 연도 키워드]` 형식 (의미론적 임시 인용)
- 확정 시: Zotero BibTeX → `[@citekey]`
- LaTeX 변환: `[@citekey]` → `\cite{citekey}`
- 참고문헌 등장순 정렬 (elsarticle-num 스타일)

### LaTeX 빌드
- 빌드 순서: pdflatex → bibtex → pdflatex → pdflatex
- natbib sort&compress 사용
- siunitx 패키지로 SI 단위
- LaTeX 변환 규칙: `[Fig:label]` → `Fig.~\ref{fig:label}`, `[Table:label]` → `Table~\ref{tab:label}`, `[Eq:label]` → `Eq.~\ref{eq:label}`, `[@citekey]` → `\cite{citekey}`

## Adding New Cases
새 케이스 생성 시 `cases/CLAUDE.md`의 템플릿을 참고할 것.
