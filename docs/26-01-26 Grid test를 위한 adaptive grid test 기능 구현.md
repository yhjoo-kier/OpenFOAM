# Adaptive Grid Independence Study 기능 구현

**작성일**: 2026-01-26
**버전**: v0.2.0
**위치**: `src/grid_study/`

---

## 1. 개요

### 1.1 Grid Independence Study란?

CFD 시뮬레이션에서 **격자 독립성 검증(Grid Independence Study)**은 계산 결과가 메시 해상도에 의존하지 않음을 증명하는 필수 절차입니다. 격자를 점진적으로 조밀하게 하면서 관심 지표(온도, 압력 등)의 변화가 허용 오차 이내로 수렴하는지 확인합니다.

### 1.2 Adaptive Mode의 필요성

기존 방식은 미리 정의된 메시 레벨들을 순차적으로 실행했습니다. 하지만 이 방식은:
- 몇 개의 레벨이 필요한지 사전에 알 수 없음
- 수렴하지 않으면 수동으로 더 조밀한 메시를 추가해야 함
- AI 에이전트가 자동화하기 어려움

**Adaptive Mode**는 수렴할 때까지 자동으로 격자를 조밀화하며, 설정된 제한 조건에 도달하면 자동 중단합니다.

---

## 2. 핵심 개념

### 2.1 수렴 판정 기준

연속된 두 메시 레벨 간 변화율을 계산합니다:

```
Δ = |값_coarse - 값_fine| / |값_fine| × 100%
```

- **Δ < threshold**: 수렴 달성 ✓
- **Δ ≥ threshold**: 더 조밀한 메시 필요

### 2.2 Richardson Extrapolation

3개 이상의 레벨이 있을 때, 격자 크기가 0으로 수렴할 때의 이론적 값을 추정합니다:

```
f_extrapolated = f_fine + (f_fine - f_medium) / (r^p - 1)
```

- `r`: 격자 비율 (h_coarse / h_fine)
- `p`: 수렴 차수 (자동 계산)

### 2.3 자동 격자 조밀화

Adaptive mode에서 새 레벨 생성 시:

```python
new_mesh_factor = previous_mesh_factor × refinement_ratio
new_bl_height = previous_bl_height × bl_refinement_ratio
```

기본값 `refinement_ratio = 0.7`은 각 레벨마다 격자 크기를 30% 감소시킵니다.

---

## 3. 설정 옵션

### 3.1 GridStudyConfig 주요 필드

| 필드 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `adaptive_mode` | bool | False | Adaptive 모드 활성화 |
| `convergence_threshold` | float | 1.0 | 수렴 기준 (%) |
| `refinement_ratio` | float | 0.7 | mesh_factor 승수 |
| `bl_refinement_ratio` | float | 0.7 | boundary layer 승수 |
| `max_cells` | int | 2,000,000 | 최대 셀 수 제한 |
| `max_levels` | int | 10 | 최대 레벨 수 제한 |
| `max_runtime_per_level` | int | None | 레벨당 최대 실행 시간 (초) |

### 3.2 종료 조건

Adaptive mode는 다음 조건 중 하나를 만족하면 종료합니다:

1. **✓ converged**: 수렴 달성 (`Δ < threshold`)
2. **⚠ max_cells_exceeded**: 셀 수가 `max_cells` 초과
3. **⚠ max_levels_reached**: 레벨 수가 `max_levels` 도달
4. **⚠ max_runtime_exceeded**: 레벨 실행 시간이 `max_runtime_per_level` 초과
5. **⚠ error**: 시뮬레이션 오류 발생

---

## 4. 사용법

### 4.1 Python API - 기본 사용

```python
import sys
sys.path.insert(0, 'src')

from pathlib import Path
from grid_study import GridStudyConfig, GridStudy, MeshLevel

# 초기 메시 레벨 정의 (선택사항)
initial_levels = [
    MeshLevel(
        name="L1_coarse",
        mesh_factor=2.0,
        bl_first_height=0.001,
        bl_growth_ratio=1.2,
        bl_num_layers=0,
    ),
    MeshLevel(
        name="L2_medium",
        mesh_factor=1.3,
        bl_first_height=0.0005,
        bl_growth_ratio=1.2,
        bl_num_layers=0,
    ),
]

# Adaptive 모드 설정
config = GridStudyConfig(
    base_case_path=Path("cases/heatsink_water_cht_steady"),
    study_name="adaptive_study",
    output_dir=Path("grid_study_output"),
    mesh_levels=initial_levels,

    # 모니터링 지표
    metric_name="T_base_avg",
    metric_patch="heat_source",
    metric_field="T",
    metric_region="solid",

    # Adaptive 모드 옵션
    adaptive_mode=True,
    convergence_threshold=1.0,  # 1% 수렴 기준
    refinement_ratio=0.75,
    max_cells=500_000,
    max_levels=8,
)

# 실행
study = GridStudy(config)
analysis = study.run()
reports = study.generate_reports(analysis)

# 결과 확인
print(f"수렴 여부: {analysis.is_converged}")
print(f"종료 이유: {analysis.stop_reason}")
print(f"추천 레벨: {analysis.recommended_level}")
```

### 4.2 편의 함수 사용

```python
from grid_study import run_grid_study

analysis = run_grid_study(
    base_case="cases/heatsink_water_cht_steady",
    adaptive=True,
    threshold=1.0,
    max_cells=1_000_000,
    max_levels=10,
)
```

### 4.3 Skill 사용 (Claude Code)

```
/grid-study cases/heatsink_water_cht_steady --adaptive --threshold 1.0 --max-cells 500000
```

---

## 5. 샘플 실행 결과

### 5.1 실행 조건

- **케이스**: `cases/heatsink_water_cht_steady` (수냉 히트싱크 CHT)
- **모니터링 지표**: 베이스 플레이트 평균 온도 (`heat_source` 패치)
- **수렴 기준**: 1.0%
- **최대 셀 수**: 500,000
- **refinement_ratio**: 0.75

### 5.2 결과 테이블

```
┌─────────────────┬────────────┬────────────────┬──────────┬────────┐
│ Level           │      Cells │ T_base_avg [K] │    Δ [%] │ Status │
├─────────────────┼────────────┼────────────────┼──────────┼────────┤
│ L1_coarse       │     28,145 │       379.4407 │        - │   -    │
│ L2_medium       │     72,921 │       355.8805 │     6.62 │  FAIL  │
│ L3_adaptive     │    158,576 │       326.7191 │     8.93 │  FAIL  │
│ L4_adaptive     │    352,157 │       311.7680 │     4.80 │  FAIL  │
│ L5_adaptive     │    742,826 │       308.7909 │     0.96 │  PASS  │
└─────────────────┴────────────┴────────────────┴──────────┴────────┘

Result: ✓ CONVERGED (threshold: 1.0%)
Recommended: L4_adaptive
Richardson extrapolated: 307.9465 K
Stop reason: max_cells_exceeded (742,826 > 500,000)
```

### 5.3 분석

1. **자동 레벨 생성**: L3, L4, L5가 자동으로 생성됨
2. **수렴 달성**: L5에서 0.96% 변화율로 1% 기준 충족
3. **제한 조건 동작**: 742k 셀이 500k 제한을 초과하여 추가 조밀화 중단
4. **Richardson 외삽값**: 307.95 K (L5의 308.79 K와 0.3% 차이)
5. **추천 레벨**: L4_adaptive (수렴 달성한 쌍 중 계산 비용이 낮은 레벨)

### 5.4 온도 수렴 경향

```
379.44 K → 355.88 K → 326.72 K → 311.77 K → 308.79 K → (307.95 K 추정)
```

격자가 조밀해질수록 열전달이 더 정확하게 계산되어 온도가 낮아지는 경향을 보입니다.

---

## 6. 출력 파일

Adaptive 모드 실행 후 생성되는 파일들:

```
grid_study_output/
├── config.json                    # 설정 파일
├── cases/                         # 각 레벨별 OpenFOAM 케이스
│   ├── study_L1_coarse/
│   ├── study_L2_medium/
│   ├── study_L3_adaptive/
│   └── ...
├── meshes/                        # Gmsh 메시 파일
│   ├── L1_coarse.msh
│   └── ...
├── logs/                          # 솔버 로그
│   ├── log.chtMultiRegionSimpleFoam_L1_coarse
│   └── ...
└── results/
    ├── study_report.txt           # 텍스트 보고서
    ├── study_report.json          # JSON 보고서
    ├── study_results.csv          # CSV 데이터
    └── study_convergence.png      # 수렴 그래프
```

---

## 7. 아키텍처

### 7.1 모듈 구조

```
src/grid_study/
├── __init__.py          # 패키지 API
├── __main__.py          # CLI 진입점
├── config.py            # 설정 클래스 (MeshLevel, GridStudyConfig)
├── mesh_generator.py    # Gmsh 메시 생성
├── case_manager.py      # OpenFOAM 케이스 관리
├── runner.py            # 시뮬레이션 실행
├── extractor.py         # 결과 추출
├── analyzer.py          # 수렴 분석
├── reporter.py          # 보고서 생성
├── study.py             # 메인 오케스트레이터
└── cli.py               # 명령줄 인터페이스
```

### 7.2 Adaptive 모드 플로우

```
┌─────────────────────────────────────────────────────────────┐
│                    GridStudy.run()                          │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
              ┌─────────────────────────┐
              │  adaptive_mode == True? │
              └─────────────────────────┘
                     │ Yes
                     ▼
         ┌─────────────────────────┐
         │   _run_adaptive()       │
         └─────────────────────────┘
                     │
                     ▼
    ┌────────────────────────────────────┐
    │  while True:                       │
    │    1. Check max_levels limit       │──────► STOP
    │    2. Get/generate mesh level      │
    │    3. Run simulation               │
    │    4. Extract metric               │
    │    5. Check max_cells limit        │──────► STOP
    │    6. Check convergence            │──────► STOP (converged)
    │    7. Check max_runtime limit      │──────► STOP
    │    8. Generate next level          │
    │       (mesh_factor × 0.7)          │
    └────────────────────────────────────┘
                     │
                     ▼
         ┌─────────────────────────┐
         │   Analyze & Report      │
         └─────────────────────────┘
```

---

## 8. 제한사항 및 개선 방향

### 8.1 현재 제한사항

- 단일 모니터링 지표만 지원
- Gmsh 기반 메시 생성만 지원 (blockMesh 미지원)
- CHT 솔버 (chtMultiRegionSimpleFoam) 기본 설정

### 8.2 향후 개선 방향

- [ ] 다중 지표 모니터링
- [ ] 병렬 시뮬레이션 지원
- [ ] 메모리 사용량 기반 제한
- [ ] 비정상 해석 지원
- [ ] 자동 초기 레벨 추정

---

## 9. 참고 문헌

1. Roache, P.J. (1994). "Perspective: A Method for Uniform Reporting of Grid Refinement Studies"
2. Celik, I.B. et al. (2008). "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications"
