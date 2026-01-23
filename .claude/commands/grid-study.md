# Grid Independence Study

OpenFOAM CHT 케이스의 격자 독립성 검증을 자동으로 수행합니다.

## 사용법

```
/grid-study [case_path] [options]
```

### 인자
- `case_path`: OpenFOAM 케이스 경로 (기본값: cases/heatsink_water_cht_steady)

### 옵션
- `--levels N`: 메시 레벨 수 (기본값: 4)
- `--threshold N`: 수렴 기준 % (기본값: 1.0)
- `--patch NAME`: 모니터링 패치 (기본값: heat_source)
- `--field NAME`: 모니터링 필드 (기본값: T)
- `--region NAME`: 영역 (기본값: solid)
- `--adaptive`: 적응형 모드 활성화 (수렴까지 자동 격자 조밀화)
- `--max-cells N`: 최대 셀 수 제한 (기본값: 2,000,000)
- `--max-levels N`: 최대 레벨 수 제한 (기본값: 10)

## 실행 단계

1. 설정된 메시 레벨별로 Gmsh 메시 생성
2. gmshToFoam으로 OpenFOAM 포맷 변환
3. splitMeshRegions로 fluid/solid 영역 분리
4. chtMultiRegionSimpleFoam 시뮬레이션 실행
5. 지정된 패치에서 평균값 추출
6. 수렴성 분석 및 GCI 계산
7. 결과 보고서 생성 (TXT, JSON, CSV, PNG)

## 적응형 모드 (Adaptive Mode)

`--adaptive` 옵션 사용 시, 수렴 조건을 만족할 때까지 자동으로 격자를 조밀화합니다.

### 동작 방식
1. 초기 메시 레벨들로 시뮬레이션 실행
2. 연속된 두 레벨 간 변화율(Δ) 계산
3. Δ < threshold이면 수렴 → 종료
4. 수렴하지 않으면 다음 레벨 자동 생성 (mesh_factor × 0.7)
5. 제한 조건 도달 시 종료

### 종료 조건
- ✓ 수렴 달성 (Δ < threshold)
- ⚠ 최대 셀 수 초과 (max_cells)
- ⚠ 최대 레벨 수 도달 (max_levels)
- ⚠ 레벨당 최대 실행 시간 초과 (max_runtime_per_level)

## 출력 예시

```
┌─────────────────┬────────────┬────────────┬──────────┬────────┐
│ Level           │      Cells │ T_avg [K]  │    Δ [%] │ Status │
├─────────────────┼────────────┼────────────┼──────────┼────────┤
│ L1_coarse       │     50,000 │   330.1234 │        - │   -    │
│ L2_medium       │    150,000 │   332.5678 │   0.74   │  PASS  │
│ L3_fine         │    400,000 │   332.8901 │   0.10   │  PASS  │
└─────────────────┴────────────┴────────────┴──────────┴────────┘

Result: ✓ CONVERGED (threshold: 1.0%)
Recommended: L2_medium
```

## 프레임워크 위치

`src/grid_study/`

## Python API

```python
from grid_study import GridStudy, GridStudyConfig
from pathlib import Path

# Standard mode (predefined levels)
config = GridStudyConfig(
    base_case_path=Path("cases/heatsink_water_cht_steady"),
    study_name="my_study",
    metric_patch="heat_source",
    convergence_threshold=1.0,
)

# Adaptive mode (auto-refine until converged)
config = GridStudyConfig(
    base_case_path=Path("cases/heatsink_water_cht_steady"),
    study_name="adaptive_study",
    adaptive_mode=True,
    convergence_threshold=1.0,
    max_cells=2_000_000,
    max_levels=10,
    refinement_ratio=0.7,
)

study = GridStudy(config)
analysis = study.run()
study.generate_reports(analysis)

# Check stop reason in adaptive mode
if analysis.stop_reason:
    print(f"Stopped: {analysis.stop_reason}")
```

$ARGUMENTS
