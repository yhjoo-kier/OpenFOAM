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

## 실행 단계

1. 설정된 메시 레벨별로 Gmsh 메시 생성
2. gmshToFoam으로 OpenFOAM 포맷 변환
3. splitMeshRegions로 fluid/solid 영역 분리
4. chtMultiRegionSimpleFoam 시뮬레이션 실행
5. 지정된 패치에서 평균값 추출
6. 수렴성 분석 및 GCI 계산
7. 결과 보고서 생성 (TXT, JSON, CSV, PNG)

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

config = GridStudyConfig(
    base_case_path=Path("cases/heatsink_water_cht_steady"),
    study_name="my_study",
    metric_patch="heat_source",
    convergence_threshold=1.0,
)

study = GridStudy(config)
analysis = study.run()
study.generate_reports(analysis)
```

$ARGUMENTS
