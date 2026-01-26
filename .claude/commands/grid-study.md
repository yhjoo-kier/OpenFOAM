# Grid Independence Study

OpenFOAM CHT 케이스의 격자 독립성 검증을 자동으로 수행합니다.

## 중요: 실행 방법

**반드시 아래 Python 코드를 실행하여 grid study를 수행하세요. 직접 구현하지 마세요.**

$ARGUMENTS 값을 파싱하여 아래 코드의 파라미터로 전달하세요:

```python
import sys
sys.path.insert(0, 'src')

from grid_study import run_grid_study

# $ARGUMENTS에서 파싱한 값 사용
# 예: /grid-study cases/my_case --adaptive --threshold 2.0 --max-cells 500000
analysis = run_grid_study(
    base_case="$CASE_PATH",      # 기본값: "cases/heatsink_water_cht_steady"
    adaptive=$ADAPTIVE,          # --adaptive 플래그 있으면 True, 없으면 False
    threshold=$THRESHOLD,        # --threshold 값 (기본값: 1.0)
    max_cells=$MAX_CELLS,        # --max-cells 값 (기본값: 2_000_000)
    max_levels=$MAX_LEVELS,      # --max-levels 값 (기본값: 10)
    metric_patch="$PATCH",       # --patch 값 (기본값: "heat_source")
    metric_field="$FIELD",       # --field 값 (기본값: "T")
    metric_region="$REGION",     # --region 값 (기본값: "solid")
)

# 결과 출력
print(f"\n수렴 여부: {analysis.is_converged}")
print(f"추천 레벨: {analysis.recommended_level}")
if analysis.stop_reason:
    print(f"종료 이유: {analysis.stop_reason}")
if analysis.extrapolated_value:
    print(f"Richardson 외삽값: {analysis.extrapolated_value:.4f}")
```

## 인자 파싱 규칙

| 인자 | 변수 | 기본값 | 예시 |
|------|------|--------|------|
| 첫 번째 위치 인자 | `$CASE_PATH` | `cases/heatsink_water_cht_steady` | `cases/my_case` |
| `--adaptive` | `$ADAPTIVE` | `False` | `True` |
| `--threshold N` | `$THRESHOLD` | `1.0` | `2.0` |
| `--max-cells N` | `$MAX_CELLS` | `2_000_000` | `500000` |
| `--max-levels N` | `$MAX_LEVELS` | `10` | `8` |
| `--patch NAME` | `$PATCH` | `heat_source` | `inlet` |
| `--field NAME` | `$FIELD` | `T` | `p` |
| `--region NAME` | `$REGION` | `solid` | `fluid` |

## 출력 예시

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

## 종료 조건 (Adaptive Mode)

| 조건 | 의미 |
|------|------|
| `converged` | ✓ 수렴 달성 (Δ < threshold) |
| `max_cells_exceeded` | ⚠ 셀 수 초과 |
| `max_levels_reached` | ⚠ 레벨 수 도달 |
| `max_runtime_exceeded` | ⚠ 시간 초과 |

$ARGUMENTS
