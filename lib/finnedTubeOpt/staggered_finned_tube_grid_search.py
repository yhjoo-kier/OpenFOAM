"""
Staggered circular finned-tube bank geometry grid search & visualization.

이 모듈은 `staggered_finned_tube_bank.StaggeredFinnedTubeBank` 모델을 이용하여
고정 펌핑 동력 조건에서 선택한 기하 변수들에 대해 grid search를 수행하고,
열전달 성능(예: UA_effective, j/f^(1/3))을 contour 그래프로 시각화한다.

기본 설정
---------
- 외형 체적 (W_bank, H_bank, L_bank)과 펌핑 동력 P_pump_target 는 고정.
- 튜브 외경 d_o, 핀 두께 t_fin, 핀 높이 h_fin 은 기본적으로 고정.
- 기본 grid search 변수: S_T (가로 피치), S_L (세로 피치), S_fin (clear fin spacing).

탐색 범위 설정 근거
-------------------
이 프로젝트의 기본 공력·열전달 모델은
`staggered_finned_tube_bank.StaggeredFinnedTubeBank` 에 구현된
**Briggs & Young (1963) 열전달 상관식 + ESDU 86022 (1986) 압력강하 상관식**
을 따르므로, grid search 범위 역시 이 두 상관식의 유효 영역 **교집합**을
반영하여 잡는다.

상관식 적용 범위 (docs/ESDU_86022_상관식_정리.md 참조)
------------------------------------------------------

- Briggs & Young (1963) high-fin circular tube correlation:
  * 사용 변수: 레이놀즈 수 Re_d 와 두 개의 기하 비율 (s/h_f), (s/t).
  * 유효 Re_d: 1,100 ~ 18,000
  * 유효 s/d_o: ~ 0.40 (문헌에 따라 0.25~0.40)

- ESDU 86022 (1986) high-fin staggered tube bank friction correlation:
  * 사용 변수: Re_d, A_increase, (S_T/d_o), (S_L/d_o).
  * 유효 Re_d (V_max 기준): 10² ~ 10⁵ (층류 말단 ~ 완전 난류)
  * 유효 S_T/d_o: 1.1 ~ 4.0
  * 유효 S_L/d_o: 1.1 ~ 3.0
  * 유효 fin density: 4 ~ 11 fpi → s/d_o ≈ 0.13 ~ 0.38 (d_o=15.88mm, t=0.25mm 기준)
  * 유효 d_o: 9.5 ~ 50.8 mm
  * 유효 h_fin: 8.5 ~ 15.9 mm
    ※ 본 코드의 h_fin = 7.0 mm는 ESDU 범위를 약간 하회하나,
       Gemini 문헌 조사에서 6~16 mm 범위를 언급하여 완화 적용함.
  * 유효 d_f/d_o: 1.2 ~ 2.4

  참고:
    - docs/ESDU 86022 압력강하 상관식 적용 범위 - Gemini.md
    - docs/ESDU_86022_상관식_정리.md

두 상관식의 교집합 범위
-----------------------
  - Re_d: 1,100 ~ 18,000 (B&Y 범위가 더 제한적, ESDU는 10²~10⁵)
  - S_fin/d_o: 0.15 ~ 0.38 (ESDU fin density 4~11 fpi 범위)
  - S_T/d_o: 1.97 ~ 2.55 (ESDU 유효 범위 1.1~4.0 내)
  - S_L/d_o: 1.70 ~ 2.58 (ESDU 유효 범위 1.1~3.0 내)

따라서 기본 grid search 범위는 다음과 같이 설정한다 (무차원 비 기준):

  - 1.97 ≤ S_T/d_o ≤ 2.55
  - 1.70 ≤ S_L/d_o ≤ 2.58
  - 0.15 ≤ S_fin/d_o ≤ 0.38  (← ESDU 86022 적용으로 확장됨)

이 값들은 사용자가 필요에 따라 조정할 수 있으며, 보다 엄밀한 적용 범위가
필요한 경우 ESDU 86022 원문이나 추가 실험 데이터를 참고하여 상·하한을
재조정하는 것이 바람직하다.

확장성
------
- 본 모듈의 핵심 grid search 루틴은 임의의 BankGeometry 필드
  (`d_o`, `S_T`, `S_L`, `S_fin`, `t_fin`, `h_fin` 등)을 설계 변수로 받을 수 있도록
  일반화되어 있다.
- 기본 예제(`__main__`)에서는 요청사항에 따라 S_T, S_L, S_fin 세 변수에 대해서만
  grid search와 contour plot을 수행한다.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 스크립트 파일이 위치한 디렉토리 (출력 파일 저장용)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# matplotlib 한글 폰트 설정 (Windows 환경)
try:
    plt.rcParams['font.family'] = 'Malgun Gothic'  # 맑은 고딕
    plt.rcParams['axes.unicode_minus'] = False  # 마이너스 기호 깨짐 방지
except:
    # 폰트 설정 실패 시 기본값 유지
    pass

# =============================================================================
# 상관식 적용 범위 상수 정의
# =============================================================================
# 참조: docs/ESDU 86022 압력강하 상관식 적용 범위 - Gemini.md
#
# [레이놀즈 수 범위]
# - Briggs & Young (1963): 1,100 ~ 18,000
# - ESDU 86022 (1986):     10² ~ 10⁵ (Gemini 문헌 조사)
# → 교집합: 1,100 ~ 18,000 (B&Y가 더 제한적)
RE_D_MIN_VALID = 1.1e3  # Briggs & Young 하한
RE_D_MAX_VALID = 1.8e4  # Briggs & Young 상한

# [ESDU 86022 튜브 피치 범위]
# - S_T/d_o: 1.1 ~ 4.0
# - S_L/d_o: 1.1 ~ 3.0
# 현재 grid search 범위(1.97~2.55, 1.70~2.58)는 이 범위 내에 있음
ESDU_ST_OVER_DO_MIN, ESDU_ST_OVER_DO_MAX = 1.1, 4.0
ESDU_SL_OVER_DO_MIN, ESDU_SL_OVER_DO_MAX = 1.1, 3.0

from staggered_finned_tube_bank import (
    BankGeometry,
    GasProperties,
    StaggeredFinnedTubeBank,
    BankPerformance,
)


@dataclass
class GridSearchResult:
    """
    grid search 결과 컨테이너.

    Attributes
    ----------
    param_names:
        설계 변수 이름 목록 (예: ["S_T", "S_L", "S_fin"]).
    param_values:
        각 설계 변수에 대한 1D numpy 배열 리스트.
        `values` 배열의 각 축과 순서가 일치한다.
    objective:
        사용한 목적 함수 이름 ("UA_effective", "PEC", "h_effective" 등).
    values:
        목적 함수 값 배열. shape = (len(param_values[0]), len(param_values[1]), ...).
    best_index:
        values 배열 상에서 최댓값 위치의 인덱스 튜플.
    best_value:
        최댓값 (목적 함수 값).
    best_geometry:
        최적 형상에 해당하는 BankGeometry 인스턴스.
    best_performance:
        최적 형상에 대한 BankPerformance 인스턴스.
    """

    param_names: List[str]
    param_values: List[np.ndarray]
    objective: str
    values: np.ndarray
    best_index: Tuple[int, ...]
    best_value: float
    best_geometry: BankGeometry
    best_performance: BankPerformance


def build_default_param_grids(
    geom: BankGeometry,
    n_ST: int = 25,
    n_SL: int = 25,
    n_Sfin: int = 15,
) -> Dict[str, np.ndarray]:
    """
    기본 grid search 용 S_T, S_L, S_fin 범위를 생성한다.

    모든 범위는 상관식/문헌의 무차원 적용 범위(P_t/D, P_l/D, s/D)를
    geom.d_o 에 맞게 dimensionalize 한 값이다.

    Parameters
    ----------
    geom:
        기준 형상 (튜브 외경 d_o 를 사용하여 무차원 범위를 dimensionalize 한다).
    n_ST, n_SL, n_Sfin:
        각 변수에 대해 생성할 grid 포인트 수.

    Returns
    -------
    Dict[str, np.ndarray]
        {"S_T": ..., "S_L": ..., "S_fin": ...} 형태의 1D 배열 사전.
    """
    d_o = geom.d_o

    # 문헌 기반 무차원 범위 (상세 근거는 모듈 상단 주석 참조).
    # 튜브 피치: ESDU 86022 유효 범위 (1.1~4.0, 1.1~3.0) 내에서 보수적 범위 설정
    # 참조: docs/ESDU 86022 압력강하 상관식 적용 범위 - Gemini.md
    ST_over_do_min, ST_over_do_max = 1.97, 2.55  # ESDU 범위 1.1~4.0 내
    SL_over_do_min, SL_over_do_max = 1.70, 2.58  # ESDU 범위 1.1~3.0 내
    # S_fin: ESDU 86022 적용으로 확장
    # - 기존 (Robinson & Briggs): 0.066 ~ 0.095 → s ≈ 1.05–1.5 mm
    # - 변경 (ESDU 86022):        0.15 ~ 0.38   → s ≈ 2.4–6.0 mm (d_o = 15.88 mm 기준)
    # ※ 4~11 fpi 범위에 해당하며, ESDU 86022 검증 범위 내
    Sfin_over_do_min, Sfin_over_do_max = 0.15, 0.38

    S_T_values = np.linspace(ST_over_do_min * d_o, ST_over_do_max * d_o, n_ST)
    S_L_values = np.linspace(SL_over_do_min * d_o, SL_over_do_max * d_o, n_SL)
    S_fin_values = np.linspace(Sfin_over_do_min * d_o, Sfin_over_do_max * d_o, n_Sfin)

    return {
        "S_T": S_T_values,
        "S_L": S_L_values,
        "S_fin": S_fin_values,
    }


def run_grid_search(
    base_geometry: BankGeometry,
    gas: GasProperties | None,
    P_pump_target: float,
    param_grids: Mapping[str, np.ndarray],
    objective: str = "UA_effective",
) -> GridSearchResult:
    """
    주어진 BankGeometry 를 기준으로 여러 기하 변수를 grid search 한다.

    Briggs & Young / ESDU 86022 상관식의 적용 범위를 존중하기 위해
    레이놀즈 수가 [RE_D_MIN_VALID, RE_D_MAX_VALID] 범위를 벗어나는 조합에 대해서는
    목적 함수 값을 NaN 으로 두고, 최적 설계 후보 선정에서는 자동으로 제외한다.

    Parameters
    ----------
    base_geometry:
        검색 기준이 되는 BankGeometry 인스턴스.
        param_grids 에 포함되지 않은 변수 값은 이 인스턴스의 값을 사용한다.
    gas:
        GasProperties 인스턴스. None 이면 기본 공기(0°C) 특성을 사용한다.
    P_pump_target:
        [W] 고정 펌핑 동력.
    param_grids:
        설계 변수 이름 → 후보 값 1D 배열 매핑.
        예: {"S_T": np.linspace(...), "S_L": ..., "S_fin": ...}
    objective:
        목적 함수 선택:
        - "UA_effective" (기본값): BankPerformance.UA_effective
        - "PEC": j/f^(1/3) = BankPerformance.j_over_f_cubert
        - "h_effective": BankPerformance.h_effective

    Returns
    -------
    GridSearchResult
    """
    if gas is None:
        gas = GasProperties()

    param_names = list(param_grids.keys())
    if not param_names:
        raise ValueError("param_grids 에 적어도 하나의 설계 변수가 포함되어야 합니다.")

    param_values: List[np.ndarray] = [
        np.asarray(param_grids[name], dtype=float) for name in param_names
    ]
    shape = tuple(len(v) for v in param_values)
    values = np.empty(shape, dtype=float)

    best_value = -np.inf
    best_index: Tuple[int, ...] = tuple(0 for _ in shape)
    best_perf: BankPerformance | None = None
    best_geom: BankGeometry | None = None

    # N차원 grid 상의 모든 조합을 탐색
    for idx in np.ndindex(*shape):
        geom_kwargs = base_geometry.__dict__.copy()
        for name, grid, i in zip(param_names, param_values, idx):
            geom_kwargs[name] = float(grid[i])

        geom = BankGeometry(**geom_kwargs)
        bank = StaggeredFinnedTubeBank(
            geometry=geom,
            gas=gas,
            P_pump_target=P_pump_target,
            fin_efficiency=None,
        )
        perf = bank.evaluate(T_s=40.0, T_g_in=10.0)

        # 레이놀즈 수 유효 범위(Briggs & Young / ESDU 86022 교집합)를 벗어나면
        # 해당 조합은 최적 설계 후보에서 제외한다.
        Re_d = perf.Re_d
        if RE_D_MIN_VALID <= Re_d <= RE_D_MAX_VALID:
            if objective == "UA_effective":
                obj_val = perf.UA_effective
            elif objective == "PEC":
                obj_val = perf.j_over_f_cubert
            elif objective == "h_effective":
                obj_val = perf.h_effective
            else:
                raise ValueError(
                    f"알 수 없는 objective '{objective}'. "
                    "['UA_effective', 'PEC', 'h_effective'] 중에서 선택하세요."
                )
        else:
            obj_val = np.nan

        values[idx] = obj_val

        if np.isfinite(obj_val) and obj_val > best_value:
            best_value = obj_val
            best_index = idx
            best_perf = perf
            best_geom = geom

    if best_perf is None or best_geom is None:
        raise RuntimeError("grid search 중 최적 값을 찾지 못했습니다.")

    return GridSearchResult(
        param_names=param_names,
        param_values=param_values,
        objective=objective,
        values=values,
        best_index=best_index,
        best_value=float(best_value),
        best_geometry=best_geom,
        best_performance=best_perf,
    )


def plot_contour_from_result(
    result: GridSearchResult,
    fixed_dim_index: int | None = None,
) -> None:
    """
    GridSearchResult 에서 2D contour 플롯을 생성한다.

    - 설계 변수 차원이 2개인 경우: 전체 domain 에 대해 contourf 플롯 생성.
    - 설계 변수 차원이 3개인 경우: 한 변수(fixed_dim_index)를 고정한 단면에서
      나머지 두 변수에 대한 contourf 플롯 생성.

    Parameters
    ----------
    result:
        GridSearchResult 인스턴스.
    fixed_dim_index:
        len(param_names) == 3 인 경우, 세 번째 축에서 사용할 인덱스.
        None 이면 중앙 인덱스를 사용.
    """
    n_dim = len(result.param_names)

    if n_dim == 1:
        # 1차원인 경우에는 단순 선 그래프
        x = result.param_values[0]
        y = result.values
        plt.figure(figsize=(6, 4))
        plt.plot(x, y, "-o")
        plt.xlabel(result.param_names[0] + " [m]")
        plt.ylabel(result.objective)
        plt.grid(True, alpha=0.3)
        return

    if n_dim == 2:
        x_vals = result.param_values[0]
        y_vals = result.param_values[1]
        X, Y = np.meshgrid(x_vals, y_vals, indexing="ij")
        Z = result.values

        plt.figure(figsize=(6, 5))
        cf = plt.contourf(X, Y, Z, levels=20, cmap="viridis")
        plt.colorbar(cf, label=result.objective)
        plt.xlabel(result.param_names[0] + " [m]")
        plt.ylabel(result.param_names[1] + " [m]")
        plt.title(f"{result.objective} contour")
        return

    if n_dim == 3:
        # 하나의 변수(보통 S_fin)를 고정한 단면에서 contour 를 그림.
        # 관례상 마지막 변수를 고정하는 것이 직관적이므로,
        # fixed_dim_index 가 None 이면 마지막 축의 중앙 인덱스를 사용한다.
        fixed_axis = 2
        axis_len = result.values.shape[fixed_axis]

        if fixed_dim_index is None:
            fixed_dim_index = axis_len // 2
        else:
            if not (0 <= fixed_dim_index < axis_len):
                raise ValueError(
                    f"fixed_dim_index={fixed_dim_index} 가 유효 범위 [0, {axis_len-1}] 밖입니다."
                )

        x_vals = result.param_values[0]
        y_vals = result.param_values[1]
        X, Y = np.meshgrid(x_vals, y_vals, indexing="ij")
        Z = result.values[:, :, fixed_dim_index]

        fixed_name = result.param_names[fixed_axis]
        fixed_val = result.param_values[fixed_axis][fixed_dim_index]

        plt.figure(figsize=(6, 5))
        cf = plt.contourf(X, Y, Z, levels=20, cmap="viridis")
        plt.colorbar(cf, label=result.objective)
        plt.xlabel(result.param_names[0] + " [m]")
        plt.ylabel(result.param_names[1] + " [m]")
        plt.title(
            f"{result.objective} contour "
            f"(고정 {fixed_name} = {fixed_val:.4f} m)"
        )
        return

    raise ValueError("plot_contour_from_result 는 최대 3개 설계 변수까지만 지원합니다.")


def _print_best_result(result: GridSearchResult, file=None) -> None:
    """grid search 결과 중 최적 설계와 성능 지표를 요약 출력."""
    import sys
    out = file if file else sys.stdout
    
    print("=== Grid search 최적 설계 (objective: {}) ===".format(result.objective), file=out)
    for axis, name in enumerate(result.param_names):
        vals = result.param_values[axis]
        idx = result.best_index[axis]
        val = vals[idx]
        print(f"  {name:<6s} = {val:.6f} [m]", file=out)

    geom = result.best_geometry
    perf = result.best_performance

    print("  --- Dimensionless / performance ---", file=out)
    print(f"  Re_d          = {perf.Re_d:.1f}", file=out)
    print(f"  j             = {perf.j:.5f}", file=out)
    print(f"  K_f (ESDU)    = {perf.f:.5f}", file=out)
    print(f"  j/f^(1/3)     = {perf.j_over_f_cubert:.5f}", file=out)
    print(f"  UA_effective  = {perf.UA_effective:.2f} [W/K]", file=out)
    print(f"  h_effective   = {perf.h_effective:.2f} [W/m^2-K]", file=out)
    print(f"  A_total       = {perf.A_total:.4f} [m^2]", file=out)
    print(f"  P_pump        = {perf.P_pump:.1f} [W]", file=out)


def run_default_example() -> None:
    """
    문서에 사용된 예제 형상과 상관식 적용 범위 기반의 기본 grid search 데모.

    - 고정:
      * W_bank, H_bank, L_bank: 보고서/예제와 동일 값
      * d_o, t_fin, h_fin: 예제 값 고정
      * GasProperties: 기본 공기 (0°C)
      * P_pump_target: 1000 W (완화된 Re 범위 확인용)

    - 설계 변수:
      * S_T, S_L, S_fin

    - 목적 함수:
      * UA_effective (고정 펌핑 동력 하 열전달량에 비례하는 유효 전열계수)

    결과로 S_T–S_L 평면에서 UA_effective contour 를 그리고,
    콘솔에는 최적 설계와 성능 지표를 요약 출력한다.
    """
    gas = GasProperties()
    geom_base = BankGeometry(
        W_bank=0.42,
        H_bank=0.42,
        L_bank=0.165,
        d_o=0.01588,   # 15.88 mm
        S_T=0.04,      # S_T/d_o ≈ 2.52 (권장 범위 1.97–2.55 내)
        S_L=0.035,     # S_L/d_o ≈ 2.20 (권장 범위 1.70–2.58 내)
        S_fin=0.004,   # s ≈ 4.0 mm (권장 범위 2.4–6.0 mm 내, ESDU 86022 기준)
        t_fin=0.00025, # 0.25 mm → δ_f/d_o ≈ 0.016
        h_fin=0.007,   # 7.0 mm → L_f/d_o ≈ 0.44 (ESDU 범위 약간 하회, 완화 적용)
    )

    param_grids = build_default_param_grids(geom_base)

    result = run_grid_search(
        base_geometry=geom_base,
        gas=gas,
        P_pump_target=1000.0,
        param_grids=param_grids,
        objective="UA_effective",
    )

    _print_best_result(result)
    plot_contour_from_result(result)
    plt.tight_layout()
    # GUI 없이 파일로만 저장 (스크립트 위치 기준)
    filename = os.path.join(_SCRIPT_DIR, f'grid_search_{objective}.png')
    plt.savefig(filename, dpi=150)
    print(f"\n그래프 저장: {filename}")
    plt.close()  # GUI 창 표시 방지


if __name__ == "__main__":
    # 기본 예시 실행:
    # 주어진 체적(W_bank, H_bank, L_bank)과 펌핑 동력 조건에서
    # UA_effective 와 PEC 각각에 대해 grid search 를 수행하여 최적 형상을 비교한다.
    
    gas = GasProperties()
    geom_base = BankGeometry(
        W_bank=0.42,
        H_bank=0.42,
        L_bank=0.165,
        d_o=0.01588,   # 15.88 mm
        S_T=0.04,      # S_T/d_o ≈ 2.52
        S_L=0.035,     # S_L/d_o ≈ 2.20
        S_fin=0.004,   # s ≈ 4.0 mm
        t_fin=0.00025, # 0.25 mm
        h_fin=0.007,   # 7.0 mm
    )
    param_grids = build_default_param_grids(geom_base)

    # 1. UA_effective 최적화
    result_file = os.path.join(_SCRIPT_DIR, "comparison_result.txt")
    with open(result_file, "w", encoding="utf-8") as f:
        print("\n" + "="*60, file=f)
        print(">>> CASE 1: Maximize UA_effective (Capacity)", file=f)
        print("="*60, file=f)
        result_ua = run_grid_search(
            base_geometry=geom_base,
            gas=gas,
            P_pump_target=500.0,
            param_grids=param_grids,
            objective="UA_effective",
        )
        _print_best_result(result_ua, file=f)
        # 최적값의 S_fin 인덱스에서 contour plot 생성
        plot_contour_from_result(result_ua, fixed_dim_index=result_ua.best_index[2])
        plt.tight_layout()
        plt.savefig(os.path.join(_SCRIPT_DIR, 'grid_search_UA_effective.png'), dpi=150)
        plt.close()

        # 2. PEC (j/f^(1/3)) 최적화
        print("\n" + "="*60, file=f)
        print(">>> CASE 2: Maximize PEC (Efficiency, j/f^(1/3))", file=f)
        print("="*60, file=f)
        result_pec = run_grid_search(
            base_geometry=geom_base,
            gas=gas,
            P_pump_target=500.0,
            param_grids=param_grids,
            objective="PEC",
        )
        _print_best_result(result_pec, file=f)
        # 최적값의 S_fin 인덱스에서 contour plot 생성
        plot_contour_from_result(result_pec, fixed_dim_index=result_pec.best_index[2])
        plt.tight_layout()
        plt.savefig(os.path.join(_SCRIPT_DIR, 'grid_search_PEC.png'), dpi=150)
        plt.close()
        
    print(f"결과가 {result_file} 에 저장되었습니다.")
