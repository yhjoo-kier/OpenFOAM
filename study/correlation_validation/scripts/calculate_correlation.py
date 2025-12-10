"""
Calculate correlation predictions for validation cases.

이 스크립트는 Briggs & Young (1963) 열전달 상관식과 ESDU 86022 (1986)
압력강하 상관식을 사용하여 검증 케이스들의 Nu와 K_f를 예측합니다.

상관식:
- Briggs & Young: Nu = 0.134 * Re_d^0.681 * Pr^(1/3) * (s/h_f)^0.2 * (s/t)^0.1134
- ESDU 86022: K_f = 4.567 * Re_d^(-0.242) * A_increase^0.504 * (S_T/d_o)^(-0.376) * (S_L/d_o)^(-0.546)
"""
import json
import math
import os
import sys

# lib 경로 추가
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib/finnedTubeOpt'))

# 형상 파라미터 (검증용)
D_O = 0.016       # 튜브 외경 [m]
S_T = 0.036       # 가로 피치 [m]
S_L = 0.034       # 세로 피치 [m]
H_FIN = 0.010     # 핀 방사 높이 [m]
T_FIN = 0.0005    # 핀 두께 [m]
S_FIN = 0.004     # clear fin spacing [m]
P_FIN = 0.0045    # 핀 피치 [m]
N_FINS = 5        # 핀 개수
N_REPEAT = 5      # REV 반복 횟수

# 유체 물성 (공기 300K)
RHO = 1.177       # 밀도 [kg/m³]
MU = 1.846e-5     # 점성 [Pa·s]
K_FLUID = 0.0263  # 열전도도 [W/m-K]
CP = 1005.0       # 비열 [J/kg-K]
PR = 0.707        # 프란틀 수

# 검증 케이스 Re 수
VALIDATION_RE = [2000, 5000, 10000, 15000]


def compute_sigma():
    """최소 유동 면적비 σ = A_min/A_fr 계산"""
    g_T = S_T - D_O
    S_d = math.sqrt((S_T / 2.0)**2 + S_L**2)
    g_D = S_d - D_O
    g_min = min(g_T, g_D)
    sigma = g_min / S_T
    return sigma


def briggs_young_nu(Re_d: float) -> float:
    """
    Briggs & Young (1963) 열전달 상관식.
    
    Nu = 0.134 * Re_d^0.681 * Pr^(1/3) * (s/h_f)^0.2 * (s/t)^0.1134
    
    Parameters
    ----------
    Re_d : float
        레이놀즈 수 (V_max, d_o 기준)
    
    Returns
    -------
    float
        누셀트 수
    """
    s = S_FIN      # clear fin spacing
    h_f = H_FIN    # 핀 방사 높이
    t = T_FIN      # 핀 두께
    
    Nu = (
        0.134
        * (Re_d ** 0.681)
        * (PR ** (1.0 / 3.0))
        * (s / h_f) ** 0.2
        * (s / t) ** 0.1134
    )
    
    return Nu


def esdu_86022_kf(Re_d: float) -> float:
    """
    ESDU 86022 (1986) 압력강하 상관식.
    
    K_f = 4.567 * Re_d^(-0.242) * A_increase^0.504 * (S_T/d_o)^(-0.376) * (S_L/d_o)^(-0.546)
    
    Parameters
    ----------
    Re_d : float
        레이놀즈 수 (V_max, d_o 기준)
    
    Returns
    -------
    float
        마찰계수 K_f (열당)
    """
    # A_increase 계산 (핀 포함 표면적 / 베어튜브 표면적)
    r_i = D_O / 2.0
    r_o = r_i + H_FIN
    
    # 단일 핀 면적
    A_fin_single = (
        2.0 * math.pi * (r_o**2 - r_i**2)  # 양면
        + 2.0 * math.pi * r_o * T_FIN       # 외주
    )
    
    # 튜브 길이당 베어튜브 면적
    L_tube = N_FINS * P_FIN
    A_bare_per_tube = math.pi * D_O * L_tube
    
    # 핀 포함 전체 면적
    A_fin_total_per_tube = N_FINS * A_fin_single
    A_total_per_tube = A_bare_per_tube + A_fin_total_per_tube
    
    A_increase = A_total_per_tube / A_bare_per_tube
    
    # K_f 계산
    K_f = (
        4.567
        * (Re_d ** -0.242)
        * (A_increase ** 0.504)
        * (S_T / D_O) ** -0.376
        * (S_L / D_O) ** -0.546
    )
    
    return K_f


def compute_j_factor(Nu: float, Re_d: float) -> float:
    """Colburn j-factor 계산: j = Nu / (Re * Pr^(1/3))"""
    return Nu / (Re_d * (PR ** (1.0 / 3.0)))


def compute_h_from_nu(Nu: float) -> float:
    """열전달계수 h 계산: h = Nu * k / d_o"""
    return Nu * K_FLUID / D_O


def compute_V_max_from_Re(Re_d: float) -> float:
    """Re_d에서 V_max 계산"""
    return Re_d * MU / (RHO * D_O)


def compute_V_fr_from_Re(Re_d: float) -> tuple:
    """Re_d에서 V_fr, V_max, sigma 계산"""
    sigma = compute_sigma()
    V_max = compute_V_max_from_Re(Re_d)
    V_fr = V_max * sigma
    return V_fr, V_max, sigma


def main():
    """검증 케이스들에 대한 상관식 예측값 계산"""
    
    print("=" * 60)
    print("Correlation Predictions for Validation Cases")
    print("=" * 60)
    
    # 형상 정보 출력
    print("\n--- 형상 파라미터 ---")
    print(f"d_o:    {D_O*1000:.1f} mm")
    print(f"S_T:    {S_T*1000:.1f} mm (S_T/d_o = {S_T/D_O:.3f})")
    print(f"S_L:    {S_L*1000:.1f} mm (S_L/d_o = {S_L/D_O:.3f})")
    print(f"h_fin:  {H_FIN*1000:.1f} mm")
    print(f"t_fin:  {T_FIN*1000:.2f} mm")
    print(f"S_fin:  {S_FIN*1000:.1f} mm (s/d_o = {S_FIN/D_O:.3f})")
    
    sigma = compute_sigma()
    print(f"\nσ (A_min/A_fr): {sigma:.4f}")
    
    # 유체 물성 출력
    print("\n--- 유체 물성 (공기 300K) ---")
    print(f"ρ:  {RHO:.3f} kg/m³")
    print(f"μ:  {MU:.3e} Pa·s")
    print(f"k:  {K_FLUID:.4f} W/m-K")
    print(f"Pr: {PR:.3f}")
    
    # 각 Re에 대해 상관식 계산
    results = []
    
    print("\n--- 상관식 예측 결과 ---")
    print(f"{'Re_d':>8} {'V_fr':>8} {'V_max':>8} {'Nu':>10} {'K_f':>10} {'j':>10} {'h':>10}")
    print(f"{'':>8} {'[m/s]':>8} {'[m/s]':>8} {'[-]':>10} {'[-]':>10} {'[-]':>10} {'[W/m²K]':>10}")
    print("-" * 78)
    
    for Re_d in VALIDATION_RE:
        V_fr, V_max, _ = compute_V_fr_from_Re(Re_d)
        Nu = briggs_young_nu(Re_d)
        K_f = esdu_86022_kf(Re_d)
        j = compute_j_factor(Nu, Re_d)
        h = compute_h_from_nu(Nu)
        
        print(f"{Re_d:>8.0f} {V_fr:>8.3f} {V_max:>8.3f} {Nu:>10.3f} {K_f:>10.5f} {j:>10.5f} {h:>10.2f}")
        
        results.append({
            "Re_d": Re_d,
            "V_fr": V_fr,
            "V_max": V_max,
            "Nu": Nu,
            "K_f": K_f,
            "j": j,
            "h": h,
        })
    
    # 결과를 JSON으로 저장
    script_dir = os.path.dirname(os.path.abspath(__file__))
    study_dir = os.path.dirname(script_dir)
    output_file = os.path.join(study_dir, "correlation_predictions.json")
    
    output_data = {
        "geometry": {
            "d_o": D_O,
            "S_T": S_T,
            "S_L": S_L,
            "h_fin": H_FIN,
            "t_fin": T_FIN,
            "S_fin": S_FIN,
            "sigma": sigma,
        },
        "fluid_properties": {
            "rho": RHO,
            "mu": MU,
            "k": K_FLUID,
            "Pr": PR,
        },
        "results": results,
    }
    
    with open(output_file, "w") as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n결과 저장: {output_file}")
    
    # 적용 범위 확인
    print("\n--- 상관식 적용 범위 확인 ---")
    print(f"Briggs & Young Re_d 범위: 1,100 ~ 18,000")
    print(f"ESDU 86022 Re_d 범위:     100 ~ 100,000")
    print(f"S_T/d_o = {S_T/D_O:.3f} (범위: 1.1 ~ 4.0)")
    print(f"S_L/d_o = {S_L/D_O:.3f} (범위: 1.1 ~ 3.0)")
    print(f"s/d_o   = {S_FIN/D_O:.3f} (범위: 0.15 ~ 0.38)")
    
    for Re_d in VALIDATION_RE:
        in_range = 1100 <= Re_d <= 18000
        status = "✓" if in_range else "✗"
        print(f"Re_d = {Re_d:>5}: {status}")


if __name__ == "__main__":
    main()

