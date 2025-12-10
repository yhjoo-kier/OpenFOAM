"""
Extract CFD results from OpenFOAM simulations.

이 스크립트는 OpenFOAM CHT 시뮬레이션 결과에서 Nu와 K_f를 계산합니다.

계산 방법:
- Nu = h * d_o / k, where h = Q / (A * ΔT_lm)
- K_f: ΔP = (K_acc + N_rows * K_f) * (ρ * V_max² / 2) 에서 역산
"""
import json
import math
import os
import re
import sys
from pathlib import Path

# 형상 파라미터 (검증용) - setup_validation_case.py와 동일
D_O = 0.016       # 튜브 외경 [m]
S_T = 0.036       # 가로 피치 [m]
S_L = 0.034       # 세로 피치 [m]
H_FIN = 0.010     # 핀 방사 높이 [m]
T_FIN = 0.0005    # 핀 두께 [m]
S_FIN = 0.004     # clear fin spacing [m]
P_FIN = 0.0045    # 핀 피치 [m]
N_FINS = 5        # 핀 개수
N_REPEAT = 5      # REV 반복 횟수
N_ROWS = 2 * N_REPEAT  # 총 튜브 열 수

# 유체 물성 (공기 300K)
RHO = 1.177       # 밀도 [kg/m³]
MU = 1.846e-5     # 점성 [Pa·s]
K_FLUID = 0.0263  # 열전도도 [W/m-K]
CP = 1005.0       # 비열 [J/kg-K]
PR = 0.707        # 프란틀 수

# 검증 케이스 목록
VALIDATION_CASES = ["case_Re2000", "case_Re5000", "case_Re10000", "case_Re15000"]


def compute_sigma():
    """최소 유동 면적비 σ = A_min/A_fr 계산"""
    g_T = S_T - D_O
    S_d = math.sqrt((S_T / 2.0)**2 + S_L**2)
    g_D = S_d - D_O
    g_min = min(g_T, g_D)
    sigma = g_min / S_T
    return sigma


def read_postprocessing_value(case_dir: Path, function_name: str, field: str,
                               latest: bool = True) -> float:
    """
    OpenFOAM postProcessing 디렉토리에서 값 읽기.
    
    Parameters
    ----------
    case_dir : Path
        케이스 디렉토리 경로
    function_name : str
        함수 이름 (예: "inletPressure")
    field : str
        필드 이름 (예: "p")
    latest : bool
        최신 시간 스텝 데이터만 읽을지 여부
    
    Returns
    -------
    float
        후처리 값
    """
    post_dir = case_dir / "postProcessing" / function_name
    if not post_dir.exists():
        raise FileNotFoundError(f"후처리 디렉토리가 없습니다: {post_dir}")
    
    # 시간 디렉토리 찾기
    time_dirs = sorted([d for d in post_dir.iterdir() if d.is_dir()],
                       key=lambda x: float(x.name) if x.name.replace('.', '').isdigit() else 0)
    
    if not time_dirs:
        raise FileNotFoundError(f"시간 디렉토리가 없습니다: {post_dir}")
    
    target_dir = time_dirs[-1] if latest else time_dirs[0]
    
    # surfaceFieldValue 출력 파일 형식
    data_file = target_dir / f"surfaceFieldValue.dat"
    if not data_file.exists():
        # 대체 파일 이름 시도
        for f in target_dir.iterdir():
            if f.suffix == ".dat":
                data_file = f
                break
    
    if not data_file.exists():
        raise FileNotFoundError(f"데이터 파일이 없습니다: {target_dir}")
    
    # 마지막 라인에서 값 읽기
    with open(data_file, 'r') as f:
        lines = f.readlines()
    
    # 헤더가 아닌 마지막 데이터 라인 찾기
    for line in reversed(lines):
        line = line.strip()
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 2:
                return float(parts[-1])  # 마지막 열 값 반환
    
    raise ValueError(f"데이터를 파싱할 수 없습니다: {data_file}")


def read_case_info(case_dir: Path) -> dict:
    """케이스 정보 파일 읽기"""
    info_file = case_dir / "case_info.txt"
    if not info_file.exists():
        raise FileNotFoundError(f"케이스 정보 파일이 없습니다: {info_file}")
    
    info = {}
    with open(info_file, 'r') as f:
        for line in f:
            if ':' in line:
                key, value = line.strip().split(':', 1)
                try:
                    info[key.strip()] = float(value.strip())
                except ValueError:
                    info[key.strip()] = value.strip()
    return info


def compute_heat_transfer_area():
    """
    총 열전달 면적 계산 (튜브 + 핀).
    
    스태거드 배열에서 REV당 튜브 수:
    - Row 1: 2개의 반쪽 튜브 (y=0, y=S_T) = 1개
    - Row 2: 1개의 전체 튜브 (y=S_T/2)
    → REV당 2개 튜브 (효과적으로)
    """
    # REV당 효과적 튜브 수
    tubes_per_REV = 2.0
    total_tubes = tubes_per_REV * N_REPEAT
    
    # 튜브당 길이
    L_tube = N_FINS * P_FIN
    
    # 핀 면적 (단일 핀)
    r_i = D_O / 2.0
    r_o = r_i + H_FIN
    A_fin_single = (
        2.0 * math.pi * (r_o**2 - r_i**2)  # 양면
        + 2.0 * math.pi * r_o * T_FIN       # 외주
    )
    
    # 튜브당 핀 면적
    A_fin_per_tube = N_FINS * A_fin_single
    
    # 튜브당 베어튜브 면적 (핀 사이 노출 부분)
    # 핀 사이 간격 = S_FIN, 핀 개수 = N_FINS
    # 실제로는 핀 사이와 양 끝의 노출 부분
    exposed_length = L_tube - N_FINS * T_FIN
    A_tube_per_tube = math.pi * D_O * exposed_length
    
    # 총 면적
    A_total = total_tubes * (A_fin_per_tube + A_tube_per_tube)
    
    return A_total


def extract_results_from_vtk(case_dir: Path) -> dict:
    """
    VTK 파일 또는 OpenFOAM 필드에서 결과 추출.
    
    실제 구현에서는 PyVista를 사용하여 VTK 파일을 읽거나,
    OpenFOAM 필드 파일을 직접 파싱합니다.
    
    여기서는 후처리 함수 출력에서 값을 추출합니다.
    """
    results = {}
    
    try:
        # 입구 압력 (평균)
        p_inlet = read_postprocessing_value(case_dir, "inletPressure", "p")
        results["p_inlet"] = p_inlet
    except Exception as e:
        print(f"  입구 압력 읽기 실패: {e}")
        results["p_inlet"] = None
    
    try:
        # 출구 압력 (평균)
        p_outlet = read_postprocessing_value(case_dir, "outletPressure", "p")
        results["p_outlet"] = p_outlet
    except Exception as e:
        print(f"  출구 압력 읽기 실패: {e}")
        results["p_outlet"] = None
    
    try:
        # 입구 온도
        T_inlet = read_postprocessing_value(case_dir, "inletTemperature", "T")
        results["T_inlet"] = T_inlet
    except Exception as e:
        print(f"  입구 온도 읽기 실패: {e}")
        results["T_inlet"] = 300.0  # 기본값
    
    try:
        # 출구 온도
        T_outlet = read_postprocessing_value(case_dir, "outletTemperature", "T")
        results["T_outlet"] = T_outlet
    except Exception as e:
        print(f"  출구 온도 읽기 실패: {e}")
        results["T_outlet"] = None
    
    return results


def compute_nu_from_cfd(case_info: dict, cfd_results: dict, T_solid: float = 350.0) -> float:
    """
    CFD 결과에서 Nu 계산.
    
    방법: 에너지 균형을 통한 열전달계수 역산
    Q = m_dot * cp * (T_out - T_in)
    Q = h * A * ΔT_lm
    → h = Q / (A * ΔT_lm)
    → Nu = h * d_o / k
    """
    if cfd_results.get("T_outlet") is None:
        return float('nan')
    
    V_fr = case_info.get("V_fr", 0)
    if V_fr <= 0:
        return float('nan')
    
    # 질량 유량
    A_fr = S_T * N_FINS * P_FIN
    m_dot = RHO * V_fr * A_fr
    
    # 열전달량
    T_in = cfd_results.get("T_inlet", 300.0)
    T_out = cfd_results["T_outlet"]
    Q = m_dot * CP * (T_out - T_in)
    
    if Q <= 0:
        return float('nan')
    
    # 로그평균 온도차
    dT_in = T_solid - T_in
    dT_out = T_solid - T_out
    
    if dT_in <= 0 or dT_out <= 0:
        return float('nan')
    
    if abs(dT_in - dT_out) < 1e-6:
        dT_lm = dT_in
    else:
        dT_lm = (dT_in - dT_out) / math.log(dT_in / dT_out)
    
    # 열전달 면적
    A_total = compute_heat_transfer_area()
    
    # 열전달계수
    h = Q / (A_total * dT_lm)
    
    # Nu
    Nu = h * D_O / K_FLUID
    
    return Nu


def compute_kf_from_cfd(case_info: dict, cfd_results: dict) -> float:
    """
    CFD 결과에서 K_f 계산.
    
    ΔP = (K_acc + N_rows * K_f) * (ρ * V_max² / 2)
    K_acc = 1 + σ²
    
    → K_f = (ΔP / (ρ * V_max² / 2) - K_acc) / N_rows
    """
    p_inlet = cfd_results.get("p_inlet")
    p_outlet = cfd_results.get("p_outlet")
    
    if p_inlet is None or p_outlet is None:
        return float('nan')
    
    delta_p = p_inlet - p_outlet
    
    V_max = case_info.get("V_max", 0)
    if V_max <= 0:
        return float('nan')
    
    sigma = compute_sigma()
    K_acc = 1.0 + sigma**2
    
    # 동압
    q = 0.5 * RHO * V_max**2
    
    # K_f 역산
    total_K = delta_p / q
    K_f = (total_K - K_acc) / N_ROWS
    
    return K_f


def process_case(case_name: str, study_dir: Path) -> dict:
    """단일 케이스 처리"""
    case_dir = study_dir / case_name
    
    print(f"\n처리 중: {case_name}")
    
    if not case_dir.exists():
        print(f"  케이스 디렉토리가 없습니다: {case_dir}")
        return None
    
    try:
        case_info = read_case_info(case_dir)
    except Exception as e:
        print(f"  케이스 정보 읽기 실패: {e}")
        return None
    
    Re_d = case_info.get("Re_d", 0)
    V_fr = case_info.get("V_fr", 0)
    V_max = case_info.get("V_max", 0)
    
    print(f"  Re_d: {Re_d:.0f}")
    print(f"  V_fr: {V_fr:.3f} m/s")
    print(f"  V_max: {V_max:.3f} m/s")
    
    # CFD 결과 추출
    cfd_results = extract_results_from_vtk(case_dir)
    
    # Nu 계산
    Nu_cfd = compute_nu_from_cfd(case_info, cfd_results)
    print(f"  Nu (CFD): {Nu_cfd:.3f}" if not math.isnan(Nu_cfd) else "  Nu (CFD): N/A")
    
    # K_f 계산
    K_f_cfd = compute_kf_from_cfd(case_info, cfd_results)
    print(f"  K_f (CFD): {K_f_cfd:.5f}" if not math.isnan(K_f_cfd) else "  K_f (CFD): N/A")
    
    return {
        "case": case_name,
        "Re_d": Re_d,
        "V_fr": V_fr,
        "V_max": V_max,
        "T_inlet": cfd_results.get("T_inlet"),
        "T_outlet": cfd_results.get("T_outlet"),
        "p_inlet": cfd_results.get("p_inlet"),
        "p_outlet": cfd_results.get("p_outlet"),
        "Nu": Nu_cfd,
        "K_f": K_f_cfd,
    }


def main():
    """모든 검증 케이스에서 CFD 결과 추출"""
    
    print("=" * 60)
    print("CFD Results Extraction for Validation Cases")
    print("=" * 60)
    
    script_dir = Path(__file__).resolve().parent
    study_dir = script_dir.parent
    
    results = []
    
    for case_name in VALIDATION_CASES:
        result = process_case(case_name, study_dir)
        if result:
            results.append(result)
    
    # 결과 출력
    print("\n" + "=" * 60)
    print("CFD 결과 요약")
    print("=" * 60)
    
    if results:
        print(f"\n{'Case':<15} {'Re_d':>8} {'Nu':>10} {'K_f':>12}")
        print("-" * 47)
        for r in results:
            Nu_str = f"{r['Nu']:.3f}" if not math.isnan(r['Nu']) else "N/A"
            Kf_str = f"{r['K_f']:.5f}" if not math.isnan(r['K_f']) else "N/A"
            print(f"{r['case']:<15} {r['Re_d']:>8.0f} {Nu_str:>10} {Kf_str:>12}")
    
    # 결과를 JSON으로 저장
    output_file = study_dir / "cfd_results.json"
    
    output_data = {
        "geometry": {
            "d_o": D_O,
            "S_T": S_T,
            "S_L": S_L,
            "h_fin": H_FIN,
            "t_fin": T_FIN,
            "S_fin": S_FIN,
            "N_rows": N_ROWS,
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
        json.dump(output_data, f, indent=2, default=lambda x: None if isinstance(x, float) and math.isnan(x) else x)
    
    print(f"\n결과 저장: {output_file}")


if __name__ == "__main__":
    main()

