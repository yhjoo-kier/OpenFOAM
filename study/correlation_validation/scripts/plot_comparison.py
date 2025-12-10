"""
Generate comparison plots: Correlation vs CFD results.

이 스크립트는 상관식 예측값과 CFD 결과를 비교하는 그래프를 생성합니다.
- 그래프 1: Re vs Nu (Correlation, CFD scatter)
- 그래프 2: Re vs K_f (Correlation, CFD scatter)
"""
import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# 한글 폰트 설정 시도
try:
    plt.rcParams['font.family'] = 'DejaVu Sans'
except:
    pass
plt.rcParams['axes.unicode_minus'] = False


def load_correlation_predictions(study_dir: Path) -> dict:
    """상관식 예측값 로드"""
    pred_file = study_dir / "correlation_predictions.json"
    if not pred_file.exists():
        raise FileNotFoundError(f"상관식 예측 파일이 없습니다: {pred_file}")
    
    with open(pred_file, 'r') as f:
        return json.load(f)


def load_cfd_results(study_dir: Path) -> dict:
    """CFD 결과 로드"""
    cfd_file = study_dir / "cfd_results.json"
    if not cfd_file.exists():
        raise FileNotFoundError(f"CFD 결과 파일이 없습니다: {cfd_file}")
    
    with open(cfd_file, 'r') as f:
        return json.load(f)


def create_nu_comparison_plot(corr_data: dict, cfd_data: dict, output_dir: Path) -> None:
    """
    Re vs Nu 비교 그래프 생성.
    
    x축: Re_d
    y축: Nu
    - Correlation: 선 + 마커
    - CFD: scatter
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # 상관식 데이터
    corr_Re = [r["Re_d"] for r in corr_data["results"]]
    corr_Nu = [r["Nu"] for r in corr_data["results"]]
    
    # CFD 데이터 (유효한 값만)
    cfd_Re = []
    cfd_Nu = []
    for r in cfd_data["results"]:
        if r["Nu"] is not None and not math.isnan(r["Nu"]):
            cfd_Re.append(r["Re_d"])
            cfd_Nu.append(r["Nu"])
    
    # 상관식 곡선 (연속 범위)
    Re_range = np.linspace(1100, 18000, 100)
    # Briggs & Young correlation
    d_o = corr_data["geometry"]["d_o"]
    s = corr_data["geometry"]["S_fin"]
    h_f = 0.010  # H_FIN
    t = corr_data["geometry"]["t_fin"]
    Pr = corr_data["fluid_properties"]["Pr"]
    
    Nu_curve = 0.134 * (Re_range ** 0.681) * (Pr ** (1/3)) * (s/h_f)**0.2 * (s/t)**0.1134
    
    # 플롯
    ax.plot(Re_range, Nu_curve, 'b-', linewidth=2, label='Correlation (Briggs & Young)', alpha=0.7)
    ax.scatter(corr_Re, corr_Nu, c='blue', s=150, marker='o', edgecolors='darkblue',
               linewidths=2, label='Correlation Points', zorder=5)
    
    if cfd_Nu:
        ax.scatter(cfd_Re, cfd_Nu, c='red', s=200, marker='s', edgecolors='darkred',
                   linewidths=2, label='CFD Results', zorder=6)
        
        # 오차 표시
        for i, (re, nu_cfd) in enumerate(zip(cfd_Re, cfd_Nu)):
            # 대응하는 상관식 값 찾기
            for j, (re_corr, nu_corr) in enumerate(zip(corr_Re, corr_Nu)):
                if abs(re - re_corr) < 100:
                    error = (nu_cfd - nu_corr) / nu_corr * 100
                    ax.annotate(f'{error:+.1f}%', 
                               (re, nu_cfd), 
                               textcoords="offset points", 
                               xytext=(10, 10),
                               fontsize=10,
                               color='red')
                    break
    
    ax.set_xlabel('Reynolds Number (Re$_d$)', fontsize=14)
    ax.set_ylabel('Nusselt Number (Nu)', fontsize=14)
    ax.set_title('Heat Transfer Correlation Validation\nBriggs & Young (1963) vs CFD', fontsize=16)
    ax.legend(loc='lower right', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 20000)
    ax.set_ylim(0, max(corr_Nu + cfd_Nu if cfd_Nu else corr_Nu) * 1.2)
    
    # 적용 범위 표시
    ax.axvline(x=1100, color='gray', linestyle='--', alpha=0.5, label='_nolegend_')
    ax.axvline(x=18000, color='gray', linestyle='--', alpha=0.5, label='_nolegend_')
    ax.text(1100, ax.get_ylim()[1]*0.95, 'Valid Range', fontsize=10, 
            color='gray', ha='left', va='top')
    
    plt.tight_layout()
    
    output_file = output_dir / "nu_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Nu 비교 그래프 저장: {output_file}")


def create_kf_comparison_plot(corr_data: dict, cfd_data: dict, output_dir: Path) -> None:
    """
    Re vs K_f 비교 그래프 생성.
    
    x축: Re_d
    y축: K_f
    - Correlation: 선 + 마커
    - CFD: scatter
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # 상관식 데이터
    corr_Re = [r["Re_d"] for r in corr_data["results"]]
    corr_Kf = [r["K_f"] for r in corr_data["results"]]
    
    # CFD 데이터 (유효한 값만)
    cfd_Re = []
    cfd_Kf = []
    for r in cfd_data["results"]:
        if r["K_f"] is not None and not math.isnan(r["K_f"]):
            cfd_Re.append(r["Re_d"])
            cfd_Kf.append(r["K_f"])
    
    # 상관식 곡선 (연속 범위)
    Re_range = np.linspace(1100, 18000, 100)
    
    # ESDU 86022 correlation parameters
    d_o = corr_data["geometry"]["d_o"]
    S_T = corr_data["geometry"]["S_T"]
    S_L = corr_data["geometry"]["S_L"]
    h_fin = 0.010
    t_fin = corr_data["geometry"]["t_fin"]
    N_fins = 5
    P_fin = 0.0045
    
    # A_increase 계산
    r_i = d_o / 2.0
    r_o = r_i + h_fin
    A_fin_single = 2.0 * math.pi * (r_o**2 - r_i**2) + 2.0 * math.pi * r_o * t_fin
    L_tube = N_fins * P_fin
    A_bare = math.pi * d_o * L_tube
    A_total = A_bare + N_fins * A_fin_single
    A_increase = A_total / A_bare
    
    Kf_curve = (4.567 * (Re_range ** -0.242) * (A_increase ** 0.504) 
                * (S_T/d_o) ** -0.376 * (S_L/d_o) ** -0.546)
    
    # 플롯
    ax.plot(Re_range, Kf_curve, 'b-', linewidth=2, label='Correlation (ESDU 86022)', alpha=0.7)
    ax.scatter(corr_Re, corr_Kf, c='blue', s=150, marker='o', edgecolors='darkblue',
               linewidths=2, label='Correlation Points', zorder=5)
    
    if cfd_Kf:
        ax.scatter(cfd_Re, cfd_Kf, c='red', s=200, marker='s', edgecolors='darkred',
                   linewidths=2, label='CFD Results', zorder=6)
        
        # 오차 표시
        for i, (re, kf_cfd) in enumerate(zip(cfd_Re, cfd_Kf)):
            # 대응하는 상관식 값 찾기
            for j, (re_corr, kf_corr) in enumerate(zip(corr_Re, corr_Kf)):
                if abs(re - re_corr) < 100:
                    error = (kf_cfd - kf_corr) / kf_corr * 100
                    ax.annotate(f'{error:+.1f}%', 
                               (re, kf_cfd), 
                               textcoords="offset points", 
                               xytext=(10, 10),
                               fontsize=10,
                               color='red')
                    break
    
    ax.set_xlabel('Reynolds Number (Re$_d$)', fontsize=14)
    ax.set_ylabel('Friction Factor (K$_f$)', fontsize=14)
    ax.set_title('Pressure Drop Correlation Validation\nESDU 86022 (1986) vs CFD', fontsize=16)
    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 20000)
    ax.set_ylim(0, max(corr_Kf + cfd_Kf if cfd_Kf else corr_Kf) * 1.5)
    
    # 적용 범위 표시
    ax.axvline(x=1100, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=18000, color='gray', linestyle='--', alpha=0.5)
    ax.text(1100, ax.get_ylim()[1]*0.95, 'Valid Range', fontsize=10, 
            color='gray', ha='left', va='top')
    
    plt.tight_layout()
    
    output_file = output_dir / "kf_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"K_f 비교 그래프 저장: {output_file}")


def create_combined_comparison_plot(corr_data: dict, cfd_data: dict, output_dir: Path) -> None:
    """
    Nu와 K_f를 함께 보여주는 2x1 서브플롯 생성.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # 상관식 데이터
    corr_Re = [r["Re_d"] for r in corr_data["results"]]
    corr_Nu = [r["Nu"] for r in corr_data["results"]]
    corr_Kf = [r["K_f"] for r in corr_data["results"]]
    
    # CFD 데이터
    cfd_Re_Nu, cfd_Nu = [], []
    cfd_Re_Kf, cfd_Kf = [], []
    
    for r in cfd_data["results"]:
        if r["Nu"] is not None and not math.isnan(r["Nu"]):
            cfd_Re_Nu.append(r["Re_d"])
            cfd_Nu.append(r["Nu"])
        if r["K_f"] is not None and not math.isnan(r["K_f"]):
            cfd_Re_Kf.append(r["Re_d"])
            cfd_Kf.append(r["K_f"])
    
    # === Nu 플롯 (왼쪽) ===
    Re_range = np.linspace(1100, 18000, 100)
    d_o = corr_data["geometry"]["d_o"]
    s = corr_data["geometry"]["S_fin"]
    h_f = 0.010
    t = corr_data["geometry"]["t_fin"]
    Pr = corr_data["fluid_properties"]["Pr"]
    
    Nu_curve = 0.134 * (Re_range ** 0.681) * (Pr ** (1/3)) * (s/h_f)**0.2 * (s/t)**0.1134
    
    ax1.plot(Re_range, Nu_curve, 'b-', linewidth=2, alpha=0.7, label='Correlation')
    ax1.scatter(corr_Re, corr_Nu, c='blue', s=120, marker='o', edgecolors='darkblue',
                linewidths=2, zorder=5)
    
    if cfd_Nu:
        ax1.scatter(cfd_Re_Nu, cfd_Nu, c='red', s=180, marker='s', edgecolors='darkred',
                    linewidths=2, label='CFD', zorder=6)
    
    ax1.set_xlabel('Re$_d$', fontsize=14)
    ax1.set_ylabel('Nu', fontsize=14)
    ax1.set_title('(a) Nusselt Number\nBriggs & Young (1963)', fontsize=14)
    ax1.legend(loc='lower right', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 20000)
    ax1.axvline(x=1100, color='gray', linestyle='--', alpha=0.4)
    ax1.axvline(x=18000, color='gray', linestyle='--', alpha=0.4)
    
    # === K_f 플롯 (오른쪽) ===
    S_T = corr_data["geometry"]["S_T"]
    S_L = corr_data["geometry"]["S_L"]
    t_fin = corr_data["geometry"]["t_fin"]
    N_fins = 5
    P_fin = 0.0045
    
    r_i = d_o / 2.0
    r_o = r_i + h_f
    A_fin_single = 2.0 * math.pi * (r_o**2 - r_i**2) + 2.0 * math.pi * r_o * t_fin
    L_tube = N_fins * P_fin
    A_bare = math.pi * d_o * L_tube
    A_total = A_bare + N_fins * A_fin_single
    A_increase = A_total / A_bare
    
    Kf_curve = (4.567 * (Re_range ** -0.242) * (A_increase ** 0.504) 
                * (S_T/d_o) ** -0.376 * (S_L/d_o) ** -0.546)
    
    ax2.plot(Re_range, Kf_curve, 'b-', linewidth=2, alpha=0.7, label='Correlation')
    ax2.scatter(corr_Re, corr_Kf, c='blue', s=120, marker='o', edgecolors='darkblue',
                linewidths=2, zorder=5)
    
    if cfd_Kf:
        ax2.scatter(cfd_Re_Kf, cfd_Kf, c='red', s=180, marker='s', edgecolors='darkred',
                    linewidths=2, label='CFD', zorder=6)
    
    ax2.set_xlabel('Re$_d$', fontsize=14)
    ax2.set_ylabel('K$_f$', fontsize=14)
    ax2.set_title('(b) Friction Factor\nESDU 86022 (1986)', fontsize=14)
    ax2.legend(loc='upper right', fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 20000)
    ax2.axvline(x=1100, color='gray', linestyle='--', alpha=0.4)
    ax2.axvline(x=18000, color='gray', linestyle='--', alpha=0.4)
    
    plt.tight_layout()
    
    output_file = output_dir / "combined_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"통합 비교 그래프 저장: {output_file}")


def create_error_summary_table(corr_data: dict, cfd_data: dict, output_dir: Path) -> None:
    """오차 요약 테이블 생성"""
    
    lines = []
    lines.append("=" * 80)
    lines.append("Correlation vs CFD Comparison Summary")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"{'Re_d':>8} {'Nu_corr':>10} {'Nu_CFD':>10} {'Nu_err%':>10} "
                 f"{'Kf_corr':>12} {'Kf_CFD':>12} {'Kf_err%':>10}")
    lines.append("-" * 80)
    
    corr_results = {r["Re_d"]: r for r in corr_data["results"]}
    
    for cfd_r in cfd_data["results"]:
        Re = cfd_r["Re_d"]
        corr_r = corr_results.get(Re, {})
        
        Nu_corr = corr_r.get("Nu", float('nan'))
        Nu_cfd = cfd_r.get("Nu", float('nan'))
        Kf_corr = corr_r.get("K_f", float('nan'))
        Kf_cfd = cfd_r.get("K_f", float('nan'))
        
        if Nu_cfd is None:
            Nu_cfd = float('nan')
        if Kf_cfd is None:
            Kf_cfd = float('nan')
        
        Nu_err = (Nu_cfd - Nu_corr) / Nu_corr * 100 if not math.isnan(Nu_cfd) and not math.isnan(Nu_corr) else float('nan')
        Kf_err = (Kf_cfd - Kf_corr) / Kf_corr * 100 if not math.isnan(Kf_cfd) and not math.isnan(Kf_corr) else float('nan')
        
        Nu_corr_s = f"{Nu_corr:.3f}" if not math.isnan(Nu_corr) else "N/A"
        Nu_cfd_s = f"{Nu_cfd:.3f}" if not math.isnan(Nu_cfd) else "N/A"
        Nu_err_s = f"{Nu_err:+.1f}" if not math.isnan(Nu_err) else "N/A"
        Kf_corr_s = f"{Kf_corr:.5f}" if not math.isnan(Kf_corr) else "N/A"
        Kf_cfd_s = f"{Kf_cfd:.5f}" if not math.isnan(Kf_cfd) else "N/A"
        Kf_err_s = f"{Kf_err:+.1f}" if not math.isnan(Kf_err) else "N/A"
        
        lines.append(f"{Re:>8.0f} {Nu_corr_s:>10} {Nu_cfd_s:>10} {Nu_err_s:>10} "
                     f"{Kf_corr_s:>12} {Kf_cfd_s:>12} {Kf_err_s:>10}")
    
    lines.append("-" * 80)
    lines.append("")
    
    output_file = output_dir / "comparison_summary.txt"
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"비교 요약 저장: {output_file}")
    
    # 콘솔에도 출력
    for line in lines:
        print(line)


def main():
    """메인 함수"""
    print("=" * 60)
    print("Generating Comparison Plots")
    print("=" * 60)
    
    script_dir = Path(__file__).resolve().parent
    study_dir = script_dir.parent
    output_dir = study_dir / "report" / "images"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 데이터 로드
    try:
        corr_data = load_correlation_predictions(study_dir)
        print("상관식 예측 데이터 로드 완료")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("먼저 calculate_correlation.py를 실행하세요.")
        return
    
    try:
        cfd_data = load_cfd_results(study_dir)
        print("CFD 결과 데이터 로드 완료")
    except FileNotFoundError as e:
        print(f"Warning: {e}")
        print("CFD 결과가 없어 상관식 데이터만으로 그래프를 생성합니다.")
        cfd_data = {"results": []}
    
    # 그래프 생성
    print("\n그래프 생성 중...")
    create_nu_comparison_plot(corr_data, cfd_data, output_dir)
    create_kf_comparison_plot(corr_data, cfd_data, output_dir)
    create_combined_comparison_plot(corr_data, cfd_data, output_dir)
    create_error_summary_table(corr_data, cfd_data, study_dir / "report")
    
    print("\n모든 그래프 생성 완료!")


if __name__ == "__main__":
    main()

