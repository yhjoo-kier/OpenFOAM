# Staggered Finned-Tube 상관식 CFD 검증 연구

## 개요

이 연구는 스태거드 원형 핀-튜브 열교환기에 대한 **Briggs & Young (1963) 열전달 상관식**과 **ESDU 86022 (1986) 압력강하 상관식**의 정확도를 OpenFOAM CFD 해석을 통해 검증합니다.

## 상관식 요약

### Briggs & Young (1963) 열전달 상관식

```
Nu = 0.134 * Re_d^0.681 * Pr^(1/3) * (s/h_f)^0.2 * (s/t)^0.1134
```

적용 범위:
- Re_d (V_max 기준): 1,100 ~ 18,000

### ESDU 86022 (1986) 압력강하 상관식

```
K_f = 4.567 * Re_d^(-0.242) * A_increase^0.504 * (S_T/d_o)^(-0.376) * (S_L/d_o)^(-0.546)
```

적용 범위:
- Re_d: 100 ~ 100,000
- S_T/d_o: 1.1 ~ 4.0
- S_L/d_o: 1.1 ~ 3.0
- Fin density: 4 ~ 11 fpi (s = 2.0 ~ 6.1 mm)
- d_o: 9.5 ~ 50.8 mm
- h_fin: 8.5 ~ 15.9 mm

## 검증 형상

상관식 적용 범위 중앙값을 사용한 기준 형상:

| Parameter | Value | 무차원비 |
|-----------|-------|----------|
| d_o (튜브 외경) | 16 mm | - |
| S_T (가로 피치) | 36 mm | S_T/d_o = 2.25 |
| S_L (세로 피치) | 34 mm | S_L/d_o = 2.125 |
| h_fin (핀 높이) | 10 mm | d_f/d_o = 2.25 |
| t_fin (핀 두께) | 0.5 mm | - |
| S_fin (핀 간격) | 4 mm | s/d_o = 0.25 |
| 핀 피치 | 4.5 mm | ~5.6 fpi |

## 검증 케이스

| Case | Re_d (target) | V_max | 목적 |
|------|---------------|-------|------|
| case_Re2000 | 2,000 | ~1.9 m/s | 저 Re 영역 |
| case_Re5000 | 5,000 | ~4.7 m/s | 중간 Re |
| case_Re10000 | 10,000 | ~9.4 m/s | 중-고 Re |
| case_Re15000 | 15,000 | ~14.1 m/s | 고 Re (상한 근처) |

## 무차원수 정의

CFD와 상관식에서 동일한 정의를 사용:

- **Re_d** = ρ * V_max * d_o / μ
- **V_max** = V_fr / σ (σ = A_min/A_fr, 최소 유동 면적비)
- **Nu** = h * d_o / k (h: 튜브+핀 표면 평균 열전달계수)
- **K_f**: ΔP = (K_acc + N_rows * K_f) * (ρ * V_max² / 2) 에서 역산

## 유체 물성 (공기, 300K)

- ρ = 1.177 kg/m³
- μ = 1.846e-5 Pa·s
- k = 0.0263 W/m-K
- Pr = 0.707
- cp = 1005 J/kg-K

## 실행 방법

### 1. 메시 생성 및 케이스 설정

```bash
cd /workspaces/OpenFOAM/study/correlation_validation

# 모든 케이스에 대해 메시 생성
python scripts/generate_validation_mesh.py --case case_Re2000
python scripts/generate_validation_mesh.py --case case_Re5000
python scripts/generate_validation_mesh.py --case case_Re10000
python scripts/generate_validation_mesh.py --case case_Re15000
```

### 2. OpenFOAM 솔버 실행

```bash
# 각 케이스 디렉토리에서
cd case_Re2000
gmshToFoam mesh.msh
splitMeshRegions -cellZones -overwrite
python ../scripts/setup_validation_case.py --Re 2000
chtMultiRegionSimpleFoam > log.solver 2>&1
```

### 3. 결과 후처리 및 비교

```bash
# 상관식 계산
python scripts/calculate_correlation.py

# CFD 결과 추출
python scripts/extract_cfd_results.py

# 비교 그래프 생성
python scripts/plot_comparison.py
```

## 디렉토리 구조

```
correlation_validation/
├── README.md                   # 이 파일
├── case_Re2000/                # Re=2000 검증 케이스
├── case_Re5000/                # Re=5000 검증 케이스
├── case_Re10000/               # Re=10000 검증 케이스
├── case_Re15000/               # Re=15000 검증 케이스
├── scripts/
│   ├── generate_validation_mesh.py
│   ├── setup_validation_case.py
│   ├── calculate_correlation.py
│   ├── extract_cfd_results.py
│   └── plot_comparison.py
├── report/
│   ├── validation_report.md
│   └── images/
└── results_summary.txt
```

## 참고 문헌

1. Briggs, D.E. & Young, E.H. (1963). "Convection heat transfer and pressure drop of air flowing across triangular pitch banks of finned tubes". Chemical Engineering Progress Symposium Series, 59(41), 1-10.

2. ESDU 86022 (1986). "Pressure loss during cross flow of fluids with heat transfer over plain tube banks without baffles". Engineering Sciences Data Unit.

## 관련 코드

- 상관식 구현: `/workspaces/OpenFOAM/lib/finnedTubeOpt/staggered_finned_tube_bank.py`
- 기존 CHT 케이스: `/workspaces/OpenFOAM/cases/staggered_finned_tube_full/`

