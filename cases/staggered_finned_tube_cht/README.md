# Staggered circular finned-tube CHT case

This case provides a streamwise-periodic conjugate heat transfer setup for a finned-tube heat exchanger representative elementary volume (REV).

## REV 도메인 정의

Staggered 배열 열교환기의 반복 단위(REV)를 정의합니다.

### 기하학적 파라미터

| 파라미터 | 기호 | 기본값 | 설명 |
|----------|------|--------|------|
| Longitudinal pitch | S_L | 0.060 m | 유동 방향 튜브 간격 |
| Transverse pitch | S_T | 0.060 m | 횡방향 튜브 간격 |
| Fin pitch | P_fin | 0.050 m | 핀 중심간 간격 |
| Tube outer radius | R_tube_out | 0.012 m | 튜브 외경/2 |
| Tube inner radius | R_tube_in | 0.010 m | 튜브 내경/2 |
| Fin radius | R_fin | 0.022 m | 핀 외경/2 |
| Fin thickness | t_fin | 0.001 m | 핀 두께 |

### 도메인 크기

| 방향 | 크기 | 범위 |
|------|------|------|
| x (유동 방향) | 2 × S_L = 0.12 m | [0, 2×S_L] |
| y (횡방향) | S_T / 2 = 0.03 m | [0, S_T/2] |
| z (튜브 축방향) | P_fin = 0.05 m | [-P_fin/2, P_fin/2] |

### REV 단면 스키메틱 (z축 방향에서 본 모습)

```
                         Flow direction (x) -->

y = S_T/2  +--------------+---------------+--------------+
           |              |               |              |
           |              |      __       |              |
           |              |    /    \     |              |
           |              |   ( half )    |              |
           |              |    \ __ /     |              |
           |              |    (S_L,      |              |
           |              |     S_T/2)    |              |
           |              |               |              |
           |   __         |               |         __   |
y = 0      | /   |        |               |        |   \ |
           ||1/4 |        |               |        |1/4 ||
           | \__ |        |               |        |__ / |
           +--------------+---------------+--------------+
          x=0           x=S_L          x=2×S_L

           (0,0)                           (2×S_L, 0)
         quarter                            quarter
          tube                               tube
```

**튜브 배치:**
- **(0, 0)**: 왼쪽 아래 꼭지점에 사분원 (1/4 튜브)
- **(S_L, S_T/2)**: 중앙 위쪽 경계에 반원 (1/2 튜브)
- **(2×S_L, 0)**: 오른쪽 아래 꼭지점에 사분원 (1/4 튜브)

이 배치로 inlet(x=0)과 outlet(x=2×S_L)의 단면 형상이 동일하여 streamwise periodic 조건을 만족합니다.

**Fin 배치:** 각 튜브 중심의 z=0 평면(센터라인)에 원형 핀 배치

### 주기 경계조건

| 경계 쌍 | 방향 | Separation Vector |
|---------|------|-------------------|
| inlet ↔ outlet | x | (0.12, 0, 0) |
| bottom ↔ top | y | (0, 0.03, 0) |
| back ↔ front | z | (0, 0, 0.05) |

---

## 실행 절차 (2025-12-06 기준)
아래 명령어들은 `cases/staggered_finned_tube_cht`에서 실행.

1) 격자 생성 및 변환  
```bash
python generate_mesh.py --output mesh.msh
gmshToFoam mesh.msh
changeDictionary -constant -dict system/changeDictionaryDict -subDict dictionaryReplacement
splitMeshRegions -cellZones -overwrite
```

2) 경계/해석 설정 파일 생성  
```bash
python setup_openfoam.py --Ubar 1.0 --Tinlet 300 --Tsolid 320
```

3) 해석 실행 (약 1000 스텝)  
```bash
chtMultiRegionSimpleFoam > log.chtMultiRegionSimpleFoam
```

4) 후처리(Δp, LMTD, Nu 등 계산)  
```bash
python post_process.py . --Umean 1.0
```

5) VTK 내보내기 및 간단 시각화(슬라이스 PNG)  
```bash
foamToVTK -latestTime
foamToVTK -latestTime -region fluid
foamToVTK -latestTime -region solid
python - <<'PY'
from pathlib import Path
import numpy as np, pyvista as pv
pv.OFF_SCREEN = True
base = Path(".")
out = base/"VTK"/"plots"; out.mkdir(parents=True, exist_ok=True)
sl = (0,0,0); n = (0,0,1)
f = pv.read(base/"VTK"/"fluid"/"staggered_finned_tube_cht_1000"/"internal.vtu").cell_data_to_point_data()
f["U_mag"] = np.linalg.norm(f["U"], axis=1)
fs = f.slice(normal=n, origin=sl)
pv.Plotter(off_screen=True).add_mesh(fs, scalars="T", cmap="plasma").view_xy().show(screenshot=str(out/"fluid_temperature_slice.png"))
pv.Plotter(off_screen=True).add_mesh(fs, scalars="U_mag", cmap="viridis").view_xy().show(screenshot=str(out/"fluid_velocity_magnitude_slice.png"))
s = pv.read(base/"VTK"/"solid"/"staggered_finned_tube_cht_1000"/"internal.vtu").cell_data_to_point_data()
ss = s.slice(normal=n, origin=sl)
pv.Plotter(off_screen=True).add_mesh(ss, scalars="T", cmap="inferno").view_xy().show(screenshot=str(out/"solid_temperature_slice.png"))
print("plots saved to", out)
PY
```

## 주요 산출물 및 경로
- 솔버 로그: `log.chtMultiRegionSimpleFoam`
- 계산 결과 타임 디렉토리: `0`, `800`, `900`, `1000` 등
- 함수 객체 결과: `postProcessing/fluid/*`, `postProcessing/solid/*`
- 후처리 스크립트 출력: `python post_process.py . --Umean 1.0`
  - delta_p=0.0, T_inlet=307.275 K, T_outlet=307.275 K, T_wall=319.969 K, Q=1.312 W, hydraulic_diameter=0.01652 m, friction_factor=0.0
- VTK 내보내기: `VTK/`, `VTK/fluid/`, `VTK/solid/`
- 슬라이스 이미지: `VTK/plots/fluid_temperature_slice.png`, `VTK/plots/fluid_velocity_magnitude_slice.png`, `VTK/plots/solid_temperature_slice.png` (EGL/DISPLAY 경고는 무시 가능; 파일은 저장됨)

### 추가 메모
- 격자는 Gmsh에서 주기 조건으로 생성되며 `changeDictionary`에서 주기 패치 번역 벡터를 지정.
- `fvOptions/meanVelocityForce`로 부피 평균 유량을 맞추고, 함수 객체로 패치 평균 압력/온도와 벽 열유속 적분을 기록함.
