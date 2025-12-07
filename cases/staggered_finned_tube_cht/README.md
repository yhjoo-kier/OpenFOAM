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
| x (유동 방향) | 2 × S_L = 0.12 m | [-S_L/2, 3×S_L/2] = [-0.03, 0.09] |
| y (횡방향) | S_T = 0.06 m | [0, S_T] |
| z (튜브 축방향) | P_fin = 0.05 m | [-P_fin/2, P_fin/2] |

**참고:** x 방향 도메인이 -S_L/2 만큼 시프트되어 inlet/outlet 면이 튜브와 교차하지 않고 깨끗한 평면을 형성합니다. 이를 통해 cyclicAMI 주기 경계조건의 면 매칭이 정확하게 이루어집니다.

**튜브 배치 (Staggered 배열):**
- **Row 1 (x=0):**
  - **(0, 0)**: y=0 경계에 반원 (1/2 튜브, 대칭 경계로 절반)
  - **(0, S_T)**: y=S_T 경계에 반원 (1/2 튜브, 대칭 경계로 절반)
- **Row 2 (x=S_L):**
  - **(S_L, S_T/2)**: 도메인 중앙에 완전한 원형 튜브 (offset by S_T/2)

이 배치로 inlet(x=-S_L/2)과 outlet(x=3×S_L/2)의 단면 형상이 동일한 평면이 되어 streamwise periodic 조건(cyclicAMI)을 정확히 만족합니다.

**Fin 배치:** 각 튜브 중심의 z=0 평면(센터라인)에 원형 핀 배치

### 경계조건

| 경계 쌍 | 방향 | 유형 | Separation Vector |
|---------|------|------|-------------------|
| inlet ↔ outlet | x | cyclicAMI | (0.12, 0, 0) |
| bottom, top | y | symmetry | - |
| back, front | z | symmetry | - |

**참고:** y, z 방향은 대칭(symmetry) 경계조건을 사용하여 계산 비용을 줄이면서도 물리적으로 타당한 결과를 얻습니다.

---

## 실행 절차 (2025-12-07 기준)

아래 명령어들은 `cases/staggered_finned_tube_cht` 디렉토리에서 실행합니다.

### 사전 준비
- OpenFOAM v1912 이상
- Gmsh (Python API 포함)
- Python 3 with: `gmsh`, `numpy`, `scipy`, `pyvista`, `matplotlib`

### 1) 격자 생성 및 변환

```bash
# Gmsh로 메시 생성 (mesh resolution 조절: --lc 옵션)
python3 generate_mesh.py --output mesh.msh

# OpenFOAM 형식으로 변환
gmshToFoam mesh.msh

# 경계조건 타입 설정 (cyclicAMI, symmetry 등)
changeDictionary

# 멀티리전 분리 (fluid + solid 영역)
# 참고: 튜브가 분리되어 있어 solid, domain0, domain1 등 여러 solid 영역 생성됨
splitMeshRegions -cellZones -overwrite
```

### 2) 초기/경계조건 설정

`splitMeshRegions` 실행 후 생성된 각 영역(fluid, solid, domain0, domain1 등)에 대해
초기조건과 물성치 파일을 설정해야 합니다.

```bash
# 자동 설정 스크립트 실행 (없는 경우 수동 설정 필요)
# python3 setup_openfoam.py --Ubar 1.0 --Tinlet 300 --Tsolid 320
```

**수동 설정 시 필요한 파일:**
- `constant/regionProperties`: fluid/solid 영역 정의
- `constant/{region}/thermophysicalProperties`: 각 영역 물성치
- `0/{region}/T`, `0/{region}/p` 등: 초기/경계조건

### 3) 해석 실행

```bash
chtMultiRegionSimpleFoam > log.chtMultiRegionSimpleFoam 2>&1 &

# 진행 상황 모니터링
tail -f log.chtMultiRegionSimpleFoam
```

기본 설정으로 약 1000 iteration 수행 (system/controlDict의 endTime 참조).

### 4) 결과 확인

```bash
# 수렴 확인
grep "Time =" log.chtMultiRegionSimpleFoam | tail -5

# 열유속 확인
grep "wallHeatFlux" log.chtMultiRegionSimpleFoam | tail -5
```

### 5) VTK 변환 및 시각화

```bash
# 각 영역별 VTK 변환
foamToVTK -latestTime -region fluid
foamToVTK -latestTime -region solid
foamToVTK -latestTime -region domain0
foamToVTK -latestTime -region domain1

# 시각화 스크립트 실행
python3 visualize_results.py
```

생성되는 이미지:
- `velocity_field.png`: 속도장 (z=0 단면)
- `temperature_field.png`: 온도장 (z=0 단면)
- `geometry_check.png`: 형상 확인용

---

## 메시 해상도 조절

`generate_mesh.py`의 `--mesh-size` 옵션으로 메시 크기를 조절할 수 있습니다:

```bash
# 기본값 (mesh-size=0.0035)
python3 generate_mesh.py --output mesh.msh

# 더 조밀한 메시 (mesh-size=0.002)
python3 generate_mesh.py --output mesh.msh --mesh-size 0.002

# 더 성긴 메시 (mesh-size=0.005)
python3 generate_mesh.py --output mesh.msh --mesh-size 0.005
```

---

## 주요 산출물

| 항목 | 경로 |
|------|------|
| 솔버 로그 | `log.chtMultiRegionSimpleFoam` |
| 계산 결과 | `100/`, `200/`, ..., `1000/` |
| VTK 파일 | `VTK/fluid/`, `VTK/solid/`, `VTK/domain0/`, `VTK/domain1/` |
| 시각화 이미지 | `velocity_field.png`, `temperature_field.png`, `geometry_check.png` |

---

## 추가 메모

- **cyclicAMI**: x 방향 streamwise periodic 경계조건으로 완전발달 유동 해석
- **symmetry**: y, z 방향 대칭 경계조건
- **다중 solid 영역**: 튜브들이 물리적으로 분리되어 있어 `splitMeshRegions`가 여러 solid 영역 생성 (solid, domain0, domain1)
- **열전달 커플링**: `turbulentTemperatureCoupledBaffleMixed` 경계조건으로 fluid-solid 인터페이스 열전달
