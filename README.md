# OpenFOAM CFD Pipeline

OpenFOAM을 사용한 히트싱크 열유동 해석 파이프라인입니다.

## 개요

이 프로젝트는 히트싱크 냉각 해석을 위한 완전한 CFD 워크플로우를 제공합니다:

1. **메시 생성**: blockMesh를 통한 구조 격자 생성
2. **유동 해석**: simpleFoam을 사용한 정상상태 유동 해석
3. **열전달 해석**: scalarTransportFoam을 사용한 온도장 해석
4. **후처리**: PyVista/Matplotlib를 사용한 결과 시각화

## 프로젝트 구조

```
OpenFOAM/
├── cases/
│   ├── heatsink_flow/     # 히트싱크 열유동 해석 (3핀)
│   └── channelFlow/       # 관내유동 예제 (icoFoam)
├── results/
│   └── heatsink_cooling_with_fins/  # 해석 결과 이미지
├── scripts/
│   ├── setup_vm.sh        # VM 환경 셋업 스크립트
│   └── env.sh             # 환경변수 설정
├── src/
│   ├── geometry/          # 형상 생성 스크립트
│   └── visualization/     # 시각화 스크립트
└── README.md
```

## 설치 방법

### VM/CODEX 환경 자동 설정 (권장)

클라우드 VM 또는 CODEX 환경에서 빠르게 셋업하려면 제공된 스크립트를 사용합니다.

```bash
# 셋업 스크립트 실행 (root 권한 필요)
chmod +x scripts/setup_vm.sh
sudo ./scripts/setup_vm.sh
```

#### 설치되는 구성요소

| 구성요소 | 버전 | 용도 |
|----------|------|------|
| OpenFOAM | v1912 | CFD 솔버 |
| Gmsh | v4.12 | 메시 생성 |
| Python 패키지 | - | pyvista, matplotlib, numpy, vtk, scipy |

#### 설치 후 확인

```bash
# OpenFOAM 명령어 확인
which blockMesh simpleFoam scalarTransportFoam

# Python 패키지 확인
python3 -c "import pyvista, matplotlib, scipy; print('OK')"

# 환경변수 로드 (새 터미널에서)
source /etc/profile.d/openfoam.sh
```

### 직접 설치 (Ubuntu 22.04/24.04)

```bash
# OpenFOAM 설치
sudo apt-get update
sudo apt-get install -y openfoam

# Python 패키지 설치
pip3 install pyvista matplotlib numpy vtk scipy

# 환경변수 설정
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
```

## 사용법

### 히트싱크 열유동 해석 예제

```bash
cd cases/heatsink_flow

# 1. 메시 생성
blockMesh

# 2. 유동 해석 (simpleFoam)
simpleFoam

# 3. 온도장 해석 (scalarTransportFoam)
scalarTransportFoam

# 4. VTK 변환 및 시각화
foamToVTK -latestTime
python3 visualize_yz_section.py
```

결과 파일은 `results/heatsink_cooling_with_fins/`에 저장됩니다.

## 해석 케이스 설명

### 히트싱크 열유동 해석 (heatsink_flow)

3개의 사각핀이 있는 채널 내 강제대류 냉각 해석입니다.

#### 형상

| 항목 | 치수 |
|------|------|
| 채널 | 120 × 40 × 60 mm (x, y, z) |
| 핀 개수 | 3개 |
| 핀 크기 | 120 × 25 × 5 mm (길이 × 높이 × 폭) |

#### 경계조건

| 경계 | 조건 |
|------|------|
| 입구 (inlet) | 속도 0.1 m/s, 온도 25°C |
| 출구 (outlet) | zeroGradient |
| 바닥벽 (bottomWall) | 온도 80°C (등온) |
| 핀벽 (finWalls) | 온도 80°C (등온) |
| 상단벽 (topWall) | 단열 (zeroGradient) |
| 전후면 (front/back) | 대칭 (symmetry) |

#### 해석 결과

| 항목 | 값 |
|------|-----|
| 메시 셀 수 | 34,560 |
| 입구 온도 | 25.0°C |
| 출구 평균 온도 | 53.8°C |
| 온도 상승 | 28.8°C |
| 목표 범위 | 30-70°C ✓ |

#### 생성되는 시각화

- `yz_temperature_contour.png`: y-z 단면 온도 컨투어
- `centerline_temperature.png`: 중심선 온도 프로파일
- `wall_temperature_profile.png`: 벽면 온도 프로파일

### 관내유동 해석 (channelFlow)

간단한 2D 채널 유동 예제입니다.

```bash
cd cases/channelFlow
blockMesh
icoFoam
foamToVTK
```

- **솔버**: icoFoam (비압축성 층류)
- **용도**: OpenFOAM 기본 동작 확인

## 문제 해결

### OpenFOAM 환경변수 오류

```bash
# "Could not find mandatory etc entry 'controlDict'" 오류 시
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
```

### PyVista 렌더링 오류

Headless 환경에서는 OSMesa가 필요합니다:

```bash
sudo apt-get install libosmesa6-dev
```

## 라이선스

이 프로젝트는 MIT 라이선스를 따릅니다.

## 참고 자료

- [OpenFOAM 문서](https://www.openfoam.com/documentation)
- [PyVista 문서](https://docs.pyvista.org/)
