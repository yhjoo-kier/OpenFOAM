# OpenFOAM CFD Pipeline

Gmsh와 OpenFOAM을 사용한 CFD 해석 자동화 파이프라인입니다.

## 개요

이 프로젝트는 형상 생성부터 해석, 시각화까지 완전한 CFD 워크플로우를 제공합니다:

1. **형상 생성**: Gmsh Python API를 사용한 파라메트릭 형상 생성
2. **메시 생성**: blockMesh 또는 Gmsh를 통한 격자 생성
3. **CFD 해석**: OpenFOAM 솔버를 사용한 열/유동 해석
4. **후처리**: PyVista/Matplotlib를 사용한 결과 시각화

## 프로젝트 구조

```
OpenFOAM/
├── src/
│   ├── geometry/          # 형상 생성 스크립트
│   │   ├── heatsink.py    # Gmsh 기반 히트싱크
│   │   └── simple_heatsink.py  # blockMesh 기반 히트싱크
│   ├── mesh/              # 메시 변환 유틸리티
│   ├── solver/            # OpenFOAM 케이스 설정
│   └── visualization/     # 결과 시각화
├── cases/                 # 해석 케이스 디렉토리
│   └── heatsink_cht/      # 히트싱크 열해석 예제
├── scripts/               # 파이프라인 실행 스크립트
├── Dockerfile             # Docker 이미지 빌드
└── README.md
```

## 설치 방법

### 방법 1: Docker 사용 (권장)

```bash
# 이미지 빌드
docker build -t openfoam-pipeline .

# 컨테이너 실행
docker run -it --rm -v $(pwd)/results:/app/results openfoam-pipeline

# 파이프라인 실행
python scripts/run_pipeline.py
```

### 방법 2: 직접 설치 (Ubuntu 22.04/24.04)

#### 1. OpenFOAM 설치

```bash
# Ubuntu 패키지로 설치
sudo apt-get update
sudo apt-get install -y openfoam

# 또는 openfoam.org에서 설치
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get update
sudo apt-get install openfoam11
```

#### 2. Gmsh 설치

```bash
sudo apt-get install -y gmsh
```

#### 3. Python 환경 설정

```bash
# uv 설치 (권장)
curl -LsSf https://astral.sh/uv/install.sh | sh

# 가상환경 생성 및 패키지 설치
uv venv .venv
source .venv/bin/activate
uv pip install gmsh pyvista matplotlib numpy vtk
```

#### 4. 환경변수 설정

```bash
# OpenFOAM 환경변수 (bashrc에 추가 권장)
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM

# 또는 source 사용 (설치 방식에 따라 경로 다름)
source /opt/openfoam11/etc/bashrc
```

## 사용법

### 히트싱크 열전도 해석 예제

```bash
# 가상환경 활성화
source .venv/bin/activate

# 1. 케이스 생성
python src/geometry/simple_heatsink.py cases/my_heatsink

# 2. 메시 생성
cd cases/my_heatsink
blockMesh

# 3. 해석 실행
laplacianFoam

# 4. VTK 변환
foamToVTK

# 5. 시각화
cd ../..
python src/visualization/visualize_heatsink.py cases/my_heatsink -o results/
```

### 파라미터 조정

```python
# src/geometry/simple_heatsink.py 에서 파라미터 수정
create_heat_conduction_case(
    case_dir="cases/my_case",
    heat_flux=10000,      # W/m² (열유속)
    ambient_temp=300      # K (주변 온도)
)
```

### 전체 파이프라인 실행

```bash
python scripts/run_pipeline.py --case-name my_heatsink \
    --heat-flux 5000 \
    --end-time 1000
```

## 해석 케이스 설명

### 히트싱크 열전도 해석

- **솔버**: laplacianFoam (정상상태 열전도)
- **형상**: 베이스 (50x50x5 mm) + 핀 (50x10x20 mm)
- **재료**: 알루미늄 (k = 205 W/m·K)
- **경계조건**:
  - 바닥: 열유속 (5000 W/m²)
  - 상단: 고정온도 (300 K)
  - 측면: 단열 (adiabatic)

### 결과 예시

| 항목 | 값 |
|------|-----|
| 최대 온도 | 711.74 K (438.6 °C) |
| 최소 온도 | 300.00 K (26.9 °C) |
| 온도 상승 | 411.74 K |

## 확장 가능한 해석

### 1. CHT (Conjugate Heat Transfer)
- **솔버**: chtMultiRegionFoam
- 유체-고체 열전달 연성 해석
- `src/solver/cht_case_setup.py` 참조

### 2. 채널 유동 해석
- **솔버**: icoFoam (층류) / simpleFoam (난류)
- 2D/3D 내부 유동 해석
- `channelFlow/` 예제 참조

## 문제 해결

### OpenFOAM 환경변수 오류
```bash
# "Could not find mandatory etc entry 'controlDict'" 오류 시
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
```

### Gmsh-OpenFOAM 변환 오류
복잡한 형상에서 인터페이스 문제 발생 시 blockMesh 사용 권장:
```bash
# blockMesh 사용
python src/geometry/simple_heatsink.py cases/mycase
```

### PyVista 3D 렌더링 오류
Headless 환경에서는 OSMesa가 필요:
```bash
sudo apt-get install libosmesa6-dev
```
또는 matplotlib 기반 2D 시각화 사용

## 라이선스

이 프로젝트는 MIT 라이선스를 따릅니다.

## 참고 자료

- [OpenFOAM 문서](https://www.openfoam.com/documentation)
- [Gmsh 문서](https://gmsh.info/doc/texinfo/gmsh.html)
- [PyVista 문서](https://docs.pyvista.org/)
