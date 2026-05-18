# Sketchfab BIW 후보 메쉬 품질 사전 점검

## 대상

- Model: `[Libre] Automotive Body Shell Bracket`
- URL: https://sketchfab.com/3d-models/libre-automotive-body-shell-bracket-7a77cc9b1c0b47ec8e8bc22491566d1d
- Provider: SHINING 3D / EinScan
- License: Creative Commons Attribution / CC BY 4.0
- Description: Half automotive body frame (2300 × 700 × 1300 mm), scanned using Libre scanner in Laser Mode at 0.8 mm resolution.

## 현재 확보 가능한 메타데이터 기반 품질 정보

Sketchfab 페이지/API 상태에서 확인한 메타데이터:

- Source file: `Car_IR0.8_simplify1.stl`
- Source format: STL / Stereolithography
- Archive size: 171,091,494 bytes
- Source file size: 171,091,334 bytes
- Triangles: 3,421,825
- Vertices: 1,883,480
- Quad/polygon/line/point count: 0
- Has normals: false
- Has tangents: false
- Has vertex color: false
- UV mapped: false
- Material count: 1
- Texture count: 1
- Animation/rigging: none
- Model scale metadata flag: false

## CFD 관점 1차 판정

### 긍정적 요소

- 원본이 STL이므로 OpenFOAM `triSurface` 입력 형식으로 전환 부담이 작다.
- 삼각형 surface mesh라 `surfaceCheck`, `surfaceFeatureExtract`, `snappyHexMesh` workflow와 맞는다.
- 실제 스캔 기반 자동차 body frame/BIW 계열로, 도장부스 차체 프레임 장애물로 의미가 있다.
- 공개 라이선스가 CC BY 4.0라 저작자 표기 조건으로 연구 활용성이 좋다.

### 위험 요소

- 3.42M triangles는 초기 snappyHexMesh 입력으로 매우 무겁다. 1차 해석용으로는 decimation/cleanup 필요 가능성이 높다.
- `Has normals: false`라 normal 재계산 및 orientation 검사가 필요하다.
- 스캔 mesh이므로 hole, disconnected component, non-manifold edge, self-intersection 가능성이 높다.
- half body frame이라 full booth 해석에는 mirror 또는 symmetry 처리 전략이 필요하다.
- `Model scale metadata flag: false`라 STL 좌표 단위가 실제 mm인지 m인지 파일 로드 후 bounding box로 검증해야 한다.

## 다운로드 상태

공식 Sketchfab API/페이지 기준 모델은 `downloadable=true`이나, anonymous session에서는 다운로드 API가 `401 Unauthorized`를 반환했다.

- `https://api.sketchfab.com/v3/models/<uid>/download` → `401 Authentication credentials were not provided.`
- 페이지 내부 상태: `mayDownloadThisModel=false`

즉, 메쉬 파일 자체 품질 검사를 하려면 Sketchfab 로그인 상태의 공식 다운로드가 필요하다.

## 준비한 품질 검사 스크립트

파일을 확보하면 아래 스크립트로 바로 검사할 수 있게 준비했다.

- `scripts/inspect_mesh_quality.py`

예상 실행 예:

```bash
cd ~/projects/OpenFOAM
mkdir -p data/geometry/raw/sketchfab_body_shell
# downloaded STL을 아래 위치에 둔다고 가정
# data/geometry/raw/sketchfab_body_shell/Car_IR0.8_simplify1.stl

docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  python scripts/inspect_mesh_quality.py \
    data/geometry/raw/sketchfab_body_shell/Car_IR0.8_simplify1.stl \
    --out docs/26-05-02_sketchfab_body_shell_mesh_quality.json
```

추가로 OpenFOAM 쪽에서는 다음 검사가 필요하다.

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'surfaceCheck data/geometry/raw/sketchfab_body_shell/Car_IR0.8_simplify1.stl'
```

## 다음 액션

1. Sketchfab 로그인 계정으로 모델을 공식 다운로드한다.
2. 원본 STL을 `data/geometry/raw/sketchfab_body_shell/` 아래에 둔다.
3. `scripts/inspect_mesh_quality.py`와 `surfaceCheck`를 실행한다.
4. 결과에 따라 decimation, repair, mirror/full-body reconstruction, scale normalization을 진행한다.
