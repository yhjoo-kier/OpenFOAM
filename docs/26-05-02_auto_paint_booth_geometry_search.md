# 자동차 도장 부스 CFD용 공개 차체/BIW 형상 후보 조사

## 목적

`~/projects/OpenFOAM`의 Docker/OpenFOAM 자산을 바탕으로 자동차 도장 공정의 도장 부스 내부 CFD 해석을 시작하기 위한 공개 3D 차체 형상 후보를 조사했다. 초점은 도장 단계의 차체 프레임/BIW(body-in-white)에 가까운 형상, 또는 초기 booth airflow 검증에 쓸 수 있는 자동차 외형 CAD/STL이다.

## 검색 범위

- Web search: body-in-white, automotive body shell, car body frame, paint booth CFD car body geometry
- Sketchfab downloadable model API/page
- TUM DrivAer 공식 배포 페이지
- Harvard Dataverse DrivAerNet++
- GitHub OpenFOAM/automotive geometry repositories
- 상용 3D marketplace: CGTrader 등

## 후보 요약

### 1. Sketchfab: [Libre] Automotive Body Shell Bracket

- URL: https://sketchfab.com/3d-models/libre-automotive-body-shell-bracket-7a77cc9b1c0b47ec8e8bc22491566d1d
- 제공자: SHINING 3D / EinScan
- 설명: `Half automotive body frame (2300 × 700 × 1300 mm), scanned using Libre scanner in Laser Mode at 0.8 mm resolution.`
- 공개/다운로드: Sketchfab Download 3D Model 버튼 존재. API상 downloadable=true.
- 라이선스: CC Attribution / CC BY 4.0. 저작자 표기 필요, 상업적 사용 허용.
- 메쉬 규모: 약 3.4M triangles, 1.9M vertices.
- CFD 적합성:
  - 장점: 실제 자동차 body frame/BIW 계열에 가장 가까운 무료 공개 형상.
  - 단점: half frame이므로 전차체 해석에는 mirror/정리 필요. 스캔 mesh라 watertight 여부, non-manifold, hole, thin/sheet 구조 문제가 있을 가능성 높음.
  - 권장 용도: 초기 paint booth 내부 차체 프레임 난류/유동 장애물 해석 후보 1순위. 전처리로 decimation, surface repair, mirroring, bounding-box 정렬 필요.

### 2. CGTrader: Generic Car Body in white 3D model

- URL: https://www.cgtrader.com/3d-models/vehicle/vehicle-part/generic-car-luxury-class-body-in-white
- 설명: Highly detailed and realistic generic body in white / car frame model. Doors, hood and trunk can be opened.
- 가격: 페이지 확인 시 약 $139.44 할인가 표시.
- 파일 형식: OBJ, C4D, FBX, LWO, MAX, MA, 3DS.
- 라이선스: Royalty Free License (no AI)로 표시. 연구/논문 공개 재배포 가능성은 약관 확인 필요.
- CFD 적합성:
  - 장점: full generic BIW에 가장 가까운 고품질 상용 후보. 문/후드/트렁크 개폐 가능은 도장 공정 설정에 유리할 수 있음.
  - 단점: 유료, polygon 약 1.0M~4.48M. CFD용 watertight/simplification 필요. 라이선스가 논문/데이터 공개와 충돌할 수 있음.
  - 권장 용도: 실제 BIW 형상성이 중요한 경우 구매 검토 후보.

### 3. TUM DrivAer 공식 CAD/STL

- Geometry page: https://www.epc.ed.tum.de/en/aer/research-groups/automotive/drivaer/geometry/
- Download page: https://www.epc.ed.tum.de/en/aer/research-groups/automotive/drivaer/download/
- 제공 형식: STEP, IGES, STL.
- 공개/다운로드: 무료 다운로드 가능하나 registration 필요. 페이지상 `.step-Dateien`, `.iges-Dateien`, `.stl-Dateien` 링크가 존재.
- CFD 적합성:
  - 장점: automotive CFD에서 표준에 가까운 공개 benchmark geometry. CAD/STEP 제공으로 전처리 품질이 좋을 가능성 높음.
  - 단점: paint-stage BIW가 아니라 외부 공력용 완성차 generic geometry. cabin/opening/sheet-metal BIW 구조 없음.
  - 권장 용도: 도장부스 전체 airflow, downdraft/recirculation, car-sized blockage 효과를 빠르게 검증하는 baseline.

### 4. DrivAerNet++ / Harvard Dataverse

- GitHub: https://github.com/Mohamedelrefaie/DrivAerNet
- Dataverse: https://dataverse.harvard.edu/dataverse/DrivAerNet
- 3D Mesh dataset DOI: https://doi.org/10.7910/DVN/OYU2FG
- 라이선스: Dataverse API 확인 기준 CC BY-NC 4.0.
- 파일: DrivAerNet++ 3D Meshes는 여러 zip으로 구성. 예: `E_S_WWC_WM.zip`, `F_D_WM_WW_*.zip` 등. 파일 규모가 매우 큼(수십 GB zip 포함).
- CFD 적합성:
  - 장점: 많은 자동차 형상 variant와 CFD 데이터가 있어 ML/benchmark 연계 가능.
  - 단점: non-commercial, 대용량, paint-stage BIW가 아니라 외부 공력용 자동차 mesh.
  - 권장 용도: 향후 AI surrogate/geometry variation 연구에는 좋지만, 단일 도장부스 CFD 시작 형상으로는 과함.

### 5. GitHub: Ahmed bluff body CFD validation

- URL: https://github.com/nathanrooy/ahmed-bluff-body-cfd
- 포함: Ahmed body CAD + STL, OpenFOAM v7 validation cases.
- Geometry files: `geometry/ahmed_25deg_m.stl`, `geometry/ahmed_35deg_m.stl`, FreeCAD files.
- 라이선스: GitHub API상 명시 license 없음.
- CFD 적합성:
  - 장점: OpenFOAM validation/case 구조가 이미 있어 booth solver/mesh workflow smoke test에 적합.
  - 단점: 자동차 BIW/도장 단계 형상이 아님.
  - 권장 용도: paint booth 본해석 전 OpenFOAM pipeline 검증용 단순 surrogate.

## 결론

현재 공개적으로 바로 확인된 것 중 도장 단계 차체 프레임/BIW에 가장 가까운 무료 후보는 Sketchfab의 `[Libre] Automotive Body Shell Bracket`이다. 다만 half scan mesh라 CFD 전처리가 필요하다.

현실적인 시작 순서는 다음과 같다.

1. `[Libre] Automotive Body Shell Bracket` 다운로드 및 mesh 품질 확인.
2. MeshLab/Blender/PyVista/Open3D로 decimation, repair, mirror/full-body reconstruction 가능성 검토.
3. OpenFOAM `snappyHexMesh`용 `constant/triSurface/*.stl`로 변환.
4. 단순 rectangular paint booth + downdraft inlet/outlet + car body wall boundary로 cold-flow RANS case 구성.
5. 전처리 난이도가 너무 높으면 DrivAer STL을 registration 후 받아 baseline geometry로 먼저 부스 유동을 구축.
6. 실제 BIW full model이 필요하면 CGTrader generic BIW 상용 모델 구매/라이선스 검토.

## 다음 액션 후보

- Sketchfab 후보를 다운로드하여 `assets/geometry/raw/` 또는 별도 ignored data directory에 보관.
- 표면 품질 리포트 작성: triangle count, bounding box, manifold/watertight 여부, connected components, hole count.
- OpenFOAM Docker 환경에서 `surfaceCheck`, `surfaceClean`, `surfaceFeatureExtract` 가능성 확인.
- 초기 도장부스 domain 치수/경계조건을 정의하는 case generator 작성.
