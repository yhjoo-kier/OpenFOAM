# OpenFOAM Solver Stabilization Compass

## 목적
이 문서는 **Gemini/사진 기반 자동 형상 → Gmsh → OpenFOAM** 파이프라인에서,
복잡 형상에 대한 수렴성 개선 작업을 진행할 때 사용하는 **실행 체크리스트**다.

앞으로의 작업은 이 문서를 나침반으로 삼아 아래 순서로 진행한다.

---

## 핵심 원칙
1. **발산은 예외가 아니라 정상 failure mode**로 취급한다.
2. **solver만 조정해서 해결하려 하지 않는다.**
   - geometry
   - mesh
   - patch/boundary quality
   - numerics
   - initialization
   을 함께 본다.
3. **photo/Gemini-derived geometry는 바로 solve-ready라고 가정하지 않는다.**
4. **checkMesh OK는 충분조건이 아니다.** 별도의 risk gating이 필요하다.
5. **자동 재시도는 정책 기반 ladder**로 관리한다.

---

## 전체 프레임워크 루프
1. scene generation
2. scene validation
3. geometry repair / simplification
4. meshing
5. mesh risk assessment
6. numerics policy selection
7. solver execution
8. failure classification
9. retry / escalation
10. visualization + report

---

# A. Geometry / Scene 단계 체크리스트

## A1. scene schema validity
- [ ] `indoor_cfd_scene_v1` validator 통과
- [ ] obstacle overlap 없음
- [ ] obstacle이 room 내부에 완전히 포함됨
- [ ] opening 개수/타입이 최소 요구 조건 충족

## A2. solver-friendly geometry rules
- [ ] inlet 바로 앞에 obstacle이 너무 가깝지 않음
- [ ] outlet 바로 앞에 recirculation trap이 심하지 않음
- [ ] 통로(aisle/corridor) 최소 폭 확보
- [ ] 지나치게 긴 slender obstacle row는 단순화 또는 분할
- [ ] wall과 obstacle 사이 최소 clear margin 확보
- [ ] opening 크기가 과도하게 크거나 작은 경우 제한

## A3. photo/Gemini-derived geometry 전용 rules
- [ ] 사진 기반 geometry는 **정밀 복원**이 아니라 **simulation-ready abstraction**으로 취급
- [ ] 사진 1장 기반 결과는 사람 검토 없이 바로 solve하지 않음
- [ ] rack/server-row 형태는 기본적으로 보수화된 v2/v3 geometry로 변환 후 사용
- [ ] floor/ceiling opening은 초기에는 피하고 west/east 또는 south/north 벽 opening으로 단순화 검토

---

# B. Meshing 단계 체크리스트

## B1. Gmsh generation
- [ ] `.geo_unrolled` 생성 성공
- [ ] `.msh` 생성 성공
- [ ] overlapping facets / self-intersection 없음

## B2. gmshToFoam import
- [ ] `gmshToFoam` 성공
- [ ] boundary patch 이름이 의도대로 매핑됨
- [ ] `defaultFaces` 존재 여부 확인

## B3. defaultFaces policy
- [ ] `defaultFaces == 0` ideally
- [ ] `defaultFaces > 0`이면 **strong warning 또는 fail-fast 후보**
- [ ] defaultFaces 처리 방식(`wall`, `noSlip`, `calculated`)이 현재 solver setup과 일관적인지 확인

---

# C. Mesh Risk Assessment 체크리스트

## C1. checkMesh must-read 지표
- [ ] max non-orthogonality
- [ ] severe non-orth face count
- [ ] max skewness
- [ ] min determinant / face weight 관련 이상 여부
- [ ] aspect ratio
- [ ] total tetra-only imported mesh 여부

## C2. risk classification
### Low risk
- [ ] maxNonOrtho < 70
- [ ] severe non-orth faces 거의 없음
- [ ] defaultFaces 없음

### Medium risk
- [ ] 70 <= maxNonOrtho <= 80
- [ ] severe non-orth faces 존재
- [ ] imported tetra mesh지만 전반 품질은 사용 가능

### High risk
- [ ] maxNonOrtho > 80
- [ ] severe non-orth faces 많음
- [ ] defaultFaces 존재
- [ ] long narrow passage + imported tetra + complex obstacles 조합

## C3. action mapping
- [ ] low risk → standard stabilized SIMPLE/SIMPLEC
- [ ] medium risk → additional non-ortho correction + conservative numerics
- [ ] high risk → geometry repair or transient warm-up 우선 고려

---

# D. Numerics / Solver 설정 체크리스트

## D1. pressure-velocity coupling
- [ ] SIMPLE만 고집하지 않고 SIMPLEC 옵션 검토
- [ ] `consistent yes;` 사용 가능성 검토
- [ ] pressure relaxation은 필요 이상으로 공격적이지 않은지 확인

## D2. non-orthogonality 대응
- [ ] `nNonOrthogonalCorrectors`를 mesh risk에 따라 자동 설정
  - [ ] < 70 → 0
  - [ ] 70~80 → 1
  - [ ] 80~85 → 2
  - [ ] > 85 → 3 이상 또는 geometry retry

## D3. linear solver choices
- [ ] pressure solver: GAMG 설정 점검
- [ ] smoother: GaussSeidel / DICGaussSeidel 비교 검토
- [ ] coarsest level controls 필요한지 검토
- [ ] relTol / tolerance가 너무 느슨하거나 과도하지 않은지 확인

## D4. under-relaxation
- [ ] steady-state는 under-relaxation 기반 안정화 사용
- [ ] 너무 작은 값으로만 버티지 않고 staged relaxation 설계
- [ ] relaxation은 `geometry risk + solver stage`에 따라 자동 선택

## D5. fvSchemes robustness
- [ ] bounded schemes 사용 여부 점검
- [ ] laplacian / snGrad non-orth handling 검토
- [ ] 가능 시 `relaxedNonOrthoLaplacian`, `relaxedSnGrad` 계열 활용 검토

---

# E. Initialization / Turbulence 체크리스트

## E1. 초기장
- [ ] `internalField U`를 inlet과 동일하게 둘지, 0으로 둘지 risk-based로 선택
- [ ] 복잡 형상은 초기에 `(0 0 0)` 쪽이 더 안전한지 확인
- [ ] pressure 초기장 일관성 확인

## E2. turbulence initialization
- [ ] `k`, `omega`, `nut` 초기값이 과도하지 않은지 확인
- [ ] 복잡 imported mesh에서 turbulence field가 지나치게 공격적으로 설정되지 않았는지 확인
- [ ] wall function compatibility 확인

## E3. wall distance / wall function
- [ ] `wallDist` method 명시 여부 확인
- [ ] 필요 시 `meshWave` vs `exactDistance` 비교 실험
- [ ] `yPlus` function object/postProcess를 진단 항목으로 추가할지 검토

---

# F. Execution Strategy 체크리스트

## F1. staged solve strategy
- [ ] 바로 final steady solve로 가지 않음
- [ ] baseline / conservative / robust ladder 정의
- [ ] 필요 시 transient warm-up → steady finish 구조 고려

## F2. recommended ladder
### Stage 0 — geometry-clean baseline
- [ ] 가장 단순한 geometry로 solve sanity check

### Stage 1 — conservative steady
- [ ] 낮은 inlet velocity
- [ ] internal U = 0
- [ ] stronger relaxation
- [ ] increased non-ortho correction

### Stage 2 — robust steady
- [ ] turbulence 더 완화
- [ ] solver tolerances / smoother 재조정
- [ ] non-orth 및 pressure stabilization 강화

### Stage 3 — transient warm-up
- [ ] `pimpleFoam` 또는 유사 transient warm start 짧게 수행
- [ ] 이후 `simpleFoam`로 이어받기

### Stage 4 — geometry/mesh retry
- [ ] obstacle 축소
- [ ] corridor widening
- [ ] opening 위치/크기 수정
- [ ] mesh coarsen/refine 조정

---

# G. Failure Classification 체크리스트

## G1. 로그에서 자동 탐지할 신호
- [ ] `Floating point exception`
- [ ] `FOAM FATAL ERROR`
- [ ] `bounding omega`
- [ ] `bounding k`
- [ ] `bounding U`
- [ ] continuity error 급증
- [ ] pressure solve iteration 비정상 폭증

## G2. failure type 분류
### pressure-driven failure
- [ ] p/GAMG residual 폭주
- [ ] continuity explosion 동반

### turbulence-driven failure
- [ ] `bounding omega`, `k`, `nut` 비정상
- [ ] wall-function/wallDist 의심

### mesh-import failure
- [ ] `defaultFaces`
- [ ] patch 누락/오분류
- [ ] severe non-orth / skewness

### geometry-driven failure
- [ ] 너무 좁은 통로
- [ ] inlet/outlet 바로 앞 obstacle
- [ ] long rack rows

## G3. failure → action mapping
- [ ] pressure failure → non-ortho correction 증가, pressure solver setting 재조정
- [ ] turbulence failure → laminar/transient warm-up 또는 turbulence 초기값 완화
- [ ] import failure → patch/defaultFaces 우선 수정
- [ ] geometry failure → repair/re-simplification 우선

---

# H. 산출물 / 기록 체크리스트

## H1. 각 run마다 남길 것
- [ ] scene JSON
- [ ] repaired scene JSON (if any)
- [ ] `.geo_unrolled`
- [ ] `.msh`
- [ ] `checkMesh` 요약
- [ ] solver preset별 로그
- [ ] stabilization summary JSON
- [ ] 2D preview
- [ ] 3D preview

## H2. 성공 run에서 추가
- [ ] final VTK
- [ ] comparison figure
- [ ] 3D pipeline figure
- [ ] Google Drive copy 여부 기록

---

# I. 당장 해야 할 실무 우선순위

## Priority 1 — pipeline gating 강화
- [ ] `defaultFaces`를 hard warning 이상으로 승격
- [ ] `checkMesh` 기반 risk score 도입
- [ ] geometry simplification rule을 solver-aware하게 강화

## Priority 2 — numerics ladder 재설계
- [ ] SIMPLEC preset 추가
- [ ] non-orth adaptive preset 추가
- [ ] bounded numerics preset 추가
- [ ] wallDist 옵션 추가

## Priority 3 — transient warm-start 도입
- [ ] short transient warm-up stage 설계
- [ ] warm-up → steady handoff 자동화

## Priority 4 — failure classifier 자동화
- [ ] 로그 parser에서 failure 유형 분류
- [ ] failure type별 retry policy 연결

## Priority 5 — 연구용 평가 프레임 구축
- [ ] solve success rate
- [ ] retry count
- [ ] runtime overhead
- [ ] geometry repair 횟수
- [ ] photo-derived scene success/failure statistics

---

# 현재까지의 판단
- **photo/Gemini 기반 geometry generation은 가능하다.**
- **mesh-ready geometry까지 가는 것도 가능하다.**
- 하지만 **stable OpenFOAM solve는 별도의 stabilization loop 없이는 신뢰할 수 없다.**
- 따라서 이 프레임워크의 핵심 기술은 **solver stabilization loop + geometry repair loop**다.

---

## 앞으로의 운영 원칙
앞으로 indoor/OpenFOAM 자동화 작업은,
**이 문서를 체크리스트로 삼아 단계별로 통과/실패를 기록하면서 진행한다.**
특히 새로운 photo-derived 또는 Gemini-derived geometry를 다룰 때는,
반드시:

1. geometry risk 점검
2. mesh risk 점검
3. numerics policy 선택
4. failure classification
5. retry strategy 실행

순서를 지킨다.
