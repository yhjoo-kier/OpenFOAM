# Paint Booth Baseline 유동해석 전략

Date: 2026-05-02 KST

## 1. 현재 기준 형상 및 case 상태

현재 기준 형상은 `cases/paint_booth_baseline/constant/triSurface/simplified_car_shell.stl`의 원안 smooth single-shell car body다.

형상 요약:

- Length: 4.5 m
- Width: 1.8 m
- Height range: z = 0.18 ~ 1.61095 m
- Vertices: 4,482
- Triangles: 8,960
- Surface: closed, single connected component, watertight
- Strategy: `single_closed_smooth_shell`

계산 영역:

- Domain: x = -1.5 ~ 6.5 m, y = -2.0 ~ 2.0 m, z = 0 ~ 3.0 m
- Booth-like domain size: 8.0 m × 4.0 m × 3.0 m
- Vehicle centerline: y = 0
- Vehicle located roughly from x = 0 to 4.5 m

현재 smoke-test boundary condition:

- x-min `inlet`: fixedValue U = (0.35 0 0) m/s
- x-max `outlet`: pressure outlet
- `ceiling`, `floorExhaust`, `sideWalls`, `carBody`: wall/no-slip

현재 case는 OpenFOAM 검증용 skeleton이다. 실제 도장부스의 downdraft ceiling inlet / floor exhaust 물리조건은 아직 반영하지 않았다.

## 2. 해석 목적 정의

초기 유동해석의 목적은 도장부스 내 차체 주변 유동장이 다음 조건을 만족하는지 평가하는 것이다.

1. 천장 필터에서 내려오는 downdraft가 차체 주변에서 얼마나 균일하게 유지되는가?
2. 차체 hood/roof/trunk 주변에서 recirculation 또는 low-velocity pocket이 발생하는가?
3. 바닥 배기 조건이 차체 주변 체류영역을 줄이는 데 충분한가?
4. overspray scalar/particle 해석으로 확장할 수 있는 안정적인 base flow를 만들 수 있는가?

따라서 1차 목표는 detailed spray가 아니라 **steady incompressible booth airflow baseline**이다.

## 3. 권장 해석 단계

### Stage 0 — 현재 skeleton 유지 검증

목적:

- 현재 geometry, snappyHexMesh, simpleFoam workflow가 계속 정상 작동하는지 확인.

이미 완료:

- `surfaceCheck`: closed / single part / consistent zone
- `blockMesh`: success
- `snappyHexMesh -overwrite`: success
- `checkMesh`: `Mesh OK`
- `simpleFoam`: time 200 정상 종료

이 단계는 향후 변경 전 regression baseline으로 유지한다.

### Stage 1 — 실제 도장부스형 downdraft boundary condition

핵심 변경:

- `ceiling`: wall → velocity inlet
- `floorExhaust`: wall → outlet 또는 pressure outlet
- 기존 x-direction `inlet/outlet`: wall 또는 symmetry-like side boundary로 재정의

권장 초기 조건:

- Ceiling inlet velocity: U = (0 0 -0.30) ~ (0 0 -0.50) m/s
- 초기 baseline: U = (0 0 -0.35) m/s
- Floor exhaust: `pressureInletOutletVelocity` + p = 0 또는 `inletOutlet`
- Side/front/rear walls: no-slip wall 또는 far-field/symmetry 대안 비교
- Car body: no-slip wall

이 단계의 핵심 질문:

- 차체 주변 vertical velocity가 충분히 아래 방향으로 유지되는가?
- roof/hood/trunk 후류에 큰 재순환이 생기는가?
- floor outlet 전체를 열었을 때 유량 밸런스가 자연스럽게 맞는가?

권장 solver:

- `simpleFoam`
- incompressible steady RANS
- turbulence: `kOmegaSST`

### Stage 2 — floor exhaust/grating 모델 구체화

Stage 1에서는 floor 전체를 outlet으로 둘 수 있지만, 실제 도장부스는 floor grating/slot/underfloor plenum을 갖는다. 따라서 두 번째 단계에서는 floor boundary를 분할한다.

권장 patch 구성:

- `floorSolid`: wall
- `floorExhaustCenter`: outlet
- 선택: `floorExhaustLeft`, `floorExhaustRight` 또는 slot array

구현 방법:

1. `blockMeshDict`의 bottom face를 여러 patch로 분할하거나,
2. 별도 thin baffle/patch surface를 만들어 `topoSet`/`createPatch`로 patch 분리하거나,
3. generator에서 multi-block domain을 생성한다.

권장 시작은 1번, 즉 generator 기반 multi-block 또는 bottom face split이다.

비교할 설계 변수:

- Exhaust width: 0.8, 1.2, 1.6 m
- Exhaust length: full booth vs car-zone only
- Exhaust placement: car centerline below vs side slots

### Stage 3 — mesh refinement/grid independence

현재 mesh는 smoke-test용으로 충분하지만, 차체 주변 정량 평가에는 grid study가 필요하다.

권장 grid levels:

- Coarse: current 수준, car refinement level (2 3)
- Medium: base cell size 0.20 m, car refinement level (3 4)
- Fine: base cell size 0.15 m 또는 local refinement 강화

관찰 지표:

- Total cells
- Max non-orthogonality
- Max skewness
- Solver convergence
- Car 주변 평균/최소 vertical velocity
- Recirculation volume fraction
- Selected probe velocity convergence

현재 repo에 grid-test script 자산이 있으므로, paint booth case에도 동일한 batch pattern을 확장할 수 있다.

### Stage 4 — 후처리 metric 정의

도장부스 airflow 관점에서 단순 residual보다 중요한 지표를 정의해야 한다.

권장 metric:

1. **Downdraft uniformity**
   - car 주변 control volume에서 Uz 평균, 표준편차, coefficient of variation
   - target: Uz가 대부분 음수이며 공간 변동이 작을수록 좋음

2. **Low-speed pocket**
   - `|U| < 0.05 m/s` 또는 `|Uz| < 0.05 m/s` 영역 비율
   - hood/roof/trunk 주변 separate reporting

3. **Reverse-flow / recirculation indicator**
   - downdraft case에서 `Uz > 0`인 cell volume fraction
   - car wake 주변 `Ux`, `Uy`, `Uz` sign 변화 확인

4. **Near-body ventilation index**
   - carBody에서 일정 거리 이내 shell/volume의 평균 velocity magnitude
   - 추후 overspray deposition risk proxy로 사용 가능

5. **Mass balance**
   - ceiling inlet integrated flow rate
   - floor outlet integrated flow rate
   - imbalance check

후처리 방식:

- `foamToVTK` 후 PyVista 분석
- 또는 OpenFOAM `sample`, `surfaces`, `probes`, `postProcess` functionObjects 사용

### Stage 5 — scalar/particle overspray 확장

Base airflow가 안정화된 뒤 다음 중 하나로 확장한다.

1. Passive scalar transport
   - `scalarTransportFoam` 또는 custom passive scalar field
   - paint mist concentration proxy
   - 빠르고 안정적

2. Lagrangian particles
   - droplet residence time, near-body deposition tendency 평가 가능
   - 모델 복잡도 증가

3. Wall deposition proxy
   - carBody/floor/sideWalls 근방 scalar flux 또는 particle hit count

초기에는 passive scalar가 더 실용적이다.

## 4. 추천 case variant 구조

현재 `cases/paint_booth_baseline/`은 generator output이므로, 물리조건 variant를 다음처럼 나누는 것이 좋다.

- `cases/paint_booth_skeleton_streamwise/`
  - 현재 smoke-test skeleton 보존
- `cases/paint_booth_downdraft_v035/`
  - ceiling inlet 0.35 m/s, floor outlet
- `cases/paint_booth_downdraft_v050/`
  - ceiling inlet 0.50 m/s
- `cases/paint_booth_floor_slot_center/`
  - center floor exhaust slot
- `cases/paint_booth_floor_slot_side/`
  - side exhaust slot variant

또는 generator CLI에 `--flow-mode streamwise|downdraft`, `--ceiling-velocity`, `--floor-mode full|center-slot|side-slots`를 추가해서 재현성을 높이는 것이 좋다.

## 5. 즉시 실행 가능한 다음 구현 단위

가장 추천하는 다음 구현은 **Stage 1 downdraft case generator option 추가**다.

구현 항목:

1. `blockMeshDict`
   - ceiling/floor patch type을 `patch`로 변경 가능한 옵션 추가
   - front/rear/side wall 이름 정리

2. `0/U`
   - `ceiling`: fixedValue `(0 0 -0.35)`
   - `floorExhaust`: pressureInletOutletVelocity 또는 inletOutlet
   - remaining walls: noSlip
   - `carBody`: noSlip

3. `0/p`
   - `floorExhaust`: fixedValue 0
   - `ceiling`: zeroGradient
   - walls: zeroGradient

4. turbulence fields
   - ceiling inlet k/omega fixedValue 추정값 설정
   - outlet/walls는 기존 OpenFOAM 관례 유지

5. 검증
   - `blockMesh`
   - `snappyHexMesh -overwrite`
   - `checkMesh`
   - `simpleFoam`
   - `foamToVTK`
   - mid-plane slice render

## 6. 리스크 및 주의점

- Ceiling을 velocity inlet으로 바꾸고 floor를 outlet으로 바꾸면 유량 밸런스가 solver 안정성에 민감할 수 있다.
- Floor 전체 outlet은 현실적이지 않지만 첫 downdraft baseline으로는 좋다.
- Floor grating을 바로 상세하게 만들면 geometry/mesh 복잡도가 커지므로, slot patch 수준부터 시작하는 것이 좋다.
- 현재 car shell은 wheel/underbody가 없어서 하부 유동은 실제와 다르다. 그러나 초기 booth-scale downdraft/recirculation 파악에는 충분하다.
- 실제 도장부스에서는 열/습도/입자/전착/도장 로봇 등 변수가 많으므로, 현재 단계에서는 airflow-only cold-flow를 명확히 baseline으로 정의한다.

## 7. 결론

현재 smooth single-shell car body 기준으로는 다음 순서가 가장 실용적이다.

1. 현재 streamwise smoke-test case는 regression baseline으로 보존.
2. 바로 downdraft ceiling inlet + floor outlet case를 만든다.
3. 차체 주변 Uz uniformity, low-speed pocket, reverse-flow fraction을 metric으로 정의한다.
4. floor exhaust를 full outlet에서 center/side slot으로 단계적으로 구체화한다.
5. airflow baseline이 안정화되면 passive scalar 또는 Lagrangian particle로 overspray proxy를 추가한다.

즉, 다음 개발 목표는 **`paint_booth_downdraft_v035` variant 생성 및 후처리 metric pipeline 추가**가 적절하다.
