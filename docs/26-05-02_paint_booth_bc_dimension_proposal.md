# Paint Booth CFD Boundary Condition 및 Dimension 제안

Date: 2026-05-02 KST

## 1. 기준 형상

현재 기준 형상은 원안 smooth single-shell car body다.

- Car length: 4.5 m
- Car width: 1.8 m
- Car z range: 0.18 ~ 1.61095 m
- Triangles: 8,960
- Surface: closed, single connected component

현재 계산영역은 다음과 같다.

- x: -1.5 ~ 6.5 m
- y: -2.0 ~ 2.0 m
- z: 0.0 ~ 3.0 m
- Domain size: 8.0 m × 4.0 m × 3.0 m

현재 차체 위치 기준 clearance:

- Upstream/front clearance: 1.5 m
- Downstream/rear clearance: 2.0 m
- Lateral clearance each side: 1.1 m
- Roof clearance above car: 약 1.389 m
- Ground clearance: 0.18 m

이 크기는 초기 paint booth cold-flow baseline으로는 적절하다. 실제 full-size booth보다는 축약형이지만, 차체 주변 downdraft/recirculation 경향을 보기에 충분한 시작점이다.

## 2. 권장 booth dimension

### Baseline domain

우선 현재 도메인을 그대로 유지한다.

- Booth length: 8.0 m
- Booth width: 4.0 m
- Booth height: 3.0 m

이유:

- 현재 mesh/solver 검증 완료
- 차체 전후/좌우/상부 clearance가 과도하게 작지 않음
- grid study와 boundary-condition 비교에 적합
- full-size industrial booth로 확장하기 전 계산비용 제어 가능

### 향후 확장 domain

후속 고정밀 단계에서는 다음 크기를 고려한다.

- Length: 9.0 ~ 10.0 m
- Width: 4.5 ~ 5.0 m
- Height: 3.2 ~ 3.5 m

확장 이유:

- downstream wake 회복길이 증가
- 측벽 영향 감소
- ceiling plenum/filter와 car roof 간 거리 확보

하지만 현 단계에서는 8.0 × 4.0 × 3.0 m를 기준으로 진행하는 것이 좋다.

## 3. Stage 1 권장 boundary condition: full-floor downdraft baseline

목적은 가장 단순하고 안정적인 도장부스 downdraft baseline을 만드는 것이다.

### Patch 구성

기존 blockMesh patch를 다음처럼 재정의한다.

- `ceilingInlet`
  - 위치: z = 3.0 m 전체 천장면
  - 면적: 8.0 × 4.0 = 32.0 m²
  - patch type: `patch`

- `floorExhaust`
  - 위치: z = 0.0 m 전체 바닥면
  - 면적: 8.0 × 4.0 = 32.0 m²
  - patch type: `patch`

- `frontWall`
  - 위치: x = -1.5 m
  - patch type: `wall`

- `rearWall`
  - 위치: x = 6.5 m
  - patch type: `wall`

- `sideWalls`
  - 위치: y = ±2.0 m
  - patch type: `wall`

- `carBody`
  - snappyHexMesh에서 생성
  - patch type: `wall`

### 0/U

- `ceilingInlet`
  - type: `fixedValue`
  - value: `uniform (0 0 -0.35)`

- `floorExhaust`
  - type: `pressureInletOutletVelocity`
  - value: `uniform (0 0 0)`

- `frontWall`, `rearWall`, `sideWalls`, `carBody`
  - type: `noSlip`

### 0/p

- `ceilingInlet`
  - type: `zeroGradient`

- `floorExhaust`
  - type: `fixedValue`
  - value: `uniform 0`

- walls and `carBody`
  - type: `zeroGradient`

### turbulence fields

For `kOmegaSST`:

- `ceilingInlet`
  - `k`: fixedValue
  - `omega`: fixedValue
- `floorExhaust`
  - `k`: inletOutlet or zeroGradient
  - `omega`: inletOutlet or zeroGradient
- walls and `carBody`
  - `k`: `kqRWallFunction`
  - `omega`: `omegaWallFunction`
  - `nut`: `nutkWallFunction`

초기 inlet turbulence estimate:

- Turbulence intensity: 5%
- Length scale: 0.1 ~ 0.3 m

처음에는 기존 값과 비슷한 안정적 fixedValue로 시작하고, 후속 단계에서 inlet turbulence를 parameter화한다.

### 유량 규모

Ceiling 전체를 inlet으로 쓰면 면적은 32.0 m²다.

- U = 0.30 m/s → Q = 9.6 m³/s
- U = 0.35 m/s → Q = 11.2 m³/s
- U = 0.50 m/s → Q = 16.0 m³/s

초기 baseline은 `0.35 m/s`를 추천한다.

## 4. Stage 2 권장 boundary condition: ceiling filter panel + floor grating

Stage 1의 full-ceiling/full-floor는 수치적으로 안정적인 첫 모델이다. 그 다음은 실제 부스에 가까운 patch 분할이 필요하다.

### Ceiling filter panel

천장 전체가 아니라 가장자리 0.4 ~ 0.5 m를 제외한 중앙 filter panel을 inlet으로 둔다.

권장 크기:

- x range: -1.0 ~ 6.0 m
- y range: -1.6 ~ 1.6 m
- z = 3.0 m
- Filter length: 7.0 m
- Filter width: 3.2 m
- Filter area: 22.4 m²

Boundary:

- `ceilingFilterInlet`
  - type: `patch`
  - U: `fixedValue uniform (0 0 -0.35)`
  - p: `zeroGradient`

- `ceilingSolid`
  - type: `wall`
  - U: `noSlip`
  - p: `zeroGradient`

유량:

- U = 0.30 m/s → Q = 6.72 m³/s
- U = 0.35 m/s → Q = 7.84 m³/s
- U = 0.50 m/s → Q = 11.2 m³/s

### Floor grating option A: center slot exhaust

권장 시작안:

- x range: -0.75 ~ 5.25 m
- y range: -0.60 ~ 0.60 m
- z = 0.0 m
- Slot length: 6.0 m
- Slot width: 1.2 m
- Area: 7.2 m²

Boundary:

- `floorExhaustCenter`
  - type: `patch`
  - U: `pressureInletOutletVelocity`
  - p: `fixedValue uniform 0`

- `floorSolid`
  - type: `wall`
  - U: `noSlip`
  - p: `zeroGradient`

주의:

- Ceiling filter에서 Q = 7.84 m³/s가 들어오면, center slot의 평균 outlet velocity scale은 약 1.09 m/s다.
- 즉 inlet downdraft 0.35 m/s보다 outlet에서 훨씬 빠른 흡입이 발생한다.
- 이는 floor grating/duct 흡입을 단순화한 모델로는 가능하지만, 차체 하부 유동을 강하게 만들 수 있다.

### Floor grating option B: side slot exhaust

권장 시작안:

- Left slot:
  - x range: -0.75 ~ 5.25 m
  - y range: -1.45 ~ -1.00 m
  - width: 0.45 m
- Right slot:
  - x range: -0.75 ~ 5.25 m
  - y range: 1.00 ~ 1.45 m
  - width: 0.45 m
- Total area: 2 × 6.0 × 0.45 = 5.4 m²

Boundary:

- `floorExhaustLeft`, `floorExhaustRight`
  - type: `patch`
  - U: `pressureInletOutletVelocity`
  - p: `fixedValue uniform 0`

주의:

- Ceiling filter Q = 7.84 m³/s 기준 outlet velocity scale은 약 1.45 m/s다.
- side exhaust는 차체 하부/측면 유동을 강하게 끌 수 있으므로, center slot과 비교하는 연구 설계로 좋다.

## 5. Stage별 case 추천

### Case A: `paint_booth_downdraft_full_floor_v035`

가장 먼저 구현할 case.

- Domain: 8.0 × 4.0 × 3.0 m
- Ceiling inlet: full ceiling, U = (0 0 -0.35) m/s
- Floor outlet: full floor, p = 0
- Front/rear/side walls: wall
- Car: wall

장점:

- 가장 안정적
- mass balance 단순
- downdraft baseline 파악 가능

단점:

- 실제 grating보다 이상화됨

### Case B: `paint_booth_filter_full_floor_v035`

- Ceiling inlet: central filter panel, 7.0 × 3.2 m
- Floor outlet: full floor
- 나머지 ceiling: wall

장점:

- 실제 ceiling filter panel에 가까움
- floor outlet은 아직 단순해서 안정적

### Case C: `paint_booth_filter_center_slot_v035`

- Ceiling inlet: central filter panel, 7.0 × 3.2 m
- Floor outlet: center slot, 6.0 × 1.2 m
- Floor solid: wall

장점:

- floor grating/duct 효과를 보기 시작할 수 있음

단점:

- outlet velocity가 커져 차체 하부 유동이 과장될 수 있음

### Case D: `paint_booth_filter_side_slots_v035`

- Ceiling inlet: central filter panel
- Floor outlet: left/right side slots

장점:

- exhaust 배치가 recirculation에 미치는 영향 비교 가능

단점:

- 초기 안정성은 Case A/B보다 낮을 수 있음

## 6. 구현 우선순위

추천 순서:

1. `paint_booth_downdraft_full_floor_v035`
2. `paint_booth_filter_full_floor_v035`
3. `paint_booth_filter_center_slot_v035`
4. `paint_booth_filter_side_slots_v035`

첫 구현에서는 Case A만 만들어도 충분하다. 이후 후처리 metric이 준비되면 B/C/D를 비교한다.

## 7. OpenFOAM patch naming 제안

기존 이름을 다음처럼 정리한다.

- `inlet` → `frontWall` 또는 삭제/변경
- `outlet` → `rearWall` 또는 삭제/변경
- `ceiling` → `ceilingInlet` 또는 `ceilingFilterInlet` + `ceilingSolid`
- `floorExhaust` → `floorExhaust` 또는 `floorExhaustCenter`/`floorExhaustLeft`/`floorExhaustRight`
- `sideWalls` 유지
- `carBody` 유지

Stage 1 full-floor case에서는 다음 patch set을 추천한다.

- `ceilingInlet`
- `floorExhaust`
- `frontWall`
- `rearWall`
- `sideWalls`
- `carBody`

Stage 2 filter/slot case에서는 다음 patch set을 추천한다.

- `ceilingFilterInlet`
- `ceilingSolid`
- `floorExhaustCenter` 또는 `floorExhaustLeft`/`floorExhaustRight`
- `floorSolid`
- `frontWall`
- `rearWall`
- `sideWalls`
- `carBody`

## 8. 후처리 dimension/region 제안

차체 주변 metric을 계산할 control region을 명확히 둔다.

### Near-car control volume

- x: -0.25 ~ 4.75 m
- y: -1.25 ~ 1.25 m
- z: 0.10 ~ 2.00 m

목적:

- car 주변 downdraft uniformity
- low-speed pocket
- reverse-flow fraction

### Roof region

- x: 1.2 ~ 3.8 m
- y: -1.0 ~ 1.0 m
- z: 1.2 ~ 2.2 m

목적:

- roof 위쪽 recirculation/low-speed pocket 확인

### Hood/front region

- x: -0.25 ~ 1.5 m
- y: -1.1 ~ 1.1 m
- z: 0.3 ~ 1.5 m

### Trunk/rear wake region

- x: 3.0 ~ 5.25 m
- y: -1.1 ~ 1.1 m
- z: 0.3 ~ 1.7 m

## 9. 결론

현재 형상 기준으로는 다음 설정이 가장 적절한 첫 유동해석 설계다.

- Domain: 8.0 × 4.0 × 3.0 m 유지
- Car: current smooth single-shell
- Solver: `simpleFoam`
- Turbulence: `kOmegaSST`
- First physical case: `paint_booth_downdraft_full_floor_v035`
- Ceiling inlet: full ceiling, fixedValue U = (0 0 -0.35) m/s
- Floor exhaust: full floor, p = 0, pressureInletOutletVelocity
- Front/rear/side walls: noSlip
- Car body: noSlip

이후 ceiling filter panel과 floor grating/slot을 단계적으로 추가한다.
