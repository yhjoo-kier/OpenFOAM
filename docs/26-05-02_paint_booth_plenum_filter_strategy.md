# Paint Booth Plenum / Filter 포함 유동해석 설계 수정안

Date: 2026-05-02 KST

## 1. 사용자 의견 반영 요약

기존 제안은 ceiling inlet에서 바로 작업공간으로 uniform downdraft를 넣는 단순화였다. 하지만 실제 자동차 도장부스에서는 상부 duct에서 들어온 유동이 곧바로 jet처럼 작업공간으로 내려오지 않는다.

일반적인 구조는 다음에 가깝다.

1. 상부 duct 또는 supply opening을 통해 공기가 상부 plenum으로 유입된다.
2. Plenum 공간에서 압력/유동이 어느 정도 완화된다.
3. Ceiling filter, diffuser, perforated plate 등 저항체를 통과한다.
4. 저항체 하부에서 비교적 균일한 downdraft가 작업공간으로 내려온다.
5. 바닥 grating/exhaust를 통해 배출된다.

따라서 CFD 모델도 단순 ceiling velocity inlet 대신, **상부 plenum + filter resistance + 작업공간**을 포함하는 방향으로 수정한다.

## 2. 수정된 domain dimension 제안

기존 작업공간은 유지하되, 상부 plenum 공간을 추가한다.

### 작업공간(work zone)

- x: -1.5 ~ 6.5 m
- y: -2.0 ~ 2.0 m
- z: 0.0 ~ 3.0 m
- Size: 8.0 m × 4.0 m × 3.0 m

현재 car 및 주변 clearance가 검증되어 있으므로, 작업공간 자체는 유지한다.

### 상부 plenum

추천 baseline:

- z: 3.0 ~ 3.8 m
- Plenum height: 0.8 m

대안:

- Compact: 0.6 m
- Baseline: 0.8 m
- More realistic/smoothing: 1.0 m

### 전체 domain

Baseline 추천:

- x: -1.5 ~ 6.5 m
- y: -2.0 ~ 2.0 m
- z: 0.0 ~ 3.8 m
- Size: 8.0 m × 4.0 m × 3.8 m

즉 기존보다 수직 방향을 0.8 m 높인다.

## 3. 권장 상부 구조

### Filter plane

Filter/diffuser는 작업공간과 plenum 사이에 둔다.

- 위치: z = 3.0 m
- Filter panel x range: -1.0 ~ 6.0 m
- Filter panel y range: -1.6 ~ 1.6 m
- Filter size: 7.0 m × 3.2 m
- Filter area: 22.4 m²

Filter panel 바깥쪽 천장/프레임 영역:

- z = 3.0 m
- filter panel 외곽 영역
- solid frame 또는 blocked/near-wall region으로 취급

### Plenum top/side supply duct

실제 부스에서 supply duct가 어디 있는지는 설계별로 다르다. 초기 모델에서는 다음 중 하나를 선택한다.

#### Option A: top supply opening

- 위치: z = 3.8 m, plenum 상단
- opening size: 1.5 m × 1.2 m
- area: 1.8 m²
- 중심: x = 2.5 m, y = 0

Filter 평균 유속 0.35 m/s 기준:

- Filter flow rate: 22.4 × 0.35 = 7.84 m³/s
- Duct velocity: 7.84 / 1.8 = 4.36 m/s

이 경우 plenum 내부에는 jet가 생기지만, filter 저항 때문에 작업공간으로는 완화된 downdraft가 내려가야 한다.

#### Option B: rear/side plenum slot supply

- 위치: plenum 후방 또는 측면
- 예: rear slot size = 0.6 m × 3.2 m
- area: 1.92 m²
- duct velocity at Q = 7.84 m³/s: 4.08 m/s

#### Option C: two side supply slots

- 좌우 각각: 0.5 m × 3.0 m
- total area: 3.0 m²
- duct velocity at Q = 7.84 m³/s: 2.61 m/s

초기 추천은 Option A 또는 C다.

- Option A: 구조가 단순하고 구현 쉬움
- Option C: plenum 내 유동이 더 대칭적이고 균일화에 유리

## 4. Filter/저항체 모델링 방식

### 권장 1차 모델: porous baffle 또는 porous layer

실제 filter를 상세 섬유 구조로 모델링하지 않고, filter plane 또는 얇은 filter layer를 압력손실 저항체로 모델링한다.

가능한 접근:

1. Porous baffle / pressure-jump face
   - filter를 내부 face/baffle로 두고 pressure jump를 준다.
   - mesh가 간단하면 좋지만 OpenFOAM patch 생성 절차가 조금 필요하다.

2. Thin porous volume layer
   - z = 2.95 ~ 3.05 m 정도의 얇은 cellZone을 만들고 fvOptions porous resistance를 준다.
   - snappy/blockMesh multi-region 또는 topoSet/cellZone 작업 필요.

3. Simplified uniform outlet/inlet pair
   - plenum 해석은 생략하고 filter 아래 boundary에 uniform velocity를 직접 부여한다.
   - 가장 쉽지만, 사용자 의견인 duct jet + filter smoothing 현상을 모사하지 못한다.

이번 수정 방향에서는 1 또는 2를 추천한다.

### 현실적 구현 우선순위

초기 구현 난이도와 안정성을 고려하면:

1. `plenum + filter outlet plane을 porous baffle로 구현`
2. 어려우면 `thin porous volume layer`로 대체
3. 그래도 복잡하면 임시로 `filter plane 아래 uniform downdraft`를 regression baseline으로 유지

## 5. Pressure drop / resistance scale

초기 filter target pressure drop은 다음 범위를 추천한다.

- Low resistance: 30 Pa at 0.35 m/s
- Baseline: 50 ~ 80 Pa at 0.35 m/s
- High resistance: 100 ~ 150 Pa at 0.35 m/s

동압 기준으로 무차원 pressure-loss coefficient를 잡으면:

\[
\Delta p = K \cdot \frac{1}{2}\rho U^2
\]

공기 밀도 ρ = 1.2 kg/m³, U = 0.35 m/s일 때 동압은 약 0.0735 Pa다. 따라서 filter의 K는 매우 크게 나온다.

- Δp = 30 Pa → K ≈ 408
- Δp = 50 Pa → K ≈ 680
- Δp = 80 Pa → K ≈ 1088
- Δp = 100 Pa → K ≈ 1361
- Δp = 150 Pa → K ≈ 2041

이 값이 큰 것은 정상이다. filter/HEPA/pre-filter는 매우 낮은 속도에서도 수십 Pa 압력손실을 만드는 저항체다.

초기 baseline 추천:

- Target filter velocity: 0.35 m/s
- Target pressure drop: 50 Pa
- Equivalent K: 약 680

단, OpenFOAM 구현 방식에 따라 Darcy/Forchheimer coefficient 입력 형태가 다르므로, 실제 coefficient는 사용하는 모델에 맞춰 변환해야 한다.

## 6. 수정된 Stage별 case 제안

### Case 0: current skeleton

- `paint_booth_skeleton_streamwise`
- 목적: geometry/snappy/simpleFoam regression
- 현재 완료 상태

### Case 1: idealized downdraft without plenum

- `paint_booth_downdraft_full_floor_v035`
- 목적: 가장 단순한 top-down baseline
- ceiling에서 바로 uniform U = (0 0 -0.35)
- floor 전체 outlet

이 case는 여전히 필요하다. 이유는 plenum/filter case의 비교 기준이 되기 때문이다.

### Case 2: plenum + filter + full floor exhaust

추천 주력 case:

- `paint_booth_plenum_filter_full_floor_v035`

Dimension:

- Total domain: 8.0 × 4.0 × 3.8 m
- Work zone: z = 0 ~ 3.0 m
- Plenum: z = 3.0 ~ 3.8 m
- Filter panel: 7.0 × 3.2 m at z = 3.0 m
- Supply duct/opening: top 1.5 × 1.2 m 또는 side slots
- Floor exhaust: full floor outlet

Boundary condition:

- `supplyInlet`
  - type: fixedValue
  - U magnitude set by desired flow rate Q = 7.84 m³/s
  - If area = 1.8 m², U ≈ 4.36 m/s into plenum

- `filterPorousBaffle` or `filterPorousZone`
  - pressure drop target: 50 Pa at 0.35 m/s
  - role: jet smoothing / flow equalization

- `floorExhaust`
  - p = 0
  - U = pressureInletOutletVelocity

- `solidCeilingFrame`, `plenumWalls`, `frontWall`, `rearWall`, `sideWalls`, `carBody`
  - noSlip wall

### Case 3: plenum + filter + floor center slot

- `paint_booth_plenum_filter_center_slot_v035`

Dimension:

- Same plenum/filter as Case 2
- Floor exhaust center slot: 6.0 × 1.2 m
- Floor solid region: wall

Purpose:

- 실제 floor grating/underfloor exhaust effect 반영
- 차체 하부/측면 recirculation 변화 평가

### Case 4: resistance sweep

- Filter Δp target: 30, 50, 80, 120 Pa
- Same geometry and Q

Purpose:

- filter resistance가 downdraft uniformity와 plenum pressure distribution에 미치는 영향 확인

## 7. OpenFOAM 구현 전략

### Geometry/blockMesh 변경

현재 single-block domain은 bottom/top/front/rear/side만 갖는다. Plenum/filter를 넣으려면 다음 중 하나가 필요하다.

#### 방법 A: multi-block blockMesh

z 방향을 다음처럼 분할한다.

- z = 0.0: floor
- z = 3.0: filter plane / work-zone top
- z = 3.8: plenum top

x/y 방향도 filter panel과 solid frame을 나누려면 split이 필요하다.

권장 split coordinates:

- x: -1.5, -1.0, 6.0, 6.5
- y: -2.0, -1.6, 1.6, 2.0
- z: 0.0, 3.0, 3.8

이렇게 하면 z=3.0 plane에서 중앙 7.0 × 3.2 m filter patch와 외곽 solid frame을 구분할 수 있다.

#### 방법 B: 기존 domain + internal baffle 생성

blockMesh는 전체 domain만 만들고, z=3.0에 faceSet을 만든 뒤 createBaffles/createPatch로 filter baffle을 만든다.

장점:

- geometry 확장에 유연

단점:

- 자동화 스크립트가 조금 복잡

#### 방법 C: thin porous volume layer

z = 2.95 ~ 3.05 m를 filter layer로 만들고 cellZone을 지정한다.

장점:

- fvOptions porous source 적용 가능

단점:

- layer 주변 mesh 제어 필요
- filter panel 외곽 solid 영역 분리 필요

### 추천 구현 경로

현 프로젝트 generator 구조에서는 방법 A가 가장 재현성이 좋다.

즉, `create_paint_booth_baseline.py`에 다음 옵션을 추가한다.

```bash
--flow-mode streamwise|downdraft|plenum-filter
--booth-height 3.8
--workzone-height 3.0
--plenum-height 0.8
--filter-length 7.0
--filter-width 3.2
--filter-dp 50
--supply-mode top|side-slots
--floor-mode full|center-slot|side-slots
```

## 8. Boundary condition 상세 제안

### Plenum-filter case: velocity and pressure

`supplyInlet`:

- If top opening normal is -z:
  - U = `(0 0 -4.36)` m/s for 1.5 × 1.2 m opening and Q = 7.84 m³/s
- If side slot normal is ±y:
  - U magnitude adjusted by total side-slot area

`floorExhaust`:

- U: `pressureInletOutletVelocity`
- p: `fixedValue uniform 0`

Walls:

- U: `noSlip`
- p: `zeroGradient`

Filter:

- porous pressure loss or porous zone source
- target Δp = 50 Pa at U_filter = 0.35 m/s

### Turbulence condition

`supplyInlet` turbulence can be set based on duct inlet velocity.

First estimate:

- Turbulence intensity: 5 ~ 10%
- Length scale: 0.05 ~ 0.2 m for duct inlet

For filter exit / work zone, turbulence will be affected by porous resistance. If porous model is not turbulence-aware, post-filter turbulence may need simple fixed/limited values.

## 9. 평가 metric 수정

기존 work-zone metric은 유지하되, plenum/filter metric을 추가한다.

### Plenum metric

- Plenum pressure uniformity at z = 3.2 ~ 3.6 m
- Plenum velocity magnitude distribution
- Jet impingement on filter top surface

### Filter metric

- Filter face-normal velocity uniformity
- Filter pressure drop distribution
- Filter flow imbalance

### Work-zone metric

- Near-car Uz uniformity
- Low-speed pocket
- Reverse-flow fraction
- Roof/hood/trunk wake
- Mass balance

핵심은 filter 바로 아래 z = 2.8 ~ 2.9 m plane에서 Uz가 얼마나 균일한지를 보는 것이다.

## 10. 결론

사용자 의견을 반영하면, 단순 ceiling velocity inlet은 실제성 측면에서 부족하다. 다음 전략이 더 타당하다.

1. 작업공간은 기존 8.0 × 4.0 × 3.0 m 유지.
2. 상부 plenum 0.8 m를 추가하여 전체 높이를 3.8 m로 확장.
3. z = 3.0 m에 7.0 × 3.2 m filter panel을 둔다.
4. duct/supply opening은 plenum 상단 또는 측면에 배치한다.
5. filter는 pressure-drop 저항체로 모델링한다.
6. 작업공간으로 들어오는 유동은 duct jet가 아니라 filter 통과 후 균일화된 downdraft가 되도록 한다.
7. 첫 plenum case는 `paint_booth_plenum_filter_full_floor_v035`로 잡고, 이후 floor slot과 filter resistance sweep으로 확장한다.
