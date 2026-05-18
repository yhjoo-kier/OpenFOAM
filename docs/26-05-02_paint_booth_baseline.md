# Simplified Car Shell 기반 Paint Booth Baseline 생성 및 검증

## 목적

도장부스 CFD 연구를 위해 과도하게 상세한 스캔 BIW 형상 대신, CFD 친화적인 절차적 simplified car shell을 생성하고 OpenFOAM baseline case skeleton이 실제로 mesh/solver까지 통과하는지 확인했다.

## 현재 채택 형상

현재 baseline은 **원안 smooth single-shell car body**로 유지한다.

한 차례 lower body + roof rail + A/B/C pillars + side sill 방식의 side-window opening prototype을 생성해 검토했지만, box형 필러가 현재 단계의 surrogate로는 시각적으로 부자연스러웠다. 따라서 해당 prototype은 generator 내부에서 실험 옵션으로 남기되, 기본 baseline은 원안 smooth shell로 되돌렸다.

## 생성 파일

- Generator: `scripts/create_paint_booth_baseline.py`
- Mesh inspector: `scripts/inspect_mesh_quality.py`
- Case directory: `cases/paint_booth_baseline/`
- Generated STL: `cases/paint_booth_baseline/constant/triSurface/simplified_car_shell.stl`
- Mesh quality JSON: `docs/26-05-02_paint_booth_simplified_car_mesh_quality.json`
- Current render image: `docs/figures/26-05-02_paint_booth_simplified_car_shell_reverted_original.png`

## 형상 파라미터

절차적 sedan-like envelope를 superellipse cross-section으로 생성했다.

- Length: 4.5 m
- Width: 1.8 m
- Ground clearance: 0.18 m
- Approx. height range: z = 0.18 ~ 1.61095 m
- Vertices: 4,482
- Triangles: 8,960
- Window openings: false
- Opening strategy: `single_closed_smooth_shell`

기존 Sketchfab scan 후보와 비교:

- Sketchfab 후보: 약 3.42M triangles
- 이번 simplified shell: 8,960 triangles
- 약 380배 이상 가벼운 surface mesh

## Surface 품질 검사

Docker image: `openfoam-pipeline-local:latest`

실행:

```bash
surfaceCheck constant/triSurface/simplified_car_shell.stl
```

주요 결과:

- Bounding Box: `(0 -0.899943 0.18)` to `(4.5 0.899943 1.61095)`
- Surface closed: yes
- All edges connected to two faces: yes
- Number of unconnected parts: 1
- Number of zones with consistent normal: 1

판정: OpenFOAM surface 기준으로 closed / single connected / consistent normal 조건을 만족한다.

## Python/PyVista 품질 검사

실행:

```bash
python scripts/inspect_mesh_quality.py \
  cases/paint_booth_baseline/constant/triSurface/simplified_car_shell.stl \
  --out docs/26-05-02_paint_booth_simplified_car_mesh_quality.json
```

주요 결과:

- Points: 4,482
- Cells/Triangles: 8,960
- Open/non-manifold edges: 0
- Connected regions: 1
- Watertight likely: true
- Single component likely: true
- Very large for snappyHexMesh: false

## OpenFOAM mesh 검증

실행 순서:

```bash
blockMesh
snappyHexMesh -overwrite
checkMesh
```

`checkMesh` 주요 결과:

- Points: 63,920
- Faces: 160,900
- Internal faces: 148,509
- Cells: 49,206
- Domain bounding box: `(-1.5 -2 0)` to `(6.5 2 3)`
- Max aspect ratio: 3.97981
- Max non-orthogonality: 53.2065
- Average non-orthogonality: 10.9944
- Max skewness: 2.71408
- Result: `Mesh OK.`

## Solver smoke test

실행:

```bash
simpleFoam
```

결과:

- Solver completed to time step 200.
- Exit code: 0
- Final cumulative continuity error: ~1.70548e-04
- Final clock time: 약 24 s
- Result: `End`

판정: 현재 원안 smooth single-shell baseline은 OpenFOAM simpleFoam smoke test까지 통과한다.

## Window-opening prototype에 대한 판단

시도한 방식:

- Smooth lower body
- Roof rail
- A/B/C pillars
- Side sill
- Front/rear header

검증 결과 자체는 성공했다.

- Surface closed: yes
- Open/non-manifold edges: 0
- checkMesh: `Mesh OK`
- simpleFoam: exit code 0

하지만 시각적으로 box형 필러가 너무 인위적이어서 현재 baseline으로 채택하지 않았다. 이후 필러/창문을 다시 넣는다면 axis-aligned box가 아니라 tapered/curved pillar 또는 CAD-like parametric surface 방식이 필요하다.

## 현재 case의 한계

현재 case는 아직 "도장부스 최종 물리 모델"이 아니라 skeleton이다.

- 현재 유동은 streamwise inlet/outlet cold-flow 형태다.
- Ceiling downdraft / floor exhaust 물리조건은 아직 실제 도장부스 형태로 분리하지 않았다.
- ceiling/floorExhaust는 wall patch로 두어 wall-function smoke test가 통과하도록 했다.
- 실제 paint booth 모델에서는 ceiling filter inlet과 floor grating outlet을 patch로 재정의해야 한다.
- Wheel, underbody, booth plenum, conveyor, porous floor/grating은 아직 생략했다.

## 다음 단계 제안

1. `ceiling`을 downdraft velocity inlet으로 변경.
2. `floorExhaust`를 pressure outlet 또는 porous outlet으로 변경.
3. side inlet/outlet 방식이 아니라 도장부스 전형 조건인 top-down flow로 변경.
4. floor grating 영역을 전체 floor가 아닌 slot/perforated/porous patch로 분리.
5. carBody 주변 속도 균일도와 recirculation zone post-processing 추가.
6. 이후 overspray particle tracking 또는 scalar transport를 추가.
7. 형상 디테일은 현재처럼 안정적인 baseline 이후, tapered/curved pillar 방식으로 재검토한다.
