# Benchmark Evaluation Aggregate Results (posthoc_uniform_longest_span_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span`
- Task count in this aggregate: **100**
- Aggregate mean structural score: **0.7813**
- Aggregate mean CFD score: **0.4537**
- Mean room bbox relative error: `Lx=0.0039`, `Ly=0.2194`, `Lz=0.1983`
- Mean room volume relative error: **0.3663**
- Room-kind match rate: **0.95**
- Opening-wall match rate: **0.7**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **2** tasks.
- Laminar fallback count: **5**

## View-level results

- **perspective**: structural `0.7824`, CFD `0.4252`, bbox-rel `(Lx=0.0043, Ly=0.3338, Lz=0.2193)`, opening-wall match `0.65`
- **birdseye**: structural `0.8`, CFD `0.4812`, bbox-rel `(Lx=0.0038, Ly=0.2123, Lz=0.2199)`, opening-wall match `0.85`
- **floorplan**: structural `0.884`, CFD `0.5405`, bbox-rel `(Lx=0.0038, Ly=0.1003, Lz=0.1938)`, opening-wall match `1.0`
- **wireframe**: structural `0.7603`, CFD `0.4613`, bbox-rel `(Lx=0.004, Ly=0.2851, Lz=0.2299)`, opening-wall match `0.75`
- **section**: structural `0.6788`, CFD `0.3578`, bbox-rel `(Lx=0.0038, Ly=0.1802, Lz=0.132)`, opening-wall match `0.25`

## Category-level results

- **A1**: structural `0.845`, CFD `0.4599`, bbox-rel `(Lx=0.0, Ly=0.172, Lz=0.1816)`, opening-wall match `0.88`
- **A2**: structural `0.707`, CFD `0.4868`, bbox-rel `(Lx=0.0153, Ly=0.1846, Lz=0.2049)`, opening-wall match `0.48`
- **A3**: structural `0.783`, CFD `0.4558`, bbox-rel `(Lx=0.0, Ly=0.2474, Lz=0.2131)`, opening-wall match `0.72`
- **A4**: structural `0.791`, CFD `0.4087`, bbox-rel `(Lx=0.0, Ly=0.2795, Lz=0.1936)`, opening-wall match `0.72`

## Room-kind split

- **rectangular**: structural `0.776`, CFD `0.4733`, bbox-rel `(Lx=0.0077, Ly=0.1783, Lz=0.1932)`, opening-wall match `0.68`
- **composite**: structural `0.7869`, CFD `0.4328`, bbox-rel `(Lx=0.0, Ly=0.2631, Lz=0.2036)`, opening-wall match `0.766`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 2 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
- **rectangular_blockage_layout_failure**: 14 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
- **dense_composite_structure_physics_gap**: 13 tasks
  - Dense composite task looked structurally acceptable but still lost CFD fidelity.
- **section_room_kind_collapse**: 2 tasks
  - Section-view task collapsed a composite room into a rectangular prediction.
- **mesh025_fallback**: 2 tasks
  - Task needed the smaller 0.25 mesh-size fallback to succeed.
- **ultra_robust_escalation**: 4 tasks
  - Task required ultra_robust solver escalation.
- **laminar_fallback_success**: 5 tasks
  - Task finished only after falling back to a laminar run.

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a4_03` (A4, composite): structural `0.881`, CFD `0.3114`, bbox-rel `(Lx=0.0, Ly=0.2256, Lz=0.3176)`
- `bench_a4_04` (A4, composite): structural `0.85`, CFD `0.3567`, bbox-rel `(Lx=0.0, Ly=0.3038, Lz=0.1694)`
- `bench_a4_02` (A4, composite): structural `0.881`, CFD `0.3592`, bbox-rel `(Lx=0.0, Ly=0.216, Lz=0.1801)`
- `bench_a3_02` (A3, composite): structural `0.775`, CFD `0.367`, bbox-rel `(Lx=0.0, Ly=0.2023, Lz=0.1892)`
- `bench_a1_05` (A1, rectangular): structural `0.85`, CFD `0.3904`, bbox-rel `(Lx=0.0, Ly=0.1327, Lz=0.135)`

### Highest mean CFD

- `bench_a1_02` (A1, rectangular): structural `0.825`, CFD `0.6309`, bbox-rel `(Lx=0.0, Ly=0.2457, Lz=0.2293)`
- `bench_a4_05` (A4, composite): structural `0.5958`, CFD `0.5916`, bbox-rel `(Lx=0.0, Ly=0.3725, Lz=0.1098)`
- `bench_a3_04` (A3, composite): structural `0.6354`, CFD `0.5765`, bbox-rel `(Lx=0.0, Ly=0.344, Lz=0.3636)`
- `bench_a3_03` (A3, composite): structural `0.9`, CFD `0.5661`, bbox-rel `(Lx=0.0, Ly=0.2533, Lz=0.1828)`
- `bench_a2_04` (A2, rectangular): structural `0.55`, CFD `0.5553`, bbox-rel `(Lx=0.0, Ly=0.1442, Lz=0.1686)`

## Artifact

- JSON summary: `26-03-22_phase2_aggregate_results.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span/*/*/evaluation_summary.json`
