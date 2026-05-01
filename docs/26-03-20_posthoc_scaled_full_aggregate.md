# Benchmark Evaluation Aggregate Results (posthoc_uniform_longest_span_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span`
- Task count in this aggregate: **100**
- Aggregate mean structural score: **0.7813**
- Aggregate mean CFD score: **0.4909**
- Mean room bbox relative error: `Lx=0.0039`, `Ly=0.2194`, `Lz=0.1983`
- Mean room volume relative error: **0.3663**
- Room-kind match rate: **0.95**
- Opening-wall match rate: **0.7**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **2** tasks.
- Laminar fallback count: **5**

## View-level results

- **perspective**: structural `0.7824`, CFD `0.4677`, bbox-rel `(Lx=0.0043, Ly=0.3338, Lz=0.2193)`, opening-wall match `0.65`
- **birdseye**: structural `0.8`, CFD `0.509`, bbox-rel `(Lx=0.0038, Ly=0.2123, Lz=0.2199)`, opening-wall match `0.85`
- **floorplan**: structural `0.884`, CFD `0.5379`, bbox-rel `(Lx=0.0038, Ly=0.1003, Lz=0.1938)`, opening-wall match `1.0`
- **wireframe**: structural `0.7603`, CFD `0.5082`, bbox-rel `(Lx=0.004, Ly=0.2851, Lz=0.2299)`, opening-wall match `0.75`
- **section**: structural `0.6788`, CFD `0.4304`, bbox-rel `(Lx=0.0038, Ly=0.1802, Lz=0.132)`, opening-wall match `0.25`

## Category-level results

- **A1**: structural `0.845`, CFD `0.5452`, bbox-rel `(Lx=0.0, Ly=0.172, Lz=0.1816)`, opening-wall match `0.88`
- **A2**: structural `0.707`, CFD `0.426`, bbox-rel `(Lx=0.0153, Ly=0.1846, Lz=0.2049)`, opening-wall match `0.48`
- **A3**: structural `0.783`, CFD `0.5336`, bbox-rel `(Lx=0.0, Ly=0.2474, Lz=0.2131)`, opening-wall match `0.72`
- **A4**: structural `0.791`, CFD `0.458`, bbox-rel `(Lx=0.0, Ly=0.2795, Lz=0.1936)`, opening-wall match `0.72`

## Room-kind split

- **rectangular**: structural `0.776`, CFD `0.4856`, bbox-rel `(Lx=0.0077, Ly=0.1783, Lz=0.1932)`, opening-wall match `0.68`
- **composite**: structural `0.7869`, CFD `0.4966`, bbox-rel `(Lx=0.0, Ly=0.2631, Lz=0.2036)`, opening-wall match `0.766`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 9 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
- **rectangular_blockage_layout_failure**: 17 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
- **dense_composite_structure_physics_gap**: 8 tasks
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

- `bench_a2_02` (A2, rectangular): structural `0.675`, CFD `0.3601`, bbox-rel `(Lx=0.0246, Ly=0.1483, Lz=0.1931)`
- `bench_a4_02` (A4, composite): structural `0.881`, CFD `0.3603`, bbox-rel `(Lx=0.0, Ly=0.216, Lz=0.1801)`
- `bench_a4_05` (A4, composite): structural `0.5958`, CFD `0.3962`, bbox-rel `(Lx=0.0, Ly=0.3725, Lz=0.1098)`
- `bench_a2_01` (A2, rectangular): structural `0.785`, CFD `0.4043`, bbox-rel `(Lx=0.0, Ly=0.1889, Lz=0.1881)`
- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4129`, bbox-rel `(Lx=0.0, Ly=0.2748, Lz=0.2902)`

### Highest mean CFD

- `bench_a1_03` (A1, rectangular): structural `0.85`, CFD `0.6947`, bbox-rel `(Lx=0.0, Ly=0.2247, Lz=0.2031)`
- `bench_a3_03` (A3, composite): structural `0.9`, CFD `0.5944`, bbox-rel `(Lx=0.0, Ly=0.2533, Lz=0.1828)`
- `bench_a1_01` (A1, rectangular): structural `0.825`, CFD `0.5648`, bbox-rel `(Lx=0.0, Ly=0.1298, Lz=0.0699)`
- `bench_a1_02` (A1, rectangular): structural `0.825`, CFD `0.5628`, bbox-rel `(Lx=0.0, Ly=0.2457, Lz=0.2293)`
- `bench_a4_03` (A4, composite): structural `0.881`, CFD `0.5583`, bbox-rel `(Lx=0.0, Ly=0.2256, Lz=0.3176)`

## Artifact

- JSON summary: `26-03-20_posthoc_scaled_full_aggregate.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span/*/*/evaluation_summary.json`
