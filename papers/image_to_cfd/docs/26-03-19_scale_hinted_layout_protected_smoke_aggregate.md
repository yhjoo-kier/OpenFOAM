# Benchmark Evaluation Aggregate Results (scale_hinted_longest_span_layout_protected_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_layout_protected`
- Task count in this aggregate: **20**
- Aggregate mean structural score: **0.7854**
- Aggregate mean CFD score: **0.5017**
- Mean room bbox relative error: `Lx=0.0004`, `Ly=0.1945`, `Lz=0.1305`
- Mean room volume relative error: **0.2908**
- Room-kind match rate: **1.0**
- Opening-wall match rate: **0.7**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **1** tasks.
- Laminar fallback count: **2**

## View-level results

- **perspective**: structural `0.7604`, CFD `0.4601`, bbox-rel `(Lx=0.0004, Ly=0.3379, Lz=0.1253)`, opening-wall match `0.75`
- **birdseye**: structural `0.8333`, CFD `0.4836`, bbox-rel `(Lx=0.0004, Ly=0.2061, Lz=0.0832)`, opening-wall match `0.75`
- **floorplan**: structural `0.9688`, CFD `0.7111`, bbox-rel `(Lx=0.0004, Ly=0.0711, Lz=0.0491)`, opening-wall match `1.0`
- **wireframe**: structural `0.7188`, CFD `0.4983`, bbox-rel `(Lx=0.0004, Ly=0.311, Lz=0.1318)`, opening-wall match `0.75`
- **section**: structural `0.6458`, CFD `0.3555`, bbox-rel `(Lx=0.0004, Ly=0.0465, Lz=0.263)`, opening-wall match `0.25`

## Category-level results

- **A1**: structural `0.725`, CFD `0.4652`, bbox-rel `(Lx=0.0006, Ly=0.1673, Lz=0.1022)`, opening-wall match `0.8`
- **A2**: structural `0.7167`, CFD `0.4546`, bbox-rel `(Lx=0.0005, Ly=0.2407, Lz=0.2416)`, opening-wall match `0.2`
- **A3**: structural `0.8`, CFD `0.6145`, bbox-rel `(Lx=0.0005, Ly=0.1568, Lz=0.0676)`, opening-wall match `0.8`
- **A4**: structural `0.9`, CFD `0.4725`, bbox-rel `(Lx=0.0001, Ly=0.2133, Lz=0.1105)`, opening-wall match `1.0`

## Room-kind split

- **rectangular**: structural `0.7208`, CFD `0.4599`, bbox-rel `(Lx=0.0005, Ly=0.204, Lz=0.1719)`, opening-wall match `0.5`
- **composite**: structural `0.85`, CFD `0.5435`, bbox-rel `(Lx=0.0003, Ly=0.185, Lz=0.0891)`, opening-wall match `0.9`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 2 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
- **rectangular_blockage_layout_failure**: 4 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
- **dense_composite_structure_physics_gap**: 3 tasks
  - Dense composite task looked structurally acceptable but still lost CFD fidelity.
- **mesh025_fallback**: 1 tasks
  - Task needed the smaller 0.25 mesh-size fallback to succeed.
- **ultra_robust_escalation**: 1 tasks
  - Task required ultra_robust solver escalation.
- **laminar_fallback_success**: 2 tasks
  - Task finished only after falling back to a laminar run.

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4546`, bbox-rel `(Lx=0.0005, Ly=0.2407, Lz=0.2416)`
- `bench_a1_01` (A1, rectangular): structural `0.725`, CFD `0.4652`, bbox-rel `(Lx=0.0006, Ly=0.1673, Lz=0.1022)`
- `bench_a4_03` (A4, composite): structural `0.9`, CFD `0.4725`, bbox-rel `(Lx=0.0001, Ly=0.2133, Lz=0.1105)`
- `bench_a3_03` (A3, composite): structural `0.8`, CFD `0.6145`, bbox-rel `(Lx=0.0005, Ly=0.1568, Lz=0.0676)`

### Highest mean CFD

- `bench_a3_03` (A3, composite): structural `0.8`, CFD `0.6145`, bbox-rel `(Lx=0.0005, Ly=0.1568, Lz=0.0676)`
- `bench_a4_03` (A4, composite): structural `0.9`, CFD `0.4725`, bbox-rel `(Lx=0.0001, Ly=0.2133, Lz=0.1105)`
- `bench_a1_01` (A1, rectangular): structural `0.725`, CFD `0.4652`, bbox-rel `(Lx=0.0006, Ly=0.1673, Lz=0.1022)`
- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4546`, bbox-rel `(Lx=0.0005, Ly=0.2407, Lz=0.2416)`

## Artifact

- JSON summary: `26-03-19_scale_hinted_layout_protected_smoke_aggregate.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_layout_protected/*/*/evaluation_summary.json`
