# Benchmark Evaluation Aggregate Results (scale_hinted_longest_span_view_guarded_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_view_guarded`
- Task count in this aggregate: **20**
- Aggregate mean structural score: **0.7933**
- Aggregate mean CFD score: **0.4926**
- Mean room bbox relative error: `Lx=0.0004`, `Ly=0.2308`, `Lz=0.1306`
- Mean room volume relative error: **0.3154**
- Room-kind match rate: **1.0**
- Opening-wall match rate: **0.7**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **0** tasks.
- Laminar fallback count: **0**

## View-level results

- **perspective**: structural `0.8125`, CFD `0.524`, bbox-rel `(Lx=0.0004, Ly=0.3314, Lz=0.1071)`, opening-wall match `0.75`
- **birdseye**: structural `0.8333`, CFD `0.5164`, bbox-rel `(Lx=0.0004, Ly=0.231, Lz=0.1336)`, opening-wall match `0.75`
- **floorplan**: structural `0.9688`, CFD `0.5822`, bbox-rel `(Lx=0.0004, Ly=0.0453, Lz=0.0934)`, opening-wall match `1.0`
- **wireframe**: structural `0.7083`, CFD `0.4518`, bbox-rel `(Lx=0.0004, Ly=0.3959, Lz=0.2298)`, opening-wall match `0.75`
- **section**: structural `0.6438`, CFD `0.3888`, bbox-rel `(Lx=0.0004, Ly=0.1503, Lz=0.0889)`, opening-wall match `0.25`

## Category-level results

- **A1**: structural `0.75`, CFD `0.4731`, bbox-rel `(Lx=0.0006, Ly=0.2304, Lz=0.1022)`, opening-wall match `0.8`
- **A2**: structural `0.7333`, CFD `0.4207`, bbox-rel `(Lx=0.0005, Ly=0.2984, Lz=0.1506)`, opening-wall match `0.2`
- **A3**: structural `0.75`, CFD `0.5696`, bbox-rel `(Lx=0.0005, Ly=0.2141, Lz=0.1589)`, opening-wall match `0.8`
- **A4**: structural `0.94`, CFD `0.5071`, bbox-rel `(Lx=0.0001, Ly=0.1803, Lz=0.1105)`, opening-wall match `1.0`

## Room-kind split

- **rectangular**: structural `0.7417`, CFD `0.4469`, bbox-rel `(Lx=0.0005, Ly=0.2644, Lz=0.1264)`, opening-wall match `0.5`
- **composite**: structural `0.845`, CFD `0.5384`, bbox-rel `(Lx=0.0003, Ly=0.1972, Lz=0.1347)`, opening-wall match `0.9`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 3 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
- **rectangular_blockage_layout_failure**: 4 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
- **dense_composite_structure_physics_gap**: 2 tasks
  - Dense composite task looked structurally acceptable but still lost CFD fidelity.
- **ultra_robust_escalation**: 3 tasks
  - Task required ultra_robust solver escalation.

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a2_03` (A2, rectangular): structural `0.7333`, CFD `0.4207`, bbox-rel `(Lx=0.0005, Ly=0.2984, Lz=0.1506)`
- `bench_a1_01` (A1, rectangular): structural `0.75`, CFD `0.4731`, bbox-rel `(Lx=0.0006, Ly=0.2304, Lz=0.1022)`
- `bench_a4_03` (A4, composite): structural `0.94`, CFD `0.5071`, bbox-rel `(Lx=0.0001, Ly=0.1803, Lz=0.1105)`
- `bench_a3_03` (A3, composite): structural `0.75`, CFD `0.5696`, bbox-rel `(Lx=0.0005, Ly=0.2141, Lz=0.1589)`

### Highest mean CFD

- `bench_a3_03` (A3, composite): structural `0.75`, CFD `0.5696`, bbox-rel `(Lx=0.0005, Ly=0.2141, Lz=0.1589)`
- `bench_a4_03` (A4, composite): structural `0.94`, CFD `0.5071`, bbox-rel `(Lx=0.0001, Ly=0.1803, Lz=0.1105)`
- `bench_a1_01` (A1, rectangular): structural `0.75`, CFD `0.4731`, bbox-rel `(Lx=0.0006, Ly=0.2304, Lz=0.1022)`
- `bench_a2_03` (A2, rectangular): structural `0.7333`, CFD `0.4207`, bbox-rel `(Lx=0.0005, Ly=0.2984, Lz=0.1506)`

## Artifact

- JSON summary: `26-03-19_scale_hinted_view_guarded_smoke_aggregate.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_view_guarded/*/*/evaluation_summary.json`
