# Benchmark Evaluation Aggregate Results (scale_hinted_longest_span_plus_height_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_span_height`
- Task count in this aggregate: **20**
- Aggregate mean structural score: **0.8063**
- Aggregate mean CFD score: **0.4801**
- Mean room bbox relative error: `Lx=0.0004`, `Ly=0.1867`, `Lz=0.001`
- Mean room volume relative error: **0.2208**
- Room-kind match rate: **1.0**
- Opening-wall match rate: **0.65**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **1** tasks.
- Laminar fallback count: **0**

## View-level results

- **perspective**: structural `0.8125`, CFD `0.5167`, bbox-rel `(Lx=0.0004, Ly=0.2187, Lz=0.001)`, opening-wall match `0.75`
- **birdseye**: structural `0.8542`, CFD `0.4959`, bbox-rel `(Lx=0.0004, Ly=0.2128, Lz=0.001)`, opening-wall match `0.75`
- **floorplan**: structural `0.9688`, CFD `0.5748`, bbox-rel `(Lx=0.0004, Ly=0.0284, Lz=0.001)`, opening-wall match `1.0`
- **wireframe**: structural `0.7604`, CFD `0.4464`, bbox-rel `(Lx=0.0004, Ly=0.3202, Lz=0.001)`, opening-wall match `0.5`
- **section**: structural `0.6354`, CFD `0.3667`, bbox-rel `(Lx=0.0004, Ly=0.1532, Lz=0.001)`, opening-wall match `0.25`

## Category-level results

- **A1**: structural `0.775`, CFD `0.524`, bbox-rel `(Lx=0.0006, Ly=0.1584, Lz=0.0004)`, opening-wall match `0.8`
- **A2**: structural `0.7167`, CFD `0.4138`, bbox-rel `(Lx=0.0005, Ly=0.264, Lz=0.0013)`, opening-wall match `0.2`
- **A3**: structural `0.8`, CFD `0.5125`, bbox-rel `(Lx=0.0005, Ly=0.1884, Lz=0.0014)`, opening-wall match `0.6`
- **A4**: structural `0.9334`, CFD `0.4701`, bbox-rel `(Lx=0.0001, Ly=0.136, Lz=0.001)`, opening-wall match `1.0`

## Room-kind split

- **rectangular**: structural `0.7458`, CFD `0.4689`, bbox-rel `(Lx=0.0005, Ly=0.2112, Lz=0.0009)`, opening-wall match `0.5`
- **composite**: structural `0.8667`, CFD `0.4913`, bbox-rel `(Lx=0.0003, Ly=0.1622, Lz=0.0012)`, opening-wall match `0.8`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 2 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
- **rectangular_blockage_layout_failure**: 4 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
- **dense_composite_structure_physics_gap**: 4 tasks
  - Dense composite task looked structurally acceptable but still lost CFD fidelity.
- **mesh025_fallback**: 1 tasks
  - Task needed the smaller 0.25 mesh-size fallback to succeed.
- **ultra_robust_escalation**: 2 tasks
  - Task required ultra_robust solver escalation.
- **nonblocking_repair_sidecar_warning**: 1 tasks
  - Original scene succeeded, but a repair-sidecar attempt still left warning traces worth preserving as robustness metadata.

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4138`, bbox-rel `(Lx=0.0005, Ly=0.264, Lz=0.0013)`
- `bench_a4_03` (A4, composite): structural `0.9334`, CFD `0.4701`, bbox-rel `(Lx=0.0001, Ly=0.136, Lz=0.001)`
- `bench_a3_03` (A3, composite): structural `0.8`, CFD `0.5125`, bbox-rel `(Lx=0.0005, Ly=0.1884, Lz=0.0014)`
- `bench_a1_01` (A1, rectangular): structural `0.775`, CFD `0.524`, bbox-rel `(Lx=0.0006, Ly=0.1584, Lz=0.0004)`

### Highest mean CFD

- `bench_a1_01` (A1, rectangular): structural `0.775`, CFD `0.524`, bbox-rel `(Lx=0.0006, Ly=0.1584, Lz=0.0004)`
- `bench_a3_03` (A3, composite): structural `0.8`, CFD `0.5125`, bbox-rel `(Lx=0.0005, Ly=0.1884, Lz=0.0014)`
- `bench_a4_03` (A4, composite): structural `0.9334`, CFD `0.4701`, bbox-rel `(Lx=0.0001, Ly=0.136, Lz=0.001)`
- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4138`, bbox-rel `(Lx=0.0005, Ly=0.264, Lz=0.0013)`

## Artifact

- JSON summary: `26-03-19_scale_hinted_span_height_smoke_aggregate.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_span_height/*/*/evaluation_summary.json`
