# Benchmark Evaluation Aggregate Results (scale_hinted_longest_horizontal_span_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted`
- Task count in this aggregate: **20**
- Aggregate mean structural score: **0.7913**
- Aggregate mean CFD score: **0.4905**
- Mean room bbox relative error: `Lx=0.0154`, `Ly=0.1919`, `Lz=0.1251`
- Mean room volume relative error: **0.3062**
- Room-kind match rate: **1.0**
- Opening-wall match rate: **0.7**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **2** tasks.
- Laminar fallback count: **1**

## View-level results

- **perspective**: structural `0.7812`, CFD `0.5164`, bbox-rel `(Lx=0.0253, Ly=0.311, Lz=0.1318)`, opening-wall match `0.75`
- **birdseye**: structural `0.8229`, CFD `0.54`, bbox-rel `(Lx=0.0004, Ly=0.1792, Lz=0.1336)`, opening-wall match `0.75`
- **floorplan**: structural `0.9062`, CFD `0.4953`, bbox-rel `(Lx=0.0503, Ly=0.0769, Lz=0.0752)`, opening-wall match `1.0`
- **wireframe**: structural `0.8021`, CFD `0.4955`, bbox-rel `(Lx=0.0004, Ly=0.2987, Lz=0.1318)`, opening-wall match `0.75`
- **section**: structural `0.6438`, CFD `0.4052`, bbox-rel `(Lx=0.0004, Ly=0.0937, Lz=0.1529)`, opening-wall match `0.25`

## Category-level results

- **A1**: structural `0.75`, CFD `0.5352`, bbox-rel `(Lx=0.0604, Ly=0.1914, Lz=0.0961)`, opening-wall match `0.8`
- **A2**: structural `0.7167`, CFD `0.4382`, bbox-rel `(Lx=0.0005, Ly=0.2703, Lz=0.1819)`, opening-wall match `0.2`
- **A3**: structural `0.775`, CFD `0.5078`, bbox-rel `(Lx=0.0005, Ly=0.1369, Lz=0.08)`, opening-wall match `0.8`
- **A4**: structural `0.9233`, CFD `0.4806`, bbox-rel `(Lx=0.0001, Ly=0.169, Lz=0.1424)`, opening-wall match `1.0`

## Room-kind split

- **rectangular**: structural `0.7333`, CFD `0.4867`, bbox-rel `(Lx=0.0304, Ly=0.2308, Lz=0.139)`, opening-wall match `0.5`
- **composite**: structural `0.8492`, CFD `0.4942`, bbox-rel `(Lx=0.0003, Ly=0.153, Lz=0.1112)`, opening-wall match `0.9`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 1 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
- **rectangular_blockage_layout_failure**: 4 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
- **dense_composite_structure_physics_gap**: 3 tasks
  - Dense composite task looked structurally acceptable but still lost CFD fidelity.
- **mesh025_fallback**: 2 tasks
  - Task needed the smaller 0.25 mesh-size fallback to succeed.
- **ultra_robust_escalation**: 2 tasks
  - Task required ultra_robust solver escalation.
- **laminar_fallback_success**: 1 tasks
  - Task finished only after falling back to a laminar run.
- **nonblocking_repair_sidecar_warning**: 1 tasks
  - Original scene succeeded, but a repair-sidecar attempt still left warning traces worth preserving as robustness metadata.

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4382`, bbox-rel `(Lx=0.0005, Ly=0.2703, Lz=0.1819)`
- `bench_a4_03` (A4, composite): structural `0.9233`, CFD `0.4806`, bbox-rel `(Lx=0.0001, Ly=0.169, Lz=0.1424)`
- `bench_a3_03` (A3, composite): structural `0.775`, CFD `0.5078`, bbox-rel `(Lx=0.0005, Ly=0.1369, Lz=0.08)`
- `bench_a1_01` (A1, rectangular): structural `0.75`, CFD `0.5352`, bbox-rel `(Lx=0.0604, Ly=0.1914, Lz=0.0961)`

### Highest mean CFD

- `bench_a1_01` (A1, rectangular): structural `0.75`, CFD `0.5352`, bbox-rel `(Lx=0.0604, Ly=0.1914, Lz=0.0961)`
- `bench_a3_03` (A3, composite): structural `0.775`, CFD `0.5078`, bbox-rel `(Lx=0.0005, Ly=0.1369, Lz=0.08)`
- `bench_a4_03` (A4, composite): structural `0.9233`, CFD `0.4806`, bbox-rel `(Lx=0.0001, Ly=0.169, Lz=0.1424)`
- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4382`, bbox-rel `(Lx=0.0005, Ly=0.2703, Lz=0.1819)`

## Artifact

- JSON summary: `26-03-19_scale_hinted_smoke_aggregate.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted/*/*/evaluation_summary.json`
