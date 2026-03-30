# Benchmark Evaluation Aggregate Results (posthoc_uniform_longest_span_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span`
- Task count in this aggregate: **20**
- Aggregate mean structural score: **0.8307**
- Aggregate mean CFD score: **0.5329**
- Mean room bbox relative error: `Lx=0.0`, `Ly=0.2209`, `Lz=0.2151`
- Mean room volume relative error: **0.3861**
- Room-kind match rate: **1.0**
- Opening-wall match rate: **0.75**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **0** tasks.
- Laminar fallback count: **1**

## View-level results

- **perspective**: structural `0.8229`, CFD `0.5522`, bbox-rel `(Lx=0.0, Ly=0.3458, Lz=0.2037)`, opening-wall match `0.75`
- **birdseye**: structural `0.7708`, CFD `0.5637`, bbox-rel `(Lx=0.0, Ly=0.2415, Lz=0.2727)`, opening-wall match `0.75`
- **floorplan**: structural `0.9479`, CFD `0.598`, bbox-rel `(Lx=0.0, Ly=0.0862, Lz=0.2726)`, opening-wall match `1.0`
- **wireframe**: structural `0.7887`, CFD `0.5422`, bbox-rel `(Lx=0.0, Ly=0.3676, Lz=0.2351)`, opening-wall match `0.75`
- **section**: structural `0.8229`, CFD `0.4086`, bbox-rel `(Lx=0.0, Ly=0.0634, Lz=0.0915)`, opening-wall match `0.5`

## Category-level results

- **A1**: structural `0.825`, CFD `0.5648`, bbox-rel `(Lx=0.0, Ly=0.1298, Lz=0.0699)`, opening-wall match `1.0`
- **A2**: structural `0.7167`, CFD `0.4141`, bbox-rel `(Lx=0.0, Ly=0.2748, Lz=0.2902)`, opening-wall match `0.2`
- **A3**: structural `0.9`, CFD `0.5944`, bbox-rel `(Lx=0.0, Ly=0.2533, Lz=0.1828)`, opening-wall match `0.8`
- **A4**: structural `0.881`, CFD `0.5583`, bbox-rel `(Lx=0.0, Ly=0.2256, Lz=0.3176)`, opening-wall match `1.0`

## Room-kind split

- **rectangular**: structural `0.7708`, CFD `0.4894`, bbox-rel `(Lx=0.0, Ly=0.2023, Lz=0.18)`, opening-wall match `0.6`
- **composite**: structural `0.8905`, CFD `0.5764`, bbox-rel `(Lx=0.0, Ly=0.2395, Lz=0.2502)`, opening-wall match `0.9`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 3 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
- **rectangular_blockage_layout_failure**: 4 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
- **laminar_fallback_success**: 1 tasks
  - Task finished only after falling back to a laminar run.

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4141`, bbox-rel `(Lx=0.0, Ly=0.2748, Lz=0.2902)`
- `bench_a4_03` (A4, composite): structural `0.881`, CFD `0.5583`, bbox-rel `(Lx=0.0, Ly=0.2256, Lz=0.3176)`
- `bench_a1_01` (A1, rectangular): structural `0.825`, CFD `0.5648`, bbox-rel `(Lx=0.0, Ly=0.1298, Lz=0.0699)`
- `bench_a3_03` (A3, composite): structural `0.9`, CFD `0.5944`, bbox-rel `(Lx=0.0, Ly=0.2533, Lz=0.1828)`

### Highest mean CFD

- `bench_a3_03` (A3, composite): structural `0.9`, CFD `0.5944`, bbox-rel `(Lx=0.0, Ly=0.2533, Lz=0.1828)`
- `bench_a1_01` (A1, rectangular): structural `0.825`, CFD `0.5648`, bbox-rel `(Lx=0.0, Ly=0.1298, Lz=0.0699)`
- `bench_a4_03` (A4, composite): structural `0.881`, CFD `0.5583`, bbox-rel `(Lx=0.0, Ly=0.2256, Lz=0.3176)`
- `bench_a2_03` (A2, rectangular): structural `0.7167`, CFD `0.4141`, bbox-rel `(Lx=0.0, Ly=0.2748, Lz=0.2902)`

## Artifact

- JSON summary: `26-03-20_posthoc_scaled_smoke_aggregate.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span/*/*/evaluation_summary.json`
