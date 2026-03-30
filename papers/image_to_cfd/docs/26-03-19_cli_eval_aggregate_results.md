# Benchmark Evaluation Aggregate Results (scale_hinted_longest_span_guard_weighted_v1)

> Date: 2026-03-19

## Headline

- Evaluation root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_guard_weighted`
- Task count in this aggregate: **6**
- Aggregate mean structural score: **0.8542**
- Aggregate mean CFD score: **0.4865**
- Mean room bbox relative error: `Lx=0.0003`, `Ly=0.1814`, `Lz=0.0734`
- Mean room volume relative error: **0.247**
- Room-kind match rate: **0.8333**
- Opening-wall match rate: **0.8333**
- Repair was needed in **0** tasks; mesh-size `0.25` fallback appeared in **0** tasks.
- Laminar fallback count: **0**

## View-level results

- **perspective**: structural `0.9584`, CFD `0.4997`, bbox-rel `(Lx=0.0003, Ly=0.2036, Lz=0.0451)`, opening-wall match `1.0`
- **floorplan**: structural `1.0`, CFD `0.4321`, bbox-rel `(Lx=0.0001, Ly=0.1963, Lz=0.0332)`, opening-wall match `1.0`
- **wireframe**: structural `0.75`, CFD `0.6214`, bbox-rel `(Lx=0.0006, Ly=0.3502, Lz=0.0598)`, opening-wall match `1.0`
- **section**: structural `0.7292`, CFD `0.4332`, bbox-rel `(Lx=0.0003, Ly=0.0674, Lz=0.1285)`, opening-wall match `0.5`

## Category-level results

- **A1**: structural `0.6875`, CFD `0.5085`, bbox-rel `(Lx=0.0006, Ly=0.1887, Lz=0.158)`, opening-wall match `0.5`
- **A3**: structural `1.0`, CFD `0.546`, bbox-rel `(Lx=0.0005, Ly=0.2282, Lz=0.0768)`, opening-wall match `1.0`
- **A4**: structural `0.9167`, CFD `0.4521`, bbox-rel `(Lx=0.0001, Ly=0.161, Lz=0.0159)`, opening-wall match `1.0`

## Room-kind split

- **rectangular**: structural `0.6875`, CFD `0.5085`, bbox-rel `(Lx=0.0006, Ly=0.1887, Lz=0.158)`, opening-wall match `0.5`
- **composite**: structural `0.9375`, CFD `0.4756`, bbox-rel `(Lx=0.0002, Ly=0.1778, Lz=0.0311)`, opening-wall match `1.0`

## Derived interpretation tags

- **dense_composite_structure_physics_gap**: 3 tasks
  - Dense composite task looked structurally acceptable but still lost CFD fidelity.
- **section_room_kind_collapse**: 1 tasks
  - Section-view task collapsed a composite room into a rectangular prediction.

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a4_03` (A4, composite): structural `0.9167`, CFD `0.4521`, bbox-rel `(Lx=0.0001, Ly=0.161, Lz=0.0159)`
- `bench_a1_01` (A1, rectangular): structural `0.6875`, CFD `0.5085`, bbox-rel `(Lx=0.0006, Ly=0.1887, Lz=0.158)`
- `bench_a3_03` (A3, composite): structural `1.0`, CFD `0.546`, bbox-rel `(Lx=0.0005, Ly=0.2282, Lz=0.0768)`

### Highest mean CFD

- `bench_a3_03` (A3, composite): structural `1.0`, CFD `0.546`, bbox-rel `(Lx=0.0005, Ly=0.2282, Lz=0.0768)`
- `bench_a1_01` (A1, rectangular): structural `0.6875`, CFD `0.5085`, bbox-rel `(Lx=0.0006, Ly=0.1887, Lz=0.158)`
- `bench_a4_03` (A4, composite): structural `0.9167`, CFD `0.4521`, bbox-rel `(Lx=0.0001, Ly=0.161, Lz=0.0159)`

## Artifact

- JSON summary: `26-03-19_cli_eval_aggregate_results.json` (or CLI-provided target)
- Generated from: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_scale_hinted_guard_weighted/*/*/evaluation_summary.json`
