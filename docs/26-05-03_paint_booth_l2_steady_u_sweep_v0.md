# Paint-booth L2 steady U-sweep v0 results

Date: 2026-05-03

## Objective

Run a moderate one-parameter `U_supply` sweep using the validated L2 paint-booth CFD baseline, then export the completed cases into ML-ready cell-centered samples.

This follows the prior 3-point QA gate and is intended as the first steady CFD dataset for neural-operator/GINO pipeline development.

## Execution

Runner command:

```bash
python3 scripts/run_paint_booth_l2_flow_qa.py \
  --root cases/paint_booth_l2_u_sweep_v0 \
  --supply 2.18 2.72 3.27 3.81 4.36 4.91 5.45 6.00 6.54
```

Execution log:

```text
logs/26-05-03_paint_booth_l2_u_sweep_v0.log
```

Summary outputs:

```text
cases/paint_booth_l2_u_sweep_v0/flow_qa_summary.json
cases/paint_booth_l2_u_sweep_v0/flow_qa_summary.csv
```

Dataset export command:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 scripts/export_paint_booth_ml_dataset.py \
    --case-root cases/paint_booth_l2_u_sweep_v0 \
    --output-root data/paint_booth_ml_dataset/l2_steady_u_sweep_v0 \
    --include "l2_flow_qa_*" \
    --force'
```

Dataset export path:

```text
data/paint_booth_ml_dataset/l2_steady_u_sweep_v0/
```

Generated dataset artifacts are ignored by Git via `/data/paint_booth_ml_dataset/`.

## Baseline configuration

- Runtime image: `openfoam-pipeline-local:latest`
- Solver: `simpleFoam`
- Turbulence model: `kOmegaSST`
- Mesh/model: L2 paint-booth central filter panel + sealed frame baseline
- Filter Forchheimer coefficient: `10000`
- Base cell size: `0.125 m`
- Filter z cells: `6`
- carBody refinement: `3/4`
- carBody prism layers: `5`
- Absolute final layer thickness: `0.004 m`

## QA result table

All 9 cases passed the basic QA gate:

- `mesh_ok == true`
- `simpleFoam_ended == true`
- `abs(relative_imbalance) < 1%`
- carBody `y+ p95 < 5`
- carBody layer coverage `> 95%`

Per-case summary:

- `U = 2.18 m/s`
  - `Q_in = 3.9240 m3/s`
  - mass imbalance: `-0.0660%`
  - filter dp proxy: `16.71 Pa`
  - work-zone `Uz_mean = -0.0932 m/s`
  - near-car `Uz_mean = -0.0656 m/s`
  - near-car reverse fraction: `24.59%`
  - carBody `y+ p95 = 0.504`, `y+ max = 9.961`
  - QA: pass

- `U = 2.72 m/s`
  - `Q_in = 4.8960 m3/s`
  - mass imbalance: `+0.0134%`
  - filter dp proxy: `26.07 Pa`
  - work-zone `Uz_mean = -0.1621 m/s`
  - near-car `Uz_mean = -0.1270 m/s`
  - near-car reverse fraction: `21.84%`
  - carBody `y+ p95 = 1.532`, `y+ max = 25.667`
  - QA: pass

- `U = 3.27 m/s`
  - `Q_in = 5.8860 m3/s`
  - mass imbalance: `-0.0653%`
  - filter dp proxy: `37.60 Pa`
  - work-zone `Uz_mean = -0.1460 m/s`
  - near-car `Uz_mean = -0.1039 m/s`
  - near-car reverse fraction: `23.76%`
  - carBody `y+ p95 = 1.040`, `y+ max = 16.369`
  - QA: pass

- `U = 3.81 m/s`
  - `Q_in = 6.8580 m3/s`
  - mass imbalance: `-0.0638%`
  - filter dp proxy: `51.04 Pa`
  - work-zone `Uz_mean = -0.1727 m/s`
  - near-car `Uz_mean = -0.1235 m/s`
  - near-car reverse fraction: `23.48%`
  - carBody `y+ p95 = 1.369`, `y+ max = 19.674`
  - QA: pass

- `U = 4.36 m/s`
  - `Q_in = 7.8480 m3/s`
  - mass imbalance: `-0.0626%`
  - filter dp proxy: `66.84 Pa`
  - work-zone `Uz_mean = -0.2000 m/s`
  - near-car `Uz_mean = -0.1435 m/s`
  - near-car reverse fraction: `23.27%`
  - carBody `y+ p95 = 1.737`, `y+ max = 23.496`
  - QA: pass

- `U = 4.91 m/s`
  - `Q_in = 8.8380 m3/s`
  - mass imbalance: `-0.0624%`
  - filter dp proxy: `84.77 Pa`
  - work-zone `Uz_mean = -0.2275 m/s`
  - near-car `Uz_mean = -0.1637 m/s`
  - near-car reverse fraction: `23.06%`
  - carBody `y+ p95 = 2.128`, `y+ max = 27.309`
  - QA: pass

- `U = 5.45 m/s`
  - `Q_in = 9.8100 m3/s`
  - mass imbalance: `-0.0643%`
  - filter dp proxy: `104.44 Pa`
  - work-zone `Uz_mean = -0.2547 m/s`
  - near-car `Uz_mean = -0.1837 m/s`
  - near-car reverse fraction: `22.94%`
  - carBody `y+ p95 = 2.526`, `y+ max = 30.656`
  - QA: pass

- `U = 6.00 m/s`
  - `Q_in = 10.8000 m3/s`
  - mass imbalance: `-0.0640%`
  - filter dp proxy: `126.58 Pa`
  - work-zone `Uz_mean = -0.2826 m/s`
  - near-car `Uz_mean = -0.2042 m/s`
  - near-car reverse fraction: `22.84%`
  - carBody `y+ p95 = 2.942`, `y+ max = 33.532`
  - QA: pass

- `U = 6.54 m/s`
  - `Q_in = 11.7720 m3/s`
  - mass imbalance: `-0.0633%`
  - filter dp proxy: `150.39 Pa`
  - work-zone `Uz_mean = -0.3103 m/s`
  - near-car `Uz_mean = -0.2246 m/s`
  - near-car reverse fraction: `22.74%`
  - carBody `y+ p95 = 3.357`, `y+ max = 36.367`
  - QA: pass

## Trend checks

The runner reported:

- Supply inflow increases with supply velocity: `true`
- Filter pressure-drop proxy increases with supply velocity: `true`
- Absolute work-zone downdraft magnitude increases monotonically: `false`
- Absolute near-car downdraft magnitude increases monotonically: `false`
- Work-zone `Uz` becomes monotonically more negative: `false`
- Near-car `Uz` becomes monotonically more negative: `false`

Interpretation:

- The mass-flow and porous-pressure-drop trends are robust and monotonic.
- The local work-zone/near-car downdraft metrics show one notable kink around `U_supply = 2.72 m/s`.
- Specifically, the `2.72 m/s` case has stronger mean downdraft than the neighboring `3.27 m/s` case while still passing mesh, solver, mass-balance, and y+ QA.
- From `3.27 m/s` through `6.54 m/s`, work-zone and near-car mean downdraft trends are smooth and monotonic again.

This means the 9-point dataset is usable for ML pipeline testing and steady surrogate prototyping, but the `2.72 m/s` sample should be flagged for additional diagnostic review before treating the full sweep as a clean monotonic response dataset.

## Dataset export verification

Post-export verification:

```text
sample_count 9
l2_u2p18 U 2.18 cells 224769 qa True shape_ok True
l2_u2p72 U 2.72 cells 224769 qa True shape_ok True
l2_u3p27 U 3.27 cells 224769 qa True shape_ok True
l2_u3p81 U 3.81 cells 224769 qa True shape_ok True
l2_u4p36 U 4.36 cells 224769 qa True shape_ok True
l2_u4p91 U 4.91 cells 224769 qa True shape_ok True
l2_u5p45 U 5.45 cells 224769 qa True shape_ok True
l2_u6p00 U 6.0 cells 224769 qa True shape_ok True
l2_u6p54 U 6.54 cells 224769 qa True shape_ok True
```

Each exported `.npz` sample contains:

- `cell_centers`: `(224769, 3)`
- `U`: `(224769, 3)`
- `p`: `(224769,)`

## Recommendation

Next actions:

1. Use this exported dataset to validate the ML data loader and normalization pipeline.
2. Before training claims, inspect the `2.72 m/s` and `3.27 m/s` cases visually/locally to determine whether the kink is physical, numerical, or metric-region sensitivity.
3. Consider either:
   - keeping all 9 samples with a QA flag for `2.72 m/s`, or
   - using the monotonic subset `3.27..6.54 m/s` for first surrogate sanity training.
4. If monotonic low-flow behavior is important, rerun/refine low-flow samples around `2.4..3.4 m/s` and compare region metrics.
