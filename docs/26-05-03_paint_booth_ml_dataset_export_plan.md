# Paint-booth L2 CFD-to-ML dataset export plan

Date: 2026-05-03

## Objective

Convert the validated L2 paint-booth CFD baseline into a first ML-ready steady dataset pipeline for neural-operator/GINO-style surrogate development.

This plan follows the 2026-05-02 three-point L2 flow-rate QA result: the current `kOmegaSST` + carBody boundary-layer mesh passed the low/nominal/high supply-velocity sanity gate, so the next step is a moderate one-parameter `U_supply` sweep with per-case QA metadata.

## Baseline to preserve

- Case generator: `scripts/create_paint_booth_plenum_filter_case.py`
- Existing QA runner: `scripts/run_paint_booth_l2_flow_qa.py`
- Runtime image: `openfoam-pipeline-local:latest`
- Baseline mesh/model settings:
  - Turbulence model: `kOmegaSST`
  - Central filter panel + sealed high-resistance frame: enabled
  - Filter Forchheimer coefficient: `10000`
  - Base cell size: `0.125 m`
  - Filter z cells: `6`
  - carBody refinement: `3/4`
  - carBody prism layers: `5`
  - Absolute final layer thickness: `0.004 m`
- QA gate reference: `docs/26-05-02_paint_booth_l2_flow_qa_for_ml_dataset.md`

## Dataset scope for first iteration

### Parameter axis

Use a one-parameter steady sweep over supply velocity only.

Recommended initial sweep:

- Range: `U_supply = 2.18 .. 6.54 m/s`
- Include the already validated QA points: `2.18`, `4.36`, `6.54 m/s`
- Add intermediate points for a moderate first dataset, e.g. `2.18, 2.72, 3.27, 3.81, 4.36, 4.91, 5.45, 6.00, 6.54 m/s`

Rationale: this keeps the geometry, porous resistance, turbulence model, and mesh fixed while generating enough monotonic response variation for the first surrogate experiment.

### Output representation

First exporter target: cell-centered graph/point-cloud samples from `foamToVTK` output.

Per case sample fields:

- Coordinates: cell-center `x, y, z`
- Inputs/conditioning:
  - `supply_velocity_mps`
  - `filter_forchheimer`
  - optional geometry/region tags if available later
- Targets:
  - `U` vector
  - `p`
  - optional `k`, `omega`, `nut` if present
- QA metadata:
  - mass imbalance
  - `checkMesh` status and cell count
  - carBody layer coverage
  - carBody `y+` median/p95/max
  - residual tail summary
  - reverse-flow fractions and downdraft metrics

Preferred file layout:

```text
data/paint_booth_ml_dataset/l2_steady_u_sweep_v0/
├── manifest.json
├── manifest.csv
├── samples/
│   ├── l2_u2p18.npz
│   ├── l2_u2p72.npz
│   └── ...
└── metadata/
    ├── l2_u2p18.json
    ├── l2_u2p72.json
    └── ...
```

The `data/` output should be treated as generated artifact by default. Commit the exporter and documentation first; commit dataset artifacts only if explicitly requested.

## Implementation tasks

### Task 1: Add a reusable exporter script

Create `scripts/export_paint_booth_ml_dataset.py`.

Requirements:

1. Read completed OpenFOAM cases from a root directory.
2. Locate latest `foamToVTK` time directory and `internal.vtu`.
3. Extract cell centers and available cell fields using PyVista.
4. Write one compressed `.npz` sample per case.
5. Write one JSON metadata file per case.
6. Write dataset-level `manifest.json` and `manifest.csv`.
7. Provide `--dry-run` for metadata-only discovery without writing large samples.

### Task 2: Export the existing three QA cases as a seed dataset

Command pattern:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 scripts/export_paint_booth_ml_dataset.py \
    --case-root cases/paint_booth_l2_flow_qa \
    --output-root data/paint_booth_ml_dataset/l2_flow_qa_seed_v0 \
    --include "l2_flow_qa_*" \
    --force'
```

Verification:

- `manifest.json` lists 3 samples.
- Each sample has the same number of cells as the corresponding mesh summary.
- Required arrays exist: `cell_centers`, `U`, `p`.
- Metadata includes `qa_pass_basic: true` for all three seed cases.

### Task 3: Add a moderate sweep runner mode

Extend or wrap `scripts/run_paint_booth_l2_flow_qa.py` so it can run the 9-point dataset sweep without changing the proven QA defaults.

Recommended command:

```bash
python3 scripts/run_paint_booth_l2_flow_qa.py \
  --root cases/paint_booth_l2_u_sweep_v0 \
  --supply 2.18 2.72 3.27 3.81 4.36 4.91 5.45 6.00 6.54
```

Because this is long-running CFD work, run it as a tracked/background job and notify on completion.

### Task 4: Export the 9-point sweep

After all CFD cases finish and pass QA, export them with the same exporter:

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

### Task 5: Document dataset card

Create a dataset card after the seed export and update it after the 9-point sweep.

Suggested path:

```text
docs/26-05-03_paint_booth_l2_steady_dataset_card.md
```

Minimum contents:

- Purpose and intended ML use
- CFD baseline configuration
- Parameter range
- Sample representation
- QA criteria and pass/fail summary
- Known limitations:
  - steady RANS only
  - not physically calibrated production booth data
  - one geometry and one filter coefficient
  - high-flow `y+` max outliers tracked separately
  - near-car reverse-flow fraction remains an important QA feature

## Acceptance criteria for v0

- Exporter can process the existing 3 QA cases without rerunning CFD.
- Seed dataset has 3 `.npz` samples and matching metadata.
- Manifest contains enough information to reproduce each sample from the original case.
- All exported seed cases preserve QA pass/fail status from the CFD postprocessing.
- No generated dataset artifacts are committed unless explicitly requested.

## Next decision after v0

Once the 9-point one-parameter sweep is exported, decide whether to expand along:

1. `U_supply` only with denser sampling,
2. `U_supply × filter_forchheimer` calibration grid,
3. geometry perturbations / vehicle pose changes,
4. transient or history-conditioned data for controller dynamics.
