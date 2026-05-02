# Paint-booth L2 steady CFD dataset card — seed v0

Date: 2026-05-03

## Dataset identity

- Dataset name: `l2_flow_qa_seed_v0`
- Dataset path: `data/paint_booth_ml_dataset/l2_flow_qa_seed_v0/`
- Source cases: `cases/paint_booth_l2_flow_qa/l2_flow_qa_*`
- Exporter: `scripts/export_paint_booth_ml_dataset.py`
- Export command:

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

## Purpose

This is the first ML-ready export of the validated L2 paint-booth CFD baseline. It is intended as a seed dataset for checking neural-operator/GINO data loading, feature normalization, graph/point-cloud batching, and simple surrogate training loops before launching a larger CFD sweep.

This dataset is not a physically calibrated production booth dataset. It is a CFD-consistent steady RANS baseline dataset.

## CFD baseline

- Solver: OpenFOAM `simpleFoam`
- Turbulence model: `kOmegaSST`
- Geometry: simplified car shell in paint-booth work zone
- Inflow model: supply jet into upper plenum
- Filter model: central filter panel + sealed high-resistance frame using porous resistance
- Filter Forchheimer coefficient: `10000`
- Base cell size: `0.125 m`
- Filter z cells: `6`
- carBody refinement: `3/4`
- carBody prism layers: `5`
- Absolute final layer thickness: `0.004 m`
- Runtime image: `openfoam-pipeline-local:latest`

## Samples

- Sample count: `3`
- Supply velocities:
  - `2.18 m/s` (`0.5x` nominal)
  - `4.36 m/s` (`1.0x` nominal)
  - `6.54 m/s` (`1.5x` nominal)
- Cells per sample: `224,769`
- Representation: cell-centered OpenFOAM `foamToVTK` point cloud / graph sample

Each `.npz` sample contains:

- `cell_centers`: shape `(224769, 3)`
- `U`: shape `(224769, 3)`
- `p`: shape `(224769,)`
- optional exported fields when present: `k`, `omega`, `nut`, `yPlus`

Metadata files are stored under:

```text
data/paint_booth_ml_dataset/l2_flow_qa_seed_v0/metadata/
```

Dataset manifest files:

```text
data/paint_booth_ml_dataset/l2_flow_qa_seed_v0/manifest.json
data/paint_booth_ml_dataset/l2_flow_qa_seed_v0/manifest.csv
```

## QA summary

All three seed samples preserve the prior L2 flow-rate QA pass state.

- `l2_u2p18`
  - `U_supply`: `2.18 m/s`
  - cells: `224,769`
  - QA pass: `true`
  - relative imbalance: `-0.0006598`
  - work-zone `Uz_mean`: `-0.0932 m/s`
  - near-car `Uz_mean`: `-0.0656 m/s`
  - carBody `y+ p95`: `0.504`
  - carBody `y+ max`: `9.961`

- `l2_u4p36`
  - `U_supply`: `4.36 m/s`
  - cells: `224,769`
  - QA pass: `true`
  - relative imbalance: `-0.0006260`
  - work-zone `Uz_mean`: `-0.2000 m/s`
  - near-car `Uz_mean`: `-0.1435 m/s`
  - carBody `y+ p95`: `1.737`
  - carBody `y+ max`: `23.496`

- `l2_u6p54`
  - `U_supply`: `6.54 m/s`
  - cells: `224,769`
  - QA pass: `true`
  - relative imbalance: `-0.0006330`
  - work-zone `Uz_mean`: `-0.3103 m/s`
  - near-car `Uz_mean`: `-0.2246 m/s`
  - carBody `y+ p95`: `3.357`
  - carBody `y+ max`: `36.367`

Verification command used after export:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 - <<"PY"
import json, numpy as np
from pathlib import Path
root=Path("data/paint_booth_ml_dataset/l2_flow_qa_seed_v0")
manifest=json.loads((root/"manifest.json").read_text())
print("sample_count", manifest["sample_count"])
for s in manifest["samples"]:
    arr=np.load(root/s["sample_path"])
    print(s["sample_id"], s["n_cells"], arr["cell_centers"].shape, arr["U"].shape, arr["p"].shape, s["qa_pass_basic"])
PY'
```

Observed verification:

```text
sample_count 3
l2_u2p18 n_cells 224769 centers (224769, 3) U (224769, 3) p (224769,) qa True
l2_u4p36 n_cells 224769 centers (224769, 3) U (224769, 3) p (224769,) qa True
l2_u6p54 n_cells 224769 centers (224769, 3) U (224769, 3) p (224769,) qa True
```

## Known limitations

- Steady RANS only; no transient controller dynamics yet.
- One vehicle geometry, one mesh topology, and one filter resistance.
- Only three supply-velocity points; enough for pipeline validation, not enough for production training.
- High-flow case has high `y+ max` outliers, although `y+ p95` passes the current QA threshold.
- Near-car reverse-flow fraction remains substantial and should be retained as a dataset QA feature.
- Dataset artifacts under `data/` are generated outputs and should not be committed unless explicitly requested.

## Recommended next step

Run a moderate 9-point `U_supply` sweep using the same baseline settings:

```bash
python3 scripts/run_paint_booth_l2_flow_qa.py \
  --root cases/paint_booth_l2_u_sweep_v0 \
  --supply 2.18 2.72 3.27 3.81 4.36 4.91 5.45 6.00 6.54
```

After completion, export it to:

```text
data/paint_booth_ml_dataset/l2_steady_u_sweep_v0/
```
