# L2 paint-booth flow-rate QA sweep for ML dataset readiness

Date: 2026-05-02

## Objective

Before generating many CFD snapshots for GINO/neural-operator training, run a small representative flow-rate QA sweep using the current L2 kOmegaSST + carBody boundary-layer paint-booth baseline.

The purpose is not final physical calibration.  The purpose is to check whether the current CFD baseline is numerically stable and produces physically consistent trends when the supply flow rate is changed.

## Case configuration

- Generator: `scripts/create_paint_booth_plenum_filter_case.py`
- Runner/aggregator: `scripts/run_paint_booth_l2_flow_qa.py`
- Case root: `cases/paint_booth_l2_flow_qa/`
- Turbulence model: `kOmegaSST`
- Central filter panel + sealed high-resistance frame: enabled
- Filter Forchheimer coefficient: `10000`
- Base cell size: `0.125 m`
- Filter z cells: `6`
- carBody refinement: `3/4`
- carBody prism layers: `5`
- Layer sizing: absolute, final thickness `0.004 m`, min thickness `0.0003 m`, expansion ratio `1.2`
- Runtime: Docker image `openfoam-pipeline-local:latest`

Flow-rate samples:
- `l2_flow_qa_u2p18`: U_supply = `2.18 m/s` (0.5x nominal), Q_in = `3.9240 m^3/s`
- `l2_flow_qa_u4p36`: U_supply = `4.36 m/s` (1.0x nominal), Q_in = `7.8480 m^3/s`
- `l2_flow_qa_u6p54`: U_supply = `6.54 m/s` (1.5x nominal), Q_in = `11.7720 m^3/s`

## QA summary

- All basic QA pass: `True`
- All cases reached `simpleFoam` final time 200 and reported `End`.
- All meshes had `checkMesh` OK with identical cell count because only boundary velocity changed.
- Mass imbalance remained below 0.1% for all three cases.
- carBody y+ p95 stayed below 5 for all three cases.

## Flow trend table

| U_supply [m/s] | Q_in [m3/s] | imbalance [%] | filter dp proxy [Pa] | work Uz mean [m/s] | near-car Uz mean [m/s] | work reverse [%] | near-car reverse [%] | y+ median | y+ p95 | y+ max | QA |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 2.18 | 3.9240 | -0.0660 | 16.71 | -0.0932 | -0.0656 | 11.69 | 24.59 | 0.135 | 0.504 | 9.961 | True |
| 4.36 | 7.8480 | -0.0626 | 66.84 | -0.2000 | -0.1435 | 9.83 | 23.27 | 0.529 | 1.737 | 23.496 | True |
| 6.54 | 11.7720 | -0.0633 | 150.39 | -0.3103 | -0.2246 | 9.04 | 22.74 | 1.194 | 3.357 | 36.367 | True |

## Monotonic trend checks

- `supply_inflow_increases_with_supply_velocity`: `True`
- `filter_dp_increases_with_supply_velocity`: `True`
- `abs_work_zone_downdraft_increases_with_supply_velocity`: `True`
- `abs_near_car_downdraft_increases_with_supply_velocity`: `True`
- `work_zone_Uz_more_negative_with_supply_velocity`: `True`
- `near_car_Uz_more_negative_with_supply_velocity`: `True`

Interpretation: supply inflow, filter pressure-drop proxy, work-zone downdraft magnitude, and near-car downdraft magnitude all changed monotonically in the expected direction over the 0.5x/1.0x/1.5x sweep.

## Important caveats

- This confirms numerical/trend sanity for a steady RANS baseline; it does not yet validate real booth physics.
- The high-flow case has carBody y+ max above 30 (`36.37`) on a very small outlier fraction, while p95 remains acceptable (`3.36`). This is not a blocker for baseline QA, but it should be tracked when expanding to higher flows.
- Reverse-flow fraction near the car decreases only mildly with flow rate and remains around 23-25%, so this metric should remain a key dataset QA feature.
- The current data are suitable for a quasi-steady flow-map surrogate first. True transient response requires transient CFD or history-conditioned surrogate data later.

## Decision

The current L2 CFD setup passes the three-point flow-rate QA gate for generating a first steady CFD dataset sweep for neural-operator/GINO development, with the caveat that the dataset should be labeled as a CFD-consistent surrogate baseline rather than physically calibrated production booth data.

Recommended next step: implement the dataset exporter and run a moderate one-parameter U_supply sweep, while writing per-case QA metadata alongside every ML sample.
