# Paint-booth transient PoC implementation plan

> **For Hermes:** Use OpenFOAM CFD workflow discipline and keep this as a small PoC before launching large transient datasets.

**Goal:** Convert the validated L2 steady paint-booth baseline into a minimal `pimpleFoam` step-response case and verify that transient execution, metrics, and future dataset export are feasible.

**Architecture:** Reuse an already-converged L2 steady case as the initial condition. Copy the latest steady fields into `0/`, switch the solver from `simpleFoam` to `pimpleFoam`, change the `supplyInlet` velocity to a new target value at `t=0`, and run a short transient smoke test before longer 60–120 s cases.

**Tech Stack:** OpenFOAM `pimpleFoam`, `kOmegaSST`, existing Docker image `openfoam-pipeline-local:latest`, Python case-prep script, existing PyVista postprocessing/export pattern.

---

## Context

Current steady CFD-to-ML state:

- Validated 9-point steady L2 `U_supply` sweep exists under `cases/paint_booth_l2_u_sweep_v0/`.
- ML-ready steady dataset exists under `data/paint_booth_ml_dataset/l2_steady_u_sweep_v0/` and Google Drive handoff package.
- All 9 steady cases pass basic QA, but low-flow `2.72 m/s` has a local metric kink.

Transient goal is not immediate production-scale data generation. The first goal is to answer:

1. Can the current mesh/model run stably with `pimpleFoam`?
2. Can a step change in `U_supply` produce meaningful `Uz(t)` / region response metrics?
3. What time scale and write interval are practical before building history-conditioned datasets?

## PoC case definition

Initial steady source case:

```text
cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u4p36
```

PoC output case:

```text
cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke
```

Step response:

```text
initial condition: steady field at U_supply = 4.36 m/s
transient boundary: supplyInlet U = (0 0 -5.45) m/s from t = 0
```

Smoke-test runtime:

```text
endTime = 5 s
initial deltaT = 0.02 s
adjustTimeStep = yes
maxCo = 1.0
maxDeltaT = 0.1 s
writeInterval = 1 s
```

Longer follow-up after smoke test:

```text
endTime = 60..120 s
writeInterval = 1..2 s
```

## Task 1: Create transient case-prep script

**Objective:** Add a script that prepares a `pimpleFoam` transient case from an existing steady case.

**Files:**

- Create: `scripts/create_paint_booth_transient_poc.py`

**Requirements:**

1. Accept `--source-case`, `--case-dir`, `--target-supply-velocity`, `--end-time`, `--delta-t`, `--max-delta-t`, `--write-interval`, `--force`.
2. Copy required OpenFOAM directories from source case:
   - `constant/`
   - `system/`
   - latest numeric time directory fields into destination `0/`
3. Exclude generated heavy outputs from the copied case root:
   - `VTK/`
   - `postProcessing/`
   - `log.*`
4. Rewrite `system/controlDict` for `pimpleFoam` transient execution.
5. Rewrite `system/fvSchemes` so `ddtSchemes` uses `Euler`, not `steadyState`.
6. Rewrite `system/fvSolution` with a `PIMPLE` block.
7. Rewrite `0/U` `supplyInlet` fixedValue to the target velocity.
8. Write metadata JSON describing source, target velocity, and run settings.

## Task 2: Validate generated dictionaries without running long CFD

**Objective:** Ensure the generated case is structurally valid.

**Commands:**

```bash
python3 scripts/create_paint_booth_transient_poc.py \
  --source-case cases/paint_booth_l2_u_sweep_v0/l2_flow_qa_u4p36 \
  --case-dir cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke \
  --target-supply-velocity 5.45 \
  --end-time 5 \
  --delta-t 0.02 \
  --max-delta-t 0.1 \
  --write-interval 1 \
  --force
```

Then inside Docker:

```bash
checkMesh -constant
pimpleFoam -help
```

Expected:

- `checkMesh` succeeds with mesh OK.
- `pimpleFoam` exists in the Docker runtime.

## Task 3: Run short smoke transient

**Objective:** Prove `pimpleFoam` can advance the case for a short physical time.

**Command:**

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD/cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke" \
  openfoam-pipeline-local:latest \
  bash -lc 'checkMesh -constant > log.checkMesh 2>&1 && pimpleFoam > log.pimpleFoam 2>&1 && foamToVTK -latestTime > log.foamToVTK 2>&1'
```

Expected:

- `log.pimpleFoam` contains `End`.
- Latest time is approximately `5`.
- Continuity errors remain bounded.
- VTK export exists for the latest time.

## Task 4: Postprocess smoke response

**Objective:** Extend existing scalar-region postprocessing to read transient latest time first, then later add time-series metrics.

For smoke stage, reuse latest-time postprocessing:

```bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  openfoam-pipeline-local:latest \
  bash -lc 'python3 scripts/postprocess_paint_booth_plenum_metrics.py cases/paint_booth_l2_transient_step_poc/u4p36_to_u5p45_smoke'
```

Expected:

- `post_plenum_metrics.json` is created.
- Work-zone and near-car metrics are finite.

## Acceptance criteria for this PoC

- Generated case can run `checkMesh -constant`.
- `pimpleFoam` reaches `End` for a 5 s smoke test.
- Latest-time VTK and scalar metrics can be created.
- The result is documented with explicit caveats.
- Large transient case outputs remain untracked unless explicitly requested.

## Next step after smoke pass

If the smoke case passes, run a 60–120 s step-response case and add time-series postprocessing over all written times:

- `supplyInlet` flow rate vs time
- `floorExhaust` flow rate vs time
- work-zone `Uz_mean(t)`
- near-car `Uz_mean(t)`
- reverse-flow fraction vs time
- response time estimate: time to 63%, 90%, and 95% of final metric change

These quantities will guide transient snapshot spacing and history-window length for future neural-operator/controller datasets.
