# Benchmark Autopilot TODO

> Topic: SurfClaw / topic 2242
> Mode: autonomous progress with user-intervention only when a real decision is needed
> Cron: `openfoam-benchmark-autopilot-topic2242`
> Job ID: `543fca81-0332-442a-beec-e8ef06e5e73d`
> Last initialized: 2026-03-18

## Operating rule

- Continue autonomously unless a genuine user decision is required.
- If a blocking decision is required, disable the cron job first, then ask the user.
- Keep status compact and actionable.

## Current objective

Build and validate the benchmark dataset pipeline for the OpenFOAM Image-to-CFD paper:
1. rule-based benchmark scene generation
2. reference CFD generation
3. result visualization
4. 3D-to-2D rendering feasibility
5. success/failure analysis and stabilization

## Status legend

- [ ] pending
- [x] completed
- [-] in progress
- [!] failed / needs fix
- [~] on hold

## Active TODO

### A. Benchmark scene generation
- [x] Implement `scripts/generate_benchmark_scenes.py`
- [x] Generate initial 8 benchmark scenes (A1/A2/A3/A4 × 2)
- [x] Validate all generated JSON scenes with `validate_indoor_scene.py`
- [x] Improve generator so simple categories can optionally guarantee at least 1 obstacle when desired
- [x] Add category/seed/geometry manifest for generated scenes
- [x] Expand the frozen scene set from 8 to 12 scenes (A1/A2/A3/A4 × 3)
- [x] Expand the frozen scene set from 12 to 16 scenes (A1/A2/A3/A4 × 4)

### B. Reference CFD pipeline
- [x] Reuse stabilized pipeline entrypoint for scene-JSON-driven benchmark runs
- [x] Implement batch wrapper: `scripts/run_benchmark_reference_batch.py`
- [x] Fix `run_indoor_stabilized.py` so `--end-time` propagates into case generation
- [x] Confirm at least one rectangular benchmark scene reaches CFD + visualization end-to-end
- [x] Confirm at least one composite benchmark scene reaches CFD + visualization end-to-end
- [x] Raise first-pass success rate across the 8-scene pilot set
- [x] Run the newly added 4 scene variants (`*_03`) through the reference CFD batch flow

### C. Composite-room compatibility
- [x] Add `room.blocks` support to validator and Gmsh path
- [x] Confirm image-based composite-room generation feasibility
- [x] Identify legacy repair-stage incompatibility with `room.blocks`
- [x] Prevent repair-stage crash from blocking composite runs in `run_indoor_stabilized.py`
- [x] Update `visualize_indoor_case.py` for composite rooms
- [x] Update `render_indoor_pipeline_3d.py` for composite rooms
- [x] Ensure successful repaired-scene runs render against the actual scene JSON used by the solver
- [-] Add proper composite-room support to `repair_indoor_scene.py` instead of fallback/skip behavior

### D2. Evaluation scaffolding
- [x] Create canonical `benchmark/evaluations/` scaffold for the frozen-12 × 5-view bundle
- [x] Materialize per-task manifests/symlinks for 60 image-conditioned evaluation tasks
- [x] Refresh the evaluation scaffold to the frozen-16 × 5-view bundle (80 tasks) while preserving prior task state
- [x] Add a runner that executes one scaffolded evaluation task end-to-end and writes `evaluation_summary.json`
- [~] Smoke-test one image-conditioned task in the current cron environment (`bench_a1_01/perspective`) — runner now classifies this correctly as an environment block (`status=blocked`) because `GEMINI_API_KEY` is not set for the API backend

### D. Failure analysis / stabilization
- [x] Batch-run the initial 8-scene pilot set
- [x] Summarize current success/failure status
- [x] Analyze failure root cause for `bench_a1_02`
- [x] Analyze failure root cause for `bench_a3_01`
- [x] Analyze failure root cause for `bench_a3_02`
- [x] Decide whether failures are best fixed at generator-level, repair-level, meshing-level, or solver-level
- [x] Re-run failed pilot scenes after fixes
- [x] Update pilot success-rate summary after re-runs

### E. Rendering / dataset feasibility
- [x] Confirm 3D comparison rendering is produced for successful benchmark runs
- [x] Define canonical rendering output structure for dataset cases
- [x] Add automated 2D rendering export path for benchmark scenes
- [x] Separate render modes needed for paper benchmark (perspective / bird's-eye / floor-plan / wireframe)
- [x] Add canonical `section` view export so the benchmark view set matches the paper-plan 5-view matrix
- [x] Confirm one case can produce both reference CFD and benchmark input image assets cleanly
- [x] Verify the integrated batch flow on the newly added `*_03` subset
- [x] Refresh the frozen-12 rendering manifests to include `section` assets for all scenes

### F. Documentation / records
- [x] Write generator progress note: `docs/26-03-18_benchmark_scene_generator.md`
- [x] Create this autopilot todo file
- [x] Write pilot benchmark status note with success/failure table
- [x] Keep manifest files updated after each major batch rerun
- [x] Write a follow-up note for the 12-scene frozen subset once the new 4 reference runs finish
- [x] Write a follow-up note for the 16-scene frozen subset and manifest-refresh fixes
- [-] Commit meaningful project changes in project repo when a stable checkpoint is reached (do not push unless explicitly requested)

## Current pilot status snapshot

Recovered 8-scene pilot summary:
- Success: `bench_a1_01`, `bench_a1_02`, `bench_a2_01`, `bench_a2_02`, `bench_a3_01`, `bench_a3_02`, `bench_a4_01`, `bench_a4_02`
- Success rate: 8/8 reference CFD

Current expansion status:
- Newly frozen additional scenes for the 12-scene tranche: `bench_a1_03`, `bench_a2_03`, `bench_a3_03`, `bench_a4_03`
- Newly frozen additional scenes for the 16-scene tranche: `bench_a1_04`, `bench_a2_04`, `bench_a3_04`, `bench_a4_04`
- Benchmark input renderings for the `*_03` and `*_04` scenes: done
- Reference CFD runs for the `*_03` scenes: done (`12/12` frozen subset recovered)
- Reference CFD runs for the `*_04` scenes: done (`16/16` frozen subset recovered)
- Notable recoveries in the 12-scene tranche: `bench_a2_03` and `bench_a4_03` converged after falling back from `mesh_size=0.35` to `0.25`
- 16-scene tranche note: all `*_04` scenes solved at the initial `0.35` mesh; `bench_a1_04` needed the `ultra_robust` preset, while `bench_a2_04`, `bench_a3_04`, and `bench_a4_04` converged with `robust`

Notes:
- The earlier pilot failures were traced to meshing/repair compatibility regressions rather than flawed benchmark categories.
- The batch reference runner now also exports benchmark input views, so the reference-CFD path and input-image path are no longer tracked separately.
- During the 16-scene expansion, a manifest bookkeeping regression was found: aggregate rendering/evaluation scaffolds were drifting because `renderings_manifest.json` and downstream scaffolding were tied to the most recent batch instead of the full artifact set. This was fixed by regenerating aggregate manifests from on-disk case manifests and by making the scaffold runner consume a generic aggregate reference-status manifest.

## Next default actions

1. Re-run one or more scaffolded evaluation tasks once the Gemini API backend environment is available in the runtime shell.
2. Extend the new geometry-aware evaluation summaries toward CFD/result-side metrics once image-conditioned runs are available.
3. Expand the frozen benchmark generator toward the next tranche (20 scenes) and watch whether the current stabilization ladder still holds.
4. Take a clean local commit checkpoint for the 16-scene + manifest-refresh state.
nes) and watch whether the current stabilization ladder still holds.
4. Take a clean local commit checkpoint for the 16-scene + manifest-refresh state.
