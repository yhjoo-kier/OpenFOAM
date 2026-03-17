# Pipeline map

## Canonical project root
- `/home/yhjoo/projects/OpenFOAM`

## Goal of the framework

This framework turns one of the following inputs:
- indoor photo
- rendering / perspective view
- text scene description
- existing scene JSON
- partially prepared mesh or case

into one or more of these outputs:
- simulation-ready indoor scene abstraction (`indoor_cfd_scene_v1` JSON)
- Gmsh geometry / mesh
- OpenFOAM case and stabilized solve
- VTK outputs
- comparison and paper-facing figures

## Canonical entrypoints

### Primary end-to-end entrypoint
- `scripts/run_indoor_stabilized.py`
- Use this first for most end-to-end or retry-aware runs.
- Current responsibilities:
  - Gemini scene generation from text or image
  - scene repair fallback
  - Gmsh meshing
  - OpenFOAM case generation
  - mesh-risk-aware preset selection
  - `simpleFoam` solve attempts
  - VTK export
  - 2D comparison figure generation
  - 3D comparison figure generation
  - stabilization summary generation

### Prototype / lower-level entrypoint
- `scripts/run_indoor_pipeline.py`
- Simpler prototype path without the full stabilization ladder.
- Useful for lower-level debugging or stage-isolated execution.

## Stage-specific scripts

### Scene generation / abstraction
- `scripts/generate_indoor_scene_with_gemini.py`
- Supports:
  - `--backend api|cli`
  - image-assisted scene generation
  - model fallback tracking
  - generation summary output

### Scene validation / repair
- `scripts/validate_indoor_scene.py`
- `scripts/repair_indoor_scene.py`

### Geometry and mesh generation
- `scripts/scene_to_gmsh.py`

### OpenFOAM case creation
- `scripts/create_indoor_openfoam_case.py`

### Visualization / figure production
- `scripts/visualize_indoor_case.py`
- `scripts/render_indoor_pipeline_3d.py`

## Canonical output locations

### Generated artifacts
- `generated/<run_name>.json`
- `generated/<run_name>_repaired.json`
- `generated/<run_name>.geo_unrolled`
- `generated/<run_name>.msh`

### Case artifacts
- `cases/<run_name>/...`

### Results artifacts
- `results/<run_name>/stabilization_summary.json`
- `results/<run_name>/comparison_1x2.png`
- `results/<run_name>/log.simpleFoam_*`
- `cases/<run_name>/VTK/...`

## Code-reading shortcut

If an agent only has a minute to orient itself, read in this order:
1. `scripts/run_indoor_stabilized.py`
2. `scripts/generate_indoor_scene_with_gemini.py`
3. `scripts/scene_to_gmsh.py`
4. `scripts/create_indoor_openfoam_case.py`
5. `docs/26-03-15_openfoam_solver_stabilization_compass.md`
